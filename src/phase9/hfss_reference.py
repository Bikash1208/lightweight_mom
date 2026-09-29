"""Sec. 6-9: HFSS reference-cell data. **No HFSS export exists anywhere in this repository**
(grep-verified: no CST/HFSS/measurement data file anywhere in the project -- the same finding
every earlier phase already made). This module defines the ingestion schema and honestly reports
absence (`NOT_AVAILABLE_FROM_CURRENT_DATA`, Rule 1) rather than fabricating a stand-in file.

CSV schema, one row per (cell, frequency), two-port (default):
    cell_id, L_mm, W_mm, frequency_GHz, S11_real, S11_imag, S21_real, S21_imag

Phase 12C: a ONE-PORT variant (grounded cells, `backing="pec"`, have no port 2 -- S21 is
undefined, CONVENTIONS.md) omits the S21_real/S21_imag columns entirely (never a present-but-zero
placeholder, matching `SParameters.grounded`'s own "undefined, not silently zero" rule). Detected
from the header (S21 columns present or absent), not from a caller-supplied flag, so a single
`load_hfss_reference` call handles either file unchanged for the existing two-port case (verified
bit-identical for every existing two-port CSV/caller).
"""

import csv
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

REQUIRED_HFSS_COLUMNS_COMMON = ("cell_id", "L_mm", "W_mm", "frequency_GHz", "S11_real", "S11_imag")
REQUIRED_HFSS_COLUMNS_TWO_PORT_ONLY = ("S21_real", "S21_imag")
REQUIRED_HFSS_COLUMNS = REQUIRED_HFSS_COLUMNS_COMMON + REQUIRED_HFSS_COLUMNS_TWO_PORT_ONLY

DEFAULT_HFSS_REFERENCE_PATH = (
    Path(__file__).resolve().parent.parent.parent / "results" / "hfss_reference" / "hfss_reference_cells.csv"
)


@dataclass(frozen=True)
class HfssRecord:
    cell_id: str
    l_mm: float
    w_mm: float
    frequency_ghz: float
    s11: complex
    s21: "complex | None" = None  # Phase 12C: None for a one-port (grounded) record -- undefined,
    # not zero; every pre-Phase-12C (two-port) record always sets it.


class HfssReferenceError(ValueError):
    pass


def load_hfss_reference(path=DEFAULT_HFSS_REFERENCE_PATH):
    """Returns a tuple of `HfssRecord` if ``path`` exists and validates, or `None` if it does not
    exist -- callers must report `NOT_AVAILABLE_FROM_CURRENT_DATA`, never substitute solver data
    or a fabricated value (Rule 1, Sec. 45's explicit prohibition on 'fabricated HFSS results').

    Two-port (S21 columns present, the pre-Phase-12C schema, unchanged): every `HfssRecord.s21` is
    a complex value, exactly as before. One-port (S21 columns absent, Phase 12C): every
    `HfssRecord.s21` is `None`."""
    path = Path(path)
    if not path.exists():
        return None
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise HfssReferenceError(f"{path}: no data rows")
    header = set(rows[0])
    missing_common = [c for c in REQUIRED_HFSS_COLUMNS_COMMON if c not in header]
    if missing_common:
        raise HfssReferenceError(f"{path}: missing required column(s) {missing_common}")
    has_s21 = all(c in header for c in REQUIRED_HFSS_COLUMNS_TWO_PORT_ONLY)
    partial_s21 = any(c in header for c in REQUIRED_HFSS_COLUMNS_TWO_PORT_ONLY) and not has_s21
    if partial_s21:
        raise HfssReferenceError(
            f"{path}: partial S21 columns present ({[c for c in REQUIRED_HFSS_COLUMNS_TWO_PORT_ONLY if c in header]}); "
            f"give both S21_real and S21_imag (two-port) or neither (one-port, grounded)"
        )
    records = []
    for lineno, row in enumerate(rows, start=2):
        try:
            s21 = complex(float(row["S21_real"]), float(row["S21_imag"])) if has_s21 else None
            records.append(
                HfssRecord(
                    cell_id=row["cell_id"],
                    l_mm=float(row["L_mm"]),
                    w_mm=float(row["W_mm"]),
                    frequency_ghz=float(row["frequency_GHz"]),
                    s11=complex(float(row["S11_real"]), float(row["S11_imag"])),
                    s21=s21,
                )
            )
        except (TypeError, ValueError) as e:
            raise HfssReferenceError(f"{path}:{lineno}: {e}") from None
    return tuple(records)


# ============================================================================
# Phase 12D: generic multi-port Touchstone (.sNp) reader -- a real HFSS export
# now exists (Floquet-port supercells/single cells at the grounded 8 mm design point). This is a
# GENERAL N-port reader (format-line agnostic: any unit prefix, MA/RI/DB, optional reference
# impedance), not hard-coded to the 4- or 36-port cases actually on disk in this project; it never
# assumes a particular line-wrapping convention -- Touchstone files wrap continuation lines at
# whatever width the writer chose (HFSS wraps at 4 complex pairs/line; other tools differ), so this
# reads every non-comment, non-header token into one flat stream and chunks it by the *known* block
# size (1 frequency + N*N complex pairs), which is robust to any wrapping.
# ============================================================================

_TOUCHSTONE_FREQ_UNITS = {"HZ": 1.0, "KHZ": 1e3, "MHZ": 1e6, "GHZ": 1e9}
_TOUCHSTONE_FORMATS = ("MA", "RI", "DB")
_PORT_COMMENT_RE = re.compile(r"^!\s*Port\[(\d+)\]\s*=\s*(\S+)")


class TouchstoneError(ValueError):
    pass


@dataclass(frozen=True)
class TouchstoneNetwork:
    """A parsed N-port Touchstone (.sNp) file: real, measured HFSS data (or any other tool's),
    never fabricated (Rule 1). ``port_labels[k]`` (1-indexed in the file, 0-indexed here) is the
    ``! Port[k] = <name>`` comment text (e.g. ``"FloquetPort1:1"``) if present, else ``None``.
    ``s`` has shape ``(n_freq, n_ports, n_ports)``, ``s[f, i, j]`` = response at port i for
    excitation at port j (Touchstone's own row-major convention: row = receiving port)."""

    path: Path
    n_ports: int
    port_labels: tuple  # length n_ports, each a str or None
    reference_impedance_ohm: "float | None"  # the header's global R, if given (None: default 50 ohm, unused here)
    frequency_hz: np.ndarray  # (n_freq,)
    s: np.ndarray  # (n_freq, n_ports, n_ports) complex

    def port_mode_groups(self):
        """Groups 0-based column/row indices by the port name before the ``:mode`` suffix, e.g.
        ``{"FloquetPort1": (0, 1, ..., 17), "FloquetPort2": (18, ..., 35)}``, in first-seen order.
        Raises if any label is missing or has no ``:`` (Floquet-port files always have one)."""
        groups = {}
        for idx, label in enumerate(self.port_labels):
            if label is None or ":" not in label:
                raise TouchstoneError(
                    f"{self.path}: port {idx + 1} has no '<PortName>:<mode>' comment "
                    f"(got {label!r}); cannot group by Floquet port"
                )
            name = label.rsplit(":", 1)[0]
            groups.setdefault(name, []).append(idx)
        return {name: tuple(idxs) for name, idxs in groups.items()}

    def component(self, i, j):
        """S[i, j] (1-indexed, matching Touchstone port numbering) at every frequency, shape (n_freq,)."""
        if not (1 <= i <= self.n_ports and 1 <= j <= self.n_ports):
            raise TouchstoneError(f"{self.path}: port index out of range 1..{self.n_ports}: ({i}, {j})")
        return self.s[:, i - 1, j - 1]


def read_touchstone(path):
    """Parses a standard N-port Touchstone file (N from the ``.sNp`` extension). Supports the
    ``# <freq_unit> S <MA|RI|DB> [R <ohm>]`` option line, ``! Port[k] = <name>`` port-naming
    comments (any other ``!`` comment, e.g. HFSS's per-frequency ``! Gamma``/``! Port Impedance``
    lines, is ignored -- this reader only needs the S-matrix itself), and any continuation-line
    wrapping width (flat-token chunking, Sec. above). Raises :class:`TouchstoneError` on a
    malformed file; never fabricates or extrapolates missing data (Rule 1)."""
    path = Path(path)
    m = re.search(r"\.s(\d+)p$", path.name, re.IGNORECASE)
    if not m:
        raise TouchstoneError(f"{path}: filename must end in .sNp (N = port count)")
    n_ports = int(m.group(1))

    option_line = None
    port_labels = [None] * n_ports
    data_lines = []
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                option_line = line
                continue
            if line.startswith("!"):
                pm = _PORT_COMMENT_RE.match(line)
                if pm:
                    k = int(pm.group(1))
                    if 1 <= k <= n_ports:
                        port_labels[k - 1] = pm.group(2)
                continue
            data_lines.append(line)

    if option_line is None:
        raise TouchstoneError(f"{path}: missing '#' option line")
    opt = option_line.lstrip("#").split()
    freq_unit = next((t.upper() for t in opt if t.upper() in _TOUCHSTONE_FREQ_UNITS), "GHZ")
    fmt = next((t.upper() for t in opt if t.upper() in _TOUCHSTONE_FORMATS), "MA")
    ref_z = None
    if "R" in [t.upper() for t in opt]:
        r_idx = [t.upper() for t in opt].index("R")
        ref_z = float(opt[r_idx + 1])

    tokens = []
    for line in data_lines:
        tokens.extend(line.split())
    block = 1 + 2 * n_ports * n_ports
    if len(tokens) % block != 0:
        raise TouchstoneError(
            f"{path}: {len(tokens)} data token(s) is not a multiple of the expected "
            f"{block}-token block (1 frequency + {n_ports}x{n_ports} complex pairs) -- "
            f"file truncated or port count wrong"
        )
    n_freq = len(tokens) // block
    arr = np.asarray(tokens, dtype=float).reshape(n_freq, block)
    frequency_hz = arr[:, 0] * _TOUCHSTONE_FREQ_UNITS[freq_unit]
    pairs = arr[:, 1:].reshape(n_freq, n_ports, n_ports, 2)

    if fmt == "MA":
        s = pairs[..., 0] * np.exp(1j * np.deg2rad(pairs[..., 1]))
    elif fmt == "RI":
        s = pairs[..., 0] + 1j * pairs[..., 1]
    else:  # "DB": magnitude in dB, angle in degrees
        s = (10.0 ** (pairs[..., 0] / 20.0)) * np.exp(1j * np.deg2rad(pairs[..., 1]))

    return TouchstoneNetwork(path, n_ports, tuple(port_labels), ref_z, frequency_hz, s)
