"""Package-layer STEP clean-up for the sub6 HFSS packages (no change under src/).

Root cause found (see results/sub6_step_fix/STEP_FIX.md): the repo's snap layer `src.transmitarray.export.snap_step_text` snaps CARTESIAN_POINT / DIRECTION / VECTOR literals
within SNAP_TOL_MM = 1e-7 of values listed in the builder's `intended()` table.  Literals it does NOT reach: (a) 2-D parameter-space B-spline control points of the face pcurves,
which the kernel approximates with ~6e-7 error ('-0.0350006', '6.E-7'), i.e. above the 1e-7 tolerance, and (b) mirrored values (-T, -T/2) that were missing from the table;
(c) COLOUR_RGB reals carry float32 noise ('0.250000000672').  `exactify_step_text` removes all three classes at the file level: every numeric literal of the DATA section (except
COLOUR_RGB and UNCERTAINTY entities) within `tol` (2e-6 mm) of a multiple of 1/den (0.001 mm) is replaced by exactly that grid value, colours are rounded to 4 decimals.
Irrational design values (circle/line intersections such as 5.411952050785) are not near the grid and are left untouched.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.transmitarray.export import step_real  # noqa: E402

_NUM = re.compile(r"(?<![\w.])-?\d+\.\d*(?:E[-+]?\d+)?")
_ENTITY = re.compile(r"#\d+\s*=\s*(.*?);(?=\s*(?:#\d+\s*=|ENDSEC))", re.S)


def grid_value(v, den=1000):
    """Nearest multiple of 1/den (default 0.001 mm) as the double a decimal literal 'n/den' parses to."""
    return round(v * den) / den


def _near_grid(v, den, tol):
    g = grid_value(v, den)
    return g if abs(v - g) < tol else None


def _split_data(text):
    i = text.index("DATA;")
    return text[: i + 5], text[i + 5 :]


def _sub_outside_quotes(s, fn):
    parts = s.split("'")
    for k in range(0, len(parts), 2):
        parts[k] = _NUM.sub(fn, parts[k])
    return "'".join(parts)


def exactify_step_text(text, den=1000, tol=2e-6, colour_decimals=4):
    """Returns (new_text, stats).  Stats: literals examined, snapped (changed), max change, colours rounded."""
    head, data = _split_data(text)
    stats = dict(n_literals=0, n_snapped=0, max_change=0.0, n_colours=0)

    def colour(m):
        stats["n_colours"] += 1
        return "COLOUR_RGB(" + _sub_outside_quotes(m.group(1), lambda t: step_real(round(float(t.group(0)), colour_decimals))) + ")"

    data = re.sub(r"COLOUR_RGB\((.*?)\)(?=;)", colour, data, flags=re.S)

    def fix(m):
        v = float(m.group(0))
        stats["n_literals"] += 1
        g = _near_grid(v, den, tol)
        if g is None:
            return m.group(0)
        if g != v:
            stats["n_snapped"] += 1
            stats["max_change"] = max(stats["max_change"], abs(g - v))
        return step_real(g)

    out = []
    pos = 0
    for m in re.finditer(r"#\d+\s*=\s*.*?;(?=\s*(?:#\d+\s*=|ENDSEC))", data, re.S):
        out.append(data[pos : m.start()])
        ent = m.group(0)
        if "UNCERTAINTY_MEASURE_WITH_UNIT" in ent or "COLOUR_RGB" in ent:
            out.append(ent)
        else:
            out.append(_sub_outside_quotes(ent, fix))
        pos = m.end()
    out.append(data[pos:])
    return head + "".join(out), stats


def audit_step_text(text, den=1000, tol=2e-6):
    """Literals (DATA section, outside UNCERTAINTY / quotes) that are within tol of a grid value but are not exactly it, plus colour reals with > 4 decimals."""
    head, data = _split_data(text)
    bad = []
    for m in re.finditer(r"#\d+\s*=\s*.*?;(?=\s*(?:#\d+\s*=|ENDSEC))", data, re.S):
        ent = m.group(0)
        if "UNCERTAINTY_MEASURE_WITH_UNIT" in ent:
            continue
        is_col = "COLOUR_RGB" in ent
        for k, part in enumerate(ent.split("'")):
            if k % 2:
                continue
            for t in _NUM.finditer(part):
                v = float(t.group(0))
                if is_col:
                    if abs(v - round(v, 4)) > 0:
                        bad.append(("COLOUR_RGB", t.group(0)))
                    continue
                g = _near_grid(v, den, tol)
                if g is not None and g != v:
                    bad.append((ent[:30].replace("\n", " "), t.group(0)))
    return bad


def literals(text):
    """All float literals of the DATA section (outside quotes), for tests."""
    head, data = _split_data(text)
    out = []
    for k, part in enumerate(data.split("'")):
        if k % 2 == 0:
            out += [float(t.group(0)) for t in _NUM.finditer(part)]
    return out
