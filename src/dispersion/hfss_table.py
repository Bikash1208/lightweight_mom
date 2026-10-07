"""Part 3, HMSA item 1: import of an HFSS eigenmode dispersion table (geometry parameter vs beta
or surface reactance X_s at the design frequency) as a library type
(`labels.DISPERSION_FROM_HFSS`) -- usable now, since this project's own MoM has no eigenmode
solver yet (`mom_scaffold.py`).

CSV schema (one row per geometry sample, at a single fixed design frequency the caller states):
    geometry_param, beta_rad_per_m, x_s_ohm
`geometry_param` is whatever the library varies (e.g. patch length in mm) -- its name/units are
the caller's responsibility to track; this module treats it as an opaque, monotonically-useful
sort key. At least one of `beta_rad_per_m`/`x_s_ohm` must be present in every row; a completely
empty table, a non-numeric value, or a duplicate `geometry_param` is a `CONFIG_INVALID` error
(no silent fallback).
"""

import csv
from bisect import bisect_left
from dataclasses import dataclass
from pathlib import Path

from ..errors import config_invalid
from ..labels import DISPERSION_FROM_HFSS

REQUIRED_COLUMNS = ("geometry_param",)
OPTIONAL_COLUMNS = ("beta_rad_per_m", "x_s_ohm")


@dataclass(frozen=True)
class DispersionSample:
    geometry_param: float
    beta_rad_per_m: "float | None"
    x_s_ohm: "float | None"


@dataclass(frozen=True)
class DispersionTable:
    source_path: str
    frequency_hz: float
    samples: tuple  # DispersionSample, sorted by geometry_param
    label: str = DISPERSION_FROM_HFSS

    @property
    def geometry_range(self):
        return self.samples[0].geometry_param, self.samples[-1].geometry_param

    def beta_at(self, geometry_param):
        """Linear interpolation in `geometry_param`; raises `errors.LwmomError`
        (`IMPEDANCE_UNREACHABLE`... no -- out-of-range geometry, not impedance) if outside the
        table's range. Raises `ValueError` if this table has no `beta_rad_per_m` column."""
        return self._interpolate(geometry_param, "beta_rad_per_m")

    def x_s_at(self, geometry_param):
        return self._interpolate(geometry_param, "x_s_ohm")

    def _interpolate(self, geometry_param, field):
        values = [getattr(s, field) for s in self.samples]
        if any(v is None for v in values):
            raise ValueError(f"DispersionTable {self.source_path!r} has no complete {field!r} column")
        lo, hi = self.geometry_range
        if not (lo <= geometry_param <= hi):
            from ..errors import impedance_unreachable

            raise impedance_unreachable(
                f"geometry_param={geometry_param!r} outside this table's range [{lo!r}, {hi!r}]"
            )
        xs = [s.geometry_param for s in self.samples]
        i = bisect_left(xs, geometry_param)
        if i == 0:
            return values[0]
        if i >= len(xs) or xs[i] != geometry_param:
            i -= 1
        if xs[i] == geometry_param or i == len(xs) - 1:
            return values[i]
        x0, x1 = xs[i], xs[i + 1]
        v0, v1 = values[i], values[i + 1]
        t = (geometry_param - x0) / (x1 - x0)
        return v0 + t * (v1 - v0)

    def geometry_for_beta(self, beta_target):
        """Inverse lookup (nearest achievable geometry_param for a target beta), for mapping a
        hologram's required beta(x,y) onto a cell dimension. Returns
        `(geometry_param, achieved_beta, within_range)` -- `within_range=False` flags a value the
        caller should label `IMPEDANCE_UNREACHABLE` rather than silently clamp."""
        return self._inverse("beta_rad_per_m", beta_target)

    def geometry_for_x_s(self, x_s_target):
        return self._inverse("x_s_ohm", x_s_target)

    def _inverse(self, field, target):
        values = [getattr(s, field) for s in self.samples]
        if any(v is None for v in values):
            raise ValueError(f"DispersionTable {self.source_path!r} has no complete {field!r} column")
        lo, hi = min(values), max(values)
        within_range = lo <= target <= hi
        best = min(self.samples, key=lambda s: abs(getattr(s, field) - target))
        return best.geometry_param, getattr(best, field), within_range


def read_dispersion_csv(path, frequency_hz):
    """Reads a dispersion CSV (see module docstring). `frequency_hz` is recorded as metadata (the
    single frequency this table is valid at) -- not read from the file, since an eigenmode table
    is normally exported at one fixed design frequency."""
    path = Path(path)
    try:
        with open(path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            fieldnames = reader.fieldnames or []
            if "geometry_param" not in fieldnames:
                raise config_invalid(message=f"{path}: missing required column 'geometry_param'")
            has_beta = "beta_rad_per_m" in fieldnames
            has_xs = "x_s_ohm" in fieldnames
            if not has_beta and not has_xs:
                raise config_invalid(message=f"{path}: needs at least one of beta_rad_per_m/x_s_ohm")

            samples = []
            seen_params = set()
            for row in reader:
                try:
                    gp = float(row["geometry_param"])
                    beta = (
                        float(row["beta_rad_per_m"])
                        if has_beta and row.get("beta_rad_per_m") not in (None, "")
                        else None
                    )
                    xs = float(row["x_s_ohm"]) if has_xs and row.get("x_s_ohm") not in (None, "") else None
                except ValueError as exc:
                    raise config_invalid(exc, message=f"{path}: non-numeric value in row {row!r}") from exc
                if gp in seen_params:
                    raise config_invalid(message=f"{path}: duplicate geometry_param={gp!r}")
                seen_params.add(gp)
                samples.append(DispersionSample(gp, beta, xs))
    except FileNotFoundError as exc:
        raise config_invalid(exc, message=f"dispersion table not found: {path}") from exc

    if not samples:
        raise config_invalid(message=f"{path}: no data rows")
    samples.sort(key=lambda s: s.geometry_param)
    return DispersionTable(source_path=str(path), frequency_hz=frequency_hz, samples=tuple(samples))
