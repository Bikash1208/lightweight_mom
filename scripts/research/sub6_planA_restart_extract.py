"""Open one finished batch project in a private AEDT instance and export Stage 1 data."""
import json
import math
import re
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from ansys.aedt.core import Hfss


def utc():
    return datetime.now(timezone.utc).isoformat()


def ffd_data(path):
    with path.open() as f:
        theta = [float(x) for x in f.readline().split()]
        phi = [float(x) for x in f.readline().split()]
        assert len(theta) == 3 and len(phi) == 3
        f.readline()
        f.readline()
        v = np.loadtxt(f)
    nt, np_ = int(theta[2]), int(phi[2])
    assert v.shape == (nt * np_, 4), (v.shape, nt, np_)
    et = (v[:, 0] + 1j * v[:, 1]).reshape((nt, np_))
    ep = (v[:, 2] + 1j * v[:, 3]).reshape((nt, np_))
    th = np.linspace(theta[0], theta[1], nt)
    ph = np.linspace(phi[0], phi[1], np_)
    return th, ph, et, ep


def point(h, expression, theta=None, phi=None):
    kwargs = {"expressions": expression, "setup_sweep_name": "Setup1 : LastAdaptive"}
    if theta is not None:
        kwargs.update(report_category="Far Fields", context="Global",
                      variations={"Freq": ["5.8GHz"], "Theta": [f"{theta}deg"],
                                  "Phi": [f"{phi}deg"]})
    data = h.post.get_solution_data(**kwargs)
    assert data, (expression, theta, phi)
    _, values = data.get_expression_data(expression=expression, formula="real")
    assert len(values) == 1 and math.isfinite(float(values[0])), (expression, values)
    return float(values[0])


def main():
    project = Path(sys.argv[1]).resolve()
    out = Path(sys.argv[2]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    report = out / f"{project.stem}.json"
    d = {"run_id": project.stem, "project": str(project), "status": "EXTRACTING",
         "started_utc": utc(), "peak_batch_ram_gb": float(sys.argv[3]),
         "batch_seconds": float(sys.argv[4])}
    report.write_text(json.dumps(d, indent=2))
    h = None
    try:
        h = Hfss(project=str(project), design="AnchorPreflight", solution_type="Modal",
                 version="2025.2", non_graphical=True, new_desktop=True,
                 close_on_exit=False)
        d["aedt_version"] = str(h.desktop_class.aedt_version_id)
        d["own_aedt_pid"] = h.desktop_class.aedt_process_id
        assert d["aedt_version"] == "2025.2"
        assert "Setup1 : LastAdaptive" in h.existing_analysis_sweeps, "NO_SOLVE_IN_PROJECT"
        d["profile_file"] = str(h.export_profile("Setup1", output_file=str(out / "solve.prof")))
        d["convergence_file"] = str(h.export_convergence("Setup1", output_file=str(out / "solve.conv")))
        far = out / "farfield"
        far.mkdir(exist_ok=True)
        ok = h.export_antenna_metadata(frequencies=[5.8e9], setup="Setup1 : LastAdaptive",
                                       sphere="Global", output_dir=str(far),
                                       export_element_pattern=True, export_objects=False,
                                       export_touchstone=True, export_power=True)
        assert ok, "FFD_EXPORT_FAILED"
        ffd = list(far.rglob("*.ffd"))
        assert len(ffd) == 1, f"FFD_COUNT_{len(ffd)}"
        d["ffd_file"] = str(ffd[0])
        d["s11_db"] = point(h, "dB(S(Port1,Port1))")
        for theta, phi in ((25, 0), (25, 180), (24, 180)):
            key = f"th{theta}_ph{phi}"
            d[key] = {"gain_rhcp_db": point(h, "dB(GainRHCP)", theta, phi),
                      "gain_lhcp_db": point(h, "dB(GainLHCP)", theta, phi)}
        th, ph, et, ep = ffd_data(ffd[0])
        assert th[0] == 0 and th[-1] == 180 and ph[0] == -180 and ph[-1] == 180
        intensity = np.abs(et)**2 + np.abs(ep)**2
        weights = np.sin(np.deg2rad(th))[:, None]
        # The final phi column duplicates the first and is omitted in the integral.
        total = float(np.sum(intensity[:, :-1] * weights))
        below = float(np.sum(intensity[th > 90, :-1] * weights[th > 90]))
        assert total > 0
        d["radiated_power_fraction_theta_gt_90"] = below / total
        rh_plus = np.abs(et + 1j * ep)**2
        rh_minus = np.abs(et - 1j * ep)**2
        i24 = int(np.argmin(np.abs(th - 24)))
        i0 = int(np.argmin(np.abs(ph - 0)))
        db_r = point(h, "dB(GainRHCP)", 24, 0)
        db_l = point(h, "dB(GainLHCP)", 24, 0)
        plus_ratio = 10 * math.log10(max(rh_plus[i24, i0], 1e-300) / max(rh_minus[i24, i0], 1e-300))
        handed = "plus" if abs(plus_ratio - (db_r-db_l)) < abs(-plus_ratio - (db_r-db_l)) else "minus"
        rh = rh_plus if handed == "plus" else rh_minus
        peak_idx = np.unravel_index(int(np.argmax(rh[:, :-1])), rh[:, :-1].shape)
        peak_th, peak_ph = float(th[peak_idx[0]]), float(ph[peak_idx[1]])
        d["peak_rhcp"] = {"theta_deg": peak_th, "phi_deg": peak_ph,
                          "gain_db": point(h, "dB(GainRHCP)", peak_th, peak_ph),
                          "ffd_handedness": handed}
        d["exported_files"] = [str(p) for p in out.rglob("*") if p.is_file() and p != report]
        d["status"] = "COMPLETE"
    except BaseException:
        d["status"] = "FAILED"
        d["traceback"] = traceback.format_exc()
        raise
    finally:
        d["ended_utc"] = utc()
        report.write_text(json.dumps(d, indent=2, default=str))
        if h is not None:
            h.release_desktop(close_projects=True, close_desktop=True)


if __name__ == "__main__":
    main()
