"""Extract per-hand full-board gains, S11 and pattern metrics from an owned solved project."""
import argparse
import json
import math
import os
import sys
import traceback
from pathlib import Path

import numpy as np
import psutil

os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
from ansys.aedt.core import Hfss


def ffd(path):
    with path.open() as stream:
        th = [float(v) for v in stream.readline().split()]
        ph = [float(v) for v in stream.readline().split()]
        stream.readline()
        stream.readline()
        matrix = np.loadtxt(stream)
    ntheta, nphi = int(th[2]), int(ph[2])
    assert matrix.shape == (ntheta * nphi, 4)
    theta = np.linspace(th[0], th[1], ntheta)
    phi = np.linspace(ph[0], ph[1], nphi)
    et = (matrix[:, 0] + 1j * matrix[:, 1]).reshape(ntheta, nphi)
    ep = (matrix[:, 2] + 1j * matrix[:, 3]).reshape(ntheta, nphi)
    return theta, phi, et, ep


def solution_point(h, expression, theta=None, phi=None):
    opts = {"expressions": expression, "setup_sweep_name": "Setup1 : LastAdaptive"}
    if theta is not None:
        opts.update(report_category="Far Fields", context="Global",
                    variations={"Freq": ["5.68GHz"], "Theta": [f"{theta}deg"], "Phi": [f"{phi}deg"]})
    data = h.post.get_solution_data(**opts)
    assert data, (expression, theta, phi)
    _, values = data.get_expression_data(expression=expression, formula="real")
    assert len(values) == 1
    return float(values[0])


def direction_separation(th, ph, th0, ph0):
    a = np.deg2rad(th)[:, None]
    b = np.deg2rad(ph)[None, :]
    a0, b0 = np.deg2rad(th0), np.deg2rad(ph0)
    cos = np.sin(a) * math.sin(a0) * np.cos(b - b0) + np.cos(a) * math.cos(a0)
    return np.rad2deg(np.arccos(np.clip(cos, -1, 1)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--state", required=True)
    parser.add_argument("--design", default=None)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    temp = out / "aedt_temp"
    temp.mkdir(exist_ok=True)
    os.environ.update(TEMP=str(temp), TMP=str(temp))
    rep = {"state": args.state, "status": "EXTRACTING"}
    report = out / "metrics.json"
    report.write_text(json.dumps(rep, indent=1))
    h = None
    try:
        h = Hfss(project=str(args.project.resolve()), design=args.design or f"GATE_{args.state}", solution_type="Modal",
                 version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
        rep["own_aedt_pid"] = h.desktop_class.aedt_process_id
        rep["version"] = str(h.desktop_class.aedt_version_id)
        assert rep["version"] == "2026.1"
        h.desktop_class.temp_directory = temp
        assert Path(h.desktop_class.temp_directory).resolve() == temp.resolve()
        assert "Setup1 : LastAdaptive" in h.existing_analysis_sweeps
        target = out / "farfield"
        target.mkdir(exist_ok=True)
        assert h.export_antenna_metadata(frequencies=[5.68e9], setup="Setup1 : LastAdaptive", sphere="Global",
                                         variations={}, output_dir=str(target), export_element_pattern=True,
                                         export_objects=False, export_touchstone=True, export_power=True), "FFD_EXPORT_FAILED"
        files = list(target.rglob("*.ffd"))
        assert len(files) == 1, f"FFD_COUNT_{len(files)}"
        rep["ffd_path"] = str(files[0])
        theta, phi, et, ep = ffd(files[0])
        rep["grid"] = {"theta": [theta[0], theta[-1], len(theta)], "phi": [phi[0], phi[-1], len(phi)]}
        rep["s11_db"] = solution_point(h, "dB(S(Port1,Port1))")
        rep["accepted_power_fraction"] = 1 - 10 ** (rep["s11_db"] / 10)
        ix_th = int(np.argmin(abs(theta - 30)))
        ix_ph = int(np.argmin(abs(phi - 0)))
        plus = np.abs(et + 1j * ep) ** 2
        minus = np.abs(et - 1j * ep) ** 2
        hfss_diff = solution_point(h, "dB(GainRHCP)", 30, 0) - solution_point(h, "dB(GainLHCP)", 30, 0)
        ffd_diff = 10 * math.log10(max(plus[ix_th, ix_ph], 1e-300) / max(minus[ix_th, ix_ph], 1e-300))
        rh = plus if abs(ffd_diff - hfss_diff) < abs(-ffd_diff - hfss_diff) else minus
        lh = minus if rh is plus else plus
        rep["ffd_rhcp_convention"] = "et_plus_j_ep" if rh is plus else "et_minus_j_ep"
        # Query HFSS at each peak to establish an absolute gain scale for the FFD.
        for hand, intensity in (("RHCP", rh), ("LHCP", lh)):
            visible = intensity[theta <= 90, :-1]
            k = np.unravel_index(int(np.argmax(visible)), visible.shape)
            peak_th, peak_ph = float(theta[k[0]]), float(phi[k[1]])
            peak_gain = solution_point(h, f"dB(Gain{hand})", peak_th, peak_ph)
            scale_db = peak_gain - 10 * math.log10(max(float(intensity[k]), 1e-300))
            gain_grid = 10 * np.log10(np.maximum(intensity, 1e-300)) + scale_db
            sep = direction_separation(theta, phi, peak_th, peak_ph)
            sll = float(np.max(gain_grid[(sep >= 10) & (theta[:, None] <= 90)])) - peak_gain
            wide = float(np.max(gain_grid[(theta[:, None] >= 60) & (theta[:, None] <= 90) &
                                          np.ones((1, len(phi)), dtype=bool)]))
            rep[hand] = {"peak_gain_db": peak_gain, "peak_theta_deg": peak_th, "peak_phi_deg": peak_ph,
                         "sll_db_10deg_exclusion": sll, "wide_angle_max_gain_db_theta60_to90": wide,
                         "at_30_phi0_db": solution_point(h, f"dB(Gain{hand})", 30, 0),
                         "at_26_phi180_db": solution_point(h, f"dB(Gain{hand})", 26, 180),
                         "at_30_phi180_db": solution_point(h, f"dB(Gain{hand})", 30, 180),
                         "at_36_phi180_db": solution_point(h, f"dB(Gain{hand})", 36, 180)}
        rep["status"] = "COMPLETE"
        report.write_text(json.dumps(rep, indent=1))
        print("METRICS", rep["RHCP"], rep["LHCP"], flush=True)
    except BaseException:
        rep["status"] = "FAILED"
        rep["traceback"] = traceback.format_exc()
        report.write_text(json.dumps(rep, indent=1))
        print(rep["traceback"], file=sys.stderr, flush=True)
        raise
    finally:
        if h is not None:
            own_pid = int(h.desktop_class.aedt_process_id)
            try:
                own = psutil.Process(own_pid)
                created = own.create_time()
            except psutil.NoSuchProcess:
                own = None
            h.release_desktop(close_projects=True, close_desktop=True)
            if own is not None:
                try:
                    own.wait(timeout=10)
                except psutil.TimeoutExpired:
                    if own.create_time() == created:
                        own.terminate()
                except psutil.NoSuchProcess:
                    pass


if __name__ == "__main__":
    main()
