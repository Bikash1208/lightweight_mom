"""Phase 29 design point B3 overnight driver (unattended): waits for the two GAs, final evaluation of both cones, Frey 2->8-feed trend, decision, then
(a pass) the HFSS package for the preferred passing cone, or (no pass) one bounded second design pass (<= 2 h) on the more promising cone and stop (no package);
then the new tests, the fast suite, STATUS / brief notes, commit and push. Progress: results/phase29/B3_progress.txt, summary: results/phase29/B3_summary.json.

Conservative readings (also in the morning report): pass = pointing <= 2 deg AND gain dip <= 3 dB over the WHOLE cone with the tuned weights; both-pass -> cone a;
a package is built only after a FIRST-pass success (a second-pass success is recorded as 'ready', no package); heavy steps are skipped after 06:30 IST; hard stop 07:30."""

import datetime as dt
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results/phase29"
PY = sys.executable
ENV = dict(os.environ, PYTHONPATH=str(ROOT))
PROG = RES / "B3_progress.txt"
ATTR = "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>\nClaude-Session: https://claude.ai/code/session_01UaAn25MEtRcvrHeagn5k3z"


def log(msg):
    with open(PROG, "a", encoding="utf-8") as f:
        f.write(f"{dt.datetime.now():%m-%d %H:%M} {msg}\n")


def before(hour, minute):
    """True while the local time is before the NEXT hour:minute (the morning cutoff)."""
    n = dt.datetime.now()
    t = n.replace(hour=hour, minute=minute, second=0, microsecond=0)
    if n.hour >= 12:
        t += dt.timedelta(days=1)
    return n < t


def run(args, logname, timeout=None):
    with open(RES / logname, "w", encoding="utf-8") as f:
        try:
            r = subprocess.run([PY, *args], cwd=ROOT, env=ENV, stdout=f, stderr=subprocess.STDOUT, timeout=timeout)
            return r.returncode
        except subprocess.TimeoutExpired:
            return -9


def step(name, args, logname, timeout=None):
    """One attempt, one retry on failure; failures are recorded, never fatal."""
    for attempt in (1, 2):
        rc = run(args, logname, timeout)
        log(f"{name}: attempt {attempt} rc={rc}")
        if rc == 0:
            return True
    return False


def tail(logname, n=3):
    try:
        return " | ".join((RES / logname).read_text(encoding="utf-8", errors="replace").strip().splitlines()[-n:])[:400]
    except OSError:
        return ""


def git(*a):
    r = subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True)
    log(f"git {a[0]} rc={r.returncode} {r.stdout.strip()[-120:]} {r.stderr.strip()[-160:]}")
    return r.returncode == 0


def summarize(cone):
    p = RES / f"pointB3_{cone}.json"
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    t, pm = d["scan_tuned"], d["ports"]
    cp = d["coupling_db"]
    return {
        "ring_over_lambda0": round(d["ring_over_lambda0"], 2),
        "cells": d["n_cells"],
        "slots": d["slots"],
        "port_gain_dbi": [round(q["peak_gain_dbi"], 1) for q in pm],
        "port_hpbw_deg": [round(q["hpbw_deg"]) for q in pm],
        "port_sll_db": [round(q["sll_db"], 1) for q in pm],
        "port_peak_theta_phi": [[round(x) for x in q["peak_theta_phi"]] for q in pm],
        "crossover_mean_db": {k: round(v["mean_db"], 1) for k, v in d["crossover"].items()},
        "sum": {
            k: (round(v, 1) if isinstance(v, float) else [round(x) for x in v]) for k, v in d["sum_diff"]["sum"].items()
        },
        "delta_null_db": {k: round(d["sum_diff"][k]["null_depth_at_sum_peak_db"], 1) for k in ("delta_x", "delta_y")},
        "scan_tuned": {
            k: (round(v, 2) if isinstance(v, float) else v)
            for k, v in t.items()
            if k not in ("method_counts", "cone_edge_between_outer_beams")
        },
        "scan_mrc": {k: round(v, 2) for k, v in d["scan_mrc"].items()},
        "cone_edge": [
            [r["phi_deg"], round(r["pointing_error_deg"], 1), round(r["peak_gain_dbi"], 1)]
            for r in t["cone_edge_between_outer_beams"]
        ],
        "coupling_raw_db_range": [
            round(min(v["raw_db"] for v in cp.values()), 1),
            round(max(v["raw_db"] for v in cp.values()), 1),
        ],
        "coupling_calibrated_db_range": [
            round(min(v["calibrated_db"] for v in cp.values()), 1),
            round(max(v["calibrated_db"] for v in cp.values()), 1),
        ],
        "ceiling_gap_db": {k: round(v, 1) for k, v in d["ceiling_gap_db"].items()},
        "pass": bool(t["pass_pointing_2deg"] and t["pass_dip_3db"]),
        "merit": t["max_pointing_error_deg"] / 2.0 + t["max_gain_dip_db"] / 3.0,
    }


def save(summary):
    (RES / "B3_summary.json").write_text(json.dumps(summary, indent=1))


def main():
    S = {"started": f"{dt.datetime.now():%Y-%m-%d %H:%M}", "decisions": []}
    log("driver started")
    # 1. wait for both GAs (launched before the driver)
    while not ((RES / "pointB3_a_ga.json").exists() and (RES / "pointB3_b_ga.json").exists()) and before(6, 30):
        time.sleep(60)
    log("GAs done" if (RES / "pointB3_b_ga.json").exists() else "GAs NOT finished before the cutoff")
    # 2. final evaluation of both cones
    for c in ("a", "b"):
        if (RES / f"pointB3_{c}_ga.json").exists() and before(6, 30):
            step(f"final {c}", ["scripts/research/phase29_pointB3.py", "final", "--cone", c], f"pointB3_{c}_final.log")
    # 3. trend sanity check
    if before(6, 30):
        step("frey trend", ["scripts/research/phase29_frey_trend.py"], "frey_trend.log")
    S["a"], S["b"] = summarize("a"), summarize("b")
    if (RES / "frey_trend.json").exists():
        tr = json.loads((RES / "frey_trend.json").read_text())
        S["frey_trend"] = {
            "checks": tr["checks"],
            "gain_frey_dbi": [round(x["gain_frey_dbi"], 1) for x in tr["trend"]],
            "gain_diff_db": [round(x["gain_diff_db"], 2) for x in tr["trend"]],
            "sll_frey_db": [round(x["sll_frey_db"], 1) for x in tr["trend"]],
        }
    save(S)
    # 4. decision
    passing = [c for c in ("a", "b") if S[c] and S[c]["pass"]]
    S["decisions"].append(
        f"first pass: passing cones {passing or 'none'} (pass = pointing <= 2 deg and dip <= 3 dB over the whole cone, tuned weights)"
    )
    if passing and before(6, 30):
        cone = passing[0]
        S["package_cone"] = cone
        ok = step("export package", ["scripts/research/phase29_B3_export.py", "--cone", cone], "B3_export.log")
        S["package_built"] = bool(ok)
        S["decisions"].append(f"HFSS package built for cone {cone} (preferred a if both pass): {ok}")
    elif S["a"] or S["b"]:
        cands = [c for c in ("a", "b") if S[c]]
        cone = min(cands, key=lambda c: S[c]["merit"])
        S["decisions"].append(
            f"no cone passes: second design pass (<= 100 min GA) on cone {cone} (lower merit {S[cone]['merit']:.2f} = pointing/2 + dip/3)"
        )
        S["pass1_" + cone] = S[cone]
        if before(5, 0):
            for f in (f"pointB3_{cone}_ga.json", f"pointB3_{cone}.json"):
                (RES / f).replace(RES / f.replace(".json", "_pass1.json"))
            (RES / f"pointB3_{cone}_ga_pass1.json").write_text((RES / f"pointB3_{cone}_ga_pass1.json").read_text())
            ga_src = json.loads((RES / f"pointB3_{cone}_ga_pass1.json").read_text())
            (RES / f"pointB3_{cone}_ga.json").write_text(json.dumps(ga_src))  # the warm start reads this
            okg = step(
                "second GA",
                [
                    "scripts/research/phase29_pointB3.py",
                    "ga",
                    "--cone",
                    cone,
                    "--pop",
                    "64",
                    "--gens",
                    "40",
                    "--max-minutes",
                    "100",
                    "--warm",
                    "--seed",
                    "2",
                ],
                f"pointB3_{cone}_ga2.log",
            )
            if okg and before(6, 30):
                step(
                    "final pass 2",
                    ["scripts/research/phase29_pointB3.py", "final", "--cone", cone],
                    f"pointB3_{cone}_final2.log",
                )
            S[cone + "_pass2"] = summarize(cone)
            S["decisions"].append(
                f"second pass done on cone {cone}; pass2 result passes={bool(S[cone + '_pass2'] and S[cone + '_pass2']['pass'])}; no package built (second-pass success is only recorded: run phase29_B3_export.py --cone {cone})"
            )
        else:
            S["decisions"].append("too late for a second pass (after 05:00): skipped")
    save(S)
    # 5. tests, fast suite
    if before(7, 0):
        rc = run(
            [
                "-m",
                "pytest",
                "tests/test_shared_aperture_frey.py",
                "tests/test_aedt_stub_phase29.py",
                "tests/test_shared_aperture_multiport.py",
                "tests/test_seven_hex_cluster.py",
                "tests/test_compact7_export.py",
                "-q",
                "-p",
                "no:cacheprovider",
            ],
            "B3_newtests.log",
        )
        S["new_tests"] = {"rc": rc, "tail": tail("B3_newtests.log", 2)}
        if before(6, 40):
            rc = run(
                ["-m", "pytest", "-m", "not slow", "-q", "-n", "40", "-p", "no:cacheprovider"],
                "B3_fast_suite.log",
                timeout=75 * 60,
            )
            S["fast_suite"] = {"rc": rc, "tail": tail("B3_fast_suite.log", 2)}
        else:
            S["fast_suite"] = {"rc": None, "tail": "skipped (time)"}
    save(S)
    # 6. docs, commit, push
    write_docs(S)
    git("add", "-A")
    msg = "Phase 29 point B3: 7-port shared 150 mm hexagon, Frey angle-dependent superposition + GA, two scan cones (MODEL_ONLY), overnight run"
    git("commit", "-q", "-m", f"{msg}\n\n{ATTR}")
    git("push", "-u", "hub", "phase29-seven-hex")
    log("driver finished")
    (RES / "B3_DONE").write_text("done\n")


def write_docs(S):
    def line(c):
        x = S.get(c)
        if not x:
            return f"cone {c}: not evaluated"
        t = x["scan_tuned"]
        return (
            f"cone {c}: gains {min(x['port_gain_dbi'])}-{max(x['port_gain_dbi'])} dBi, HPBW {min(x['port_hpbw_deg'])}-{max(x['port_hpbw_deg'])} deg, scan worst pointing {t['max_pointing_error_deg']} deg "
            f"(target 2), dip {t['max_gain_dip_db']} dB (target 3), pass={x['pass']}"
        )

    row = (
        f"\n## Phase 29 design point B3 (shared 150 mm hexagon, 7 ports, Frey angle-dependent superposition + GA) - MODEL_ONLY, overnight run {S['started']}\n\n"
        f"{line('a')}; {line('b')}. Decisions: {' | '.join(S['decisions'])}. HFSS package: {'`hfss_exports/phase29_B3/` (cone ' + S['package_cone'] + ')' if S.get('package_built') else 'not built'}. "
        "Numbers: `results/phase29/B3_summary.json`, `pointB3_<cone>.json`.\n\n"
    )
    p = ROOT / "docs/STATUS.md"
    s = p.read_text(encoding="latin-1")
    i = s.index("## Open next steps")
    p.write_text(s[:i] + row + s[i:], encoding="latin-1", newline="")
    b = ROOT / "docs/briefs/PHASE29_SEVEN_HEX_CLUSTER.md"
    b.write_text(
        b.read_text(encoding="utf-8")
        + "\n## Addendum 2026-10-27 - design point B3 (replaces B2; user asleep, unattended)\n"
        "One shared 150 mm hexagon, 7 monopoles (centre + ring 1.35-1.5 lambda0), Frey et al. angle-dependent superposition + GA; cones a (+-20 deg, ~13 deg beams, tilt 13) and b (+-25 deg, ~17 deg, tilt 17). "
        "Targets: pointing <= 2 deg, gain dip <= 3 dB over the whole cone. Driver: `scripts/research/phase29_B3_driver.py` (conservative readings in its docstring). Package only after a first-pass success.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
