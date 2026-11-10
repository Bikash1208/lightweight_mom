"""FAILED-row report: by geometry and state, with a flag for a possible PHYSICS / SETUP issue (not just a software retry).

Rules (user, 2026-11-10):
  * overall failure rate > 5 % of the rows attempted  -> FLAG
  * failures clustering on ONE state / direction / phase / geometry family / project variant -> FLAG
    (cluster = category holds >= 60 % of all failures, at least 3 failures, and its own failure rate is >= 2x the rate of the other rows)

python failed_rows_report.py --runs screen_bar90_nostub_run screen_bar90_trace_run screen_bar90_full_run --out screen_failed_report
"""
import argparse
import csv
import json
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
RID = re.compile(r"^(?P<stem>.+?)_(?P<geom>.+)_(?P<phase>\d+)_(?P<state>OFF|ON)_(?P<dir>across|along)$")


def parse(run_id, stem):
    rest = run_id[len(stem) + 1:] if run_id.startswith(stem + "_") else run_id
    m = re.match(r"^(?P<geom>.+)_(?P<phase>\d+)_(?P<state>OFF|ON)_(?P<dir>across|along)$", rest)
    if not m:
        return {"geometry": rest, "phase": "?", "state": "?", "direction": "?", "family": "?"}
    geom = m.group("geom")
    family = "NOMINAL" if geom.startswith("NOMINAL") else re.sub(r"_-?[\d.]+(?=_|$)", "", geom)
    return {"geometry": geom, "phase": m.group("phase"), "state": m.group("state"), "direction": m.group("dir"), "family": family}


def variant_of(stem):
    return "nostub" if "nostub" in stem else "trace" if "trace" in stem else "full"


def collect(run_dirs):
    """Returns (all_rows, failed) where each item is a dict with category fields."""
    all_rows, failed = [], []
    for run in run_dirs:
        run = run if run.is_absolute() else HERE / run
        rep_path = run / "parallel_report.json"
        if not rep_path.exists():
            continue
        files = {r["file"] for r in json.loads(rep_path.read_text())}
        for fname in sorted(files):
            stem = Path(fname).stem
            csv_path = HERE / fname
            if not csv_path.exists():
                continue
            for row in csv.DictReader(csv_path.open(newline="", encoding="utf-8")):
                d = parse(row["*"], stem); d.update(run_id=row["*"], file=stem, variant=variant_of(stem)); all_rows.append(d)
            ff = run / f"FAILED_ROWS_{stem}.json"
            if ff.exists():
                for x in json.loads(ff.read_text()):
                    d = parse(x["run_id"], stem); d.update(run_id=x["run_id"], file=stem, variant=variant_of(stem),
                                                          error=(x.get("errors") or ["?"])[-1][:90]); failed.append(d)
    return all_rows, failed


def analyse(all_rows, failed):
    n, nf = len(all_rows), len(failed)
    rate = nf / n if n else 0.0
    out = {"rows_attempted": n, "rows_failed": nf, "failure_rate_pct": round(100 * rate, 2), "flags": [], "by": {}}
    fid = {x["run_id"] for x in failed}
    for dim in ("state", "direction", "phase", "family", "variant", "geometry"):
        tot, bad = {}, {}
        for r in all_rows:
            tot[r[dim]] = tot.get(r[dim], 0) + 1
            if r["run_id"] in fid:
                bad[r[dim]] = bad.get(r[dim], 0) + 1
        out["by"][dim] = {k: {"failed": bad.get(k, 0), "total": tot[k], "rate_pct": round(100 * bad.get(k, 0) / tot[k], 1)}
                          for k in sorted(tot) if bad.get(k, 0)}
        for k, b in bad.items():
            rest_tot, rest_bad = n - tot[k], nf - b
            rest_rate = rest_bad / rest_tot if rest_tot else 0.0
            if b >= 3 and nf and b / nf >= 0.6 and (b / tot[k]) >= 2 * rest_rate and len(tot) > 1:
                out["flags"].append(f"CLUSTER on {dim} = {k}: {b} of {nf} failures ({100*b/tot[k]:.0f}% of that category's rows vs "
                                    f"{100*rest_rate:.1f}% elsewhere) -> possible PHYSICS or SETUP issue, not only a software retry")
    if rate > 0.05:
        out["flags"].insert(0, f"FAILURE RATE {100*rate:.1f}% > 5% of {n} rows -> possible PHYSICS or SETUP issue, not only a software retry")
    errs = {}
    for x in failed:
        errs[x["error"]] = errs.get(x["error"], 0) + 1
    out["error_texts"] = dict(sorted(errs.items(), key=lambda kv: -kv[1]))
    return out


def render(out, failed):
    L = [f"# FAILED rows: {out['rows_failed']} of {out['rows_attempted']} attempted ({out['failure_rate_pct']} %)", ""]
    L += ["**FLAGS:** " + ("; ".join(out["flags"]) if out["flags"] else "none (rate <= 5 % and no clustering)"), ""]
    for dim in ("state", "direction", "phase", "variant", "family"):
        if out["by"][dim]:
            L.append(f"## by {dim}"); L.append("| value | failed | of | rate |"); L.append("|---|---|---|---|")
            for k, v in out["by"][dim].items():
                L.append(f"| {k} | {v['failed']} | {v['total']} | {v['rate_pct']} % |")
            L.append("")
    L += ["## each failed row (geometry / state / phase / direction)", "| variant | geometry | state | phase | dir | error |", "|---|---|---|---|---|---|"]
    for x in sorted(failed, key=lambda r: (r["variant"], r["geometry"], r["state"], r["phase"])):
        L.append(f"| {x['variant']} | {x['geometry']} | {x['state']} | {x['phase']} | {x['direction']} | {x['error']} |")
    L += ["", "## error texts", *[f"- {c} x  {t}" for t, c in out["error_texts"].items()]]
    return "\n".join(L) + "\n"


def write_report(run_dirs, out_dir):
    all_rows, failed = collect(run_dirs)
    out = analyse(all_rows, failed)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "FAILED_ROWS_REPORT.json").write_text(json.dumps({"summary": out, "failed": failed}, indent=2), encoding="utf-8")
    (out_dir / "FAILED_ROWS_REPORT.md").write_text(render(out, failed), encoding="utf-8")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", type=Path, nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    o = write_report(a.runs, a.out if a.out.is_absolute() else HERE / a.out)
    print(json.dumps({k: o[k] for k in ("rows_attempted", "rows_failed", "failure_rate_pct", "flags")}, indent=1))
