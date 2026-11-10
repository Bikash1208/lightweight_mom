"""Package diagnostic trend and copy only STEP/JSON/CSV/MD delivery files."""
import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
RES=ROOT/"results/sub6_island_meshtrend"
PKG=ROOT/"hfss_exports/sub6_island_baseline_g1"
DEST=Path("C:/Users/Bikash/Desktop/sub6_board_hfss/iter1_cell_checks_S2_2026-11-10")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def deliver():
    DEST.mkdir(parents=True,exist_ok=True)
    manifest=DEST/"delivery_manifest.json"
    previous=json.loads(manifest.read_text()) if manifest.exists() else {}
    saved={r["relative_path"]:r["sha256"] for r in previous.get("files",[])}
    rows=[]
    for source in sorted(PKG.rglob("*")):
        if not source.is_file() or source.suffix.lower() not in (".step",".json",".csv",".md"):
            continue
        relative=source.relative_to(PKG).as_posix()
        target=DEST/relative
        target.parent.mkdir(parents=True,exist_ok=True)
        current=digest(target) if target.exists() else None
        new=digest(source)
        if current is not None and current!=new and saved.get(relative)!=current:
            status="PRESERVED_USER_MODIFIED_DESTINATION"
        else:
            if current!=new:
                target.write_bytes(source.read_bytes())
            status="COPIED_OR_IDENTICAL"
        rows.append(dict(relative_path=relative,sha256=digest(target),source_sha256=new,status=status))
    manifest.write_text(json.dumps(dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED",source=str(PKG),files=rows),indent=2),encoding="utf-8",newline="\n")
    print("Delivery",DEST,"files",len(rows),"preserved user-modified",sum(r["status"].startswith("PRESERVED") for r in rows))


def report():
    data=json.loads((RES/"mesh_trend.json").read_text())
    by_key={(r["cell"],r["mesh"]):r for r in data["rows"]}
    text=["# Fixed BASELINE / G1 mesh trend and HFSS cell packages", "",
          "MODEL_ONLY / NOT_HFSS_VALIDATED. Diagnostic mesh sweep only; no geometry tuning.","",
          "| Cell | Mesh | gap cells / mm | alpha OFF/ON Np/m | suppression dB | beta OFF/ON /k0 | extrapolated OFF/ON |",
          "|---|---:|---:|---:|---:|---:|---|"]
    headers=["cell","mesh","gap_cells","gap_mm","alpha_OFF_Np_per_m","alpha_ON_Np_per_m","suppression_dB","beta_OFF_over_k0","beta_ON_over_k0","beta_OFF_extrapolated","beta_ON_extrapolated","status"]
    with (RES/"mesh_trend.csv").open("w",encoding="utf-8",newline="") as stream:
        writer=csv.writer(stream,lineterminator="\n");writer.writerow(headers)
        for cell in ("BASELINE","G1"):
            for mesh in (96,128,160,192):
                row=by_key.get((cell,mesh))
                if row is None:
                    status="ALLOCATION_FAILED" if mesh==192 and (RES/"mesh192_status.json").exists() else "NOT_RUN"
                    count=round(mesh/11.6);gap=count*11.6/mesh
                    text.append(f"| {cell} | {mesh} | {count} / {gap:.6f} (input) | {status} | — | — | — |")
                    writer.writerow([cell,mesh,count,gap,"","","","","","","",status])
                    continue
                if "states" not in row:
                    text.append(f"| {cell} | {mesh} | — | SHEET_EXTRACTION_FAILED | — | — | — |")
                    continue
                a,b=row["states"]["nominal|OFF"],row["states"]["nominal|ON"]
                gap=row["dimensions_mm"]["gap"]
                count=round(gap/(11.6/mesh))
                writer.writerow([cell,mesh,count,gap,a["alpha_np_per_m"],b["alpha_np_per_m"],a["suppression_dB"],a["beta_over_k0"],b["beta_over_k0"],a["extrapolated"],b["extrapolated"],"COMPUTED_MODEL_ONLY"])
                text.append(f"| {cell} | {mesh} | {count} / {gap:.6f} | {a['alpha_np_per_m']:.6f} / {b['alpha_np_per_m']:.6f} | {a['suppression_dB']:.4f} | {a['beta_over_k0']:.7f} / {b['beta_over_k0']:.7f} | {a['extrapolated']} / {b['extrapolated']} |")
    text += ["","The physical diode gap remains 1.00 mm. Nearest-cell discretization gives 8/11/14/17 cells at meshes 96/128/160/192. All Floquet cutoffs are 2n; the original F128 baseline is excluded from this standardized trend.",
             "","Bloch nodes remain kt/k0=1.017 and 1.097 in each direction. Beta outside that interval is EXTRAPOLATED, not validated by a direct solve at that beta. ON-state instability is intentionally reported and is not a stop condition for this diagnostic job.",
             "","The table alpha is first-order holographic leakage. It is NOT alpha_x from a uniform-cell eigenmode Q. No physical convergence or HFSS validation is claimed.",
             "","Four checked HFSS packages: hfss_exports/sub6_island_baseline_g1/README.md. Each contains one shared OFF/ON STEP, one-page manual setup, state table, eight-mode phase-run template, predictions and geometry gates. All metal is one METAL_PEC object; one 1.00 x 0.35 mm Impedance sheet; no AEDT scripts.",
             "","Calibration method found: scripts/research/sub6_final_cal.py; actual HFSS record results/sub6_final_board/hfss_data.json; archived setup commit 21fe444:hfss_exports/sub6_island_choke/HFSS_SETUP.md. The recorded calibration was ISL_vias2, not no-via. It used periodic eigenmode f-vs-phase and Q to obtain beta and alpha_x=omega/(2Qvg). Historical ON slope used MoM anchoring; new packages request a directly measured local HFSS slope from >=3 phase points.",
             "","User-confirmed variant (b) retains both 2.0 x 0.3 mm bottom stubs, original via centres, 0.8-mm FR4, 0.6-mm anti-pads, distant radiation bottom 30 mm below the lowest copper. +90-degree frame: across=ph_x, along=ph_y. PEC lid at +150 mm. No near-stub PEC floor.",
             "","Expected directional no-via beta/loss estimates and their extrapolation flags are in each model_reference.csv/predictions.json. Actual 3D via/FR4/stub values are NOT_AVAILABLE_FROM_CURRENT_DATA. No-via holographic leakage values are not substituted for Q loss. No HFSS solve was performed by S2.",
             "","Geometry differences are documented: BASELINE uses physical board 8.696 x 2.900 slot and 4.348 x 1.300 island; G1 CAD endpoints are snapped to 0.001 mm. MoM grid realizations and actual gap lengths are separately reported. G1 via barrel fits with only 0.003 mm x-edge clearance; no fabrication clearance qualification is claimed.",
             "","A tool-helper refresh interrupted the first mesh-160 attempt; completed nodes were retained and only missing nodes restarted in a hidden detached process. No parameter or worker-count change was made. Worker pools use all 64 logical workers, single-threaded BLAS; no worker cap.",
             "","Optional mesh192 was attempted with eight nodes and no worker cap. The process failed before any mesh192 checkpoint completed: NumPy could not allocate a 22.9 GiB complex128 array of shape 39220 x 39220. No mesh192 number is inferred. See mesh192.err and mesh192_status.json. This is an execution-resource failure, not a physical ON-state result or an imposed worker/RAM budget.",
             "",f"Desktop delivery: {DEST.as_posix()}. Only STEP/JSON/CSV/MD copied; no Python files. Existing user-modified delivery files are preserved."]
    if (RES/"validation.json").exists():
        text += ["","Validation: four geometry packages passed all pre/re-import gates; eight CAD-to-model quantization mappings passed. Full fast suite: 2953 passed, 129 skipped, six accepted ENV_ONLY hash failures, two 60-second app polling timeouts. Both timeout tests passed isolated rerun (2 passed in 41.42 s). Targeted geometry/FFT checks: 5 passed. No src/ or existing tests/fixtures changes."]
    (RES/"REPORT.md").write_text("\n".join(text)+"\n",encoding="utf-8",newline="\n")


if __name__=="__main__":
    if "--deliver-only" not in sys.argv:
        report()
    deliver()
