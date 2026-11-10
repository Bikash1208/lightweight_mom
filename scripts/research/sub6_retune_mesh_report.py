"""Summarize computed candidate mesh/sensitivity gates; no solves."""
import json
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parent))
import sub6_iter1_island_retune as R
from sub6_retune_analyze import read_stage

selected = json.loads((R.OUT/"selected.json").read_text())
rows,rejected = read_stage("refine")
by_key = {(R.key(r["spec"]),r["mesh"]):r for r in rows}
result = dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED",candidates=[],rejected=rejected,
              resonance_gate="NOT_RUN_PENDING_MESH_GATE")
text = ["# BAR50-02V island retune (S2, 2026-11-10)","",
        "MODEL_ONLY / NOT_HFSS_VALIDATED. One diode; physical gap 1.00 mm.","",
        "All numbers below are computed with the confirmed series-L circuit.","",
        "| ID | mesh | alpha OFF/ON Np/m | suppression dB | beta OFF/ON /k0 | delta beta/k0 | OFF-corner suppression min/max dB |",
        "|---|---:|---:|---:|---:|---:|---:|"]
maximum = 0.
complete = True
for index,chosen in enumerate(selected,1):
    ident = R.key(chosen["spec"])
    intent = chosen["spec"].copy()
    intent["island_w"] = intent["slot_w"]-intent["clear"]-1.0
    intent["gap"] = 1.0
    item = dict(id=f"G{index}",spec=chosen["spec"],physical_intent_mm=intent,meshes={})
    for mesh in (96,128):
        row = by_key.get((ident,mesh))
        if row is None:
            complete = False
            continue
        states = row["states"]
        off,on = states["nominal|OFF"],states["nominal|ON"]
        corners = [v for k,v in states.items() if k.endswith("|OFF") and not k.startswith("nominal")]
        suppression = [v["suppression_dB"] for v in corners]
        item["meshes"][str(mesh)] = dict(dimensions_mm=row["dimensions_mm"],states=states,
            corner_suppression_min_dB=min(suppression),corner_suppression_max_dB=max(suppression),
            corner_alpha_OFF_min_np_per_m=min(v["alpha_np_per_m"] for v in corners),
            corner_alpha_OFF_max_np_per_m=max(v["alpha_np_per_m"] for v in corners))
        text.append(f"| G{index} | {mesh} | {off['alpha_np_per_m']:.6f} / {on['alpha_np_per_m']:.6f} | {off['suppression_dB']:.4f} | {off['beta_over_k0']:.7f} / {on['beta_over_k0']:.7f} | {off['beta_ON_minus_OFF_over_k0']:.7f} | {min(suppression):.4f} / {max(suppression):.4f} |")
    if len(item["meshes"]) == 2:
        changes = {}
        for case in item["meshes"]["96"]["states"]:
            a = item["meshes"]["96"]["states"][case]
            b = item["meshes"]["128"]["states"][case]
            changes[case] = {field:abs(b[field]-a[field])/abs(a[field])
                for field in ("alpha_np_per_m","beta_over_k0")}
            da,db = complex(*a["Delta_mS"]),complex(*b["Delta_mS"])
            changes[case]["relative_Delta_difference"] = abs(db-da)/abs(da)
        item["mesh_changes"] = changes
        worst = max(v["alpha_np_per_m"] for v in changes.values())
        item["max_relative_alpha_change"] = worst
        maximum = max(maximum,worst)
        item["mesh_gate"] = "HARD_STOP_ALPHA_CHANGE_GT_5_PERCENT" if worst > .05 else "ALPHA_MESH_CHANGE_LE_5_PERCENT"
    result["candidates"].append(item)
result["complete_mesh_sensitivity_data"] = complete
result["max_relative_alpha_change"] = maximum
result["status"] = ("INCOMPLETE" if not complete else
    "HARD_STOP_MESH_NOT_CONVERGED" if maximum > .05 else "MESH_ALPHA_GATE_PASS")
if result["status"] == "HARD_STOP_MESH_NOT_CONVERGED":
    result["resonance_gate"] = "NOT_RUN_HARD_STOP_MESH_NOT_CONVERGED"
text += ["",f"Status: **{result['status']}**.","",
         "Dimensions below are actual mesh-96 dimensions; mesh-128 realizations are in mesh_report.json.","",
         "| ID | slot L x W mm | island L x W mm | bottom/end clearance mm | diode x mm |",
         "|---|---:|---:|---:|---:|"]
mesh_notes = []
for item in result["candidates"]:
    if "96" in item["meshes"]:
        d = item["meshes"]["96"]["dimensions_mm"]
        text.append(f"| {item['id']} | {d['slot_l']:.6f} x {d['slot_w']:.6f} | {d['island_l']:.6f} x {d['island_w']:.6f} | {d['clearance']:.6f} / {d['end_clearance']:.6f} | {d['diode_x']:.6f} |")
        if "max_relative_alpha_change" in item:
            mesh_notes.append(f"G{item['id'][1:]} maximum relative alpha change: {item['max_relative_alpha_change']*100:.3f}%; {item['mesh_gate']}.")
text += ["",*mesh_notes,"",
         "Physical intent uses gap exactly 1.00 mm and island width = slot width - bottom clearance - 1.00 mm; these widths differ from the grid-realized widths above. See physical_intent_mm in mesh_report.json. No STEP is exported from these mesh masks."]
text += ["","Sensitivity: C={0.105,0.12} pF crossed with Rp={3000,5000} ohm; ON is unchanged.",
         "Baseline nominal C=0.11 pF, Rp=4000 ohm. Both states include L=0.6 nH; ON R=3 ohm.",
         "All cases use Floquet order 2n; accepted gap realization 8 cells (96), 11 cells (128).",
         "Sheet extraction uses the existing two-node self-consistent mean-sheet/first-order leakage model; raw tensors and extrapolation flags are preserved.",
         "","Resonance/stop-band clearance is NOT_ESTABLISHED: frequency screening was not launched after the mandated >5% mesh-change hard stop.",
         "","No candidate is qualified for a board or HFSS export. The table contains raw two-node model outputs, not validated propagation constants.",
         "","Post-stop sheet audit (existing tensors only): nominal fixed-point residuals are <=1e-5, but all ON states extrapolate outside the 1.017..1.097 kt/k0 node interval. Mesh-128 ON beta/k0=1.265..1.430 is far outside that interval; direct Bloch validation at those values has not been run. See sheet_iteration_audit.json.",
         "","Validation: production baseline alpha/beta reproduced within 2.12e-12 relative; existing FFT assembly audit <=1.03e-15 relative sampled matrix error; exact load-update/direct-solve tensor audit <=1.08e-12. No src/ solver changes.",
         "","Fast suite in previous S2 env: 2954 passed, 129 skipped, six accepted ENV_ONLY S11 hash failures, one 300-second app timeout. App test rerun passed in 68.55 seconds; four new geometry/circuit checks passed. Original regression tests and fixtures unchanged. The initially reported four ENV_ONLY failures are documented in ENV_DIAGNOSIS.md.",
         "","HFSS_SETUP.md and PACKAGE_9CELL_PLAN.md corrected: OFF sheet impedance 5.424-j79.320 ohm/sq includes series 0.6 nH; former -86.97 omitted it. ON remains 1.05+j7.65. No prior HFSS run is relabelled.",
         "","Original C:/lwmom_final modified CSV: only sheet_z_thick column moved after sheet_centre_y; all 240 rows' named values unchanged. Original artifacts preserved."]
(R.OUT/"mesh_report.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
(R.OUT/"REPORT.md").write_text("\n".join(text)+"\n",encoding="utf-8")
print(result["status"],"maximum relative alpha change",maximum)
print("\n".join(text))
