"""Expected readouts from existing tensors; no new MoM solves or tuning.

Uniform-cell eigenmode Q measures loss, while first-order tensor leakage
describes a modulated hologram. Keep them separate. Optional directional
sheet-loss estimates use the transverse-resonance perturbation -G/dBdbeta;
all roots/extrapolation and invalid signs are explicitly retained.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(Path(__file__).resolve().parent))
import sub6_island_meshtrend as M
import sub6_iter1_island_retune as R
from src.dispersion.sheet_model import _y_air,_y_slab

OUT = ROOT/"hfss_exports/sub6_island_baseline_g1"
F = 5.8e9
K0 = R.free_space_wavenumber(F)
SUB = R.S.make_sub(R.SUB)


def directional(nodes,state,direction):
    ks = sorted(nodes)
    component = 0 if direction=="x" else 1
    values = [nodes[k][direction][state][component,component] for k in ks]
    slope = (values[-1]-values[0])/(ks[-1]-ks[0])
    def total(q):
        return values[0]+slope*(q-ks[0])+_y_air(q*K0,2*math.pi*F,K0)+_y_slab(q*K0,2*math.pi*F,SUB.eps_complex,K0,SUB.thickness)
    qs = np.linspace(1.0001,math.sqrt(2.2)-.0001,4000)
    vals = np.imag(total(qs))
    roots = []
    for i in np.flatnonzero(np.diff(np.sign(vals))):
        q = float(brentq(lambda t:total(t).imag,float(qs[i]),float(qs[i+1])))
        h = 1e-6
        derivative = (total(q+h).imag-total(q-h).imag)/(2*h*K0)
        conductance = float(total(q).real)
        alpha = -conductance/derivative
        valid = alpha>=0 and abs(alpha)/(q*K0)<.01
        roots.append(dict(beta_over_k0=q,beta_rad_per_m=q*K0,phase_deg=q*K0*.0116*180/math.pi,
            estimated_air_field_amplitude_at_150mm_lid=math.exp(-K0*math.sqrt(q*q-1)*.150),
            extrapolated=bool(q<ks[0] or q>ks[-1]),alpha_loss_perturbative_raw_np_per_m=alpha,
            loss_estimate_gate="POSITIVE_WEAK_LOSS_ESTIMATE" if valid else "PERTURBATIVE_LOSS_UNRELIABLE",
            conductance_S=conductance,d_total_susceptance_d_beta=derivative,
            alpha_x_estimate_np_per_m=alpha if valid else None))
    return dict(method="directional linear spatial sheet interpolation; all transverse-resonance roots; alpha_loss=-Re(Ytotal)/dIm(Ytotal)/dBeta",
                validation="MODEL_ONLY / NOT_HFSS_VALIDATED; ON roots mesh-unreliable",roots=roots)


def main():
    groups = {}
    for line in (M.OUT/"nodes.jsonl").read_text().splitlines():
        row = json.loads(line)
        groups.setdefault((R.key(row["spec"]),row["mesh"]),[]).append(row)
    trends = json.loads((M.OUT/"mesh_trend.json").read_text())["rows"]
    summary = {(r["cell"],r["mesh"]):r for r in trends}
    for cell,spec in M.SPECS.items():
        model = []
        for (ident,n),raw in groups.items():
            if ident!=R.key(spec) or len(raw)!=4:
                continue
            nodes = {}
            for r in raw:
                nodes.setdefault(r["kt"],{})[r["direction"]] = {}
                for st,y in r["tensors"].items():
                    a=np.asarray(y);nodes[r["kt"]][r["direction"]][st]=a[:,:,0]+1j*a[:,:,1]
            for state in ("OFF","ON"):
                mean = summary[(cell,n)]["states"]["nominal|"+state]
                model.append(dict(state=state,mesh=n,mean_sheet_beta_over_k0=mean["beta_over_k0"],mean_beta_extrapolated=mean["extrapolated"],
                    alpha_holographic_leakage_np_per_m=mean["alpha_np_per_m"],
                    beta_and_loss_by_direction={d:directional(nodes,"nominal|"+state,d) for d in "xy"}))
        for variant in ("novia","vias_fr4_stubs"):
            folder=OUT/(cell+"_"+variant)
            result=dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED",cell=cell,variant=variant,
                readout_required="per state and direction: f(phi), Q(phi), top-TM mode identification; beta58, alpha_x58=omega/(2Qvg)",
                actual_3D_variant_beta_prediction="NOT_AVAILABLE_FROM_CURRENT_DATA" if variant!="novia" else "directional model roots below; extrapolation flags apply",
                actual_3D_variant_alpha_x_prediction="NOT_AVAILABLE_FROM_CURRENT_DATA" if variant!="novia" else "perturbative directional sheet loss estimates below; not HFSS validated",
                no_via_model_references=sorted(model,key=lambda r:(r["mesh"],r["state"])),
                alpha_warning="Do not equate first-order holographic leakage alpha with uniform-cell eigenmode Q-derived loss alpha_x. No vias/FR4/stubs in MoM. Directional loss estimates are perturbative, not independently validated.",
                frame_mapping="MoM x=along slot=HFSS y; MoM y=across slot=HFSS x; +90deg package rotation",
                diode_boundary="one Impedance sheet; specified values at 5.8GHz; record frozen vs frequency-dependent boundary when extracting vg",
                HFSS_reference_method=dict(script="scripts/research/sub6_final_cal.py",data="results/sub6_final_board/hfss_data.json",setup_commit="21fe444",setup_path="hfss_exports/sub6_island_choke/HFSS_SETUP.md"))
            (folder/"predictions.json").write_text(json.dumps(result,indent=2),encoding="utf-8",newline="\n")
            with (folder/"model_reference.csv").open("w",encoding="utf-8",newline="") as stream:
                writer=csv.writer(stream,lineterminator="\n")
                writer.writerow(["state","mesh","mean_beta_over_k0","mean_beta_extrapolated","alpha_leak_Np_per_m_NOT_alpha_x","direction_MoM","direction_HFSS","root_beta_over_k0","root_extrapolated","alpha_x_sheet_estimate_Np_per_m","loss_gate","variant_3D_prediction_available"])
                for row in model:
                    for d,data in row["beta_and_loss_by_direction"].items():
                        for root in data["roots"]:
                            writer.writerow([row["state"],row["mesh"],row["mean_sheet_beta_over_k0"],row["mean_beta_extrapolated"],row["alpha_holographic_leakage_np_per_m"],d,"along_y" if d=="x" else "across_x",root["beta_over_k0"],root["extrapolated"],root["alpha_x_estimate_np_per_m"],root["loss_estimate_gate"],variant=="novia"])
            def cell_value(state,field):
                row=next(r for r in model if r["mesh"]==96 and r["state"]==state)
                values=[]
                for d in "xy":
                    roots=row["beta_and_loss_by_direction"][d]["roots"]
                    if len(roots)!=1:
                        values.append("multiple/no root")
                    elif roots[0][field] is None:
                        values.append("unreliable")
                    else:
                        values.append(f"{roots[0][field]:.5f}"+(" E" if roots[0]["extrapolated"] else ""))
                return " / ".join(values)
            marker="<!-- MODEL_READOUT_REFERENCE -->"
            manual=(folder/"HFSS_SETUP.md").read_text(encoding="utf-8").split(marker)[0].rstrip()
            table = f"\n\n{marker}\n**Expected no-via sheet-model readouts, mesh 96, MODEL_ONLY (along / across):**\n\n| Quantity | OFF | ON |\n|---|---:|---:|\n| beta/k0 | {cell_value('OFF','beta_over_k0')} | {cell_value('ON','beta_over_k0')} |\n| alpha_x loss estimate, Np/m | {cell_value('OFF','alpha_x_estimate_np_per_m')} | {cell_value('ON','alpha_x_estimate_np_per_m')} |\n\nE = extrapolated; not acceptance limits. Mesh trends and all roots: model_reference.csv. {'Actual via/FR4/stub alpha_x and beta: NOT_AVAILABLE_FROM_CURRENT_DATA; table is no-via reference only.' if variant!='novia' else 'ON is mesh-unreliable; determine its branch and slope from HFSS.'}\n"
            table+="Check Ez is small at the 150-mm lid; otherwise flag AIR_HEIGHT_NOT_CONVERGED and verify a higher lid before interpreting the mode.\n"
            (folder/"HFSS_SETUP.md").write_text(manual+table,encoding="utf-8",newline="\n")
            print(cell,variant,"prediction rows",len(model),flush=True)


if __name__ == "__main__":
    main()
