"""Put corrected diode full-board FFD metrics into the existing v3.1 summary shape."""

import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
from full_board_gate_v31_extract import ffd
from sidequest_sll_def import analyze

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"results/sub6_meander_model/full_board_gate/v31_diode_v2"
REF=ROOT/"results/sub6_meander_model/full_board_gate/v31/metrics_summary.json"


def main():
    both=json.loads(REF.read_text(encoding="utf-8"))["states"]["BOTH"]
    states={"BOTH":both}
    for state in ("ON_A","OFF_D"):
        path=OUT/f"{state}_extract/metrics.json"
        if not path.exists():continue
        metrics=json.loads(path.read_text(encoding="utf-8"))
        assert metrics["status"]=="COMPLETE"
        theta,phi,et,ep=ffd(Path(metrics["ffd_path"]))
        plus=np.abs(et+1j*ep)**2
        minus=np.abs(et-1j*ep)**2
        rh,lh=(plus,minus) if metrics["ffd_rhcp_convention"]=="et_plus_j_ep" else (minus,plus)
        for hand,power in (("RHCP",rh),("LHCP",lh)):
            row=metrics[hand]
            pattern=analyze(theta,phi,power,
                            peak=(row["peak_theta_deg"],row["peak_phi_deg"]))
            row["sll_first_null_db"]=pattern["SLL_db"]
            row["sll_direction_deg"]=pattern["SLL_dir"]
            row["hpbw_mean_deg"]=pattern["HPBW_mean_deg"]
            row["wide_angle_relative_db"]=(row["wide_angle_max_gain_db_theta60_to90"]-
                                            row["peak_gain_db"])
            row["peak_change_vs_both_db"]=row["peak_gain_db"]-both[hand]["peak_gain_db"]
        solve=json.loads((OUT/f"{state}_solve_status.json").read_text(encoding="utf-8"))
        assert solve["status"]=="COMPLETE"
        metrics["adaptive_passes"]=solve["profile"]["adaptive_passes"]
        metrics["delta_s_per_pass"]=solve["profile"]["delta_s_per_pass"]
        metrics["solve_wall_s"]=solve["solve_wall_s"]
        metrics["peak_memory_gb"]=solve["profile"]["peak_memory_gb"]
        metrics["max_tets"]=solve["profile"]["max_tets"]
        states[state]=metrics
    assert "ON_A" in states
    out={"label":"HFSS_2026_R1_V31_DIODE_V2","frequency_ghz":5.68,"states":states,
         "sll_definition":"Highest same-hand lobe outside direction-dependent first null; if no null, exclude 2x mean HPBW. Above horizon only.",
         "wide_angle_definition":"Maximum same-hand gain for 60 <= theta <= 90 deg."}
    (OUT/"metrics_summary.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print("METRICS_SUMMARY",list(states),flush=True)


if __name__=="__main__":main()
