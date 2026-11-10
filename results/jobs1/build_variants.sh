#!/bin/bash
cd /c/lwmom_final/scripts/research
PY=/c/lightweight_mom/.venv313/Scripts/python.exe
export RF_OUT_DIR=/c/lwmom_final/hfss_exports/sub6_jobs1_variants
run() { name=$1; shift; env "$@" RF_OUT_NAME=$name $PY sub6_final_rf_std.py test9 > /c/lwmom_final/results/jobs1/build_$name.log 2>&1; echo "$name rc=$?"; }
run BOARD_V0_baseline.step RF_VIA_MODE=base
run BOARD_V1_fr4_0p2.step RF_VIA_MODE=base RF_FR4_T=0.2
run BOARD_V2_antipad_0p6.step RF_VIA_MODE=base RF_R_ANTI=0.3
run BOARD_V4_shielded.step RF_VIA_MODE=shield
run BOARD_V5_novias.step RF_VIA_MODE=none
echo BUILD_DONE
