#!/bin/bash
# staged Freerouting: stage 1 = link nets N*, stage 2 = + drive leads S*, stage 3 = + RETURN.  Previous stage tracks are locked.
cd /c/lwmom_final/results/sub6_final_board/mesh96_13pad/kicad
PY=/c/lightweight_mom/.venv313/Scripts/python.exe
KP=/c/kicad_tools/kicad/bin/python.exe
J=/c/kicad_tools/jdk25/jdk-25.0.5.0.2+1/bin/java.exe
SC=/c/lwmom_final/scripts/research
KD="C:/lwmom_final/results/sub6_final_board/mesh96_13pad/kicad"
PASSES=${PASSES:-25}
prev=""
for st in "N" "NS" "NSR"; do
  echo "=== stage $st"
  $PY $SC/sub6_final_kicad_build.py $st $prev > /dev/null || exit 1
  cp sub6_board.kicad_pro sub6_board_route.kicad_pro
  $KP $SC/sub6_kicad_dsn.py $KD | tail -1
  rm -f sub6_board_route.ses
  $J -Xmx24g -jar /c/kicad_tools/dl/fr.jar --gui.enabled=false -de sub6_board_route.dsn -do sub6_board_route.ses -mp $PASSES > fr_$st.log 2>&1
  grep "Auto-routing pass #" fr_$st.log | tail -1 | grep -o "pass #[0-9]*\|([0-9]* unrouted" | tr '\n' ' '; echo
  $KP $SC/sub6_kicad_stage_post.py $KD tracks_$st.json | tail -2
  prev=tracks_$st.json
done
echo STAGES_DONE
