PYV=C:/Users/Bikash/Desktop/lwmom_labB/.venv_lab/Scripts/python.exe
B=C:/Users/Bikash/Desktop/sub6_board_hfss/board_v2/crop
kill_own() { pid=$(grep -o "OWN_AEDT_PID [0-9]*" $1 | tail -1 | cut -d' ' -f2); powershell -NoProfile -Command "Stop-Process -Id $pid -Force -ErrorAction SilentlyContinue; Start-Sleep 2; Remove-Item $2 -Force -ErrorAction SilentlyContinue"; }
for V in A0 A90 B0 B90; do
  P=$B/hfss_br$V/CROP_BR$V.aedt; mkdir -p $B/hfss_br$V
  timeout 1800 $PYV hfss_build_crop_br.py $B/out_br$V/BOARD_crop_$V.step $P $B/hfss_br$V/readback.json 70 100 0 21 sw > $B/hfss_br$V/build.log 2>&1; kill_own $B/hfss_br$V/build.log $P.lock
  if grep -q "^DONE" $B/hfss_br$V/build.log && ! grep -q Traceback $B/hfss_br$V/build.log; then
    timeout 900 $PYV mk_short.py $P > $B/hfss_br$V/mk_short.log 2>&1; kill_own $B/hfss_br$V/mk_short.log $P.lock
    timeout 3000 $PYV solve_crop.py $P SHORT_IDEAL $B/res_br$V > $B/res_br$V.log 2>&1; kill_own $B/res_br$V.log $P.lock
    timeout 900 $PYV post_crop.py $P SHORT_IDEAL "Rs='0.001ohm' Xs='0ohm'" $B/res_br$V/br$V > $B/post_br$V.log 2>&1; kill_own $B/post_br$V.log $P.lock
  fi
  echo "$V done" >> $B/br.progress
done
echo BR_DONE > $B/br.done
