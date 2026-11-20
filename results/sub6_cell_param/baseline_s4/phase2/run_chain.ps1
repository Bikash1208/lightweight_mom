param([string]$Tag,[string]$Step,[string]$State)
$env:PYTHONDONTWRITEBYTECODE=1
$env:ANSYSLMD_LICENSE_FILE='1055@10.100.105.13'
$py='C:\work\lwmom_s4\.venv\Scripts\python.exe'; $sc='C:\work\lwmom_s4\scripts\research'
$out="C:\work\s4_work\$Tag"; New-Item -ItemType Directory -Force $out | Out-Null
$t=Get-Date
& $py "$sc\full_board_gate_v31_hfss.py" --step $Step --out $out --state $State *> "$out\hfss.log"
"hfss exit $LASTEXITCODE  $((Get-Date)-$t)" | Out-File "$out\chain.log"
if($LASTEXITCODE -ne 0){ "FAILED_HFSS" | Out-File -Append "$out\chain.log"; exit 1 }
$proj="$out\BOARD_R150_GATE_V31_$State.aedt"
& $py "$sc\full_board_gate_v31_extract.py" --project $proj --state $State --out "$out\extract" *> "$out\extract.log"
"extract exit $LASTEXITCODE" | Out-File -Append "$out\chain.log"
& $py "$sc\full_board_gate_v31_profiles.py" --state $State --project $proj --out "$out\profile" *> "$out\profile.log"
"profile exit $LASTEXITCODE" | Out-File -Append "$out\chain.log"
"DONE" | Out-File -Append "$out\chain.log"
