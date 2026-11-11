#!/bin/bash
# one worker per pending row (no worker cap), 10 s stagger, unique gRPC ports
cd "$(dirname "$0")"
PY=/c/Users/Bikash/sub6_runtime/pyaedt312/Scripts/python.exe
i=0
for id in $($PY pending.py); do
  i=$((i+1))
  SUB6_AEDT_VERSION=2025.1 nohup $PY gate_hfss.py worker --run-id $id --port $((50400+i)) --cores 3 > gate_work/logs/$id.out 2> gate_work/logs/$id.err &
  sleep 10
done
wait
