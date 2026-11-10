"""Fixed BASELINE/G1 mesh trend, no optimization; MODEL_ONLY."""
import argparse
import hashlib
import itertools
import json
import os
import platform
import sys
from datetime import datetime,timezone
from multiprocessing import Pool
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["RETUNE_FFT"] = "1"
os.environ["RETUNE_WOODBURY"] = "1"
sys.path.insert(0,str(Path(__file__).resolve().parent))
import sub6_iter1_island_retune as R

OUT = R.ROOT/"results/sub6_island_meshtrend"
OLD = R.ROOT/"results/sub6_iter1_island_retune"
G1 = json.loads((OLD/"selected.json").read_text())[0]["spec"]
SPECS = {"BASELINE":R.BASE,"G1":G1}


def report(rows):
    result = []
    groups = {}
    for row in rows:
        groups.setdefault((R.key(row["spec"]),row["mesh"]),[]).append(row)
    for (ident,n),nodes in groups.items():
        if len({(r["kt"],r["direction"]) for r in nodes}) < 4:
            continue
        name = next(k for k,v in SPECS.items() if R.key(v)==ident)
        try:
            summary = R.summaries(nodes)[0]
            summary["cell"] = name
            result.append(summary)
        except Exception as error:
            result.append(dict(cell=name,mesh=n,status="SHEET_EXTRACTION_FAILED",error=str(error)))
    data = dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED; alpha is first-order leakage, not eigenmode-Q loss",
                no_optimization=True,rows=sorted(result,key=lambda r:(r["cell"],r["mesh"])))
    (OUT/"mesh_trend.json").write_text(json.dumps(data,indent=2),encoding="utf-8",newline="\n")
    for r in data["rows"]:
        if "states" in r:
            a,b = r["states"]["nominal|OFF"],r["states"]["nominal|ON"]
            print(r["cell"],r["mesh"],"alpha",a["alpha_np_per_m"],b["alpha_np_per_m"],"supp",a["suppression_dB"],"beta",a["beta_over_k0"],b["beta_over_k0"],"extrap",a["extrapolated"],b["extrapolated"],flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stage",choices=("initial","192","report"))
    args = parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    path = OUT/"nodes.jsonl"
    rows = [json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []
    if not rows:
        for file in ("search_nodes.jsonl","refine_nodes.jsonl"):
            for line in (OLD/file).read_text().splitlines():
                row = json.loads(line)
                if row["floquet_order"]!=2*row["mesh"] or R.key(row["spec"]) not in {R.key(v) for v in SPECS.values()}:
                    continue
                ident = (R.key(row["spec"]),row["mesh"],row["kt"],row["direction"])
                if ident in {(R.key(r["spec"]),r["mesh"],r["kt"],r["direction"]) for r in rows}:
                    continue
                row["tensors"] = {k:v for k,v in row["tensors"].items() if k in ("nominal|OFF","nominal|ON")}
                row["reused_from"] = str((OLD/file).relative_to(R.ROOT))
                rows.append(row)
        path.write_text("".join(json.dumps(r)+"\n" for r in rows),encoding="utf-8",newline="\n")
    if args.stage == "report":
        report(rows)
        return
    meshes = [96,128,160] if args.stage == "initial" else [192]
    done = {(R.key(r["spec"]),r["mesh"],r["kt"],r["direction"]) for r in rows}
    tasks = [(s,n,5.8e9,kt,d,2*n,False) for s,n,kt,d in itertools.product(SPECS.values(),meshes,(1.017,1.097),"xy")
             if (R.key(s),n,kt,d) not in done]
    provenance = dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED",stage=args.stage,
        utc=datetime.now(timezone.utc).isoformat(),machine=platform.node(),python=sys.version,
        executable=sys.executable,numpy=R.np.__version__,workers=os.cpu_count(),blas_threads=1,
        mesh_gap_cells={str(n):round(n/11.6) for n in meshes},parameters=SPECS,
        node_interval=[1.017,1.097],floquet_order="2n",tasks=len(tasks),
        driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/(args.stage+"_provenance.json")).write_text(json.dumps(provenance,indent=2),encoding="utf-8",newline="\n")
    print(json.dumps(provenance),flush=True)
    with path.open("a",encoding="utf-8",newline="\n") as stream,Pool() as pool:
        for row in pool.imap_unordered(R.node,tasks):
            rows.append(row)
            stream.write(json.dumps(row)+"\n");stream.flush()
            print("done",row["mesh"],row["kt"],row["direction"],round(row["seconds"],1),flush=True)
    report(rows)


if __name__ == "__main__":
    main()
