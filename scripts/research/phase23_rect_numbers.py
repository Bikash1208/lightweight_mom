"""Phase 23 rectangle gate, numeric version: the NUMBERS behind `phase23_rect_snapshot.py`'s hashes.

`phase23_rect_snapshot.py` hashes the complete `/synthesize_ports` job result byte for byte; a hash can only say "same" or "different". This tool dumps every value that hash covers
(same traversal and the same two excluded keys: `patch_layout_png_url`, `achieved_phase_deg`) into a JSON file, so two runs can be compared NUMBER BY NUMBER: exact match, or the max
absolute / relative difference and where it is. Floats are written with `repr` (round-trip exact). Use at 1 BLAS thread (the result is reproducible at a fixed thread count; it is also
thread-count independent on this server).

    python scripts/research/phase23_rect_numbers.py dump    <out.json> [--provenance <out.provenance.json>]
    python scripts/research/phase23_rect_numbers.py compare <reference.json> <new.json>      # exit 0 identical, 1 numeric differences, 2 different structure

Gate (from 2026-10-24, server 1): `results/phase23/rect_snapshot_server1.json` is the reference (with `.provenance.json`). The Mac file `rect_snapshot_after_C.json` stored hashes only, so the
cross-machine difference of this synthesis cannot be measured from it (the unit-cell MoM equivalence Mac vs server was INTERCHANGEABLE, 4.3e-10 relative).
"""

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

SKIP_KEYS = ("patch_layout_png_url", "achieved_phase_deg")


def _flat(obj, path, out):
    """Append (path, value) leaves; value is a float, an `{"nd": ...}` array record, `{"re","im"}` or a plain JSON-able object."""
    if isinstance(obj, np.ndarray):
        out[path] = {"nd": str(obj.dtype), "shape": list(obj.shape), "data": obj.ravel().tolist()}
    elif isinstance(obj, dict):
        for k in sorted(obj, key=str):
            if k in SKIP_KEYS:
                continue
            _flat(obj[k], f"{path}/{k}", out)
    elif isinstance(obj, (list, tuple)):
        out[path + "#len"] = len(obj)
        for i, v in enumerate(obj):
            _flat(v, f"{path}[{i}]", out)
    elif isinstance(obj, (bool, np.bool_)):
        out[path] = bool(obj)
    elif isinstance(obj, (float, np.floating)):
        out[path] = float(obj)
    elif isinstance(obj, complex):
        out[path] = {"re": obj.real, "im": obj.imag}
    elif isinstance(obj, (int, np.integer)):
        out[path] = int(obj)
    else:
        out[path] = repr(obj)  # str, None, anything else: compared by repr, exactly as the hash does


def flatten(result):
    out = {}
    _flat(result, "", out)
    return out


def _num_pairs(a, b):
    """Yield (x, y) float pairs from two leaf values of the same kind, or None when the kinds differ."""
    if isinstance(a, dict) and "nd" in a:
        if not (isinstance(b, dict) and "nd" in b) or a["shape"] != b["shape"] or a["nd"] != b["nd"]:
            return None
        return list(zip(a["data"], b["data"], strict=True))
    if isinstance(a, dict) and "re" in a:
        if not (isinstance(b, dict) and "re" in b):
            return None
        return [(a["re"], b["re"]), (a["im"], b["im"])]
    if isinstance(a, float) and isinstance(b, float):
        return [(a, b)]
    return None


def compare_numbers(ref, new):
    """Compare two `{case: {path: value}}` dumps. Returns dict(identical, n_values, n_float, n_diff, max_abs, max_rel, where_abs, where_rel, structure)."""
    info = {
        "identical": True,
        "n_values": 0,
        "n_float": 0,
        "n_diff": 0,
        "max_abs": 0.0,
        "max_rel": 0.0,
        "where_abs": None,
        "where_rel": None,
        "structure": [],
    }
    if sorted(ref) != sorted(new):
        info["structure"].append(f"cases differ: {sorted(set(ref) ^ set(new))}")
    for case in sorted(set(ref) & set(new)):
        a, b = ref[case], new[case]
        if sorted(a) != sorted(b):
            info["structure"].append(f"{case}: paths differ: {sorted(set(a) ^ set(b))[:5]}")
        for path in sorted(set(a) & set(b)):
            info["n_values"] += 1
            x, y = a[path], b[path]
            pairs = _num_pairs(x, y)
            if pairs is None:
                if x != y:
                    info["n_diff"] += 1
                    info["structure"].append(f"{case}{path}: {str(x)[:40]} != {str(y)[:40]}")
                continue
            for u, v in pairs:
                info["n_float"] += 1
                if (u == v) or (isinstance(u, float) and isinstance(v, float) and math.isnan(u) and math.isnan(v)):
                    continue
                d = abs(u - v)
                r = d / max(abs(u), abs(v), 1e-300)
                info["n_diff"] += 1
                if d > info["max_abs"]:
                    info["max_abs"], info["where_abs"] = d, f"{case}{path}"
                if r > info["max_rel"]:
                    info["max_rel"], info["where_rel"] = r, f"{case}{path}"
    info["identical"] = info["n_diff"] == 0 and not info["structure"]
    return info


def dump(out_path, prov_path=None):
    import hashlib

    import phase23_rect_snapshot as snap

    from src.app.routers.hmsa import _run_synthesize_ports_job
    from src.export.exact_board import exact_board_from_ports_result

    cases = {}
    hashes = {}
    for name, body in snap.CASES.items():
        if body is None:
            continue
        res = _run_synthesize_ports_job(snap._Job(f"p23snap_{name}"), body)
        h = hashlib.sha256()
        snap._digest(res, h)
        eb = exact_board_from_ports_result(res, body["period_mm"], body["port_thickness_mm"])
        h2 = hashlib.sha256()
        for p in eb.patches:
            h2.update(f"{p.port},{p.i},{p.j},{p.x0},{p.y0},{p.x1},{p.y1};".encode())
        h2.update(repr((eb.feeds_mm, eb.port_rects_mm, eb.x_min, eb.x_max, eb.y_min, eb.y_max)).encode())
        cases[name] = flatten(res)
        cases[name]["#exact_board_sha256"] = h2.hexdigest()
        cases[name]["#n_patches"] = len(eb.patches)
        hashes[name] = h.hexdigest()
        print(name, len(cases[name]), "leaves", hashes[name][:16], flush=True)
    doc = {"label": "PHASE23_RECT_NUMBERS", "result_sha256": hashes, "cases": cases}
    Path(out_path).write_text(json.dumps(doc, sort_keys=True, indent=0, allow_nan=True), encoding="utf-8", newline="\n")
    if prov_path:
        import provenance as prov

        rec = prov.finish(
            prov.provenance(
                {"tool": "phase23_rect_numbers.py dump", "out": str(out_path), "n_cases": len(cases)}, workers=1
            )
        )
        Path(prov_path).write_text(json.dumps(rec, indent=1, sort_keys=True), encoding="utf-8", newline="\n")


def compare(ref_path, new_path):
    ref = json.loads(Path(ref_path).read_text(encoding="utf-8"))
    new = json.loads(Path(new_path).read_text(encoding="utf-8"))
    info = compare_numbers(ref["cases"], new["cases"])
    hash_same = ref["result_sha256"] == new["result_sha256"]
    if info["identical"]:
        print(
            f"IDENTICAL: {info['n_values']} values ({info['n_float']} floats) match exactly; result hashes {'equal' if hash_same else 'DIFFER (unexpected)'}"
        )
        return 0
    if info["structure"]:
        print("DIFFERENT STRUCTURE / non-float values:")
        for s in info["structure"][:20]:
            print("  ", s)
    print(
        f"NUMERIC DIFFERENCES: {info['n_diff']} of {info['n_float']} floats differ; max abs {info['max_abs']:.3e} at {info['where_abs']}; max rel {info['max_rel']:.3e} at {info['where_rel']}"
    )
    return 2 if info["structure"] else 1


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "dump":
        prov = sys.argv[sys.argv.index("--provenance") + 1] if "--provenance" in sys.argv else None
        dump(sys.argv[2], prov)
    elif cmd == "compare":
        sys.exit(compare(sys.argv[2], sys.argv[3]))
    else:
        raise SystemExit(__doc__)
