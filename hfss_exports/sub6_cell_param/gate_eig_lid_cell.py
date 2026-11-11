"""SD gate lid-height extension in a newly launched local AEDT instance."""
import argparse
import json
from pathlib import Path

import cell_param_pyaedt as cp
import cell2_cell as base


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", choices=("SD_w2p9_ylines",), required=True)
    ap.add_argument("--project", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--csv", type=Path)
    ap.add_argument("--cores", type=int, default=8)
    ap.add_argument("--lid-height-mm", type=float, choices=(40, 150, 300), required=True)
    ap.add_argument("--project-variant")
    args = ap.parse_args()
    base.configure()
    g = base.geometry_named(args.geometry)
    args.output.mkdir(parents=True, exist_ok=True)
    if args.build and args.project.exists():
        raise RuntimeError("Refusing to overwrite project")
    if not args.build and not args.project.exists():
        raise RuntimeError("Built project missing")
    args.project.parent.mkdir(parents=True, exist_ok=True)
    h = cp.connect(args.project)
    if args.build:
        try:
            print(json.dumps(base.build(h, g, args.output, args.lid_height_mm, "PEC"), indent=2), flush=True)
        finally:
            h.release_desktop()
    elif args.csv:
        base.solve_csv(h, g, args.csv, args.output, args.cores, args.project,
                       args.lid_height_mm, args.project_variant)
    else:
        h.release_desktop()
        raise RuntimeError("Choose --build or --csv")


if __name__ == "__main__":
    main()
