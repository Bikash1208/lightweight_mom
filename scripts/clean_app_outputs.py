"""Remove old job folders from the app's output directory (`make clean-app-outputs`).

Folders not modified for `--days` days are removed, EXCEPT any folder that
`app_outputs/library_registry.json` points at (a registered library lives inside its job folder).
`--dry-run` only reports. Data dir: `$LWMOM_APP_DATA_DIR`, default `~/.lwmom`."""

import argparse
import json
import os
import shutil
import time
from pathlib import Path


def _size(path):
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def referenced_job_ids(root):
    reg = root / "library_registry.json"
    if not reg.is_file():
        return set()
    text = reg.read_text()
    ids = set()
    for d in root.iterdir():
        if d.is_dir() and d.name in text:
            ids.add(d.name)
    return ids


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--days", type=float, default=14.0)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    data = Path(os.environ.get("LWMOM_APP_DATA_DIR") or Path.home() / ".lwmom")
    root = data / "app_outputs"
    if not root.is_dir():
        print(f"nothing to do: {root} does not exist")
        return
    keep = referenced_job_ids(root)
    cutoff = time.time() - args.days * 86400.0
    n = total = kept = 0
    for d in sorted(root.iterdir()):
        if not d.is_dir() or d.stat().st_mtime >= cutoff:
            continue
        if d.name in keep:
            kept += 1
            continue
        total += _size(d)
        n += 1
        if not args.dry_run:
            shutil.rmtree(d)
    verb = "would remove" if args.dry_run else "removed"
    print(f"{verb} {n} folder(s), {total / 1e6:.0f} MB; kept {kept} old folder(s) referenced by library_registry.json")


if __name__ == "__main__":
    main()
