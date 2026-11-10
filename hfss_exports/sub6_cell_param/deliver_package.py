"""Copy reviewable outputs; preserve destination files modified by the user."""
import hashlib
import json
from pathlib import Path
import shutil

SOURCE=Path(__file__).resolve().parent
DEST=Path(__import__("os").environ.get("SUB6_DELIVER_DEST", str(Path(__file__).resolve().parent.parent / "sub6_cell_param_deliver")))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    DEST.mkdir(parents=True,exist_ok=True)
    manifest=DEST/"delivery_manifest.json"
    previous=json.loads(manifest.read_text()) if manifest.exists() else {}
    updated=previous.copy()
    files=[p for p in SOURCE.iterdir() if p.is_file() and p.name not in
           {".gitignore","collected_tests.log","fast_suite.log"}]
    files+=list((SOURCE/"runs").glob("*.json"))
    conflicts=[]
    for p in files:
        rel=p.relative_to(SOURCE).as_posix()
        target=DEST/rel
        if target.exists() and digest(target)!=digest(p) and previous.get(rel)!=digest(target):
            conflicts.append(rel)
            continue
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(p,target)
        updated[rel]=digest(target)
    manifest.write_text(json.dumps(updated,indent=2)+"\n")
    print(json.dumps({"files":len(updated),"preserved_user_modified_conflicts":conflicts},indent=2))
    if conflicts:
        raise SystemExit(1)


if __name__=="__main__":
    main()
