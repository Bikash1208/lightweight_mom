"""Phase 13 Part 2: the run manifest -- one JSON record per job, saved alongside its outputs (and
bundled into the app's "download diagnostic bundle" action with logs and inputs). Captures enough
to reproduce a result and to file a useful bug report: inputs, resolved config, code version/git
hash, dependency versions, timings, labels, warnings, output file hashes.
"""

import hashlib
import importlib.metadata
import json
import platform
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent

_TRACKED_DEPENDENCIES = ("numpy", "scipy", "PyYAML", "fastapi", "uvicorn", "cadquery", "matplotlib")


def git_hash():
    """The repo's current commit hash, or `None` if not a git checkout / git unavailable (never
    raises -- a manifest without provenance is still useful, an app crash over it is not)."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=_REPO_ROOT, capture_output=True, text=True, timeout=5,
        )
        return out.stdout.strip() if out.returncode == 0 else None
    except Exception:
        return None


def git_dirty():
    """`True` if the working tree has uncommitted changes relative to `git_hash()`, `None` if
    that can't be determined."""
    try:
        out = subprocess.run(
            ["git", "status", "--porcelain"], cwd=_REPO_ROOT, capture_output=True, text=True, timeout=5,
        )
        return bool(out.stdout.strip()) if out.returncode == 0 else None
    except Exception:
        return None


def dependency_versions(names=_TRACKED_DEPENDENCIES):
    versions = {}
    for name in names:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def hash_file(path):
    path = Path(path)
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass
class RunManifest:
    run_id: str
    job_kind: str                    # e.g. "unit_cell_solve", "reflectarray_design", "hmsa_design"
    inputs: dict = field(default_factory=dict)
    resolved_config: dict = field(default_factory=dict)
    labels: tuple = ()
    warnings: tuple = ()
    timings_s: dict = field(default_factory=dict)
    output_files: dict = field(default_factory=dict)   # {relative path: sha256 hash}
    code_git_hash: "str | None" = field(default=None)
    code_git_dirty: "bool | None" = field(default=None)
    dependency_versions: dict = field(default_factory=dict)
    python_version: str = field(default_factory=platform.python_version)
    platform: str = field(default_factory=platform.platform)

    def __post_init__(self):
        if self.code_git_hash is None:
            self.code_git_hash = git_hash()
        if self.code_git_dirty is None:
            self.code_git_dirty = git_dirty()
        if not self.dependency_versions:
            self.dependency_versions = dependency_versions()

    def add_output_file(self, relative_path, absolute_path):
        self.output_files[str(relative_path)] = hash_file(absolute_path)

    def to_dict(self):
        return {
            "run_id": self.run_id, "job_kind": self.job_kind, "inputs": self.inputs,
            "resolved_config": self.resolved_config, "labels": list(self.labels),
            "warnings": list(self.warnings), "timings_s": self.timings_s,
            "output_files": self.output_files, "code_git_hash": self.code_git_hash,
            "code_git_dirty": self.code_git_dirty, "dependency_versions": self.dependency_versions,
            "python_version": self.python_version, "platform": self.platform,
        }

    def write(self, path):
        Path(path).write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True, default=str))
        return Path(path)


def build_diagnostic_bundle(manifest: RunManifest, log_file_path, input_files, output_dir):
    """Zips the manifest, the relevant slice of the log file (lines matching `manifest.run_id`),
    and any `input_files` into one archive under `output_dir` -- the "download diagnostic bundle"
    action. Returns the archive's path."""
    import zipfile

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    bundle_path = output_dir / f"diagnostic_bundle_{manifest.run_id}.zip"
    with zipfile.ZipFile(bundle_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("manifest.json", json.dumps(manifest.to_dict(), indent=2, sort_keys=True, default=str))
        log_file_path = Path(log_file_path)
        if log_file_path.exists():
            lines = [
                line for line in log_file_path.read_text(errors="replace").splitlines()
                if manifest.run_id in line
            ]
            zf.writestr("log_excerpt.jsonl", "\n".join(lines))
        for input_file in input_files:
            input_file = Path(input_file)
            if input_file.exists():
                zf.write(input_file, arcname=f"inputs/{input_file.name}")
    return bundle_path
