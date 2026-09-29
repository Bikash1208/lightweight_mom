"""Canonical serialisation and hashing, so identifiers do not depend on dict order, the
platform or PYTHONHASHSEED. Floats are written with ``repr`` (shortest round-trip form)."""

import hashlib
import json


def canonical_json(obj):
    """Deterministic compact JSON: sorted keys, no NaN/inf (they have no JSON spelling)."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=False, ensure_ascii=True)


def short_hash(obj, length=16):
    """First ``length`` hex digits of the SHA-256 of ``canonical_json(obj)``."""
    return hashlib.sha256(canonical_json(obj).encode("ascii")).hexdigest()[:length]


def cell_id(context, parameters):
    """Identifier of one unit cell: everything that defines the solved problem except the frequency."""
    return short_hash({"context": context, "parameters": parameters})
