"""Platform-aware comparison against committed numeric references.

The committed references (hfss_exports/, results/phase26/, tests/data/r1_theta0_reference.json) were produced on macOS arm64. On any other CPU/BLAS the same code
differs in the last digits (~1e-15), so a bit-for-bit compare can only hold on the reference platform. Rule: on the reference platform the comparison is exact; elsewhere
every float must agree to `CROSS_PLATFORM_RTOL` (relative) and everything non-numeric (structure, keys, strings, integers, text between the numbers) must still match exactly.
"""

import math
import platform
import re
import sys

ON_REFERENCE_PLATFORM = sys.platform == "darwin" and platform.machine() == "arm64"

# Why 1e-12: BLAS/LAPACK last-digit differences between platforms are ~1e-15 relative (seen on the Windows server vs the macOS references: 30.757828525351485 vs ...282, 28.47791131215888 vs ...874); 1e-12 leaves
# three orders of margin and is still far below any physical or modelling change (a real change is >= 1e-9 relative, which tests/test_platform_tolerance.py checks must fail).
CROSS_PLATFORM_RTOL = 1e-12

# Entries far below the largest value of the same array/file (e.g. the ~1e-16 off-diagonal S-parameters of a symmetric cell, which are rounding noise whose value is arbitrary
# at the last digit) have no meaningful relative difference; they are held to RTOL * (largest magnitude) instead.
NOISE_FLOOR_REL = 1e-9

_NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def _is_float_token(tok):
    return any(c in tok for c in ".eE")


def _check_pairs(pairs, ctx, rtol):
    scale = max((abs(r) for _, r in pairs if math.isfinite(r)), default=0.0)
    worst = 0.0
    for got, ref in pairs:
        if got == ref or (math.isnan(got) and math.isnan(ref)):
            continue
        diff = abs(got - ref)
        mag = max(abs(got), abs(ref))
        allowed = rtol * (mag if mag >= NOISE_FLOOR_REL * scale else scale)
        worst = max(worst, diff / allowed if allowed > 0 else math.inf)
        assert diff <= allowed, f"{ctx}: {got!r} vs reference {ref!r} (|diff| {diff:.3e}, allowed {allowed:.3e})"
    return worst


def _pairs(got, ref, ctx, out):
    if isinstance(ref, dict):
        assert isinstance(got, dict) and got.keys() == ref.keys(), f"{ctx}: keys differ"
        for k in ref:
            _pairs(got[k], ref[k], f"{ctx}/{k}", out)
    elif isinstance(ref, (list, tuple)):
        assert isinstance(got, (list, tuple)) and len(got) == len(ref), f"{ctx}: length differs"
        for i, (g, r) in enumerate(zip(got, ref, strict=True)):
            _pairs(g, r, f"{ctx}[{i}]", out)
    elif isinstance(ref, float):
        assert isinstance(got, float), f"{ctx}: {got!r} is not a float"
        out.append((got, ref))
    else:  # str, int, bool, None: exact (bool/int are never loosened)
        assert type(got) is type(ref) and got == ref, f"{ctx}: {got!r} != {ref!r}"


def assert_matches_reference(got, ref, ctx="", *, reference_platform=ON_REFERENCE_PLATFORM, rtol=CROSS_PLATFORM_RTOL):
    """`got` (nested dict / list / float / str / int ...) against the committed `ref`: exact on the reference platform, within `rtol` elsewhere."""
    if reference_platform:
        assert got == ref, f"{ctx}: differs from the committed reference (reference platform: exact compare)"
        return
    pairs = []
    _pairs(got, ref, ctx, pairs)
    _check_pairs(pairs, ctx, rtol)


def assert_text_matches_reference(
    got, ref, ctx="", *, reference_platform=ON_REFERENCE_PLATFORM, rtol=CROSS_PLATFORM_RTOL
):
    """Byte-compare of a generated text file (bytes or str): exact on the reference platform; elsewhere the numbers are parsed and compared within `rtol`, the text between
    them and every integer token must match exactly."""
    if isinstance(got, bytes):
        got = got.decode("utf-8", "surrogateescape")
    if isinstance(ref, bytes):
        ref = ref.decode("utf-8", "surrogateescape")
    if reference_platform:
        assert got == ref, f"{ctx}: differs from the committed reference (reference platform: exact compare)"
        return
    assert _NUMBER.sub("#", got) == _NUMBER.sub("#", ref), f"{ctx}: text between the numbers differs"
    gt, rt = _NUMBER.findall(got), _NUMBER.findall(ref)
    assert len(gt) == len(rt), f"{ctx}: number of numeric tokens differs"
    pairs = []
    for g, r in zip(gt, rt, strict=True):
        if _is_float_token(r) or _is_float_token(g):
            pairs.append((float(g), float(r)))
        else:
            assert g == r, f"{ctx}: integer token {g} != {r}"
    _check_pairs(pairs, ctx, rtol)
