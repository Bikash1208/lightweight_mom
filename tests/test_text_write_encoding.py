"""Windows portability: every text write under src/ names its encoding.

Frozen report-writer packages (phase8b ... phase13) are excluded; they are a known Windows limitation listed in docs/STATUS.md.
"""

import ast
import re
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"
FROZEN = re.compile(r"phase(8b|8c|8d_rev|8e|9|10|11|12[a-d]|13)$")


def _offenders():
    bad = []
    for p in sorted(SRC.rglob("*.py")):
        rel = p.relative_to(SRC)
        if rel.parts and FROZEN.match(rel.parts[0]):
            continue
        for n in ast.walk(ast.parse(p.read_text(encoding="utf-8"))):
            if not isinstance(n, ast.Call):
                continue
            name = n.func.attr if isinstance(n.func, ast.Attribute) else getattr(n.func, "id", "")
            kw = {k.arg: k.value for k in n.keywords}
            if name == "write_text" and "encoding" not in kw:
                bad.append(f"{rel}:{n.lineno} write_text without encoding")
            elif name == "open":
                mode = n.args[1].value if len(n.args) > 1 and isinstance(n.args[1], ast.Constant) else getattr(kw.get("mode"), "value", None)
                if isinstance(mode, str) and ("w" in mode or "a" in mode) and "b" not in mode and "encoding" not in kw:
                    bad.append(f"{rel}:{n.lineno} open({mode!r}) without encoding")
    return bad


def test_no_unencoded_text_writes_in_src():
    assert _offenders() == []
