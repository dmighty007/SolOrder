"""Run the code of the beginner's guide and compare what it prints with the page.

Each page's ```python blocks run in order in one namespace, from
docs/guide/data. When a ```text block follows a python block directly, the
printed output must match it line by line; a line holding only "..." in the
expected text stands for any number of output lines.
"""

import contextlib
import io
import os
import re
from pathlib import Path

import pytest

GUIDE = Path(__file__).resolve().parents[1]/"docs/guide"
PAGES = sorted(p.name for p in GUIDE.glob("*.md"))
FENCE = re.compile(r"```(\w+)\n(.*?)```", re.S)
# Blocks that illustrate other file formats or optional backends, not the example data.
SKIP = ('selection="type 1"', "pip install")


def blocks(page):
    """(code, expected output or None) for every python block of a page."""
    fences = [(m.group(1), m.group(2), m.start(), m.end()) for m in FENCE.finditer(page)]
    out = []
    for i, (lang, body, _, end) in enumerate(fences):
        if lang != "python":
            continue
        following = fences[i+1] if i+1 < len(fences) else None
        expected = None
        if following and following[0] == "text" and not page[end:following[2]].strip():
            expected = following[1]
        out.append((body, expected))
    return out


def matches(output, expected):
    """Line-by-line equality, where an expected "..." line skips any number of lines."""
    got = [line.rstrip() for line in output.strip("\n").splitlines()]
    want = [line.rstrip() for line in expected.strip("\n").splitlines()]
    pattern = "".join(r"(?:[^\n]*\n)*?" if line.strip() == "..." else re.escape(line) + r"\n" for line in want)
    return re.fullmatch(pattern, "".join(line + "\n" for line in got)) is not None


@pytest.mark.parametrize("name", PAGES)
def test_guide_page(name, monkeypatch):
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    monkeypatch.setattr(plt, "show", lambda *args, **kwargs: plt.close("all"))
    monkeypatch.chdir(GUIDE/"data")
    namespace = {"__name__": "__guide__"}
    try:
        import ovito  # noqa: F401
        has_ptm = True
    except ImportError:
        has_ptm = False
    for code, expected in blocks((GUIDE/name).read_text()):
        if any(marker in code for marker in SKIP) or ("template_matching" in code and not has_ptm):
            continue
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            exec(compile(code, f"{name}", "exec"), namespace)
        if expected is not None:
            assert matches(buffer.getvalue(), expected), (
                f"{name}: output differs\n--- code\n{code}--- expected\n{expected}--- got\n{buffer.getvalue()}")
    for cache in (GUIDE/"data").glob(".*_offsets.*"):
        os.remove(cache)
