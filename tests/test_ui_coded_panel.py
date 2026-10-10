"""Real-browser checks of the Coded panel (Playwright). Skipped when Playwright or its Chromium build is
not installed; the static-HTML test below always runs."""

import re
import socket
import threading
import time
from pathlib import Path

import pytest

INDEX = Path(__file__).resolve().parent.parent / "src" / "app" / "static" / "index.html"


def test_advanced_section_is_not_hidden_by_bit_count_in_the_page_source():
    """Static guard: only the 'fix the square's side' block depends on the bit count, never the whole Advanced section."""
    html = INDEX.read_text()
    assert 'id="coded-auto-fixed-wrap"' in html
    assert re.search(r'getElementById\("coded-auto-advanced"\)\.style\.display', html) is None
    assert re.search(r'getElementById\("coded-auto-fixed-wrap"\)\.style\.display\s*=\s*oneBit', html)
    assert html.index('id="coded-auto-min-gap"') > html.index('id="coded-auto-advanced"')
    # the gap field sits outside the one-bit-only block
    wrap = html[html.index('id="coded-auto-fixed-wrap"') : html.index('id="coded-auto-min-gap"')]
    assert "</div>" in wrap


@pytest.fixture(scope="module")
def live_server(tmp_path_factory):
    uvicorn = pytest.importorskip("uvicorn")
    import os

    os.environ["LWMOM_APP_DATA_DIR"] = str(tmp_path_factory.mktemp("ui_appdata"))
    from src.app.main import create_app

    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]
    server = uvicorn.Server(uvicorn.Config(create_app(), host="127.0.0.1", port=port, log_level="error"))
    thread = threading.Thread(target=server.run, daemon=True)
    thread.start()
    deadline = time.time() + 20
    while not server.started and time.time() < deadline:
        time.sleep(0.05)
    yield f"http://127.0.0.1:{port}/"
    server.should_exit = True
    thread.join(timeout=5)


@pytest.fixture()
def page(live_server):
    sync_api = pytest.importorskip("playwright.sync_api")
    with sync_api.sync_playwright() as p:
        try:
            browser = p.chromium.launch()
        except Exception as exc:  # the matching build is not installed: fall back to any cached Chromium
            cached = sorted(Path.home().glob("Library/Caches/ms-playwright/chromium-*/chrome-mac*/Google Chrome for Testing.app/Contents/MacOS/Google Chrome for Testing")) \
                + sorted(Path.home().glob(".cache/ms-playwright/chromium-*/chrome-linux*/chrome"))
            if not cached:
                pytest.skip(f"Chromium not available for Playwright: {exc}")
            browser = p.chromium.launch(executable_path=str(cached[-1]))
        pg = browser.new_page(viewport={"width": 1300, "height": 1600})
        pg.goto(live_server)
        pg.click("#skip-to-multiport-btn")
        pg.click('#method-control button[data-method="coded"]')
        yield pg
        browser.close()


def test_min_gap_is_visible_for_every_bit_count_and_fixed_side_only_for_one_bit(page):
    for bits, fixed_visible in ((1, True), (2, False), (3, False), (1, True)):
        page.click(f'#coded-bits-control button[data-bits="{bits}"]')
        page.click("#coded-auto-advanced > summary") if not page.get_attribute("#coded-auto-advanced", "open") is not None else None
        assert page.is_visible("#coded-auto-advanced"), f"{bits}-bit: Advanced section hidden"
        page.evaluate("document.getElementById('coded-auto-advanced').open = true")
        assert page.is_visible("#coded-auto-min-gap"), f"{bits}-bit: minimum gap field hidden"
        assert page.is_visible("#coded-auto-fixed-L") == fixed_visible, f"{bits}-bit: fixed-side visibility"
