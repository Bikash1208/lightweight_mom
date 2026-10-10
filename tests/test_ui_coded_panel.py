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


def test_every_two_state_button_has_text_in_1_2_and_3_bit_modes(page):
    """Regression: data-label-1 is not camel-cased by `dataset` (it read undefined and blanked the button)."""
    for bits in (1, 2, 3, 1):
        page.click(f'#coded-bits-control button[data-bits="{bits}"]')
        texts = page.eval_on_selector_all(
            "#coded-two-state-control button", "els => els.map(e => [e.textContent.trim(), e.style.display])"
        )
        assert texts and all(t for t, _ in texts), f"{bits}-bit: blank button in {texts}"
        visible = [t for t, display in texts if display != "none"]
        assert len(visible) == (4 if bits == 1 else 2)


def test_compute_cell_states_needs_cell_size_and_frequency_and_sends_nothing(page):
    requests = []
    page.on("request", lambda r: requests.append(r.url) if "/api/hmsa/auto_states" in r.url else None)
    page.fill("#cell-size-mm", "")
    page.fill("#synthesis-frequency-ghz", "")
    page.click("#coded-auto-compute-btn")
    assert page.inner_text("#coded-auto-status") == "Set Cell size and Frequency in the Board section first."
    assert page.evaluate("document.activeElement.id") == "cell-size-mm"
    page.fill("#cell-size-mm", "6.4")
    page.click("#coded-auto-compute-btn")
    assert page.evaluate("document.activeElement.id") == "synthesis-frequency-ghz"
    assert requests == []


def test_422_details_are_formatted_as_readable_text(page):
    text = page.evaluate(
        """LWMOM.errorText([{loc: ["body", "period_mm"], msg: "Input should be a valid number", type: "float_type"},
                           {loc: ["body", "frequency_hz"], msg: "Field required", type: "missing"}])"""
    )
    assert text == "Cell size: must be a number; Frequency: is required"
    assert page.evaluate('LWMOM.errorText("plain")') == "plain"
    assert page.evaluate('LWMOM.errorText({message: "boom"})') == "boom"
    # a real 422 from the API, through the shared postJson
    message = page.evaluate(
        """async () => { try { await LWMOM.postJson("/api/hmsa/auto_states",
              {eps_r: 2.2, tan_delta: 0.0009, thickness_mm: 1.575, period_mm: null, frequency_hz: null, n_bits: 1});
              return "no error"; } catch (e) { return e.message; } }"""
    )
    assert "[object" not in message
    assert "Cell size: must be a number" in message and "Frequency: must be a number" in message


def test_auto_table_formats_sides_with_two_decimals_and_flags_a_weak_fixed_side(page):
    auto = {
        "n_bits": 1, "fixed_L_mm": 4.0, "warnings": ["Weak pair: 4.00 mm gives 12° separation (below 150°); leave the fixed side blank to search."],
        "solve_seconds": 0.0, "n_solves": 0, "n_cached": 14,
        "states": [
            {"level": 0, "label": "empty", "geometry": {"kind": "empty"}, "magnitude_db": 0.0, "phase_deg_patch_plane": 138.1,
             "target_phase_deg": 138.1, "separation_deg": None, "sensitivity_deg_per_0p01mm": None, "mesh_uncertainty_deg": None, "labels": []},
            {"level": 1, "label": "square", "geometry": {"kind": "rectangle", "L_mm": 6.0, "W_mm": 6.0}, "magnitude_db": -0.02,
             "phase_deg_patch_plane": 150.0, "target_phase_deg": -41.9, "separation_deg": 11.9, "sensitivity_deg_per_0p01mm": 0.1,
             "mesh_uncertainty_deg": 5.0, "labels": []},
        ],
    }
    html = page.evaluate("(a) => autoStatesTableHtml(a, true)", auto)
    assert "6.00 × 6.00 mm square" in html and "6 × 6 mm" not in html
    assert "Separation with the fixed 4.00 mm square: 11.9°" in html
    assert "Weak pair: 4.00 mm gives 12° separation" in html
