.PHONY: setup test-fast test-full lint format app clean distclean clean-app-outputs

VENV := .venv/bin

setup:
	python3 -m venv .venv
	$(VENV)/pip install --upgrade pip
	$(VENV)/pip install -e ".[dev,app,cad]"
	$(VENV)/pre-commit install

test-fast:
	$(VENV)/python -m pytest -m "not slow" -q

test-full:
	$(VENV)/python -m pytest -q

lint:
	$(VENV)/ruff check src tests examples scripts
	$(VENV)/ruff format --check src tests examples scripts

format:
	$(VENV)/ruff check --fix src tests examples scripts
	$(VENV)/ruff format src tests examples scripts

app:
	@echo "Restart this command after pulling backend (src/) changes -- the frontend is static files," \
	      "reloaded on refresh, but the FastAPI process itself is not."
	$(VENV)/python -m src.app.main

clean:
	find . -type d -name "__pycache__" -not -path "./.venv/*" -exec rm -rf {} +
	rm -rf .pytest_cache .ruff_cache

# Everything `clean` removes plus the build/venv leftovers and the HEAVY regenerable data. The unit-cell
# library solve (results/unit_cell_library/_solve|database|geometry, ~550 MB, hours to rebuild with
# examples/build_unit_cell_library.py) and the FSS coefficient cache (~1 h) are only removed with
# CONFIRM=1 (docs/DEVELOPING.md "Cleaning"). Never touches tracked files, hfss_exports/ or ~/.lwmom.
distclean: clean
	rm -rf lwmom.egg-info build dist
	@if [ "$(CONFIRM)" = "1" ]; then \
	  rm -rf results/unit_cell_library/_solve results/unit_cell_library/database \
	         results/unit_cell_library/geometry scripts/adhoc/fss/coef; \
	  echo "removed the unit-cell library solve and the FSS coefficient cache"; \
	else \
	  echo "kept results/unit_cell_library/{_solve,database,geometry} and scripts/adhoc/fss/coef" \
	       "(slow to regenerate); run 'make distclean CONFIRM=1' to remove them too"; \
	fi

# The app's own job outputs under ~/.lwmom/app_outputs (one folder per job). Removes folders not touched
# for DAYS days (default 14); the downloads of a recent run are kept.
DAYS ?= 14
clean-app-outputs:
	find "$${LWMOM_APP_DATA_DIR:-$$HOME/.lwmom}/app_outputs" -mindepth 1 -maxdepth 1 -type d -mtime +$(DAYS) -exec rm -rf {} +
