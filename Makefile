.PHONY: setup test-fast test-fast-par test-full lint format app clean distclean clean-app-outputs

# POSIX default. Windows compute server: `make VENV=.venv313/Scripts <target>` (CPython 3.13.17 venv, docs/STATUS.md).
VENV ?= .venv/bin

setup:
	python3 -m venv .venv
	$(VENV)/pip install --upgrade pip
	$(VENV)/pip install -e ".[dev,app,cad]"
	$(VENV)/pre-commit install

test-fast:
	$(VENV)/python -m pytest -m "not slow" -q

# All logical cores (needs pytest-xdist; one BLAS thread per worker). Known: none of the app polling tests fail on the 64-thread server.
test-fast-par:
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 $(VENV)/python -m pytest -m "not slow" -q -n logical

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

# `clean` plus build/ and dist/. With CONFIRM=1 it also removes the HEAVY regenerable data (the unit-cell
# library solve results/unit_cell_library/{_solve,database,geometry}, ~550 MB, hours to rebuild with
# examples/build_unit_cell_library.py, and the FSS coefficient cache, ~1 h): it prints the paths and the
# total size and asks y/N first. Never touches tracked files, hfss_exports/ or ~/.lwmom.
HEAVY = results/unit_cell_library/_solve results/unit_cell_library/database \
        results/unit_cell_library/geometry scripts/adhoc/fss/coef
distclean: clean
	rm -rf build dist
ifeq ($(CONFIRM),1)
	@du -sh $(wildcard $(HEAVY)); du -ch $(wildcard $(HEAVY)) | tail -1
	@printf "Delete the paths above? [y/N] "; read ans; \
	  if [ "$$ans" = "y" ] || [ "$$ans" = "Y" ]; then rm -rf $(HEAVY); echo removed; else echo kept; fi
else
	@echo "kept $(HEAVY) (slow to regenerate); 'make distclean CONFIRM=1' lists them and asks before deleting"
endif

# Job folders in ~/.lwmom/app_outputs not touched for DAYS days (default 14), except those a registered
# library lives in (library_registry.json). `make clean-app-outputs DRY=--dry-run` only reports.
DAYS ?= 14
clean-app-outputs:
	$(VENV)/python scripts/clean_app_outputs.py --days $(DAYS) $(DRY)
