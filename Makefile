.PHONY: setup test-fast test-full lint format app clean

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
	$(VENV)/ruff check src tests examples
	$(VENV)/ruff format --check src tests examples

format:
	$(VENV)/ruff check --fix src tests examples
	$(VENV)/ruff format src tests examples

app:
	$(VENV)/python -m src.app.main

clean:
	find . -type d -name "__pycache__" -not -path "./.venv/*" -exec rm -rf {} +
	rm -rf .pytest_cache .ruff_cache
