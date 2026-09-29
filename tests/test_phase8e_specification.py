"""Sec. 24: output writers round-trip correctly."""

import csv
import json

import numpy as np

from src.phase8b.object_wave import transverse_wavenumber_deg
from src.phase8c.scattering_analysis import ScatteringRow
from src.phase8e.candidate_library_period import PeriodRequirementAssessment
from src.phase8e.specification import (
    write_library_response_map_csv,
    write_object_spectrum_csv,
    write_reference_spectrum_csv,
    write_required_library_specification_json,
    write_required_response_map_csv,
)
from src.phase8e.spectral import ObjectWaveSpectralPoint, ObjectWaveSpectrum


def test_write_reference_spectrum_csv_row_count(tmp_path):
    kx, ky = np.meshgrid([1.0, 2.0], [3.0, 4.0])
    power = np.ones_like(kx)
    path = write_reference_spectrum_csv(kx, ky, power, tmp_path / "ref.csv")
    with open(path) as fh:
        rows = list(csv.reader(fh))
    assert len(rows) == 1 + kx.size


def test_write_object_spectrum_csv_one_row_per_beam(tmp_path):
    kx0, ky0 = transverse_wavenumber_deg(12.5e9, 10.0, 0.0)
    points = (ObjectWaveSpectralPoint(10.0, 0.0, kx0, ky0, 100.0 + 0j, True),)
    spectrum = ObjectWaveSpectrum(beams=points, k_required_max_per_m=None, all_beams_propagating=True)
    path = write_object_spectrum_csv(spectrum, tmp_path / "obj.csv")
    with open(path) as fh:
        rows = list(csv.reader(fh))
    assert len(rows) == 2  # header + 1 beam


def test_write_required_response_map_csv_combines_multiple_formulations(tmp_path):
    x, y = np.meshgrid([0.0, 1e-3], [0.0, 1e-3])
    path = write_required_response_map_csv(
        [("B", x, y, np.array([[0.0, 90.0], [180.0, -90.0]])), ("C", x, y, np.array([[0.0, 90.0], [180.0, -90.0]]))],
        tmp_path / "req.csv",
    )
    with open(path) as fh:
        rows = list(csv.reader(fh))
    assert len(rows) == 1 + 8  # header + 4 points x 2 formulations
    assert {r[2] for r in rows[1:]} == {"B", "C"}


def test_write_library_response_map_csv_combines_components(tmp_path):
    rows_s11 = [ScatteringRow(1.0, 1.0, 0.1, 0.2, 0.22, -13.0, 63.0)]
    rows_s21 = [ScatteringRow(1.0, 1.0, 0.9, 0.1, 0.91, -0.8, 6.3)]
    path = write_library_response_map_csv([(rows_s11, "S11_xx"), (rows_s21, "S21_xx")], tmp_path / "lib.csv")
    with open(path) as fh:
        rows = list(csv.reader(fh))
    assert len(rows) == 3  # header + 2 rows
    assert {r[2] for r in rows[1:]} == {"S11_xx", "S21_xx"}


def test_write_required_library_specification_json_round_trips(tmp_path):
    assessment = PeriodRequirementAssessment(
        registered_periods_mm=((4.8, 4.8),),
        candidate_periods_mm=(5.0, 4.0),
        registered_period_matches_candidates=True,
        status="CURRENT_DATASET_LIMITATION",
        required_library_specification=None,
    )
    path = write_required_library_specification_json(assessment, tmp_path / "spec.json")
    doc = json.loads(path.read_text())
    assert doc["status"] == "CURRENT_DATASET_LIMITATION"
    assert doc["registered_periods_mm"] == [[4.8, 4.8]]
