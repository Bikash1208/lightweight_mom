"""Phase 13 Part 2: logging setup, run-ID propagation, JSON file formatting."""

import json
import logging

from src import logging_setup as ls


def test_configure_logging_is_idempotent(tmp_path):
    logger1 = ls.configure_logging(level="DEBUG", log_dir=tmp_path)
    n_handlers = len(logger1.handlers)
    logger2 = ls.configure_logging(level="WARNING", log_dir=tmp_path)
    assert logger1 is logger2
    assert len(logger2.handlers) == n_handlers  # no duplicate handlers added
    assert logger2.level == logging.WARNING


def test_run_id_propagates_into_log_records(tmp_path, caplog):
    ls.configure_logging(level="DEBUG", log_dir=tmp_path)
    logger = ls.get_logger("test")
    token = ls.set_run_id("abc123")
    try:
        assert ls.current_run_id() == "abc123"
        with caplog.at_level(logging.DEBUG, logger=logger.name):
            logger.debug("hello")
    finally:
        ls.reset_run_id(token)
    assert ls.current_run_id() is None


def test_json_formatter_produces_valid_json_with_run_id():
    formatter = ls.JsonFormatter()
    record = logging.LogRecord("lwmom.test", logging.INFO, __file__, 1, "hello %s", ("world",), None)
    record.run_id = "xyz"
    payload = json.loads(formatter.format(record))
    assert payload["message"] == "hello world"
    assert payload["run_id"] == "xyz"
    assert payload["level"] == "INFO"


def test_new_run_id_is_unique_hex():
    a, b = ls.new_run_id(), ls.new_run_id()
    assert a != b
    assert len(a) == 32
    int(a, 16)  # valid hex
