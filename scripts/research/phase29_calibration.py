"""SUPERSEDED 2026-11-02 (Phase 30): the 7-port HFSS run showed the raw forward model matches HFSS gain; do NOT use this margin in any reporting. Kept for history.

Phase 29 model-to-HFSS margin (two-feed board, HFSS_CHECKED, 10.54 GHz): HFSS P1 17.5 dBi / -8.9 dB, P2 16.6 dBi / -8.1 dB against the model. Applied as a margin to every Phase 29 multi-feed prediction."""

GAIN_OFFSET_DB = -2.5  # calibrated gain = model gain + GAIN_OFFSET_DB
SLL_OFFSET_DB = 3.5  # calibrated SLL (dB rel. peak) = model SLL + SLL_OFFSET_DB (worse)
LABEL = "HFSS_CHECKED_TWO_FEED_MARGIN"


def calibrate(gain_dbi, sll_db):
    return gain_dbi + GAIN_OFFSET_DB, sll_db + SLL_OFFSET_DB
