# SPDX-License-Identifier: EUPL-1.2
"""
Calibration configuration for the SCHON FT-ICR-MS pipeline.

This file is meant to be:
- human-readable,
- easy to tweak directly,
- and later overridable from a GUI.

All defaults live here.
"""

from __future__ import annotations

from dataclasses import dataclass

# Nominal CH2 mass used for homologous-series detection (Da)
DEFAULT_CH2_MASS: float = 14.01565


@dataclass
class CalibrationConfig:
    """
    Configuration options for internal calibration.

    Attributes
    ----------
    sn_threshold :
        Minimum signal-to-noise ratio for peaks to be considered for
        formula assignment / calibration.

    poly_order :
        Polynomial order for m/z calibration: mz_true = f(mz_observed).

    only_CHO :
        If True, use only CHO-type formulas as calibrants
        (no N, S, Na, 13C).

    ppm_tolerance_link :
        Maximum ppm distance when linking calibrant peaks back to the
        full peak list for reporting m/z error.

    min_series_length :
        Minimum length of a homologous series (e.g. CH2) that a peak
        must belong to in order to be accepted as a calibrant.

    series_step_mass :
        Nominal mass difference for the homologous series (default CH2).

    series_ppm_tolerance :
        PPM tolerance when deciding if two peaks differ by an integer
        number of series_step_mass units.
    """

    sn_threshold: float = 20.0
    poly_order: int = 2
    only_CHO: bool = True
    ppm_tolerance_link: float = 0.5

    min_series_length: int = 2
    series_step_mass: float = DEFAULT_CH2_MASS
    series_ppm_tolerance: float = 5.0


# A default instance that you can import in CLI scripts or tests
DEFAULT_CALIBRATION_CONFIG = CalibrationConfig()