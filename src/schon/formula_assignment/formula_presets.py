"""Shared presets and constants for SCHON formula logic.

This module is the single place where we define physical constants and
common parameter presets used across:

* CoreMS / peak-processing glue code
* Calibration pipeline
* Formula assignment logic

By importing from here we avoid duplicating values (like element masses
or default ppm tolerance) in multiple files.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Tuple

NumberRange = Tuple[int, int]

@dataclass
class RatioLimits:
    h_c_min: float
    h_c_max: float
    o_c_min: float = 0.0
    o_c_max: float = 3.0
    n_c_min: float = 0.0
    n_c_max: float = 1.3
    s_c_min: float = 0.0
    s_c_max: float = 0.8

@dataclass
class FormulaPreset:
    name: str
    ratio: RatioLimits
    chnos_limits: Dict[str, NumberRange]


# ----------------- mass constants -----------------

# Monoisotopic exact masses in Dalton
ELEMENT_MASS: Dict[str, float] = {
    "H": 1.00782503223,
    "C": 12.00000000000,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "S": 31.9720711744,
    "P": 30.97376199842,
}

# Mass of a proton (approx) for [M-H]- / [M+H]+ conversions
PROTON_MASS: float = 1.007276466812

# CH2 exact mass used for homologous spacing and Kendrick transforms
CH2_MASS: float = ELEMENT_MASS["C"] + 2 * ELEMENT_MASS["H"]


# ----------------- general presets -----------------

# Default ppm tolerance for matching theoretical and observed masses.
# Other modules (e.g. main.py, calibration_pipeline.py) should import
# this instead of defining their own value.
PPM_TOLERANCE_DEFAULT: float = 0.5


# Optional: shared element limits that can be reused by both
# calibration and formula assignment (GUI can override).
DEFAULT_ELEMENT_LIMITS: Dict[str, Tuple[int, int]] = {
    "C": (1, 100),
    "H": (1, 200),
    "O": (0, 40),
    "N": (0, 4),
    "S": (0, 2),
    "P": (0, 2),
}
# ---- PRESETS ----

CRUDE_OIL = FormulaPreset(
    name="crude_oil",
    ratio=RatioLimits(
        h_c_min=0.2,
        h_c_max=2.2,
        o_c_min=0.0, o_c_max=0.6,
        n_c_min=0.0, n_c_max=0.2,
        s_c_min=0.0, s_c_max=0.2,
    ),
    chnos_limits={
        'C':  (0, 100),
        'H':  (0, 200),
        'N':  (0, 4),
        'O':  (0, 8),
        'S':  (0, 4),
        '*C': (0, 2),
        'Na': (0, 1),
        '34S': (0, 1),
    },
)

SEDIMENTARY_ROCK = FormulaPreset(
    name="sedimentary_rock",
    ratio=RatioLimits(
        h_c_min=0.2,
        h_c_max=2.2,
        o_c_min=0.0, o_c_max=0.6,
        n_c_min=0.0, n_c_max=0.3,
        s_c_min=0.0, s_c_max=0.3,
    ),
    chnos_limits={
        'C': (0, 100),
        'H': (0, 200),
        'N': (0, 4),
        'O': (0, 6),
        'S': (0, 4),
        '*C': (0, 2),
        'Na': (0, 1),
        '34S': (0, 1),
        '35Cl': (0, 1),
        '37Cl': (0, 1),
        '54Fe': (0, 1),
        '56Fe': (0, 1),
        '57Fe': (0, 1),
        '58Fe': (0, 1),
    },
)

COAL = FormulaPreset(
    name="coal",
    ratio=RatioLimits(
        h_c_min=0.2,
        h_c_max=2.2,
        o_c_min=0.0, o_c_max=0.6,
        n_c_min=0.0, n_c_max=0.3,
        s_c_min=0.0, s_c_max=0.3,
    ),
    chnos_limits={
        'C':  (0, 100),
        'H':  (0, 200),
        'N':  (0, 2),
        'O':  (0, 12),
        'S':  (0, 2),
        '*C': (0, 2),
        'Na': (0, 1),
    },
)

NATURAL_WATER = FormulaPreset(
    name="natural_water",
    ratio=RatioLimits(
        h_c_min=0.2,
        h_c_max=2.2,
        o_c_min=0.1, o_c_max=1.0,
        n_c_min=0.0, n_c_max=0.1,
        s_c_min=0.0, s_c_max=0.1,
    ),
    chnos_limits={
        'C':  (0, 80),
        'H':  (0, 160),
        'N':  (0, 3),
        'O':  (0, 35),
        'S':  (0, 2),
        '*C': (0, 2),
        'Na': (0, 1),
        'P':  (0, 1),
    },
)

GENERIC_ESI_NEG = FormulaPreset(
    name="generic_esi_neg",
    ratio=RatioLimits(
        h_c_min=0.4,
        h_c_max=2.8,
        # keep broad default OC/NC/SC windows as in your Rule 5
    ),
    chnos_limits={
        'C':  (7, 120),
        'H':  (0, 220),
        'N':  (0, 3),
        'O':  (0, 8),
        'S':  (0, 2),
        '*C': (0, 2),
        'Na': (0, 2),
    },
)

PRESETS = {
    "crude_oil": CRUDE_OIL,
    "sedimentary_rock": SEDIMENTARY_ROCK,
    "coal": COAL,
    "natural_water": NATURAL_WATER,
    "generic_esi_neg": GENERIC_ESI_NEG,
}


def get_preset(sample_type) -> FormulaPreset:
    """Return a FormulaPreset for the given sample_type.

    Parameters
    ----------
    sample_type : str | int | None
        - If str, accepts names like "CRUDE_OIL" or "crude_oil".
        - If int, uses a legacy mapping (0..N) to named presets.
        - If None or unknown, falls back to GENERIC_ESI_NEG.
    """
    # Handle None explicitly
    if sample_type is None:
        return GENERIC_ESI_NEG

    # If a string is provided, normalize and look up directly
    if isinstance(sample_type, str):
        key = sample_type.lower()
        return PRESETS.get(key, GENERIC_ESI_NEG)

    # Backwards compatibility: integer codes
    # You can adjust this mapping if you used different numeric codes before.
    if isinstance(sample_type, int):
        legacy_map = {
            0: "generic_esi_neg",
            1: "crude_oil",
            2: "sedimentary_rock",
            3: "coal",
            4: "natural_water",
        }
        key = legacy_map.get(sample_type, "generic_esi_neg")
        return PRESETS.get(key, GENERIC_ESI_NEG)

    # Fallback for any other unexpected type
    return GENERIC_ESI_NEG