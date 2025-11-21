# SCHON/3_formula_assignment/formula_filters.py
import pandas as pd
from itertools import product
from typing import Dict, Tuple, Iterable

from schon.formula_assignment.formula_presets import get_preset, FormulaPreset

# -------- core filter function --------

def build_formula_filter(preset: FormulaPreset):
    """
    Given a FormulaPreset, return a function
    formula_filter(c, h, o, n, s, na) -> bool
    that applies all your rules (H/C, O/C, N/C, S/C, DBE, etc.)
    for [M-H]- and [M+Na-2H]-.
    """
    ratio = preset.ratio

    def formula_filter(c: int, h: int, o: int, n: int, s: int, na: int) -> bool:
        if c <= 0:
            return False

        # neutral H (for [M-H]- / [M+Na-2H]-)
        h_neutral = h + 1
        if h_neutral <= 0:
            return False

        h_c = h_neutral / c
        o_c = o / c
        n_c = n / c
        s_c = s / c

        # ratio constraints
        if not (ratio.h_c_min <= h_c <= ratio.h_c_max):
            return False
        if not (ratio.o_c_min <= o_c <= ratio.o_c_max):
            return False
        if not (ratio.n_c_min <= n_c <= ratio.n_c_max):
            return False
        if not (ratio.s_c_min <= s_c <= ratio.s_c_max):
            return False

        # DBE (same as you had)
        DBE = 1 + (2 * c + 2 + n - h) / 2  # na does not change the formula here
        if DBE < 0 or DBE % 1 != 0:
            return False

        return True

    return formula_filter


# -------- optional formula generation (precomputation) --------

ATOMIC_MASSES = {
    "12C": 12.000000000,
    "13C": 13.0033548378,
    "H":   1.007825032,
    "N":   14.003074005,
    "O":   15.994914620,
    "S":   31.972071174,
    "Na":  22.98977,
}

ELECTRON_MASS = 0.000548579909


def generate_formulas(sample_type: str) -> pd.DataFrame:
    """
    Optional helper if you ever *do* want pre-generated formulas.
    Uses the preset associated with `sample_type`.
    """
    preset = get_preset(sample_type)
    chnos = preset.chnos_limits
    filter_fn = build_formula_filter(preset)

    formulas = []

    for c in range(*chnos["C"]):
        for h in range(*chnos["H"]):
            for n in range(*chnos["N"]):
                for o in range(*chnos["O"]):
                    for s in range(*chnos["S"]):
                        for c13 in range(*chnos.get("*C", (0, 1))):
                            for na in range(*chnos.get("Na", (0, 1))):
                                if c == 0 or c13 > c:
                                    continue

                                if not filter_fn(c, h, o, n, s, na):
                                    continue

                                neutral_mass = (
                                    c * ATOMIC_MASSES["12C"]
                                    + c13 * ATOMIC_MASSES["13C"]
                                    + h * ATOMIC_MASSES["H"]
                                    + n * ATOMIC_MASSES["N"]
                                    + o * ATOMIC_MASSES["O"]
                                    + s * ATOMIC_MASSES["S"]
                                    + na * ATOMIC_MASSES["Na"]
                                )

                                if na == 0:
                                    ion_mass = neutral_mass - ATOMIC_MASSES["H"] + ELECTRON_MASS
                                else:
                                    ion_mass = neutral_mass - 2 * ATOMIC_MASSES["H"] + ELECTRON_MASS

                                h_deprot = h - 1
                                formula_str = (
                                    f"C{c}"
                                    + (f"*C{c13}" if c13 > 0 else "")
                                    + (f"H{h_deprot}" if h_deprot > 0 else "")
                                    + (f"N{n}" if n > 0 else "")
                                    + (f"O{o}" if o > 0 else "")
                                    + (f"S{s}" if s > 0 else "")
                                    + (f"Na{na}" if na > 0 else "")
                                )

                                formulas.append(
                                    (
                                        ion_mass, formula_str,
                                        c, h, n, o, s, c13, na
                                    )
                                )

    return pd.DataFrame(
        formulas,
        columns=["Mass", "Formula", "C", "H", "N", "O", "S", "*C", "Na"],
    )