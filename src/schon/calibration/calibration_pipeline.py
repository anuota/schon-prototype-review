# SPDX-License-Identifier: EUPL-1.2
"""Calibration pipeline for FT-ICR-MS data using SCHON formula logic.

This module is intended to run *after* the CoreMS pipeline has exported
peak-picked spectra into a pandas DataFrame (typically as a pickle file).

High-level idea
---------------
1. Load a CoreMS mass spectrum DataFrame (one spectrum per file).
2. Filter peaks by a signal-to-noise threshold.
3. For each candidate peak generate molecular formulas *on the fly*
   using configurable element limits (no global pre-generated library).
4. Keep only uniquely assigned peaks which are part of at least one
   homologous series (by ~CH2 mass spacing).
5. Use those peaks as internal calibrants to fit a polynomial
   m/z calibration (linear or quadratic).
6. Apply the calibration to all peaks and return / save the result.

The design keeps the element list flexible so the GUI can configure it.
We still reuse information from ``formula_filters`` (element limits and
filter ranges) and stay compatible with column naming in ``main.py``.

This file deliberately does *not* depend on the global formula database
(``formulas_data`` / ``formula_tree``) from main.py so that calibration
remains fast and independent of any pre-generated libraries.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import sys
import numpy as np
import pandas as pd

from schon.corems_pipeline.generate_peaks_from_fid import OUTPUT_PEAKS_DIR

# We only import *data* from formula_filters, not the heavy
# generate_formulas() helper, to keep calibration independent
# of the global pre-generated library.
try:
    from formula_filters import CHNOS_LIMITS, CHO_FILTERS
except ImportError:  # pragma: no cover - makes module more robust
    CHNOS_LIMITS = {}
    CHO_FILTERS = {}

# Try to reuse default ppm tolerance from main.py if available
try:  # pragma: no cover
    from main import PPM_TOLERANCE_DEFAULT as DEFAULT_PPM_TOLERANCE
except Exception:  # pragma: no cover
    DEFAULT_PPM_TOLERANCE = 0.5  # ppm


# --------------------- basic element mass table ---------------------

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
PROTON_MASS = 1.007276466812

# CH2 exact mass used for homologous spacing
CH2_MASS = ELEMENT_MASS["C"] + 2 * ELEMENT_MASS["H"]


# -------------------------- data structures -------------------------

@dataclass
class CalibrationConfig:
    """Configuration for the calibration pipeline.

    Parameters
    ----------
    snr_min:
        Minimal signal-to-noise ratio for a peak to be considered.
    ppm_tolerance:
        Mass window for matching formulas to peaks.
    poly_order:
        Order of polynomial calibration in m/z (1=linear, 2=quadratic).
    min_series_length:
        Minimal number of peaks in a homologous series to be kept.
    prefer_CHO_only:
        If True, restrict calibrant formulas to C/H/O only.
    use_negative_mode:
        If True, treat peaks as [M-H]- (negative ESI).
    """

    snr_min: float = 10.0
    ppm_tolerance: float = DEFAULT_PPM_TOLERANCE
    poly_order: int = 1
    min_series_length: int = 2
    prefer_CHO_only: bool = True
    use_negative_mode: bool = True


@dataclass
class CalibrantPeak:
    """Container for a calibrant peak used in the fit."""

    index: int
    mz_observed: float
    mz_theoretical: float
    formula: str
    snr: float


# -------------------------- helpers ---------------------------------


def get_default_element_limits(prefer_CHO_only: bool = True) -> Dict[str, Tuple[int, int]]:
    """Return element limits, defaulting to CHO subset of CHNOS_LIMITS.

    The GUI can later override these limits and pass a custom dictionary
    into the high-level pipeline.
    """

    if not CHNOS_LIMITS:
        # Reasonable fall-back limits if the import failed
        base_limits = {
            "C": (1, 100),
            "H": (1, 200),
            "O": (0, 40),
            "N": (0, 4),
            "S": (0, 2),
            "P": (0, 2),
        }
    else:
        base_limits = {
            el: (cfg["min"], cfg["max"]) for el, cfg in CHNOS_LIMITS.items()
        }

    if prefer_CHO_only:
        return {el: lim for el, lim in base_limits.items() if el in {"C", "H", "O"}}
    return base_limits


def neutral_mass_from_mz(mz: float, charge: int) -> float:
    """Convert measured m/z into a neutral mass.

    For negative mode we assume [M-H]- so that:

        m/z = (M - H) / |z|
        M   = m/z * |z| + H

    For positive mode [M+H]+:

        m/z = (M + H) / |z|
        M   = m/z * |z| - H
    """

    z = abs(charge)
    if z == 0:
        raise ValueError("Charge cannot be zero when converting m/z to neutral mass")

    if charge < 0:
        return mz * z + PROTON_MASS
    else:
        return mz * z - PROTON_MASS


def ppm_error(theoretical: float, observed: float) -> float:
    return 1e6 * (observed - theoretical) / theoretical


# ------------------ dynamic formula generation ----------------------


def generate_CHO_formulas_for_mass(
    mz: float,
    charge: int,
    elem_limits: Dict[str, Tuple[int, int]],
    ppm_tolerance: float,
) -> List[Tuple[str, float]]:
    """Generate CHO formulas for a given m/z on the fly.

    This is intentionally simple and is only meant for selecting a
    *small* set of robust calibrant peaks. It searches over C and O
    within the given limits and derives H from the mass balance.

    Parameters
    ----------
    mz:
        Observed m/z value.
    charge:
        Ion charge, usually -1 for negative ESI.
    elem_limits:
        Dict mapping element symbol to (min, max) counts.
    ppm_tolerance:
        Allowed mass error between neutral and formula mass.

    Returns
    -------
    List of (formula_str, theoretical_mz).
    """

    neutral_mass = neutral_mass_from_mz(mz, charge)

    c_min, c_max = elem_limits.get("C", (1, 100))
    h_min, h_max = elem_limits.get("H", (1, 200))
    o_min, o_max = elem_limits.get("O", (0, 40))

    results: List[Tuple[str, float]] = []

    for c in range(c_min, c_max + 1):
        remaining_after_c = neutral_mass - c * ELEMENT_MASS["C"]
        if remaining_after_c <= 0:
            break

        # Max O limited by both element ranges and remaining mass
        o_max_mass_limited = int(remaining_after_c // ELEMENT_MASS["O"])
        o_upper = min(o_max, o_max_mass_limited)

        for o in range(o_min, o_upper + 1):
            mass_CO = c * ELEMENT_MASS["C"] + o * ELEMENT_MASS["O"]
            remaining = neutral_mass - mass_CO
            if remaining <= 0:
                break

            h = int(round(remaining / ELEMENT_MASS["H"]))
            if h < h_min or h > h_max:
                continue

            mass_formula = (
                c * ELEMENT_MASS["C"]
                + h * ELEMENT_MASS["H"]
                + o * ELEMENT_MASS["O"]
            )

            err_ppm = abs(ppm_error(mass_formula, neutral_mass))
            if err_ppm <= ppm_tolerance:
                # convert back to expected m/z for this charge state
                if charge < 0:
                    theo_mz = (mass_formula - PROTON_MASS) / abs(charge)
                else:
                    theo_mz = (mass_formula + PROTON_MASS) / abs(charge)

                formula = f"C{c}H{h}O{o}"
                results.append((formula, theo_mz))

    return results


# ------------------ homologous series detection ---------------------


def detect_homologue_series(
    mz_values: np.ndarray,
    min_series_length: int = 2,
    step_mass: float = CH2_MASS,
    ppm_tolerance: float = 5.0,
) -> np.ndarray:
    """Return a boolean mask for peaks that belong to a CH2-like series.

    The algorithm is intentionally simple: we build an undirected graph
    where two peaks are connected if their m/z difference corresponds to
    an integer number of CH2 units within a ppm tolerance. Peaks that
    belong to a connected component of size >= ``min_series_length`` are
    marked as True.
    """

    n = len(mz_values)
    if n == 0:
        return np.zeros(0, dtype=bool)

    order = np.argsort(mz_values)
    mz_sorted = mz_values[order]

    neighbors: List[List[int]] = [[] for _ in range(n)]

    # brute-force O(n^2) is acceptable here because the number of
    # calibrant candidates is expected to be modest
    for i in range(n):
        for j in range(i + 1, n):
            diff = mz_sorted[j] - mz_sorted[i]
            if diff <= 0:
                continue

            steps = round(diff / step_mass)
            if steps <= 0:
                continue

            theo_diff = steps * step_mass
            err_ppm = abs(ppm_error(theo_diff, diff))
            if err_ppm <= ppm_tolerance:
                neighbors[i].append(j)
                neighbors[j].append(i)

    # BFS / DFS to get connected components
    visited = np.zeros(n, dtype=bool)
    in_series = np.zeros(n, dtype=bool)

    for i in range(n):
        if visited[i]:
            continue
        stack = [i]
        component: List[int] = []
        while stack:
            k = stack.pop()
            if visited[k]:
                continue
            visited[k] = True
            component.append(k)
            for nb in neighbors[k]:
                if not visited[nb]:
                    stack.append(nb)

        if len(component) >= min_series_length:
            for k in component:
                in_series[k] = True

    # map back to original order
    mask = np.zeros(n, dtype=bool)
    mask[order] = in_series
    return mask


# ---------------------- core pipeline steps -------------------------


def assign_formulas_for_calibration(
    df: pd.DataFrame,
    config: CalibrationConfig,
    elem_limits: Optional[Dict[str, Tuple[int, int]]] = None,
    mz_col: str = "m/z",
    snr_col: str = "S/N",
    charge: int = -1,
) -> pd.DataFrame:
    """Assign simple CHO formulas to peaks for calibration.

    This function:
      * filters peaks by S/N,
      * generates CHO formulas on the fly,
      * keeps only peaks with exactly one candidate formula,
      * adds columns ``candidate_formulas`` and ``Molecular Formula``
        (for unique assignments only), plus ``Calculated m/z``.

    It returns a *copy* of the DataFrame to avoid side effects.
    """

    if elem_limits is None:
        elem_limits = get_default_element_limits(prefer_CHO_only=config.prefer_CHO_only)

    df = df.copy()

    if snr_col not in df.columns:
        raise KeyError(f"S/N column '{snr_col}' not found in DataFrame")
    if mz_col not in df.columns:
        raise KeyError(f"m/z column '{mz_col}' not found in DataFrame")

    candidate_formulas: List[Optional[List[Tuple[str, float]]]] = []
    chosen_formula: List[Optional[str]] = []
    chosen_theo_mz: List[Optional[float]] = []

    for _, row in df.iterrows():
        snr = float(row[snr_col]) if not pd.isna(row[snr_col]) else 0.0
        mz = float(row[mz_col])

        if snr < config.snr_min:
            candidate_formulas.append(None)
            chosen_formula.append(None)
            chosen_theo_mz.append(None)
            continue

        cands = generate_CHO_formulas_for_mass(
            mz,
            charge=charge,
            elem_limits=elem_limits,
            ppm_tolerance=config.ppm_tolerance,
        )

        candidate_formulas.append(cands or None)

        if len(cands) == 1:
            f_str, theo_mz = cands[0]
            chosen_formula.append(f_str)
            chosen_theo_mz.append(theo_mz)
        else:
            chosen_formula.append(None)
            chosen_theo_mz.append(None)

    df["candidate_formulas"] = candidate_formulas
    df["Molecular Formula"] = chosen_formula
    df["Calculated m/z"] = chosen_theo_mz

    # m/z error will be filled after calibration
    return df


def select_calibrant_peaks(
    df: pd.DataFrame,
    config: CalibrationConfig,
    mz_col: str = "m/z",
    snr_col: str = "S/N",
) -> List[CalibrantPeak]:
    """Select robust calibrant peaks from a DataFrame.

    Criteria:
      * unique formula assignment (Molecular Formula not null),
      * high S/N (already enforced in formula assignment step),
      * peak belongs to a homologous series.
    """

    mask_unique = df["Molecular Formula"].notna() & df["Calculated m/z"].notna()
    df_unique = df[mask_unique].copy()

    if df_unique.empty:
        return []

    mz_vals = df_unique[mz_col].to_numpy(dtype=float)
    in_series_mask = detect_homologue_series(
        mz_vals,
        min_series_length=config.min_series_length,
    )

    df_series = df_unique[in_series_mask].copy()
    if df_series.empty:
        return []

    calibrants: List[CalibrantPeak] = []
    for idx, row in df_series.iterrows():
        mz_obs = float(row[mz_col])
        mz_theo = float(row["Calculated m/z"])
        snr = float(row[snr_col]) if not pd.isna(row[snr_col]) else 0.0
        formula = str(row["Molecular Formula"])
        calibrants.append(
            CalibrantPeak(
                index=idx,
                mz_observed=mz_obs,
                mz_theoretical=mz_theo,
                formula=formula,
                snr=snr,
            )
        )

    return calibrants


def fit_mz_calibration(
    calibrants: Sequence[CalibrantPeak],
    order: int = 1,
) -> np.ndarray:
    """Fit polynomial coefficients mapping observed m/z -> theoretical m/z."""

    if not calibrants:
        raise ValueError("No calibrant peaks provided for calibration fit")

    x = np.array([c.mz_observed for c in calibrants], dtype=float)
    y = np.array([c.mz_theoretical for c in calibrants], dtype=float)

    if order not in (1, 2):
        raise ValueError("Only linear (1) and quadratic (2) calibration are supported")

    coeffs = np.polyfit(x, y, order)
    # numpy returns highest power first; keep as-is
    return coeffs


def apply_mz_calibration(
    df: pd.DataFrame,
    coeffs: np.ndarray,
    mz_col: str = "m/z",
    out_col: str = "Calibrated m/z (SCHON)",
) -> pd.DataFrame:
    """Apply m/z calibration polynomial to a DataFrame.

    The returned DataFrame is a *copy* with an additional column
    containing the calibrated m/z values.
    """

    df = df.copy()
    if mz_col not in df.columns:
        raise KeyError(f"m/z column '{mz_col}' not found in DataFrame")

    x = df[mz_col].to_numpy(dtype=float)
    y = np.polyval(coeffs, x)

    df[out_col] = y

    if "Calculated m/z" in df.columns and "m/z Error (ppm)" in df.columns:
        # recompute error based on calibrated m/z for assigned peaks
        mask = df["Calculated m/z"].notna()
        theo = df.loc[mask, "Calculated m/z"].to_numpy(dtype=float)
        obs = df.loc[mask, out_col].to_numpy(dtype=float)
        df.loc[mask, "m/z Error (ppm)"] = ppm_error(theo, obs)

    return df


def run_calibration_pipeline(
    ms_df: pd.DataFrame,
    config: Optional[CalibrationConfig] = None,
    elem_limits: Optional[Dict[str, Tuple[int, int]]] = None,
    mz_col: str = "m/z",
    snr_col: str = "S/N",
    charge: int = -1,
) -> Tuple[pd.DataFrame, np.ndarray, List[CalibrantPeak]]:
    """High-level convenience function to run the full calibration.

    Parameters
    ----------
    ms_df:
        DataFrame with at least ``m/z`` and ``S/N`` columns produced
        by the CoreMS pipeline.
    config:
        CalibrationConfig object. If None, defaults are used.
    elem_limits:
        Optional custom element limits from the GUI.
    mz_col, snr_col:
        Column names for m/z and signal-to-noise.
    charge:
        Charge state of the ions (usually -1 for ESI negative mode).

    Returns
    -------
    calibrated_df, coeffs, calibrants
        Calibrated DataFrame, calibration polynomial coefficients and
        list of calibrant peaks used in the fit.
    """

    if config is None:
        config = CalibrationConfig()

    # 1. Assign formulas for potential calibrant peaks
    df_with_formulas = assign_formulas_for_calibration(
        ms_df,
        config=config,
        elem_limits=elem_limits,
        mz_col=mz_col,
        snr_col=snr_col,
        charge=charge,
    )

    # 2. Select robust calibrants (unique + in a homologous series)
    calibrants = select_calibrant_peaks(
        df_with_formulas,
        config=config,
        mz_col=mz_col,
        snr_col=snr_col,
    )

    if not calibrants:
        raise RuntimeError("Calibration failed: no suitable calibrant peaks found")

    # 3. Fit calibration polynomial
    coeffs = fit_mz_calibration(calibrants, order=config.poly_order)

    # 4. Apply calibration to all peaks
    calibrated_df = apply_mz_calibration(df_with_formulas, coeffs, mz_col=mz_col)

    return calibrated_df, coeffs, calibrants


# ------------------------------ CLI ---------------------------------


def _default_input_path() -> Path:
    """Default path for local testing / docker development.

    This matches the file you mentioned in the conversation. In the
    docker image it will typically live under ``/app`` instead of the
    absolute macOS path.
    """

    # In the SCHON package layout this file is in SCHON/calibration/.
    # We go two levels up and then into output_peaks/.
    root = Path(__file__).resolve().parents[2]
    return root / "output_peaks" / "ESI_neg_G017736-0_100_200sc_000001.d_full_ms.pkl"


def main(input_path: Optional[str] = None, output_path: Optional[str] = None) -> None:
    """Simple CLI entry point for manual testing.

    Examples
    --------
    Within the docker container you can run:

    ``python -m SCHON.calibration.calibration_pipeline``

    and it will try to load the default pickle file, run calibration,
    and save a calibrated pickle next to it with ``_calibrated`` suffix.
    """

    if input_path is None:
        in_path = OUTPUT_PEAKS_DIR #_default_input_path()
    else:
        in_path = Path(input_path)

    if not in_path.is_file():
        raise FileNotFoundError(f"Input mass spectrum DataFrame not found: {in_path}")

    ms_df = pd.read_pickle(in_path)

    calibrated_df, coeffs, calibrants = run_calibration_pipeline(ms_df)

    # Default output path: add _calibrated suffix
    if output_path is None:
        out_path = in_path.with_name(in_path.stem + "_calibrated.pkl")
    else:
        out_path = Path(output_path)

    calibrated_df.to_pickle(out_path)

    # Print a tiny summary so the user sees something useful in the log
    print(f"Calibration complete. Saved calibrated DataFrame to: {out_path}")
    print("Polynomial coefficients (highest power first):", coeffs)
    print(f"Number of calibrant peaks used: {len(calibrants)}")


if __name__ == "__main__":  # pragma: no cover
    main()