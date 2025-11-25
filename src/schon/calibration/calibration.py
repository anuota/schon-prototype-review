"""
Calibration pipeline for FT-ICR-MS peaks using internally assigned formulas.

Steps
-----
1. Load peaks from `<sample>_peaks.csv`.
   Expected columns:
       signalNumber, observedExactMass_ion, observedAbundance_ion, signalNoise_ratio

2. Filter peaks by S/N >= threshold (default 20).

3. Run your formula assignment on the S/N-filtered peaks to produce a formulas CSV.

4. From that formulas CSV, select robust internal calibrants:
   - isotopolog == 'Q'  (monoisotopic)
   - no "Alternative Formula"
   - optionally CHO-only (N == S == Na == 0, *C == 0)
   - must belong to a homologous series (e.g. CH2) of length >= 2

5. Fit a polynomial calibration.

6. Apply calibration and save calibrated results.

This module deliberately keeps things:
- configurable (thresholds, file paths, polynomial order),
- independent of CoreMS internals (uses numpy.polyfit),
- and easy to extend later with a more sophisticated calibration (e.g. Kozhinov–Savory).
"""

from __future__ import annotations

from dataclasses import dataclass

# Import CalibrationConfig from the config module
from schon.calibration.calibration_config import CalibrationConfig
from pathlib import Path
from typing import Optional, Sequence, Tuple, List, Union

import numpy as np
import pandas as pd

# Directory for debug / test calibration outputs
TEST_CALIB_DIR: Path = Path("/app/results/test_calibration")


# ---------------------------------------------------------------------------
# Small numerical helpers for homologue-series detection
# ---------------------------------------------------------------------------

# Nominal CH2 mass used for homologous-series detection (Da)
CH2_MASS: float = 14.01565


def ppm_error(theoretical: float, observed: float) -> float:
    """Return mass error in ppm: (observed - theoretical) / theoretical * 1e6."""
    return (observed - theoretical) / theoretical * 1e6


from typing import List  # ensure List is available for type hints



def detect_homologue_series(
    mz_values: np.ndarray,
    min_series_length: int = 2,
    step_mass: float = CH2_MASS,
    ppm_tolerance: float = 5.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (mask, series_lengths) for peaks that belong to a CH2-like series.

    mask[i] is True if peak i belongs to a homologous series of length
    >= min_series_length. series_lengths[i] gives the size of the
    connected component that peak i belongs to.

    The algorithm is intentionally simple: we build an undirected graph
    where two peaks are connected if their m/z difference corresponds to
    an integer number of CH2 units within a ppm tolerance.
    """

    n = len(mz_values)
    if n == 0:
        return np.zeros(0, dtype=bool), np.zeros(0, dtype=int)

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
    component_sizes = np.zeros(n, dtype=int)

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

        size = len(component)
        if size >= min_series_length:
            for k in component:
                in_series[k] = True
                component_sizes[k] = size

    # map back to original order
    mask = np.zeros(n, dtype=bool)
    sizes = np.zeros(n, dtype=int)
    mask[order] = in_series
    sizes[order] = component_sizes
    return mask, sizes


# ---------------------------------------------------------------------------
# Configuration dataclasses
# ---------------------------------------------------------------------------

@dataclass
class PeakColumns:
    """Column names in the raw peaks CSV produced by CoreMS."""
    mz: str = "observedExactMass_ion"
    intensity: str = "observedAbundance_ion"
    sn: str = "signalNoise_ratio"


@dataclass
class FormulaColumns:
    """Column names in the formulas CSV produced by your formula-assignment step."""
    mz: str = "m/z"
    intensity: str = "Intensity"
    isotopolog: str = "isotopolog"
    formula: str = "Formula"
    calc_mass: str = "Calculated Mass"
    c: str = "C"
    h: str = "H"
    n: str = "N"
    o: str = "O"
    s: str = "S"
    c13: str = "*C"
    na: str = "Na"
    mass_error_ppm: str = "Mass Error (ppm)"
    alt_formula: str = "Alternative Formula"
    alt_mass_error_ppm: str = "Alternative Mass Error (ppm)"




# ---------------------------------------------------------------------------
# Step 1: Load and filter peaks by S/N
# ---------------------------------------------------------------------------

def load_peaks(peaks_path: Path, cols: PeakColumns = PeakColumns()) -> pd.DataFrame:
    """Load CoreMS peaks CSV."""
    df = pd.read_csv(peaks_path)
    # sanity check
    for required in (cols.mz, cols.intensity, cols.sn):
        if required not in df.columns:
            raise ValueError(
                f"Required column '{required}' not found in {peaks_path}.\n"
                f"Available columns: {list(df.columns)}"
            )
    return df


def filter_peaks_by_snr(
    df: pd.DataFrame,
    cfg: CalibrationConfig,
    cols: PeakColumns = PeakColumns(),
) -> pd.DataFrame:
    """Keep only peaks with S/N >= threshold."""
    mask = df[cols.sn] >= cfg.sn_threshold
    return df.loc[mask].copy()


# ---------------------------------------------------------------------------
# Step 2 (external): formula assignment
# ---------------------------------------------------------------------------

def write_filtered_peaks_for_formula_assignment(
    df_filtered: pd.DataFrame,
    out_path: Path,
    peak_cols: PeakColumns = PeakColumns(),
) -> Path:
    """
    Prepare a minimal CSV for the formula-assignment step.

    We write a file with 'm/z' and 'Intensity' columns so that your
    `schon.formula_assignment.main` script can read it. You can adapt this
    function if your main expects a different format.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)

    df_out = pd.DataFrame(
        {
            "m/z": df_filtered[peak_cols.mz],
            "Intensity": df_filtered[peak_cols.intensity],
        }
    )
    df_out.to_csv(out_path, index=False)
    return out_path





def run_formula_assignment(
    peaks_for_formulas_df: pd.DataFrame,
    sample_type: str | None = None,
    ppm_tolerance: float | None = None,
    n_processes: Optional[int] = None,
) -> pd.DataFrame:
    """Run formula assignment for a (typically filtered) peak list in memory.

    This is a thin wrapper around
    ``schon.formula_assignment.formula_assignment.assign_formulas_df`` and is
    intended to be used inside the calibration flow.

    Parameters
    ----------
    peaks_for_formulas_df : pd.DataFrame
        DataFrame with at least the columns ``m/z`` and ``Intensity``,
        usually produced by ``write_filtered_peaks_for_formula_assignment``-like
        logic (but kept in memory here).
    sample_type : str, optional
        Name of the sample type preset (e.g. "crude_oil", "natural_water").
        Passed through to the underlying ``assign_formulas_df``. If ``None``,
        the module default preset is used.
    ppm_tolerance : float, optional
        Mass-accuracy window for assignment. If ``None``, the module default
        (``DEFAULT_PPM_TOLERANCE``) is used.
    n_processes : int, optional
        Number of worker processes. If ``None``, the module decides.

    Returns
    -------
    pd.DataFrame
        DataFrame of assigned formulas (same structure as the CSV that
        used to be written by the CSV-based interface).
    """
    from schon.formula_assignment import formula_assignment as formula_main

    kwargs: dict = {}
    if sample_type is not None:
        kwargs["sample_type"] = sample_type
    if ppm_tolerance is not None:
        kwargs["ppm_tolerance"] = ppm_tolerance
    if n_processes is not None:
        kwargs["n_processes"] = n_processes

    # We assume formula_main exposes a DataFrame-based API.
    return formula_main.assign_formulas_df(
        peaks_for_formulas_df,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# Step 3: Choose robust calibrant peaks from formulas
# ---------------------------------------------------------------------------

def load_formulas(formulas_path: Path, cols: FormulaColumns = FormulaColumns()) -> pd.DataFrame:
    df = pd.read_csv(formulas_path)
    # We assume your formula assignment created the standard columns you showed.
    return df



def select_calibrant_peaks(
    formulas_df: pd.DataFrame,
    cfg: CalibrationConfig,
    cols: FormulaColumns = FormulaColumns(),
) -> pd.DataFrame:
    """
    Select robust internal calibrant peaks.

    Criteria:
    - isotopolog == 'Q' (monoisotopic peaks)
    - no alternative formula (unique assignment)
    - optionally CHO-only (no N, no S, no Na, no 13C)
    - must belong to a homologous series (e.g. CH2) of length >= 2

    Returns
    -------
    pd.DataFrame
        DataFrame with at least the columns:
        - mz_obs        (observed m/z)
        - mz_true       (theoretical m/z from the formula)
        - formula       (molecular formula string)
        - Intensity     (peak intensity as in formulas table)
        - series_length (size of the homologous series component)
    """
    df = formulas_df.copy()

    mask = df[cols.isotopolog].astype(str).str.upper().eq("Q")

    # unique assignment: no alternative formula
    if cols.alt_formula in df.columns:
        mask &= df[cols.alt_formula].isna() | (df[cols.alt_formula].astype(str).str.strip() == "")

    if cfg.only_CHO:
        for col in (cols.n, cols.s, cols.na, cols.c13):
            if col in df.columns:
                mask &= df[col].fillna(0).astype(float).eq(0)

    # keep observed m/z, theoretical m/z, formula, and intensity
    calibrants = df.loc[
        mask,
        [cols.mz, cols.calc_mass, cols.formula, cols.intensity],
    ].copy()

    calibrants.rename(
        columns={
            cols.mz: "mz_obs",
            cols.calc_mass: "mz_true",
            cols.formula: "formula",
            cols.intensity: "Intensity",
        },
        inplace=True,
    )

    # drop rows with missing values
    calibrants = calibrants.dropna(subset=["mz_obs", "mz_true"])

    if calibrants.empty:
        raise RuntimeError("No calibrant peaks selected – check your filters / formulas.")

    # Additional filter: keep only peaks that are part of a homologous
    # series (e.g. CH2) of length >= 2, using the theoretical masses.
    series_mask, series_sizes = detect_homologue_series(
        calibrants["mz_true"].to_numpy(dtype=float),
        min_series_length=cfg.min_series_length,
        step_mass=cfg.series_step_mass,
        ppm_tolerance=cfg.series_ppm_tolerance,
    )

    calibrants = calibrants.loc[series_mask].copy()
    if calibrants.empty:
        raise RuntimeError(
            "No calibrant peaks remain after enforcing the homologous-series criterion."
        )

    # store series length for each calibrant
    calibrants["series_length"] = series_sizes[series_mask]

    return calibrants.sort_values("mz_obs").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Step 4: Fit polynomial calibration
# ---------------------------------------------------------------------------

def fit_calibration_polynomial(
    calibrants: pd.DataFrame,
    cfg: CalibrationConfig,
) -> np.ndarray:
    """
    Fit a polynomial m/z_true = f(mz_obs).

    Returns: numpy array of coefficients (highest power first),
    suitable for np.poly1d.
    """
    x = calibrants["mz_obs"].to_numpy(dtype=float)
    y = calibrants["mz_true"].to_numpy(dtype=float)

    if x.size <= cfg.poly_order:
        raise RuntimeError(
            f"Not enough calibrant points ({x.size}) to fit a polynomial of "
            f"order {cfg.poly_order}."
        )

    coeffs = np.polyfit(x, y, cfg.poly_order)
    return coeffs


# ---------------------------------------------------------------------------
# Step 5: Apply calibration to full peak list and write results
# ---------------------------------------------------------------------------

def apply_calibration_to_peaks(
    peaks_df: pd.DataFrame,
    coeffs: np.ndarray,
    peak_cols: PeakColumns = PeakColumns(),
) -> pd.DataFrame:
    """
    Add a 'Calibrated m/z' column to the peaks dataframe using the
    polynomial defined by coeffs.
    """
    poly = np.poly1d(coeffs)
    mz_obs = peaks_df[peak_cols.mz].to_numpy(dtype=float)
    peaks_df = peaks_df.copy()
    peaks_df["Calibrated m/z"] = poly(mz_obs)
    return peaks_df


def build_full_ms_like_output(
    calibrated_df: pd.DataFrame,
    calibrants: pd.DataFrame,
    cfg: CalibrationConfig,
    peak_cols: PeakColumns = PeakColumns(),
) -> pd.DataFrame:
    """
    Create an output dataframe with "full_ms" style columns.

    We map:
        observedExactMass_ion -> "m/z"
        observedAbundance_ion -> "Peak Height"
        signalNoise_ratio     -> "S/N"

    For calibrant peaks we also fill:
        "Calculated m/z", "m/z Error (ppm)"
    """
    df = calibrated_df.copy()

    # Basic structure
    out = pd.DataFrame()
    out["Index"] = np.arange(len(df))
    out["m/z"] = df[peak_cols.mz]
    out["Calibrated m/z"] = df["Calibrated m/z"]
    # Preserve the original (uncalibrated) m/z explicitly for downstream logic
    out["m/z_raw"] = out["m/z"]

    # PPM shift introduced by calibration (Calibrated vs original m/z)
    out["Calibration Error (ppm)"] = ppm_error(
        out["m/z"].to_numpy(dtype=float),
        out["Calibrated m/z"].to_numpy(dtype=float),
    )

    # Calculated m/z and error, to be filled only for calibrants
    out["Calculated m/z"] = np.nan

    # map obvious columns
    out["Peak Height"] = df[peak_cols.intensity]
    # If you have Peak Area & Resolving Power from CoreMS, you can add them here;
    # for now we leave them NaN.
    out["Peak Area"] = np.nan
    out["Resolving Power"] = np.nan

    out["S/N"] = df[peak_cols.sn]
    out["Ion Charge"] = -1  # negative mode default; adjust if needed

    # m/z error in ppm (only for calibrants)
    out["m/z Error (ppm)"] = np.nan

    # Add a boolean column indicating whether each peak was used for calibration
    out["Used For Calibration"] = False

    # placeholders for other columns expected in your GUI
    out["m/z Error Score"] = np.nan
    out["Isotopologue Similarity"] = np.nan
    out["Confidence Score"] = np.nan
    out["DBE"] = np.nan
    out["O/C"] = np.nan
    out["H/C"] = np.nan
    out["Heteroatom Class"] = "unassigned"
    out["Ion Type"] = ""
    out["Adduct"] = ""
    out["Is Isotopologue"] = False
    out["Mono Isotopic Index"] = np.nan
    out["Molecular Formula"] = ""

    # ------------------------------------------------------------------
    # Link calibrant list back to the full list and fill calc m/z + error
    # ------------------------------------------------------------------
    if not calibrants.empty:
        # Merge by closest m/z within tolerance (ppm)
        full = out[["Index", "m/z"]].copy()
        full_sorted = full.sort_values("m/z")
        calib_sorted = calibrants.sort_values("mz_obs")

        # use merge_asof with absolute ppm tolerance
        merged = pd.merge_asof(
            calib_sorted.rename(columns={"mz_obs": "mz_calib"}),
            full_sorted.rename(columns={"m/z": "mz"}),
            left_on="mz_calib",
            right_on="mz",
            direction="nearest",
        )

        # filter by ppm tolerance
        delta_ppm = (merged["mz_calib"] - merged["mz"]).abs() / merged["mz"] * 1e6
        ok = delta_ppm <= cfg.ppm_tolerance_link
        merged = merged.loc[ok]

        # fill in Calculated m/z and error for matched indices
        idx = merged["Index"].to_numpy(dtype=int)
        out.loc[idx, "Calculated m/z"] = merged["mz_true"].to_numpy()
        out.loc[idx, "m/z Error (ppm)"] = (
            (out.loc[idx, "Calibrated m/z"] - out.loc[idx, "Calculated m/z"]) /
            out.loc[idx, "Calculated m/z"] * 1e6
        )

        # -------- Mark peaks used for calibration --------
        # Here we mark peaks whose observed m/z matches a calibrant mz_obs.
        calibrant_mz = set(calibrants["mz_obs"].astype(float).tolist())
        # Use the raw observed m/z column ("m/z" in this full_ms table)
        out["Used For Calibration"] = out["m/z"].astype(float).isin(calibrant_mz)
    return out


def save_calibrated_results(
    full_ms_like_df: pd.DataFrame,
    sample_name: str,
    results_dir: Path = Path("/app/results"),
) -> Tuple[Path, Path]:
    """
    Save calibrated results as CSV and pickle.

    Filenames:
        <results_dir>/<sample>_calibrated_peaks.csv
        <results_dir>/<sample>_calibrated_peaks.pkl
    """
    results_dir.mkdir(parents=True, exist_ok=True)
    csv_path = results_dir / f"{sample_name}_calibrated_peaks.csv"
    pkl_path = results_dir / f"{sample_name}_calibrated_peaks.pkl"

    full_ms_like_df.to_csv(csv_path, index=False)
    full_ms_like_df.to_pickle(pkl_path)

    return csv_path, pkl_path


# ---------------------------------------------------------------------------
# End-to-end runner
# ---------------------------------------------------------------------------


def run_calibration(
    peaks: Union[Path, pd.DataFrame],
    cfg: CalibrationConfig = CalibrationConfig(),
    peak_cols: PeakColumns = PeakColumns(),
    formula_cols: FormulaColumns = FormulaColumns(),
    results_dir: Path = Path("/app/results"),
    sample_type: str | None = None,
    ppm_tolerance: float | None = None,
    n_processes: Optional[int] = None,
    sample_name: Optional[str] = None,
) :
    """End-to-end calibration starting from peaks.

    This function can be called either with:
    - a pandas DataFrame of peaks (recommended for pipeline use), or
    - a Path to a CoreMS peaks CSV (for convenience / CLI use).
    """
    # Determine peaks dataframe and sample name
    if isinstance(peaks, Path):
        # Backwards-compatible: load from CSV
        peaks_df = load_peaks(peaks, cols=peak_cols)
        if sample_name is None:
            sample_name = peaks.stem.replace("_peaks", "")
    else:
        # Already a DataFrame
        peaks_df = peaks
        if sample_name is None:
            sample_name = "sample"

    # Safety check
    if sample_name is None:
        sample_name = "sample"

    # 1) Filter by S/N
    filtered_peaks_df = filter_peaks_by_snr(peaks_df, cfg, cols=peak_cols)

    # Ensure test-calibration directory exists
    TEST_CALIB_DIR.mkdir(parents=True, exist_ok=True)

    # Save initial mass list (after S/N filtering) for inspection:
    # m/z, Intensity, S/N
    # initial_mass_list = pd.DataFrame(
    #     {
    #         "m/z": filtered_peaks_df[peak_cols.mz].to_numpy(dtype=float),
    #         "Intensity": filtered_peaks_df[peak_cols.intensity].to_numpy(dtype=float),
    #         "S/N": filtered_peaks_df[peak_cols.sn].to_numpy(dtype=float),
    #     }
    # )
    # initial_mass_list_path = TEST_CALIB_DIR / f"{sample_name}_initial_mass_list.csv"
    # # initial_mass_list.to_csv(initial_mass_list_path, index=False)

    # 2) Prepare peaks DataFrame for formula assignment (in memory)
    peaks_for_formulas_df = pd.DataFrame(
        {
            "m/z": filtered_peaks_df[peak_cols.mz],
            "Intensity": filtered_peaks_df[peak_cols.intensity],
        }
    )

    # 3) Run formula assignment on the filtered peaks (in memory)
    formulas_df = run_formula_assignment(
        peaks_for_formulas_df=peaks_for_formulas_df,
        sample_type=sample_type,
        ppm_tolerance=ppm_tolerance,
        n_processes=n_processes,
    )

    # 4) Select calibrant peaks from the formulas
    calibrants = select_calibrant_peaks(formulas_df, cfg, cols=formula_cols)

    # 5) Fit polynomial calibration
    coeffs = fit_calibration_polynomial(calibrants, cfg)

    # 6) Apply calibration to the full peak list
    calibrated_peaks = apply_calibration_to_peaks(peaks_df, coeffs, peak_cols=peak_cols)

    # 7) Build "full_ms"-style dataframe
    full_ms_like = build_full_ms_like_output(
        calibrated_peaks,
        calibrants,
        cfg,
        peak_cols=peak_cols,
    )

    # ------------------------------------------------------------------
    # Debug / test outputs in TEST_CALIB_DIR:
    #
    # (a) Selected calibrant mass list with:
    #     initial m/z, calibrated m/z, ppm error, formula, intensity,
    #     S/N, homologous series length
    # (b) Final calibrated mass list with:
    #     m/z, calibrated m/z, ppm error, intensity, S/N
    # ------------------------------------------------------------------

    # (a) Build calibrant debug table
    poly = np.poly1d(coeffs)
    calib_debug = calibrants.copy()
    calib_debug["Calibrated m/z"] = poly(calib_debug["mz_obs"].to_numpy(dtype=float))
    calib_debug["ppm_error"] = ppm_error(
        calib_debug["mz_true"].to_numpy(dtype=float),
        calib_debug["Calibrated m/z"].to_numpy(dtype=float),
    )

    # Add S/N by joining back to filtered_peaks_df on observed m/z
    sn_map = filtered_peaks_df[[peak_cols.mz, peak_cols.sn]].rename(
        columns={
            peak_cols.mz: "mz_obs",
            peak_cols.sn: "S/N",
        }
    )
    calib_debug = calib_debug.merge(sn_map, on="mz_obs", how="left")

    calibrant_mass_list = calib_debug[
        ["mz_obs", "Calibrated m/z", "ppm_error", "formula", "Intensity", "S/N", "series_length"]
    ].rename(columns={"mz_obs": "m/z"})

    calibrant_mass_list_path = TEST_CALIB_DIR / f"{sample_name}_calibrants_mass_list.csv"
    calibrant_mass_list.to_csv(calibrant_mass_list_path, index=False)

    # (b) Final calibrated mass list
    final_mass_list = pd.DataFrame(
        {
            "m/z": full_ms_like["m/z"].to_numpy(dtype=float),
            "Calibrated m/z": full_ms_like["Calibrated m/z"].to_numpy(dtype=float),
            "m/z Error (ppm)": full_ms_like["m/z Error (ppm)"].to_numpy(dtype=float),
            "Intensity": full_ms_like["Peak Height"].to_numpy(dtype=float),
            "S/N": full_ms_like["S/N"].to_numpy(dtype=float),
        }
    )
    final_mass_list_path = TEST_CALIB_DIR / f"{sample_name}_final_calibrated_mass_list.csv"
    final_mass_list.to_csv(final_mass_list_path, index=False)

    # Optionally annotate the full output with the sample_type used for formula presets
    if sample_type is not None:
        # Store a human-readable sample type label for each peak
        full_ms_like["Sample type"] = sample_type

    # Save full calibrated results (CSV + PKL) for downstream GUI use
    save_calibrated_results(full_ms_like, sample_name, results_dir=results_dir)
    print(f"Calibrated peaks saved into:\n  CSV and PKL: {results_dir}\n  ")

    # Collect debug info that higher-level code (or a GUI) may want to inspect
    debug_info = {
        "sample_name": sample_name,
        "coeffs": coeffs,
        "calibrants": calibrants,
        "calibrant_mass_list": calibrant_mass_list,
        "final_mass_list": final_mass_list,
    }

    # Main return value is the calibrated “full ms” style DataFrame;
    # debug_info contains extra diagnostics for advanced use.
    return full_ms_like, debug_info


if __name__ == "__main__":
    # Example CLI entry point (adjust paths for your Docker mount)
    default_peaks = Path("/app/results/csv_raw/ESI_neg_G017736-0_100_200sc_000001.d_peaks.csv")

    cfg = CalibrationConfig(sn_threshold=20.0, poly_order=2, only_CHO=True)
    run_calibration(peaks=default_peaks, cfg=cfg)
    # csv_out, pkl_out = run_calibration(
    #     peaks=default_peaks,
    #     cfg=cfg,
    # )
    