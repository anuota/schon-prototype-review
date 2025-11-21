"""
Calibration pipeline for FT-ICR-MS peaks using internally assigned formulas.

Steps
-----
1. Load peaks from `<sample>_peaks.csv`.
   Expected columns:
       signalNumber, observedExactMass_ion, observedAbundance_ion, signalNoise_ratio

2. Filter peaks by S/N >= threshold (default 20).

3. (External step – not called here yet)
   Run your formula assignment pipeline on the filtered peaks to produce a CSV with
   columns like:
       m/z, Intensity, isotopolog, Formula, Calculated Mass, C, H, N, O, S, *C, Na,
       Mass Error (ppm), Alternative Formula, Alternative Mass Error (ppm)

4. From that formulas CSV, select robust internal calibrants:
   - isotopolog == 'Q'  (monoisotopic)
   - no "Alternative Formula"
   - optionally CHO-only (N == S == Na == 0, *C == 0)
   - must belong to a homologous series (e.g. CH2) of length >= 2

5. Fit a polynomial calibration m/z_true = f(m/z_observed).

6. Apply this calibration to the *full* peak list and save:
   `<sample>_calibrated_peaks.csv` and `<sample>_calibrated_peaks.pkl` in /app/results
   with columns including:
       Index, m/z, Calibrated m/z, Calculated m/z, Peak Height, Peak Area, Resolving Power,
       S/N, Ion Charge, m/z Error (ppm), ...

This module deliberately keeps things:
- configurable (thresholds, file paths, polynomial order),
- independent of CoreMS internals (uses numpy.polyfit),
- and easy to extend later with a more sophisticated calibration (e.g. Kozhinov–Savory).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence, Tuple, List

import numpy as np
import pandas as pd


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
) -> np.ndarray:
    """Return a boolean mask for peaks that belong to a CH2-like series.

    The algorithm is intentionally simple: we build an undirected graph
    where two peaks are connected if their m/z difference corresponds to
    an integer number of CH2 units within a ppm tolerance. Peaks that
    belong to a connected component of size >= min_series_length are
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


@dataclass
class CalibrationConfig:
    """Config for the calibration process."""
    sn_threshold: float = 20.0
    poly_order: int = 2
    only_CHO: bool = True              # use only CHO-type formulas as calibrants
    ppm_tolerance_link: float = 0.5    # max |m/z - m/z_calibrant| (ppm) to link back


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


# NOTE: This is intentionally left as a hook rather than a hard dependency
# on your current `main.py`, because that file still contains hard-coded
# paths and pre-generated formulas. Once you refactor `main.py` into a
# clean function (e.g. `assign_formulas(input_csv: Path, output_csv: Path)`),
# you can call it from here.
#
def run_formula_assignment(input_csv: Path, output_csv: Path) -> None:
    from schon.formula_assignment import formula_assignment as formula_main
    formula_main.assign_formulas(input_csv=input_csv, output_csv=output_csv)


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

    calibrants = df.loc[mask, [cols.mz, cols.calc_mass]].copy()
    calibrants.rename(
        columns={cols.mz: "mz_obs", cols.calc_mass: "mz_true"},
        inplace=True,
    )

    # drop rows with missing values
    calibrants = calibrants.dropna(subset=["mz_obs", "mz_true"])

    if calibrants.empty:
        raise RuntimeError("No calibrant peaks selected – check your filters / formulas.")

    # ------------------------------------------------------------------
    # Additional filter: keep only peaks that are part of a homologous
    # series (e.g. CH2) of length >= 2. We use the theoretical masses
    # (mz_true) for this check.
    # ------------------------------------------------------------------
    series_mask = detect_homologue_series(
        calibrants["mz_true"].to_numpy(dtype=float),
        min_series_length=2,
        step_mass=CH2_MASS,
        ppm_tolerance=5.0,
    )

    calibrants = calibrants.loc[series_mask].copy()
    if calibrants.empty:
        raise RuntimeError(
            "No calibrant peaks remain after enforcing the homologous-series criterion."
        )

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

        # And we can add formula information later if we want

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
    peaks_path: Path,
    formulas_path: Path,
    cfg: CalibrationConfig = CalibrationConfig(),
    peak_cols: PeakColumns = PeakColumns(),
    formula_cols: FormulaColumns = FormulaColumns(),
    results_dir: Path = Path("/app/results"),
) -> Tuple[Path, Path]:
    """
    End-to-end calibration, starting from:
        - peaks_path:  ..._peaks.csv
        - formulas_path: formulas for (high S/N) peaks with columns described above.

    Returns:
        (csv_path, pkl_path) of calibrated output.
    """
    sample_name = peaks_path.stem.replace("_peaks", "")

    # 1) load peaks
    peaks_df = load_peaks(peaks_path, cols=peak_cols)

    # 2) filter by S/N (this is independent of formulas_path; you *should*
    #    generate formulas based on these filtered peaks, but here we assume
    #    formulas_path already corresponds to that subset).
    filtered_peaks_df = filter_peaks_by_snr(peaks_df, cfg, cols=peak_cols)

    # 3) load formulas and select calibrants
    formulas_df = load_formulas(formulas_path, cols=formula_cols)
    calibrants = select_calibrant_peaks(formulas_df, cfg, cols=formula_cols)

    # 4) fit polynomial
    coeffs = fit_calibration_polynomial(calibrants, cfg)

    # 5) apply calibration to full peaks (not just filtered ones)
    calibrated_peaks = apply_calibration_to_peaks(peaks_df, coeffs, peak_cols=peak_cols)

    # 6) build "full_ms" like dataframe and save
    full_ms_like = build_full_ms_like_output(calibrated_peaks, calibrants, cfg, peak_cols=peak_cols)
    return save_calibrated_results(full_ms_like, sample_name, results_dir=results_dir)


if __name__ == "__main__":
    # Example CLI entry point (adjust paths for your Docker mount)
    default_peaks = Path("/app/data/output_peaks/ESI_neg_G017736-0_100_200sc_000001.d_peaks.csv")
    default_formulas = Path("/app/data/output_formulas/ESI_neg_G017736-0_100_200sc_000001.d_formulas.csv")

    cfg = CalibrationConfig(sn_threshold=20.0, poly_order=2, only_CHO=True)

    csv_out, pkl_out = run_calibration(
        peaks_path=default_peaks,
        formulas_path=default_formulas,
        cfg=cfg,
    )
    print(f"Calibrated peaks saved to:\n  CSV: {csv_out}\n  PKL: {pkl_out}")