import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.spatial import KDTree
from multiprocessing import Pool, cpu_count
from schon.formula_assignment.formula_filters import generate_formulas
from pathlib import Path

# ----------------- Default paths and constants -----------------
DEFAULT_HOME_DIR = Path(os.environ.get("SCHON_HOME_DIR", "/app")).resolve()
DEFAULT_PEAKS_DIR = DEFAULT_HOME_DIR / "data" / "input_peaks"
DEFAULT_FORMULAS_DIR = DEFAULT_HOME_DIR / "results" / "output_formulas"
DEFAULT_FORMULAS_DIR.mkdir(parents=True, exist_ok=True)

# Numeric defaults (can be overridden via init_formula_search)
DEFAULT_PPM_TOLERANCE: float = 3.0
ISOTOPE_DIFF: float = 1.0033548378
DEFAULT_SAMPLE_TYPE: str = "generic_esi_neg"  # default preset key used in formula_presets

# Globals initialised lazily so the module can be imported without heavy work
formulas_data = None
formula_tree: KDTree | None = None
PPM_TOLERANCE: float = DEFAULT_PPM_TOLERANCE
SAMPLE_TYPE: str = DEFAULT_SAMPLE_TYPE

# ----------------- Formula search initialisation -----------------
def init_formula_search(sample_type: str = DEFAULT_SAMPLE_TYPE,
                        ppm_tolerance: float = DEFAULT_PPM_TOLERANCE,
                        force_recompute: bool = False) -> None:
    """Initialise global formula table and KDTree for fast mass lookups.

    Parameters
    ----------
    sample_type : str
        String key passed to ``generate_formulas`` (e.g. "crude_oil", "generic_esi_neg").
    ppm_tolerance : float
        Default ppm tolerance used during formula matching.
    force_recompute : bool
        If True, regenerate the formula table even if it already exists.
    """
    global formulas_data, formula_tree, PPM_TOLERANCE, SAMPLE_TYPE

    if formulas_data is not None and formula_tree is not None and not force_recompute:
        # Already initialised with some configuration; only update tolerance
        PPM_TOLERANCE = ppm_tolerance
        SAMPLE_TYPE = sample_type
        return

    print("Precomputing molecular formulas...")
    formulas = generate_formulas(sample_type)
    tree = KDTree(formulas["Mass"].values.reshape(-1, 1))
    print(f"✅ Loaded {len(formulas)} molecular formulas (including Na adducts).")

    formulas_data = formulas
    formula_tree = tree
    PPM_TOLERANCE = ppm_tolerance
    SAMPLE_TYPE = sample_type


# ---------- DYNAMIC TOLERANCE FUNCTION ----------
def dynamic_isotope_tolerance(mz):
    ppm_low = 0.5
    ppm_high = 1
    if mz <= 100:
        ppm = ppm_low
    elif mz >= 900:
        ppm = ppm_high
    else:
        ppm = ppm_low + (ppm_high - ppm_low) * ((mz - 100) / (900 - 100))
    return (ppm / 1e6) * mz


# ----------- DEISOTOPING FUNCTION (Tracking M0-M1 pairs) -------------
def deisotope_peaks(df):
    df = df.sort_values(by='m/z').reset_index(drop=True)
    df['isotopolog'] = 'Q'
    m0_m1_pairs = {}

    for i, row in df.iterrows():
        if df.loc[i, 'isotopolog'] != 'Q':
            continue
        mz_0, intensity_0 = row['m/z'], row['Intensity']
        mz_1 = mz_0 + ISOTOPE_DIFF
        tol1 = dynamic_isotope_tolerance(mz_1)
        candidate_1 = df[(np.abs(df['m/z'] - mz_1) <= tol1)]

        if not candidate_1.empty:
            idx_1 = candidate_1.index[0]
            df.at[i, 'isotopolog'] = 'M0'
            df.at[idx_1, 'isotopolog'] = 'M1'
            m0_m1_pairs[i] = idx_1

            # Candidate for M2
            mz_2 = df.loc[idx_1, 'm/z'] + ISOTOPE_DIFF
            tol2 = dynamic_isotope_tolerance(mz_2)
            candidate_2 = df[(np.abs(df['m/z'] - mz_2) <= tol2)]
            if not candidate_2.empty:
                df.at[candidate_2.index[0], 'isotopolog'] = 'M2'

    print(f"Isotopolog categories:\n{df['isotopolog'].value_counts()}")
    return df, m0_m1_pairs


# ----------- FORMULA ASSIGNMENT (Only for M0 and Q) -------------
def assign_formula(mz, intensity, isotopolog):
    global formulas_data, formula_tree
    if formulas_data is None or formula_tree is None:
        raise RuntimeError("Formula search not initialised. Call init_formula_search() first.")

    if not isinstance(mz, (int, float)):
        return (np.nan,) * 12

    mass_tol = mz * PPM_TOLERANCE / 1e6
    idx = formula_tree.query_ball_point([mz], mass_tol)
    if not idx:
        return (np.nan,) * 12

    best_match, best_error = None, float("inf")
    for i in idx:
        formula = formulas_data.iloc[i]
        c13 = formula['*C']
        na = formula['Na']
        # Q peaks cannot have 13C
        if isotopolog == 'Q' and c13 != 0:
            continue
        # Check *C consistency for M0 and M2 only
        if (isotopolog == 'M0' and c13 != 0) or (isotopolog == 'M2' and c13 != 2):
            continue
        error = abs(mz - formula['Mass']) / mz * 1e6
        if error < best_error:
            best_match, best_error = formula, error

    if best_match is not None:
        # If Na is present, subtract 2 hydrogens for [M + Na - 2H]-, otherwise subtract 1 for [M - H]-
        h_deprotonated = best_match['H'] - (2 if best_match['Na'] > 0 else 1)
        h_deprotonated = max(h_deprotonated, 0)

        # Replace hydrogen count in formula string
        formula_str = best_match['Formula'].replace(f"H{best_match['H']}", f"H{h_deprotonated}")

        return (
            formula_str, best_match['Mass'], best_match['C'], h_deprotonated,
            best_match['N'], best_match['O'], best_match['S'], best_match['*C'],
            best_match['Na'], abs(mz - best_match['Mass']) / mz * 1e6, np.nan, np.nan
        )
    return (np.nan,) * 12



# ----------------- In-memory DataFrame formula assignment -----------------
def assign_formulas_df(
    df: pd.DataFrame,
    sample_type: str = DEFAULT_SAMPLE_TYPE,
    ppm_tolerance: float = DEFAULT_PPM_TOLERANCE,
    n_processes: int | None = None,
) -> pd.DataFrame:
    """Assign molecular formulas directly on a peak list DataFrame.

    This is the in-memory counterpart of :func:`assign_formulas`, used by
    the calibration pipeline. The input ``df`` must contain at least the
    columns ``"m/z"`` and ``"Intensity"``. The function returns a *new*
    DataFrame with the following additional columns:

    - ``isotopolog`` (Q/M0/M1/M2)
    - ``Formula``
    - ``Calculated Mass``
    - ``C, H, N, O, S, *C, Na``
    - ``Mass Error (ppm)``
    - ``Alternative Formula``
    - ``Alternative Mass Error (ppm)``

    Parameters
    ----------
    df : pandas.DataFrame
        Input peak table (must have at least ``m/z`` and ``Intensity``).
    sample_type : str
        Sample-type preset name (string) forwarded to :func:`generate_formulas` via
        :func:`init_formula_search`.
    ppm_tolerance : float
        Default ppm tolerance for mass matching in :func:`assign_formula`.
    n_processes : int, optional
        Number of worker processes to use for formula assignment. If
        ``None``, uses ``max(cpu_count() - 1, 1)``.

    Returns
    -------
    pandas.DataFrame
        A copy of the input DataFrame with formula-related columns added.
    """

    # Ensure formula search is initialised
    init_formula_search(sample_type=sample_type, ppm_tolerance=ppm_tolerance)

    # Work on a copy to avoid mutating the caller's DataFrame in-place
    df = df.copy()

    # Deisotope and get mapping of M0->M1 pairs
    df, m0_m1_pairs = deisotope_peaks(df)

    # Assign formulas only for M0 and Q
    mask_assign = df["isotopolog"].isin(["M0", "Q"])
    assign_data = df.loc[mask_assign]

    if n_processes is None:
        n_processes = max(cpu_count() - 1, 1)

    if len(assign_data) > 0:
        with Pool(n_processes) as pool:
            results = list(
                tqdm(
                    pool.starmap(
                        assign_formula,
                        zip(
                            assign_data["m/z"],
                            assign_data["Intensity"],
                            assign_data["isotopolog"],
                        ),
                    ),
                    total=len(assign_data),
                    desc="Assigning Formulas",
                )
            )

        # Build a DataFrame with the results, indexed like assign_data
        result_cols = [
            "Formula",
            "Calculated Mass",
            "C",
            "H",
            "N",
            "O",
            "S",
            "*C",
            "Na",
            "Mass Error (ppm)",
            "Alternative Formula",
            "Alternative Mass Error (ppm)",
        ]
        assign_df = pd.DataFrame(results, index=assign_data.index, columns=result_cols)

        # Ensure all target columns exist in df, then fill them for assigned rows
        for col in result_cols:
            if col not in df.columns:
                df[col] = np.nan
        df.loc[assign_df.index, result_cols] = assign_df
    else:
        # No candidates; create empty columns if they don't exist
        for col in [
            "Formula",
            "Calculated Mass",
            "C",
            "H",
            "N",
            "O",
            "S",
            "*C",
            "Na",
            "Mass Error (ppm)",
            "Alternative Formula",
            "Alternative Mass Error (ppm)",
        ]:
            if col not in df.columns:
                df[col] = np.nan

    # Ensure M1 inherits M0's formula with one less 12C and one more 13C
    for m0_idx, m1_idx in m0_m1_pairs.items():
        if m0_idx not in df.index or m1_idx not in df.index:
            continue
        if pd.notna(df.loc[m0_idx, "Formula"]):
            c = df.loc[m0_idx, "C"] - 1
            c13 = df.loc[m0_idx, "*C"] + 1
            h = df.loc[m0_idx, "H"]
            n = df.loc[m0_idx, "N"]
            o = df.loc[m0_idx, "O"]
            s = df.loc[m0_idx, "S"]
            na = df.loc[m0_idx, "Na"]

            formula_parts = [f"C{int(c)}", f"*C{int(c13)}"]
            if h > 0:
                formula_parts.append(f"H{int(h)}")
            if n > 0:
                formula_parts.append(f"N{int(n)}")
            if o > 0:
                formula_parts.append(f"O{int(o)}")
            if s > 0:
                formula_parts.append(f"S{int(s)}")
            if na > 0:
                formula_parts.append(f"Na{int(na)}")
            formula_str = "".join(formula_parts)

            calc_mass_m1 = df.loc[m0_idx, "Calculated Mass"] + ISOTOPE_DIFF
            mz_m1 = df.loc[m1_idx, "m/z"]
            mass_error_m1 = abs(mz_m1 - calc_mass_m1) / mz_m1 * 1e6

            df.loc[m1_idx, [
                "Formula",
                "Calculated Mass",
                "C",
                "*C",
                "H",
                "N",
                "O",
                "S",
                "Na",
                "Mass Error (ppm)",
            ]] = [
                formula_str,
                calc_mass_m1,
                c,
                c13,
                h,
                n,
                o,
                s,
                na,
                mass_error_m1,
            ]

    # --- Ensure CoreMS-like output columns and ordering ---
    # Sample type column
    if sample_type is not None:
        df["Sample type"] = sample_type
    elif "Sample type" not in df.columns:
        df["Sample type"] = np.nan

    # Ion Charge: default to -1 for negative mode if missing
    if "Ion Charge" not in df.columns:
        df["Ion Charge"] = -1

    # S/N column: if missing, create as NaN
    if "S/N" not in df.columns:
        df["S/N"] = np.nan

    # Optional formula-derived / classification columns
    for col, default in [
        ("DBE", np.nan),
        ("O/C", np.nan),
        ("H/C", np.nan),
        ("Heteroatom Class", ""),
        ("Adduct", ""),
    ]:
        if col not in df.columns:
            df[col] = default

    # Calibrated m/z: if missing, fall back to raw m/z
    if "Calibrated m/z" not in df.columns:
        df["Calibrated m/z"] = df["m/z"]

    # Is Isotopologue: derived from isotopolog (True for non-Q)
    if "isotopolog" in df.columns:
        df["Is Isotopologue"] = df["isotopolog"].astype(str).str.upper() != "Q"
    else:
        df["Is Isotopologue"] = False

    # Reset the pandas index, but preserve any existing "Index" column
    df = df.reset_index(drop=True)
    if "Index" not in df.columns:
        # Create a stable Index column only if it does not exist yet
        df.insert(0, "Index", np.arange(len(df)))

    # Reorder columns to a CoreMS-like layout
    desired_order = [
        "Index",
        "m/z",
        "Calibrated m/z",
        "Intensity",
        "S/N",
        "Ion Charge",
        "Formula",
        "Calculated Mass",
        "Mass Error (ppm)",
        "C",
        "H",
        "N",
        "O",
        "S",
        "*C",
        "Na",
        "Alternative Formula",
        "Alternative Mass Error (ppm)",
        "DBE",
        "O/C",
        "H/C",
        "Heteroatom Class",
        "Adduct",
        "Is Isotopologue",
        "Sample type",
    ]

    existing_cols = [c for c in df.columns if c not in desired_order]
    df = df[desired_order + existing_cols]

    return df


def assign_formulas(input_csv: str | Path,
                    output_csv: str | Path | None = None,
                    sample_type: str = DEFAULT_SAMPLE_TYPE,
                    ppm_tolerance: float = DEFAULT_PPM_TOLERANCE,
                    n_processes: int | None = None) -> Path:
    """Assign molecular formulas to a single peak list CSV.

    Parameters
    ----------
    input_csv : str or Path
        Path to the input peaks CSV. Must contain at least the columns
        ``m/z`` and ``Intensity``.
    output_csv : str or Path, optional
        Where to save the annotated table. If ``None``, the file will be
        written into ``DEFAULT_FORMULAS_DIR`` with ``"_formulas"`` appended
        to the stem (and the ``"_peaks"`` suffix removed if present).
    sample_type : str
        Sample-type preset key (string) forwarded to :func:`generate_formulas` via
        :func:`init_formula_search`.
    ppm_tolerance : float
        Default ppm tolerance for mass matching in :func:`assign_formula`.
    n_processes : int, optional
        Number of worker processes to use for formula assignment. If ``None``,
        uses ``max(cpu_count() - 1, 1)``.

    Returns
    -------
    Path
        Path to the written output CSV.
    """

    # Ensure output path is set
    input_csv = Path(input_csv)
    if output_csv is None:
        stem = input_csv.stem
        if stem.endswith("_peaks"):
            stem = stem[:-6]
        output_csv = DEFAULT_FORMULAS_DIR / f"{stem}_formulas.csv"
    else:
        output_csv = Path(output_csv)
        output_csv.parent.mkdir(parents=True, exist_ok=True)

    print(f"🔍 Processing peak file: {input_csv}")
    df_in = pd.read_csv(input_csv)

    # Run in-memory formula assignment
    df_out = assign_formulas_df(
        df_in,
        sample_type=sample_type,
        ppm_tolerance=ppm_tolerance,
        n_processes=n_processes,
    )

    output_csv = output_csv.resolve()
    df_out.to_csv(output_csv, index=False)
    print(f"✅ Formula assignment complete and saved to {output_csv}")
    return output_csv


if __name__ == "__main__":
    # Simple CLI: process all *_peaks.csv files in the default peaks folder
    if not DEFAULT_PEAKS_DIR.exists():
        raise SystemExit(f"Input peaks directory does not exist: {DEFAULT_PEAKS_DIR}")

    for file in sorted(DEFAULT_PEAKS_DIR.glob("*_peaks.csv")):
        assign_formulas(file)