import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.spatial import KDTree
from multiprocessing import Pool, cpu_count
from schon.formula_assignment.formula_filters import generate_formulas

# Home directory paths
# HOME_DIR = "/Users/anya/Coding/ML_FT-ICR-MS/SCHON/"
HOME_DIR = "/app/"
PEAKS_DIR = os.path.join(HOME_DIR, "data/input_peaks")
FORMULAS_DIR = os.path.join(HOME_DIR, "results/output_formulas")
os.makedirs(FORMULAS_DIR, exist_ok=True)

# Constants
PPM_TOLERANCE = 5
ISOTOPE_DIFF = 1.0033548378

SAMPLE_TYPE = 0 #'crude_oil', 'sedimentary_rock', 'coal', 'natural_water', and other

# ---------- Precompute and load formulas globally ----------
print("Precomputing molecular formulas...")
#formulas_data = generate_formulas()
formulas_data = generate_formulas(SAMPLE_TYPE)
formula_tree = KDTree(formulas_data['Mass'].values.reshape(-1, 1))
print(f"✅ Loaded {len(formulas_data)} molecular formulas (including Na adducts).")


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


# ----------- MAIN FILE PROCESSING (Handle M1 inheritance properly) -------------
def process_peak_file(file_path):
    print(f"🔍 Processing peak file: {file_path}")
    df = pd.read_csv(file_path)
    #print(df.head())
    df, m0_m1_pairs = deisotope_peaks(df)

    # Assign formulas only for M0 and Q
    assign_data = df[df['isotopolog'].isin(['M0', 'Q'])]
    with Pool(cpu_count() - 1) as pool:
        #print(cpu_count)
        results = list(tqdm(
            pool.starmap(assign_formula, zip(assign_data['m/z'], assign_data['Intensity'], assign_data['isotopolog'])),
            total=len(assign_data),
            desc="Assigning Formulas"
        ))

    assign_df = assign_data[['m/z']].copy()
    assign_df[['Formula', 'Calculated Mass', 'C', 'H', 'N', 'O', 'S', '*C', 'Na', 'Mass Error (ppm)', 'Alternative Formula', 'Alternative Mass Error (ppm)']] = results
    df = pd.merge(df, assign_df, on='m/z', how='left')

    # Ensure M1 inherits M0's formula but adjusted: one less C, one more *C
    for m0_idx, m1_idx in m0_m1_pairs.items():
        if pd.notna(df.loc[m0_idx, 'Formula']):
            c = df.loc[m0_idx, 'C'] - 1
            c13 = df.loc[m0_idx, '*C'] + 1
            h = df.loc[m0_idx, 'H']
            n = df.loc[m0_idx, 'N']
            o = df.loc[m0_idx, 'O']
            s = df.loc[m0_idx, 'S']
            na = df.loc[m0_idx, 'Na']

            # Regenerate correct formula string for M1
            formula_parts = [f"C{c}", f"*C{c13}"]
            if h > 0:
                formula_parts.append(f"H{h}")
            if n > 0:
                formula_parts.append(f"N{n}")
            if o > 0:
                formula_parts.append(f"O{o}")
            if s > 0:
                formula_parts.append(f"S{s}")
            if na > 0:
                formula_parts.append(f"Na{na}")
            formula_str = "".join(formula_parts)

            # Calculated mass shift for M1
            calc_mass_m1 = df.loc[m0_idx, 'Calculated Mass'] + ISOTOPE_DIFF
            mz_m1 = df.loc[m1_idx, 'm/z']
            mass_error_m1 = abs(mz_m1 - calc_mass_m1) / mz_m1 * 1e6

            # Assign updated counts and recalculated values to M1
            df.loc[m1_idx, ['Formula', 'Calculated Mass', 'C', '*C', 'H', 'N', 'O', 'S', 'Na', 'Mass Error (ppm)']] = \
                [formula_str, calc_mass_m1, c, c13, h, n, o, s, na, mass_error_m1]

    # Save final dataframe
    output_file = os.path.join(FORMULAS_DIR, os.path.basename(file_path).replace("_peaks.csv", "_formulas.csv"))
    df.to_csv(output_file, index=False)
    print(f"✅ Formula assignment complete and saved to {output_file}")

# ----------- RUN ----------
if __name__ == "__main__":
    for file in os.listdir(PEAKS_DIR):
        if file.endswith("_peaks.csv"):
            process_peak_file(os.path.join(PEAKS_DIR, file))