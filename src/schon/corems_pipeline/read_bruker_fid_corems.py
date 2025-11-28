# SPDX-License-Identifier: EUPL-1.2
from pathlib import Path
import pandas as pd
from corems.transient.input import brukerSolarix as bruker_module


class PatchedReadBrukerSolarix(bruker_module.ReadBrukerSolarix):
    """
    Subclass of the CoreMS BrukerSolarix reader that avoids crashing when
    EXC_Freq_High / EXC_Freq_Low are missing from the Bruker parameters.
    """

    def fix_freq_limits(self, d_parameters):
        high = d_parameters.get("EXC_Freq_High")
        low = d_parameters.get("EXC_Freq_Low")

        if high is None or low is None:
            print(
                "Warning: EXC_Freq_High or EXC_Freq_Low is missing; "
                "skipping fix_freq_limits and using default frequency limits."
            )
            return

        return super().fix_freq_limits(d_parameters)


def load_mass_spectrum_from_d_folder(d_folder: str):
    """
    Read a Bruker Solarix .d folder using CoreMS and return a MassSpectrum object.
    """
    d_path = Path(d_folder)
    if not d_path.is_dir():
        raise NotADirectoryError(f"{d_folder} is not a directory")

    # This matches the CoreMS example: reader -> transient -> mass spectrum
    bruker_reader = PatchedReadBrukerSolarix(str(d_path))
    transient_obj = bruker_reader.get_transient()
    mass_spectrum_obj = transient_obj.get_mass_spectrum(
        plot_result=False,
        auto_process=True,
    )
    return mass_spectrum_obj


def spectrum_to_peaks_dataframe(mass_spectrum_obj) -> pd.DataFrame:
    """
    Convert a CoreMS MassSpectrum object to a pandas.DataFrame with
    headers expected by the downstream pipeline:

        - signalNumber
        - observedExactMass_ion   (m/z)
        - observedAbundance_ion   (Intensity)
        - signalNoise_ratio       (S/N)

    The function tries to map these from the columns returned by
    MassSpectrum.to_dataframe(), handling several possible name variants.
    """
    df = mass_spectrum_obj.to_dataframe()

    # --- Helper to find a column by possible aliases ---
    def find_col(possible_names):
        for name in possible_names:
            if name in df.columns:
                return name
        # also try case-insensitive match
        lower_map = {c.lower(): c for c in df.columns}
        for name in possible_names:
            if name.lower() in lower_map:
                return lower_map[name.lower()]
        return None

    # Try to locate m/z, intensity, and S/N columns in the CoreMS dataframe
    mz_col = find_col(["observedExactMass_ion", "m/z", "mz"])
    intensity_col = find_col(
        ["observedAbundance_ion", "Intensity", "abundance", "Abundance", "I","Peak Height"]
    )
    sn_col = find_col(
        ["signalNoise_ratio", "SNR", "S/N", "sn", "signal_to_noise"]
    )
    signal_number_col = find_col(["signalNumber", "signal_number", "peak_id"])

    # Build the output dataframe with the exact headers the user requested.
    out = pd.DataFrame()

    # signalNumber: use existing column if present, otherwise 1..N
    if signal_number_col is not None:
        out["signalNumber"] = df[signal_number_col].values
    else:
        out["signalNumber"] = range(1, len(df) + 1)

    # observedExactMass_ion (m/z)
    if mz_col is not None:
        out["observedExactMass_ion"] = df[mz_col].values
    else:
        out["observedExactMass_ion"] = pd.NA

    # observedAbundance_ion (Intensity)
    if intensity_col is not None:
        out["observedAbundance_ion"] = df[intensity_col].values
    else:
        out["observedAbundance_ion"] = pd.NA

    # signalNoise_ratio (S/N)
    if sn_col is not None:
        out["signalNoise_ratio"] = df[sn_col].values
    else:
        out["signalNoise_ratio"] = pd.NA

    # Sort by m/z if we have it
    if "observedExactMass_ion" in out.columns and out["observedExactMass_ion"].notna().any():
        out = out.sort_values("observedExactMass_ion").reset_index(drop=True)

    return out