#!/usr/bin/env python3
# SPDX-License-Identifier: EUPL-1.2
from pathlib import Path
import pandas as pd

# Import the patched reader and the peak extraction helper
from SCHON.corems_pipeline.read_bruker_fid_corems import (
    PatchedReadBrukerSolarix,
    spectrum_to_peaks_dataframe,
)


def load_mass_spectrum_from_d_folder(d_folder: str):
    """
    Load a Bruker .d directory using CoreMS and return
    a MassSpectrum object.

    Uses PatchedReadBrukerSolarix, which safely handles
    missing EXC_Freq_High / EXC_Freq_Low.
    """
    d_path = Path(d_folder)
    if not d_path.is_dir():
        raise NotADirectoryError(f"{d_folder} is not a directory")

    reader = PatchedReadBrukerSolarix(str(d_path))
    transient = reader.get_transient()
    spectrum = transient.get_mass_spectrum(
        plot_result=False,
        auto_process=True,
    )
    return spectrum


def export_corems_formats(d_folder: str, out_dir: str):
    """
    For a single .d directory:
        - reads the mass spectrum with CoreMS
        - saves:
            * full mass spectrum: CSV, HDF5, pandas pickle
            * peak table (m/z, intensity, S/N): CSV
    """
    d_path = Path(d_folder)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    base_name = d_path.name  # e.g. "ESI_neg_G017736-0_100_200sc_000001.d"

    print(f"🔹 Processing .d folder: {d_path}")

    spectrum = load_mass_spectrum_from_d_folder(str(d_path))

    # --- Save full mass spectrum in CoreMS output formats ---
    csv_file = out_dir / f"{base_name}_full_ms.csv"
    hdf_file = out_dir / f"{base_name}_full_ms.h5"
    pkl_file = out_dir / f"{base_name}_full_ms.pkl"

    spectrum.to_csv(str(csv_file))
    spectrum.to_hdf(str(hdf_file))
    spectrum.to_pandas(str(pkl_file))

    print(f"   ✅ Saved full spectrum CSV   → {csv_file}")
    print(f"   ✅ Saved full spectrum HDF5  → {hdf_file}")
    print(f"   ✅ Saved full spectrum PANDAS→ {pkl_file}")

    # --- Save simplified peak table (m/z, intensity, S/N) ---
    peaks_df = spectrum_to_peaks_dataframe(spectrum)
    peaks_csv = out_dir / f"{base_name}_peaks.csv"
    peaks_df.to_csv(peaks_csv, index=False)

    print(f"   ✅ Saved peaks CSV → {peaks_csv}")


if __name__ == "__main__":
    # Default example path used for testing inside Docker
    d_path = "/app/Measurement_d-files/calibrated_d-files/ESI_neg_G017736-0_100_200sc_000001.d"
    out_dir = "/app/SCHON/output_peaks"

    export_corems_formats(d_path, out_dir)