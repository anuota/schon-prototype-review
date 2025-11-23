from pathlib import Path
from typing import Optional

from schon.corems_pipeline.read_bruker_fid_corems import (
    load_mass_spectrum_from_d_folder,
    spectrum_to_peaks_dataframe,
)

# ---------- Defaults (can be overridden by GUI later) ----------
PROJECT_ROOT = Path("/app")
DATA_FOLDER = PROJECT_ROOT / "data"
DEFAULT_D_FOLDER = DATA_FOLDER / "Measurement_d-files" / "noncalibrated_d-files" / "ESI_neg_G017736-0_100_200sc_000001.d"
RESULTS_FOLDER = PROJECT_ROOT / "results"
CSV_RAW_DIR = RESULTS_FOLDER / "csv_raw"

CSV_RAW_DIR.mkdir(parents=True, exist_ok=True)


def extract_peaks_from_d_folder(
    d_path: Path,
    save_csv: bool = True,
    output_dir: Optional[Path] = None,
):
    """
    Load a Bruker .d folder using CoreMS, extract its peaks as a pandas DataFrame,
    optionally save the DataFrame into CSV, and return the DataFrame.

    Parameters
    ----------
    d_path : Path
        Path to the Bruker .d folder.
    save_csv : bool, optional
        Whether to write the resulting peaks into a CSV file. Default: True.
    output_dir : Path, optional
        If provided, CSV will be written to this directory.
        If None, defaults to CSV_RAW_DIR.

    Returns
    -------
    peaks_df : pandas.DataFrame
        DataFrame containing uncalibrated peaks.
    """

    if not d_path.exists():
        raise FileNotFoundError(f".d folder does not exist: {d_path}")
    if not d_path.is_dir() or d_path.suffix != ".d":
        raise ValueError(f"Path must be a .d folder: {d_path}")

    print(f"🔹 Loading mass spectrum from: {d_path}")
    ms_obj = load_mass_spectrum_from_d_folder(str(d_path))
    peaks_df = spectrum_to_peaks_dataframe(ms_obj)

    if save_csv:
        out_dir = output_dir or CSV_RAW_DIR
        out_dir.mkdir(parents=True, exist_ok=True)

        basename = d_path.name.rstrip("/")
        out_path = out_dir / f"{basename}_peaks.csv"

        peaks_df.to_csv(out_path, index=False)
        print(f"✅ Saved raw peaks CSV: {out_path}")

    return peaks_df


def extract_many_d_folder_or_one(
    input_path: Path = DEFAULT_D_FOLDER,
    save_csv: bool = True,
    output_dir: Optional[Path] = None,
):
    """
    Higher‑level wrapper to process:
        • a single .d folder
        • a directory containing multiple .d folders

    Returns a dict mapping each processed folder → DataFrame.
    """

    if not input_path.exists():
        raise FileNotFoundError(f"Input path does not exist: {input_path}")

    # Case 1: directly a .d folder
    if input_path.is_dir() and input_path.suffix == ".d":
        return {
            input_path.name: extract_peaks_from_d_folder(
                input_path, save_csv=save_csv, output_dir=output_dir
            )
        }

    # Case 2: a directory containing *.d folders
    d_folders = sorted(
        p for p in input_path.iterdir() if p.is_dir() and p.suffix == ".d"
    )

    if not d_folders:
        print(f"⚠ No .d folders found in directory: {input_path}")
        return {}

    results = {}
    for d in d_folders:
        results[d.name] = extract_peaks_from_d_folder(
            d, save_csv=save_csv, output_dir=output_dir
        )

    return results


def main():
    """CLI entry point."""
    extract_many_d_folder_or_one(DEFAULT_D_FOLDER, save_csv=True, output_dir=CSV_RAW_DIR)


if __name__ == "__main__":
    main()