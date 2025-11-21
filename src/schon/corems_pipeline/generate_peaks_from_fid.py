from pathlib import Path

from schon.corems_pipeline.read_bruker_fid_corems import (
    load_mass_spectrum_from_d_folder,
    spectrum_to_peaks_dataframe,
)

PROJECT_ROOT = Path("/app")  # mounted path in Docker
SCHON_DIR = PROJECT_ROOT / "SCHON"
DATA_FOLDER = Path('/app/data')
RESULTS_FOLDER = Path('/app/results/')

# DATA_ROOT = PROJECT_ROOT / "Measurement_d-files" / "noncalibrated_d-files"
 # To analyze a single .d folder temporarily:
DATA_ROOT = DATA_FOLDER / "Measurement_d-files" / "noncalibrated_d-files" / "ESI_neg_G017736-0_100_200sc_000001.d"

# To analyze all .d folders in calibrated_d-files later, switch to:
# DATA_ROOT = PROJECT_ROOT / "Measurement_d-files" / "calibrated-d-files"
# When running in Docker with "-v /Users/anya/Coding/CoreMS/tests/tests_data:/testdata",
# the CoreMS test data will be available under /testdata inside the container.
#DATA_ROOT = Path("/testdata/ftms")
OUTPUT_PEAKS_DIR = RESULTS_FOLDER / "output_peaks1"

OUTPUT_PEAKS_DIR.mkdir(parents=True, exist_ok=True)


def process_d_folder(d_path: Path):
    print(f"🔹 Processing .d folder: {d_path}")
    ms_obj = load_mass_spectrum_from_d_folder(str(d_path))
    peaks_df = spectrum_to_peaks_dataframe(ms_obj)

    basename = d_path.name.rstrip("/")  # e.g. ESI_neg_G015548-0_100_200sc_000001.d
    out_name = f"{basename}_peaks.csv"
    out_path = OUTPUT_PEAKS_DIR / out_name

    peaks_df.to_csv(out_path, index=False)
    print(f"✅ Saved peaks to {out_path}")


def main():
    if not DATA_ROOT.exists():
        raise FileNotFoundError(f"DATA_ROOT does not exist: {DATA_ROOT}")

    # If DATA_ROOT points directly to a .d folder, wrap it into a list.
    if DATA_ROOT.is_dir() and DATA_ROOT.suffix == ".d":
        d_folders = [DATA_ROOT]
    else:
        d_folders = sorted(
            p for p in DATA_ROOT.iterdir()
            if p.is_dir() and p.suffix == ".d"
        )
    if not d_folders:
        print(f"⚠ No .d folders found in {DATA_ROOT}")
        return

    for d in d_folders:
        process_d_folder(d)


if __name__ == "__main__":
    main()