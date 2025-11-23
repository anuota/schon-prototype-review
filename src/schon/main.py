from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Literal

import pandas as pd

# --- CoreMS peak extraction ---
from schon.corems_pipeline.generate_peaks_from_fid import (
    extract_many_d_folder_or_one,  # assumes your helper returns {d_path: df}
    DEFAULT_D_FOLDER,
)

from schon.calibration.calibration import run_calibration, run_formula_assignment
from schon.calibration.calibration_config import CalibrationConfig

# --------- Default paths (Docker-friendly) ---------

# Root folder where .d data lives inside the container
DEFAULT_DATA_ROOT = Path("/app/data")

# Default folder where raw peak CSVs from CoreMS go
DEFAULT_PEAKS_OUTPUT_DIR = Path("/app/results/csv_raw")

# Default folder where calibration results go
DEFAULT_CALIBRATION_RESULTS_DIR = Path("/app/results/calibration")


# --------- High-level pipeline config ---------

@dataclass
class PipelineConfig:
    """
    Configuration for the full SCHON pipeline.

    A future GUI should create one of these and pass it to `run_pipeline`.
    """

    # Where the .d folders live (mainly for convenience / defaults)
    data_root: Path = DEFAULT_DATA_ROOT

    # Where to put CoreMS peak CSVs (optional debug / archive)
    peaks_output_dir: Path = DEFAULT_PEAKS_OUTPUT_DIR

    # Where to put calibration outputs (debug tables, calibrated peaks, etc.)
    calibration_results_dir: Path = DEFAULT_CALIBRATION_RESULTS_DIR

    # Calibration configuration (S/N threshold, polynomial degree, etc.)
    calibration: CalibrationConfig = CalibrationConfig()

    # VERY IMPORTANT: sample type controls which element presets / rules
    # are used during formula assignment inside calibration
    # (e.g. 0 = CHO, 1 = CHNO, etc., whatever your formula_filters expects).
    sample_type: Optional[int] = None

    # Whether to save the raw peak CSVs from CoreMS
    save_raw_peaks_csv: bool = True

    # Whether to save calibrated peaks as CSV / PKL
    save_calibrated_csv: bool = True
    save_calibrated_pkl: bool = False


def assign_formulas_to_calibrated(
    calibrated_df: pd.DataFrame,
    cfg: PipelineConfig,
    sample_name: str,
) -> pd.DataFrame:
    """
    Assign molecular formulas to calibrated peaks and save the result as CSV.

    This uses the *calibrated* m/z and peak height to build a minimal
    two-column DataFrame (`m/z`, `Intensity`) that is passed to
    `run_formula_assignment`. The full formulas table is returned and also
    written as:

        <calibration_results_dir>/<sample_name>_calibrated_with_formulas.csv
    """
    # Prepare calibrated peaks for formula assignment: use calibrated m/z + intensity
    peaks_for_formulas_df = pd.DataFrame(
        {
            "m/z": calibrated_df["Calibrated m/z"],
            "Intensity": calibrated_df["Peak Height"],
        }
    )

    # Use defaults for ppm_tolerance and n_processes unless configured otherwise.
    # (CalibrationConfig does not carry a ppm_tolerance for formula assignment,
    # so we simply let the formula module use its own default.)
    calibrated_formulas_df = run_formula_assignment(
        peaks_for_formulas_df=peaks_for_formulas_df,
        sample_type=cfg.sample_type,
        ppm_tolerance=None,
        n_processes=None,
    )

    out_path = cfg.calibration_results_dir / f"{sample_name}_calibrated_with_formulas.csv"
    cfg.calibration_results_dir.mkdir(parents=True, exist_ok=True)
    calibrated_formulas_df.to_csv(out_path, index=False)
    print(f"✓ Assigned formulas to calibrated peaks → {out_path}")

    return calibrated_formulas_df

# ========= Step 1: peak extraction (CoreMS) =========

def run_peak_extraction_from_d(
    input_path: Path,
    cfg: PipelineConfig,
) -> Dict[Path, pd.DataFrame]:
    """
    Run CoreMS peak extraction on a .d folder or a directory of .d folders.

    Parameters
    ----------
    input_path :
        Either a single `.d` folder or a directory containing many `.d` folders.
    cfg :
        PipelineConfig with output directories and flags.

    Returns
    -------
    dict
        Mapping `{d_folder_path: peaks_df}` for each processed .d.
    """
    peaks_output_dir = cfg.peaks_output_dir
    peaks_output_dir.mkdir(parents=True, exist_ok=True)

    # Helper is expected to return {d_path: df}, but the keys may be strings.
    raw_d_to_df = extract_many_d_folder_or_one(
        input_path=input_path,
        save_csv=cfg.save_raw_peaks_csv,
        output_dir=peaks_output_dir,
    )

    # Normalize keys to Path objects so downstream code can reliably use .name
    d_to_df: Dict[Path, pd.DataFrame] = {
        Path(d_path): df for d_path, df in raw_d_to_df.items()
    }

    # If save_raw_peaks_csv is False, we can still *optionally* create CSVs here
    # (e.g. for debugging). For now we trust extract_many_d_folder_or_one to
    # handle CSV writing when requested.
    return d_to_df


# ========= Step 2: calibration (df-based) =========

def run_calibration_for_d_dict(
    d_to_peaks_df: Dict[Path, pd.DataFrame],
    cfg: PipelineConfig,
) -> None:
    """
    Run internal calibration on each peaks DataFrame.

    This assumes `run_calibration` works on DataFrames and does NOT
    perform its own CSV reading.

    Parameters
    ----------
    d_to_peaks_df :
        Mapping from `.d` folder paths to their corresponding peaks DataFrames.
    cfg :
        PipelineConfig containing CalibrationConfig and output dirs.
    """
    results_root = cfg.calibration_results_dir
    results_root.mkdir(parents=True, exist_ok=True)

    for d_path, peaks_df in d_to_peaks_df.items():
        if isinstance(d_path, Path):
            sample_name = d_path.name.rstrip("/")
        else:
            sample_name = str(d_path).rstrip("/")
        print(f"\n=== Calibrating {sample_name} ===")

        calibrated_df, debug_info = run_calibration(
            peaks=peaks_df,
            cfg=cfg.calibration,
            results_dir=results_root,
            sample_type=cfg.sample_type,
            sample_name=sample_name,
        )

        # --- Additional Step: Formula assignment on calibrated data ---
        assign_formulas_to_calibrated(
            calibrated_df=calibrated_df,
            cfg=cfg,
            sample_name=sample_name,
        )


def run_calibration_for_single_csv(
    csv_path: Path,
    cfg: PipelineConfig,
    calibrated_input: bool = False,
) -> None:
    """
    Pipeline when the input is a CSV file (either uncalibrated or calibrated).

    Parameters
    ----------
    csv_path :
        Path to CSV file. For uncalibrated CSV, it must contain columns
        expected by the calibration step (m/z, intensity, S/N).
    cfg :
        PipelineConfig.
    calibrated_input :
        If True, treat this CSV as *already calibrated* and skip the
        calibration step (placeholder for future processing).
    """
    csv_path = csv_path.resolve()
    print(f"\n▶ Loading CSV: {csv_path}")
    df = pd.read_csv(csv_path)

    sample_name = csv_path.stem
    results_root = cfg.calibration_results_dir
    results_root.mkdir(parents=True, exist_ok=True)

    if calibrated_input:
        # For now we simply (optionally) copy it into the results folder.
        # Later you can plug in other steps (formula assignment only,
        # visualization, Kendrick filters, etc.).
        print("Input CSV is marked as *already calibrated*; skipping calibration.")
        if cfg.save_calibrated_csv or cfg.save_calibrated_pkl:
            out_base = results_root / f"{sample_name}_calibrated_input"
            if cfg.save_calibrated_csv:
                df.to_csv(out_base.with_suffix(".csv"), index=False)
            if cfg.save_calibrated_pkl:
                df.to_pickle(out_base.with_suffix(".pkl"))
        return

    # Uncalibrated CSV: run the same df-based calibration as for .d inputs
    print(f"=== Calibrating CSV sample {sample_name} ===")
    calibrated_df, debug_info = run_calibration(
        peaks=df,
        cfg=cfg.calibration,
        results_dir=results_root,
        sample_type=cfg.sample_type,
        sample_name=sample_name,
    )

    # --- Additional Step: Formula assignment on calibrated data ---
    assign_formulas_to_calibrated(
        calibrated_df=calibrated_df,
        cfg=cfg,
        sample_name=sample_name,
    )


# ========= Full pipelines & dispatcher =========

InputKind = Literal["auto", "d", "uncalibrated_csv", "calibrated_csv"]


def run_pipeline(
    input_path: Optional[Path] = None,
    pipeline_cfg: Optional[PipelineConfig] = None,
    input_kind: InputKind = "auto",
) -> None:
    """
    Run the SCHON pipeline, dispatching depending on the input type.

    Parameters
    ----------
    input_path :
        - `.d` folder, or directory containing `.d` folders
        - CSV of uncalibrated peaks
        - CSV of already calibrated peaks
        If None, `DEFAULT_D_FOLDER` from CoreMS pipeline is used.
    pipeline_cfg :
        PipelineConfig instance. If None, default values are used.
    input_kind :
        One of:
        - "auto"            : guess from extension & existence
        - "d"               : treat input as .d (or folder of .d)
        - "uncalibrated_csv": treat input as raw peak CSV
        - "calibrated_csv"  : treat input as already calibrated CSV
    """
    if pipeline_cfg is None:
        pipeline_cfg = PipelineConfig()

    if input_path is None:
        # fallback to whatever default you use for testing
        input_path = DEFAULT_D_FOLDER

    input_path = input_path.resolve()
    print(f"▶ Running SCHON pipeline on: {input_path}")

    # --- auto-detect input kind if requested ---
    if input_kind == "auto":
        if input_path.is_dir() and input_path.suffix.lower() == ".d":
            kind: InputKind = "d"
        elif input_path.is_dir():
            # Directory: assume it contains .d folders
            kind = "d"
        elif input_path.suffix.lower() == ".csv":
            # We cannot automatically know if it's calibrated;
            # default: treat as *uncalibrated* unless user says otherwise.
            kind = "uncalibrated_csv"
        else:
            raise ValueError(
                f"Cannot auto-detect input kind for: {input_path}. "
                "Please specify --input-kind explicitly."
            )
    else:
        kind = input_kind

    # --- dispatch to appropriate pipeline ---
    if kind == "d":
        # Step 1: CoreMS peak extraction
        d_to_peaks_df = run_peak_extraction_from_d(input_path, pipeline_cfg)
        # Step 2: calibration (includes formula assignment internally)
        run_calibration_for_d_dict(d_to_peaks_df, pipeline_cfg)

    elif kind == "uncalibrated_csv":
        run_calibration_for_single_csv(
            csv_path=input_path,
            cfg=pipeline_cfg,
            calibrated_input=False,
        )

    elif kind == "calibrated_csv":
        run_calibration_for_single_csv(
            csv_path=input_path,
            cfg=pipeline_cfg,
            calibrated_input=True,
        )

    else:
        raise ValueError(f"Unsupported input_kind: {kind}")

    print("\n✅ SCHON pipeline finished.")


# ========= CLI entrypoint =========

def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="SCHON pipeline: CoreMS peak extraction + calibration."
    )
    parser.add_argument(
        "--input",
        type=str,
        help="Path to a .d folder, a directory containing .d folders, "
             "or a CSV file (uncalibrated or calibrated). "
             "Default: use DEFAULT_D_FOLDER from CoreMS pipeline.",
        default=None,
    )
    parser.add_argument(
        "--input-kind",
        type=str,
        choices=["auto", "d", "uncalibrated_csv", "calibrated_csv"],
        default="auto",
        help=(
            "How to interpret the input path. "
            "'auto' guesses by extension; "
            "'d' for .d folders; "
            "'uncalibrated_csv' for raw peaks; "
            "'calibrated_csv' for already calibrated peaks."
        ),
    )
    parser.add_argument(
        "--snr-min",
        type=float,
        help="Minimum S/N for peaks to be considered for calibration.",
        default=None,
    )
    parser.add_argument(
        "--poly-order",
        type=int,
        help="Polynomial order for m/z calibration (0=offset, 1=linear, 2=quadratic, ...).",
        default=None,
    )
    parser.add_argument(
        "--peaks-out",
        type=str,
        help=f"Directory to store raw peaks CSVs (default: {DEFAULT_PEAKS_OUTPUT_DIR})",
        default=None,
    )
    parser.add_argument(
        "--calib-out",
        type=str,
        help=f"Directory to store calibration results (default: {DEFAULT_CALIBRATION_RESULTS_DIR})",
        default=None,
    )
    parser.add_argument(
        "--no-save-calibrated-csv",
        action="store_true",
        help="Do not save calibrated peaks as CSV.",
    )
    parser.add_argument(
        "--save-calibrated-pkl",
        action="store_true",
        help="Also save calibrated peaks as pickle (.pkl).",
    )
    parser.add_argument(
        "--sample-type",
        type=int,
        help=(
            "Integer code for sample type / formula presets "
            "(used inside formula assignment during calibration)."
        ),
        default=None,
    )
    return parser


def main() -> None:
    parser = _build_arg_parser()
    args = parser.parse_args()

    # Base calibration config with defaults
    calib_cfg = CalibrationConfig()
    if args.snr_min is not None:
        calib_cfg.sn_threshold = args.snr_min
    if args.poly_order is not None:
        calib_cfg.poly_order = args.poly_order

    peaks_out = Path(args.peaks_out) if args.peaks_out else DEFAULT_PEAKS_OUTPUT_DIR
    calib_out = Path(args.calib_out) if args.calib_out else DEFAULT_CALIBRATION_RESULTS_DIR

    pipeline_cfg = PipelineConfig(
        peaks_output_dir=peaks_out,
        calibration_results_dir=calib_out,
        calibration=calib_cfg,
        sample_type=args.sample_type,   # <— here!
        save_calibrated_csv=not args.no_save_calibrated_csv,
        save_calibrated_pkl=args.save_calibrated_pkl,
    )

    input_path: Optional[Path] = Path(args.input) if args.input else None

    run_pipeline(
        input_path=input_path,
        pipeline_cfg=pipeline_cfg,
        input_kind=args.input_kind,
    )


if __name__ == "__main__":
    main()