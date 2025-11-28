# SPDX-License-Identifier: EUPL-1.2
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional, Literal, Dict, Any
import json
from datetime import datetime

import pandas as pd

# Core pipeline building blocks
from schon.corems_pipeline.generate_peaks_from_fid import (
    extract_many_d_folder_or_one,
)
from schon.calibration.calibration import (
    run_calibration,
    run_formula_assignment,
)
from schon.calibration.calibration_config import CalibrationConfig


InputKind = Literal["d", "uncalibrated_csv", "calibrated_csv"]


def _ensure_index_column(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure that a stable 'Index' column exists.

    If 'Index' does not exist, create it once from the current row order.
    If it does exist, leave it alone and only reset the *pandas* index
    (so later operations that call reset_index(drop=True) won't break).
    """
    df = df.copy()
    if "Index" not in df.columns:
        df = df.reset_index(drop=True)
        df.insert(0, "Index", df.index)
    else:
        # Just make the underlying pandas index clean; keep 'Index' column
        df = df.reset_index(drop=True)
    return df


@dataclass
class SampleRunConfig:
    input_kind: InputKind = "d"
    sample_type: Optional[str] = None
    calib_cfg: CalibrationConfig = field(default_factory=CalibrationConfig)

    formula_ppm_tolerance: Optional[float] = None
    formula_n_processes: Optional[int] = None

    data_root: Path = Path("/app/data")
    peaks_output_dir: Path = Path("/app/results/csv_raw")
    calibration_results_dir: Path = Path("/app/results/calibration")
    runs_dir: Path = Path("/app/results/runs")

    save_raw_peaks_csv: bool = True

    # --- NEW: serialization helpers ---

    def to_dict(self) -> dict:
        return {
            "input_kind": self.input_kind,
            "sample_type": self.sample_type,
            "calib_cfg": asdict(self.calib_cfg),
            "formula_ppm_tolerance": self.formula_ppm_tolerance,
            "formula_n_processes": self.formula_n_processes,
            "data_root": str(self.data_root),
            "peaks_output_dir": str(self.peaks_output_dir),
            "calibration_results_dir": str(self.calibration_results_dir),
            "save_raw_peaks_csv": self.save_raw_peaks_csv,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "SampleRunConfig":
        calib_cfg = CalibrationConfig(**d["calib_cfg"])

        # Ensure sample_type is always treated as a string (for presets & .lower())
        sample_type = d.get("sample_type")
        if sample_type is not None:
            sample_type = str(sample_type)

        return cls(
            input_kind=d["input_kind"],
            sample_type=sample_type,
            calib_cfg=calib_cfg,
            formula_ppm_tolerance=d.get("formula_ppm_tolerance"),
            formula_n_processes=d.get("formula_n_processes"),
            data_root=Path(d["data_root"]),
            peaks_output_dir=Path(d["peaks_output_dir"]),
            calibration_results_dir=Path(d["calibration_results_dir"]),
            save_raw_peaks_csv=d.get("save_raw_peaks_csv", True),
        )


@dataclass
class FormulaRun:
    """
    Record of a single formula-assignment run on (typically) calibrated peaks.
    """
    sample_type: Optional[str]
    ppm_tolerance: Optional[float]
    n_processes: Optional[int]
    df: pd.DataFrame


@dataclass
class SampleRun:
    """
    One logical run for a given sample.

    A "run" can start from a .d folder or from a CSV (calibrated or not),
    be calibrated with a chosen CalibrationConfig, and then have *several*
    formula-assignment passes with different sample_types / ppm tolerances.

    Attributes
    ----------
    sample_name :
        Human-readable label for the sample (usually derived from the file/folder name).
    input_path :
        Path to the .d folder or CSV used as starting point.
    cfg :
        SampleRunConfig controlling calibration + formula assignment.
    """

    sample_name: str
    input_path: Path
    cfg: SampleRunConfig

    # Dataframes for different stages
    raw_peaks_df: Optional[pd.DataFrame] = None
    calibrated_df: Optional[pd.DataFrame] = None

    # Multiple formula runs (e.g. different filters / tolerances)
    formula_runs: list[FormulaRun] = field(default_factory=list)

    # Diagnostics from calibration
    calibration_debug: Dict[str, Any] = field(default_factory=dict)

    # If the CSV input is already calibrated
    calibrated_input: bool = False

    # --- NEW: per-run identity / metadata ---
    run_id: str = field(
        default_factory=lambda: datetime.utcnow().strftime("%Y%m%d_%H%M%S_%f")
    )
    created_at: str = field(
        default_factory=lambda: datetime.utcnow().isoformat()
    )

    # Where to store metadata for this run (defaults into calibration_results_dir/runs/)
    def _runs_root(self) -> Path:
        return self.cfg.runs_dir / self.sample_name / self.run_id

    def to_metadata_dict(self) -> dict:
        """Minimal JSON-serializable snapshot of this run."""
        return {
            "sample_name": self.sample_name,
            "input_path": str(self.input_path),
            "run_id": self.run_id,
            "created_at": self.created_at,
            "calibrated_input": self.calibrated_input,
            "cfg": self.cfg.to_dict(),
            # You can add more here later: e.g. paths to result CSVs, notes, etc.
        }

    def save_metadata(self) -> Path:
        root = self._runs_root()
        root.mkdir(parents=True, exist_ok=True)
        meta = self.to_metadata_dict()

        # 1) Save per-run config
        meta_path = root / "run_config.json"
        with meta_path.open("w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2)

        # 2) Update global index of all runs
        index_path = self.cfg.runs_dir / "runs_index.json"

        try:
            if index_path.exists():
                with index_path.open("r", encoding="utf-8") as f:
                    index_data = json.load(f)
                # keep it as a list of entries
                if not isinstance(index_data, list):
                    index_data = []
            else:
                index_data = []
        except Exception:
            # if something is wrong, start a new index
            index_data = []

        # Minimal summary entry for this run
        summary = {
            "sample_name": self.sample_name,
            "run_id": self.run_id,
            "created_at": self.created_at,
            "sample_type": self.cfg.sample_type,
            "input_path": str(self.input_path),
            "calibrated_input": self.calibrated_input,
            # optional: where main outputs are
            "outputs": {
                "calibration_dir": str(self.cfg.calibration_results_dir),
                # you can add more paths here later if you want
            },
        }

        index_data.append(summary)

        with index_path.open("w", encoding="utf-8") as f:
            json.dump(index_data, f, indent=2)

        return meta_path
        
    @classmethod
    def from_metadata(cls, meta_path: Path) -> "SampleRun":
        with meta_path.open("r", encoding="utf-8") as f:
            data = json.load(f)
        cfg = SampleRunConfig.from_dict(data["cfg"])
        return cls(
            sample_name=data["sample_name"],
            input_path=Path(data["input_path"]),
            cfg=cfg,
            run_id=data["run_id"],
            created_at=data["created_at"],
            calibrated_input=data.get("calibrated_input", False),
        )
    # ------------------ Stage 1: load / extract peaks ------------------

    def load_input(self) -> None:
        """
        Load or extract peaks depending on cfg.input_kind.

        - 'd'               → run CoreMS peak extraction on the .d folder
        - 'uncalibrated_csv'→ read CSV as raw peaks
        - 'calibrated_csv'  → read CSV as *already calibrated* peaks
        """
        kind = self.cfg.input_kind
        path = self.input_path

        if kind == "d":
            # For the OOP object we assume: one run = one .d folder.
            result = extract_many_d_folder_or_one(
                input_path=path,
                save_csv=self.cfg.save_raw_peaks_csv,
                output_dir=self.cfg.peaks_output_dir,
            )
            if len(result) != 1:
                raise RuntimeError(
                    f"Expected exactly one .d folder for SampleRun, got {len(result)}"
                )
            # result keys can be strings or Paths; we just take the single df
            (_, df), = result.items()
            df = _ensure_index_column(df)
            self.raw_peaks_df = df

        elif kind in ("uncalibrated_csv", "calibrated_csv"):
            df = pd.read_csv(path)
            df = _ensure_index_column(df)
            if kind == "calibrated_csv":
                # We treat this as already calibrated: no calibration step.
                self.calibrated_input = True
                self.calibrated_df = df
            else:
                self.raw_peaks_df = df
        else:
            raise ValueError(f"Unsupported input_kind: {kind}")

    # ------------------ Stage 2: calibration ------------------

    def run_calibration(self) -> None:
        """
        Run internal calibration on raw_peaks_df (unless input is already calibrated).

        Uses the existing `schon.calibration.calibration.run_calibration`.
        The function:
        - Works on DataFrames.
        - Writes CSV/PKL side-effects into cfg.calibration_results_dir.
        - Returns (full_ms_like_df, debug_info).

        We keep the returned DataFrame in `self.calibrated_df` and diagnostics
        in `self.calibration_debug`.
        """
        if self.calibrated_input:
            # Input was already calibrated; nothing to do here.
            if self.calibrated_df is None:
                raise RuntimeError(
                    "calibrated_input=True but calibrated_df is None; "
                    "did you forget to call load_input()?"
                )
            return

        if self.raw_peaks_df is None:
            raise RuntimeError("No raw_peaks_df; call load_input() first.")

        full_ms_like, debug_info = run_calibration(
            peaks=self.raw_peaks_df,
            cfg=self.cfg.calib_cfg,
            results_dir=self.cfg.calibration_results_dir,
            sample_type=self.cfg.sample_type,
            sample_name=self.sample_name,
        )

        # Ensure 'Index' is stable and present
        full_ms_like = _ensure_index_column(full_ms_like)

        self.calibrated_df = full_ms_like
        self.calibration_debug = debug_info

    # ------------------ Stage 3: formula assignment on calibrated data ------------------

    def assign_formulas_on_calibrated(
        self,
        sample_type: Optional[str] = None,
        ppm_tolerance: Optional[float] = None,
        n_processes: Optional[int] = None,
    ) -> pd.DataFrame:
        """
        Run formula assignment on the *calibrated* peaks and record this run.
        Returns a DataFrame with canonical columns.
        """
        if self.calibrated_df is None:
            raise RuntimeError("No calibrated_df; call run_calibration() first.")

        # Effective parameters for this run
        eff_sample_type = sample_type if sample_type is not None else self.cfg.sample_type
        eff_ppm = ppm_tolerance if ppm_tolerance is not None else self.cfg.formula_ppm_tolerance
        eff_nproc = n_processes if n_processes is not None else self.cfg.formula_n_processes

        # Prepare a copy for formula assignment:
        df_for_assignment = self.calibrated_df.copy()
        # Keep original m/z in "m/z_raw" if not present
        if "m/z_raw" not in df_for_assignment.columns:
            if "m/z" in df_for_assignment.columns:
                df_for_assignment["m/z_raw"] = df_for_assignment["m/z"]
            elif "m/z_raw" in df_for_assignment.columns:
                pass
        # Use "Calibrated m/z" as "m/z" for matching
        if "Calibrated m/z" in df_for_assignment.columns:
            df_for_assignment["m/z"] = df_for_assignment["Calibrated m/z"]

        # Ensure intensity column exists
        if "Intensity" not in df_for_assignment.columns and "Peak Height" in df_for_assignment.columns:
            df_for_assignment["Intensity"] = df_for_assignment["Peak Height"]

        # Ensure the Index column is present exactly once
        df_for_assignment = _ensure_index_column(df_for_assignment)

        # Use the same in-memory wrapper as calibration uses
        formulas_df = run_formula_assignment(
            peaks_for_formulas_df=df_for_assignment,
            sample_type=eff_sample_type,
            ppm_tolerance=eff_ppm,
            n_processes=eff_nproc,
        )

        # --- Build canonical output directly from formulas_df ---

        # Мы предполагаем, что assign_formulas_df возвращает тот же DataFrame,
        # который мы ему передали (df_for_assignment), плюс новые колонки
        # (Formula, Calculated Mass, Mass Error (ppm), C, H, ...).
        formulas_df = _ensure_index_column(formulas_df)

        # Строковое представление типа образца (чтобы совпадало с CLI-логикой)
        sample_type_str = eff_sample_type if eff_sample_type is not None else ""

        # Векторно посчитаем Calibration Error (ppm), если есть сырое m/z и калиброванное
        if "Calibrated m/z" in formulas_df.columns and "m/z_raw" in formulas_df.columns:
            mz_raw = formulas_df["m/z_raw"].astype(float)
            mz_cal = formulas_df["Calibrated m/z"].astype(float)
            with pd.option_context("mode.use_inf_as_na", True):
                calib_err_ppm = (mz_cal - mz_raw) / mz_raw * 1e6
        else:
            calib_err_ppm = pd.Series([pd.NA] * len(formulas_df), index=formulas_df.index)

        # Колонка "Used For Calibration" могла быть проставлена калибровкой.
        if "Used For Calibration" in formulas_df.columns:
            used_for_cal = formulas_df["Used For Calibration"].astype(bool)
        else:
            used_for_cal = pd.Series([False] * len(formulas_df), index=formulas_df.index)

        # Собираем DataFrame в "CLI-совместимом" формате:
        canonical_cols = [
            "Index", "Calibrated m/z", "Intensity", "S/N", "Ion Charge",
            "Formula", "isotopolog", "Calculated Mass", "Mass Error (ppm)",
            "C", "H", "N", "O", "S", "*C", "Na",
            "Alternative Formula", "Alternative Mass Error (ppm)",
            "Used For Calibration",
            "Sample type", "m/z_raw", "Calibration Error (ppm)"
        ]

        out = pd.DataFrame(index=formulas_df.index)

        # Index
        out["Index"] = formulas_df["Index"].astype("Int64")

        # М/z калиброванный
        out["Calibrated m/z"] = (
            formulas_df["Calibrated m/z"] if "Calibrated m/z" in formulas_df.columns else pd.NA
        )

        # Интенсивность
        if "Intensity" in formulas_df.columns:
            out["Intensity"] = formulas_df["Intensity"]
        elif "Peak Height" in formulas_df.columns:
            out["Intensity"] = formulas_df["Peak Height"]
        else:
            out["Intensity"] = pd.NA

        # S/N
        out["S/N"] = formulas_df["S/N"] if "S/N" in formulas_df.columns else pd.NA

        # Заряд
        if "Ion Charge" in formulas_df.columns:
            out["Ion Charge"] = formulas_df["Ion Charge"]
        else:
            out["Ion Charge"] = -1

        # Формулы и атрибуты
        for col in [
            "Formula",
            "isotopolog",
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
        ]:
            out[col] = formulas_df[col] if col in formulas_df.columns else pd.NA

        # Тип образца
        out["Sample type"] = sample_type_str

        # Сырое m/z
        if "m/z_raw" in formulas_df.columns:
            out["m/z_raw"] = formulas_df["m/z_raw"]
        else:
            # на всякий случай можно использовать исходную колонку "m/z_cal" / "m/z"
            if "m/z" in formulas_df.columns:
                out["m/z_raw"] = formulas_df["m/z"]
            else:
                out["m/z_raw"] = pd.NA

        # Ошибка калибровки
        out["Calibration Error (ppm)"] = calib_err_ppm

        # Использован ли пик для калибровки
        out["Used For Calibration"] = used_for_cal

        # Гарантируем порядок колонок
        out = out[canonical_cols]

        # Сохраняем run в память
        self.formula_runs.append(
            FormulaRun(
                sample_type=eff_sample_type,
                ppm_tolerance=eff_ppm,
                n_processes=eff_nproc,
                df=out,
            )
        )

        # Пишем на диск
        out_dir = self.cfg.runs_dir / self.sample_name / self.run_id
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{self.sample_name}_calibrated_with_formulas_sampletype_{sample_type_str}.csv"
        out.to_csv(out_path, index=False)
        print(f"✓ Assigned formulas to calibrated peaks → {out_path}")

        return out

    # ------------------ Convenience ------------------

    @property
    def latest_formulas_df(self) -> Optional[pd.DataFrame]:
        """Return the DataFrame from the last formula run, if any."""
        if not self.formula_runs:
            return None
        return self.formula_runs[-1].df

    @property
    def current_df(self) -> pd.DataFrame:
        """
        Return the "most advanced" DataFrame:

        1. if formulas were assigned → last formulas df
        2. else if calibration done   → calibrated_df
        3. else                       → raw_peaks_df

        Raises if nothing was loaded.
        """
        if self.latest_formulas_df is not None:
            return self.latest_formulas_df
        if self.calibrated_df is not None:
            return self.calibrated_df
        if self.raw_peaks_df is not None:
            return self.raw_peaks_df
        raise RuntimeError("No data loaded; call load_input() first.")