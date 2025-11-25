from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

def _pick_intensity_column(df: pd.DataFrame) -> str:
    """Choose a reasonable intensity column name from the dataframe."""
    if "Intensity" in df.columns:
        return "Intensity"
    if "Peak Height" in df.columns:
        return "Peak Height"
    # fallback: any numeric column that's not m/z-like
    numeric_cols = df.select_dtypes("number").columns.tolist()
    for col in ["Calibrated m/z", "m/z", "m/z_raw"]:
        if col in numeric_cols:
            numeric_cols.remove(col)
    return numeric_cols[0] if numeric_cols else "Intensity"


def _prepare_mz_and_intensity(
    df: pd.DataFrame,
    mz_col_preference: Tuple[str, ...] = ("Calibrated m/z", "m/z", "m/z_raw"),
) -> Tuple[np.ndarray, np.ndarray, str, str]:
    """
    Return (mz, intensity, mz_col, intensity_col) from df.

    mz_col_preference is the ordered list of columns we try in df for m/z.
    """
    mz_col = None
    for col in mz_col_preference:
        if col in df.columns:
            mz_col = col
            break
    if mz_col is None:
        raise ValueError("No m/z column found in dataframe.")

    intensity_col = _pick_intensity_column(df)

    mz = pd.to_numeric(df[mz_col], errors="coerce").to_numpy(dtype=float)
    intensity = pd.to_numeric(df[intensity_col], errors="coerce").to_numpy(dtype=float)

    return mz, intensity, mz_col, intensity_col


def _select_labeled_indices(
    df: pd.DataFrame,
    intensity_col: str,
    base_mask: Optional[pd.Series] = None,
    max_labels: int = 30,
) -> pd.Index:
    """
    Choose which peaks to label (for formulas) based on intensity.

    As the number of peaks in view decreases, the max labels effectively
    increases (because we always pick the top N by intensity).
    """
    if base_mask is None:
        mask = pd.Series(True, index=df.index)
    else:
        mask = base_mask

    if mask.sum() == 0:
        return df.index[0:0]  # empty index

    df_sub = df.loc[mask, [intensity_col]].copy()
    df_sub = df_sub.sort_values(intensity_col, ascending=False)

    n = len(df_sub)
    # Make max_labels scale with number of peaks in view
    dynamic_max = min(max_labels, max(5, n // 20))
    dynamic_max = max(5, dynamic_max)  # at least a few

    return df_sub.head(dynamic_max).index


# -------------------------------------------------------------------
# 1) Spectrum with calibrants and formulas
# -------------------------------------------------------------------

def plot_spectrum_with_calibrants(
    df: pd.DataFrame,
    used_for_calib_col: str = "Used For Calibration",
    formula_col: str = "Formula",
    output_path: Optional[Path] = None,
) -> plt.Figure:
    """
    Plot spectrum with peaks as grey sticks and calibrant peaks highlighted.

    - All peaks: grey vertical lines
    - Calibrant peaks: red points on top of peaks
    - Formulas for calibrant peaks: vertical labels (subsampled by intensity)
    - If 'series_length' is present and >1, marker size scales with series length.
    """

    fig, ax = plt.subplots(figsize=(10, 4))

    if df.empty:
        ax.text(0.5, 0.5, "No data to plot", ha="center", va="center")
        ax.set_axis_off()
        if output_path:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_path, dpi=200, bbox_inches="tight")
        return fig

    mz, intensity, mz_col, intensity_col = _prepare_mz_and_intensity(df)

    # All peaks: grey sticks
    ax.vlines(mz, 0, intensity, linewidth=0.6, alpha=0.3)

    # Calibrants
    if used_for_calib_col in df.columns:
        mask_calib = df[used_for_calib_col].fillna(False).astype(bool)
    else:
        mask_calib = pd.Series(False, index=df.index)

    mz_calib = mz[mask_calib.to_numpy()]
    int_calib = intensity[mask_calib.to_numpy()]

    # Marker size: optionally depend on series_length
    if "series_length" in df.columns:
        series_len = df["series_length"].fillna(0).to_numpy()
        series_len_calib = series_len[mask_calib.to_numpy()]
        sizes = 10 + 5 * np.clip(series_len_calib, 0, 10)
    else:
        sizes = np.full_like(int_calib, 20.0, dtype=float)

    if len(mz_calib) > 0:
        ax.scatter(
            mz_calib,
            int_calib,
            s=sizes,
            c="red",
            marker="o",
            alpha=0.9,
            zorder=3,
            label="Calibrants",
        )

    # Label calibrant formulas (subsampled by intensity)
    if formula_col in df.columns and mask_calib.any():
        label_indices = _select_labeled_indices(
            df, intensity_col=intensity_col, base_mask=mask_calib
        )
        for idx in label_indices:
            row = df.loc[idx]
            f = str(row[formula_col]) if pd.notna(row[formula_col]) else ""
            if not f.strip():
                continue
            x = float(row[mz_col])
            y = float(row[intensity_col])
            ax.text(
                x,
                y * 1.05,
                f,
                rotation=90,
                ha="center",
                va="bottom",
                fontsize=7,
                alpha=0.9,
            )

    ax.set_xlabel(mz_col)
    ax.set_ylabel(intensity_col)
    ax.set_title("Spectrum with calibrants (red) and formulas")
    ax.set_ylim(bottom=0)

    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

    return fig


# -------------------------------------------------------------------
# 2) Raw vs calibrated spectrum
# -------------------------------------------------------------------

def plot_raw_vs_calibrated(
    df: pd.DataFrame,
    output_path: Optional[Path] = None,
) -> plt.Figure:
    """
    Plot original m/z in red and calibrated m/z in black on top.

    Requires:
    - 'm/z_raw' for original masses
    - 'Calibrated m/z' for calibrated masses
    """
    fig, ax = plt.subplots(figsize=(10, 4))

    if df.empty:
        ax.text(0.5, 0.5, "No data to plot", ha="center", va="center")
        ax.set_axis_off()
        if output_path:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_path, dpi=200, bbox_inches="tight")
        return fig

    if "m/z_raw" not in df.columns or "Calibrated m/z" not in df.columns:
        ax.text(
            0.5,
            0.5,
            "m/z_raw or Calibrated m/z missing",
            ha="center",
            va="center",
        )
        ax.set_axis_off()
        if output_path:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_path, dpi=200, bbox_inches="tight")
        return fig

    intensity_col = _pick_intensity_column(df)

    mz_raw = pd.to_numeric(df["m/z_raw"], errors="coerce").to_numpy(dtype=float)
    mz_cal = pd.to_numeric(df["Calibrated m/z"], errors="coerce").to_numpy(dtype=float)
    intensity = pd.to_numeric(df[intensity_col], errors="coerce").to_numpy(dtype=float)

    # Original spectrum in red, slightly thinner & more transparent
    ax.vlines(mz_raw, 0, intensity, linewidth=0.5, alpha=0.4, color="red", label="Raw")

    # Calibrated in black
    ax.vlines(
        mz_cal,
        0,
        intensity,
        linewidth=0.8,
        alpha=0.8,
        color="black",
        label="Calibrated",
    )

    ax.set_xlabel("m/z")
    ax.set_ylabel(intensity_col)
    ax.set_title("Raw (red) vs calibrated (black) spectrum")
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper right", fontsize=8)

    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

    return fig


# -------------------------------------------------------------------
# 3) Spectrum with assigned formulas after formula assignment
# -------------------------------------------------------------------

def plot_spectrum_with_assigned_formulas(
    df: pd.DataFrame,
    formula_col: str = "Formula",
    used_for_calib_col: str = "Used For Calibration",
    output_path: Optional[Path] = None,
) -> plt.Figure:
    """
    Plot spectrum after formula assignment:

    - All peaks: grey sticks
    - Peaks with assigned formulas: orange dots on top
    - Formulas: vertical labels for strongest assigned peaks in view
    - Calibrant peaks (if Used For Calibration is present): darker border
      or larger markers (if we want to distinguish them)
    """
    fig, ax = plt.subplots(figsize=(10, 4))

    if df.empty:
        ax.text(0.5, 0.5, "No data to plot", ha="center", va="center")
        ax.set_axis_off()
        if output_path:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_path, dpi=200, bbox_inches="tight")
        return fig

    mz, intensity, mz_col, intensity_col = _prepare_mz_and_intensity(df)

    # All peaks
    ax.vlines(mz, 0, intensity, linewidth=0.6, alpha=0.3)

    # Assigned peaks = non-empty formulas
    if formula_col in df.columns:
        mask_assigned = df[formula_col].astype(str).str.strip() != ""
    else:
        mask_assigned = pd.Series(False, index=df.index)

    mz_assigned = mz[mask_assigned.to_numpy()]
    int_assigned = intensity[mask_assigned.to_numpy()]

    if len(mz_assigned) > 0:
        ax.scatter(
            mz_assigned,
            int_assigned,
            s=15,
            c="orange",
            marker="o",
            alpha=0.9,
            zorder=3,
            label="Assigned formulas",
        )

    # Label strongest assigned peaks
    if mask_assigned.any() and formula_col in df.columns:
        label_indices = _select_labeled_indices(
            df, intensity_col=intensity_col, base_mask=mask_assigned
        )
        for idx in label_indices:
            row = df.loc[idx]
            f = str(row[formula_col]) if pd.notna(row[formula_col]) else ""
            if not f.strip():
                continue
            x = float(row[mz_col])
            y = float(row[intensity_col])
            ax.text(
                x,
                y * 1.05,
                f,
                rotation=90,
                ha="center",
                va="bottom",
                fontsize=7,
                alpha=0.9,
            )

    # Optionally highlight calibrants inside assigned peaks
    if used_for_calib_col in df.columns:
        mask_calib = df[used_for_calib_col].fillna(False).astype(bool) & mask_assigned
        mz_calib = mz[mask_calib.to_numpy()]
        int_calib = intensity[mask_calib.to_numpy()]
        if len(mz_calib) > 0:
            ax.scatter(
                mz_calib,
                int_calib,
                s=30,
                facecolors="none",
                edgecolors="red",
                linewidths=1.0,
                alpha=0.9,
                zorder=4,
                label="Calibrant (series) members",
            )

    ax.set_xlabel(mz_col)
    ax.set_ylabel(intensity_col)
    ax.set_title("Spectrum with assigned formulas (orange) and calibrant series")
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper right", fontsize=8)

    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

    return fig