from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import plotly.graph_objs as go


def _guess_mz_and_intensity(df: pd.DataFrame) -> tuple[str, str]:
    """
    Heuristic to pick m/z and intensity columns from a SCHON dataframe.

    Priority for m/z:
        1) "Calibrated m/z"
        2) "m/z"
        3) "m/z_raw"

    Priority for intensity:
        1) "Intensity"
        2) "Peak Height"
        3) first numeric column we find
    """
    mz_candidates = ["Calibrated m/z", "m/z", "m/z_raw"]
    mz_col = None
    for c in mz_candidates:
        if c in df.columns:
            mz_col = c
            break
    if mz_col is None:
        raise ValueError("Could not find an m/z column in dataframe.")

    if "Intensity" in df.columns:
        inten_col = "Intensity"
    elif "Peak Height" in df.columns:
        inten_col = "Peak Height"
    else:
        inten_col = None
        for c in df.columns:
            if pd.api.types.is_numeric_dtype(df[c]):
                inten_col = c
                break
        if inten_col is None:
            raise ValueError("Could not find an intensity column in dataframe.")

    return mz_col, inten_col


def _base_spectrum_figure(
    df: pd.DataFrame,
    title: str,
    showlegend: bool = True,
) -> tuple[go.Figure, str, str]:
    """
    Build a base Plotly figure with the full spectrum line.

    Returns the figure and the chosen (mz_col, inten_col).
    """
    mz_col, inten_col = _guess_mz_and_intensity(df)

    x = pd.to_numeric(df[mz_col], errors="coerce")
    y = pd.to_numeric(df[inten_col], errors="coerce")
    mask = x.notna() & y.notna()
    x = x[mask]
    y = y[mask]

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="lines",
            name="Spectrum",
            line=dict(width=1),
        )
    )

    fig.update_layout(
        title=title,
        xaxis_title=mz_col,
        yaxis_title=inten_col,
        showlegend=showlegend,
        dragmode="zoom",  # rectangular zoom in both x and y
        hovermode="x unified",
    )

    # Allow independent zoom in x and y
    fig.update_xaxes(rangeslider=dict(visible=False))
    fig.update_yaxes(scaleanchor=None)

    return fig, mz_col, inten_col


def _maybe_save_figure(fig: go.Figure, output_path: Optional[Path]) -> None:
    """
    Try to save the figure as a static PNG if a path is provided.
    If kaleido or other static exporter is not available, fail silently.
    """
    if output_path is None:
        return

    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.write_image(str(output_path))
    except Exception:
        # Do not crash the pipeline if static export is not available
        pass


def plot_spectrum_with_calibrants(
    df: pd.DataFrame,
    output_path: Optional[Path] = None,
) -> go.Figure:
    """
    Plot the spectrum with calibrant peaks highlighted.

    Expects:
        - main spectrum columns (m/z, Intensity)
        - optional boolean column "Used For Calibration"
        - optional "Formula" column for label text
    """
    fig, mz_col, inten_col = _base_spectrum_figure(
        df, title="Spectrum with calibrants", showlegend=True
    )

    if "Used For Calibration" in df.columns:
        calib_df = df[df["Used For Calibration"].astype(bool)].copy()
    else:
        calib_df = pd.DataFrame([])

    if not calib_df.empty:
        x = pd.to_numeric(calib_df[mz_col], errors="coerce")
        y = pd.to_numeric(calib_df[inten_col], errors="coerce")
        mask = x.notna() & y.notna()
        x = x[mask]
        y = y[mask]
        formulas = (
            calib_df.loc[mask, "Formula"].astype(str).replace("nan", "")
            if "Formula" in calib_df.columns
            else None
        )

        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                mode="markers+text",
                name="Calibrants",
                marker=dict(size=7, color="red"),
                text=formulas,
                textposition="top center",
                textfont=dict(size=10),
                hovertemplate="<b>m/z</b>: %{x:.6f}<br>"
                "<b>Intensity</b>: %{y:.2e}<br>"
                "<b>Formula</b>: %{text}<extra></extra>",
            )
        )

    _maybe_save_figure(fig, output_path)
    return fig


def plot_raw_vs_calibrated(
    df: pd.DataFrame,
    output_path: Optional[Path] = None,
) -> go.Figure:
    """
    Plot original (raw) m/z vs calibrated m/z on the same intensity scale.

    Expects:
        - "m/z_raw" for uncalibrated m/z
        - "Calibrated m/z" for calibrated m/z
        - intensity column ("Intensity" or similar)
    """
    fig, mz_col, inten_col = _base_spectrum_figure(
        df, title="Raw vs calibrated spectrum", showlegend=True
    )

    # Clear base spectrum trace; we'll add raw + calibrated explicitly
    fig.data = []

    # Raw spectrum
    x_raw = pd.to_numeric(df["m/z_raw"], errors="coerce") if "m/z_raw" in df.columns else None
    if x_raw is not None:
        y = pd.to_numeric(df[inten_col], errors="coerce")
        mask = x_raw.notna() & y.notna()
        x_raw = x_raw[mask]
        y_raw = y[mask]
        fig.add_trace(
            go.Scatter(
                x=x_raw,
                y=y_raw,
                mode="lines",
                name="Raw m/z",
                line=dict(width=1, color="red"),
            )
        )

    # Calibrated spectrum
    x_cal = pd.to_numeric(df["Calibrated m/z"], errors="coerce") if "Calibrated m/z" in df.columns else None
    if x_cal is not None:
        y = pd.to_numeric(df[inten_col], errors="coerce")
        mask = x_cal.notna() & y.notna()
        x_cal = x_cal[mask]
        y_cal = y[mask]
        fig.add_trace(
            go.Scatter(
                x=x_cal,
                y=y_cal,
                mode="lines",
                name="Calibrated m/z",
                line=dict(width=1, color="black"),
            )
        )

    fig.update_layout(
        xaxis_title="m/z",
        yaxis_title=inten_col,
        dragmode="zoom",
        hovermode="x unified",
    )
    fig.update_xaxes(rangeslider=dict(visible=False))
    fig.update_yaxes(scaleanchor=None)

    _maybe_save_figure(fig, output_path)
    return fig


def plot_spectrum_with_assigned_formulas(
    df: pd.DataFrame,
    output_path: Optional[Path] = None,
    max_labels: int = 50,
) -> go.Figure:
    """
    Plot spectrum and highlight peaks with assigned formulas.

    Expects:
        - spectrum columns (m/z, Intensity)
        - "Formula" column for labels (optional)
        - "isotopolog" column (optional, used only in hover)

    To avoid overcrowding, only the `max_labels` most intense labelled peaks
    are shown with text at once; users can zoom further with rectangular zoom.
    """
    fig, mz_col, inten_col = _base_spectrum_figure(
        df, title="Spectrum with assigned formulas", showlegend=True
    )

    if "Formula" not in df.columns:
        _maybe_save_figure(fig, output_path)
        return fig

    df_form = df[df["Formula"].notna()].copy()
    if df_form.empty:
        _maybe_save_figure(fig, output_path)
        return fig

    inten = pd.to_numeric(df_form[inten_col], errors="coerce")
    df_form = df_form[inten.notna()].copy()
    df_form["_Intensity_for_sort"] = inten[inten.notna()]
    df_form = df_form.sort_values("_Intensity_for_sort", ascending=False).head(max_labels)

    x = pd.to_numeric(df_form[mz_col], errors="coerce")
    y = pd.to_numeric(df_form[inten_col], errors="coerce")
    formulas = df_form["Formula"].astype(str)
    isotopolog = df_form["isotopolog"].astype(str) if "isotopolog" in df_form.columns else None

    hovertemplate = "<b>m/z</b>: %{x:.6f}<br><b>Intensity</b>: %{y:.2e}<br>"
    hovertemplate += "<b>Formula</b>: %{text}"
    if isotopolog is not None:
        hovertemplate += "<br><b>Isotopolog</b>: %{customdata}"
    hovertemplate += "<extra></extra>"

    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="markers+text",
            name="Assigned formulas",
            marker=dict(size=7, color="orange"),
            text=formulas,
            textposition="top center",
            textfont=dict(size=10),
            customdata=isotopolog,
            hovertemplate=hovertemplate,
        )
    )

    fig.update_layout(
        dragmode="zoom",
        hovermode="x unified",
    )
    fig.update_xaxes(rangeslider=dict(visible=False))
    fig.update_yaxes(scaleanchor=None)

    _maybe_save_figure(fig, output_path)
    return fig