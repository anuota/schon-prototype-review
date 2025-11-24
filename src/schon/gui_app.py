from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st

from schon.sample_run import SampleRun, SampleRunConfig
from schon.calibration.calibration_config import CalibrationConfig
from schon.formula_assignment.formula_presets import DEFAULT_ELEMENT_LIMITS, PRESETS

try:
    # Placeholder for future ESI(+) presets
    from schon.formula_assignment import formula_presets_esi_pos
except ImportError:
    formula_presets_esi_pos = None

try:
    # Placeholder for future APPI presets
    from schon.formula_assignment import formula_presets_appi
except ImportError:
    formula_presets_appi = None 

# Defaults consistent with your Docker layout
DEFAULT_DATA_ROOT = Path("/app/data")
DEFAULT_D_SUBFOLDER = DEFAULT_DATA_ROOT / "Measurement_d-files" / "noncalibrated_d-files"
DEFAULT_RESULTS_DIR = Path("/app/results")
DEFAULT_PEAKS_OUTPUT_DIR = DEFAULT_RESULTS_DIR / "csv_raw"
DEFAULT_CALIB_OUTPUT_DIR = DEFAULT_RESULTS_DIR / "calibration"


def _build_calib_cfg(
    sn_threshold: float,
    poly_order: int,
    only_CHO: bool,
    min_series_length: int,
    series_ppm_tolerance: float,
    link_ppm_tolerance: float,
) -> CalibrationConfig:
    """Helper to construct a CalibrationConfig from GUI parameters."""
    cfg = CalibrationConfig()
    cfg.sn_threshold = sn_threshold
    cfg.poly_order = poly_order
    cfg.only_CHO = only_CHO
    cfg.min_series_length = min_series_length
    cfg.series_ppm_tolerance = series_ppm_tolerance
    cfg.ppm_tolerance_link = link_ppm_tolerance
    return cfg


def _spectrum_plot(df: pd.DataFrame, mz_col: str = "Calibrated m/z") -> plt.Figure:
    """
    Build a simple spectrum plot with all peaks in light grey and
    peaks with assigned formulas highlighted in a different colour.
    """
    fig, ax = plt.subplots(figsize=(10, 4))

    if df.empty:
        ax.text(0.5, 0.5, "No data to plot", ha="center", va="center")
        ax.set_axis_off()
        return fig

    # Use calibrated m/z by default; fall back to raw m/z if needed
    if mz_col not in df.columns:
        mz_col = "m/z"

    # Intensity column: "Intensity" (from formula assignment) or "Peak Height"
    if "Intensity" in df.columns:
        y_col = "Intensity"
    elif "Peak Height" in df.columns:
        y_col = "Peak Height"
    else:
        # fallback: first numeric column after m/z
        numeric_cols = df.select_dtypes("number").columns.tolist()
        numeric_cols = [c for c in numeric_cols if c != mz_col]
        y_col = numeric_cols[0] if numeric_cols else mz_col

    mz = df[mz_col].astype(float)
    intensity = df[y_col].astype(float)

    # All peaks
    ax.vlines(mz, 0, intensity, linewidth=0.8, alpha=0.3)

    # Assigned peaks: where Formula is not NaN / empty
    if "Formula" in df.columns:
        mask_assigned = df["Formula"].astype(str).str.strip() != ""
        mz_assigned = mz[mask_assigned]
        int_assigned = intensity[mask_assigned]
        ax.vlines(
            mz_assigned,
            0,
            int_assigned,
            linewidth=1.0,
            alpha=0.9,
        )

    ax.set_xlabel(mz_col)
    ax.set_ylabel(y_col)
    ax.set_title("Spectrum (assigned peaks highlighted)")
    ax.set_ylim(bottom=0)
    return fig


def main() -> None:
    st.set_page_config(
        page_title="SCHON – FT-ICR-MS GUI",
        layout="wide",
    )

    st.title("SCHON – FT-ICR-MS processing")

    # ================= Sidebar: parameters & controls ==================
    with st.sidebar:
        st.header("Input")

        # ---- Choose input source: path inside container or upload ----
        input_source = st.radio(
            "Input source",
            ["Path inside Docker", "Upload file"],
            index=0,
        )

        input_path: Optional[Path] = None

        if input_source == "Path inside Docker":
            default_path = str(DEFAULT_D_SUBFOLDER)
            input_path_str = st.text_input(
                "Path (inside container)",
                value=default_path,
                help=(
                    "E.g. /app/data/Measurement_d-files/noncalibrated-d-files/ESI_neg_... .d "
                    "or a CSV path."
                ),
            )
            if input_path_str:
                input_path = Path(input_path_str)
        else:
            uploaded = st.file_uploader(
                "Upload CSV (for now)",
                type=["csv"],
                help="For now, use CSV for peaks; .d handling via path is recommended.",
            )
            if uploaded is not None:
                upload_dir = DEFAULT_DATA_ROOT / "gui_uploads"
                upload_dir.mkdir(parents=True, exist_ok=True)
                target_path = upload_dir / uploaded.name
                with target_path.open("wb") as f:
                    f.write(uploaded.getbuffer())
                st.info(f"Saved uploaded file to {target_path}")
                input_path = target_path

        # ---- Main parameters ----
        st.markdown("---")
        st.subheader("Main parameters")

        input_kind = st.selectbox(
            "Input type",
            options=["d", "uncalibrated_csv", "calibrated_csv"],
            format_func=lambda x: {
                "d": ".d folder (Bruker)",
                "uncalibrated_csv": "Uncalibrated peaks CSV",
                "calibrated_csv": "Already calibrated CSV",
            }[x],
        )

        ion_mode = st.selectbox(
            "Ionisation mode",
            options=["ESI-NEG", "ESI-POS", "APPI"],
            index=0,
            help="Used to choose appropriate formula presets and future mode-specific logic.",
        )

        # Sample type depends on ionisation mode
        if ion_mode == "ESI-NEG":
            sample_type_options = [
                "CRUDE_OIL",
                "SEDIMENTARY_ROCK",
                "COAL",
                "NATURAL_WATER",
                "GENERIC_ESI_NEG",
            ]
            default_sample_idx = 4  # GENERIC_ESI_NEG
        elif ion_mode == "ESI-POS":
            # Placeholder; will hook to formula_presets_esi_pos later
            sample_type_options = ["GENERIC_ESI_POS"]
            default_sample_idx = 0
        else:  # APPI
            # Placeholder; will hook to formula_presets_appi later
            sample_type_options = ["GENERIC_APPI"]
            default_sample_idx = 0

        sample_type_label = st.selectbox(
            "Sample type",
            options=sample_type_options,
            index=default_sample_idx,
            help="Sample category; used to choose default formula filters (element limits).",
        )

        # Numeric code for backend (keeps compatibility with existing sample_type usage)
        sample_type_code = sample_type_options.index(sample_type_label)

        # ---- Element limits (formula filters) directly under sample type ----
        element_limits: dict[str, tuple[int, int]] = {}
        preset_key: Optional[str] = None

        with st.expander("Element limits (formula filters)", expanded=False):
            st.caption("Min / max atom counts per element for candidate formulas.")

            # Decide which preset dict and limits to use
            if ion_mode == "ESI-NEG":
                # Map GUI sample type label to preset key
                preset_key_map = {
                    "CRUDE_OIL": "crude_oil",
                    "SEDIMENTARY_ROCK": "sedimentary_rock",
                    "COAL": "coal",
                    "NATURAL_WATER": "natural_water",
                    "GENERIC_ESI_NEG": "generic_esi_neg",
                }
                preset_key = preset_key_map.get(sample_type_label, "generic_esi_neg")
                preset = PRESETS.get(preset_key, PRESETS["generic_esi_neg"])
                default_limits = preset.chnos_limits

            elif ion_mode == "ESI-POS":
                # Placeholder: fall back to default limits unless a POS preset module exists
                preset_key = "generic_esi_pos"
                default_limits = DEFAULT_ELEMENT_LIMITS
                if (
                    formula_presets_esi_pos is not None
                    and hasattr(formula_presets_esi_pos, "PRESETS")
                ):
                    presets_pos = formula_presets_esi_pos.PRESETS
                    if presets_pos:
                        preset_key = next(iter(presets_pos.keys()))
                        default_limits = presets_pos[preset_key].chnos_limits

            else:  # APPI
                preset_key = "generic_appi"
                default_limits = DEFAULT_ELEMENT_LIMITS
                if (
                    formula_presets_appi is not None
                    and hasattr(formula_presets_appi, "PRESETS")
                ):
                    presets_appi = formula_presets_appi.PRESETS
                    if presets_appi:
                        preset_key = next(iter(presets_appi.keys()))
                        default_limits = presets_appi[preset_key].chnos_limits

            # Build editable min/max per element.
            # Keys include ion_mode and sample_type_label so that changing sample type
            # creates a new widget state with appropriate defaults.
            for elem, (emin, emax) in default_limits.items():
                col_min, col_max = st.columns(2)
                key_min = f"{elem}_min_{ion_mode}_{sample_type_label}"
                key_max = f"{elem}_max_{ion_mode}_{sample_type_label}"

                with col_min:
                    new_min = st.number_input(
                        f"{elem} min",
                        min_value=0,
                        max_value=1000,
                        value=int(emin),
                        step=1,
                        key=key_min,
                    )
                with col_max:
                    new_max = st.number_input(
                        f"{elem} max",
                        min_value=0,
                        max_value=1000,
                        value=int(emax),
                        step=1,
                        key=key_max,
                    )

                # IMPORTANT: real < operator, not HTML-escaped
                if new_max < new_min:
                    new_max = new_min

                element_limits[elem] = (int(new_min), int(new_max))

        # ---- Advanced sections ----
        st.markdown("---")
        st.subheader("Advanced settings")

        with st.expander("Calibration settings", expanded=False):
            sn_threshold = st.number_input(
                "Min S/N for calibration",
                min_value=0.0,
                value=20.0,
                step=1.0,
                help="Peaks below this S/N are not used as internal calibrants.",
            )
            poly_order = st.selectbox(
                "Calibration polynomial order",
                options=[1, 2, 3],
                index=1,
                help="1 = linear, 2 = quadratic, 3 = cubic…",
            )
            only_CHO = st.checkbox(
                "Restrict calibrants to CHO",
                value=True,
                help="If checked, only CHO formulas (no N, S, Na, 13C) are used as calibrants.",
            )
            min_series_length = st.number_input(
                "Min homologous series length",
                min_value=2,
                value=2,
                step=1,
                help="Calibrants must belong to a CH2-like series with at least this many members.",
            )
            series_ppm_tolerance = st.number_input(
                "Series matching tolerance (ppm)",
                min_value=0.1,
                value=5.0,
                step=0.1,
                help="Tolerance for detecting CH2-like series during calibrant selection.",
            )
            link_ppm_tolerance = st.number_input(
                "Link tolerance to full list (ppm)",
                min_value=0.1,
                value=2.0,
                step=0.1,
                help="Tolerance when mapping calibrants back to the full peak list.",
            )

        with st.expander("Formula assignment settings", expanded=False):
            formula_ppm_tolerance = st.number_input(
                "Formula search tolerance (ppm)",
                min_value=0.1,
                value=3.0,
                step=0.1,
                help="Default ppm window for formula matching.",
            )
            formula_n_proc = st.number_input(
                "Number of processes for formula assignment",
                min_value=1,
                max_value=32,
                value=4,
                step=1,
                help="Parallel workers for formula assignment (multiprocessing).",
            )

        st.markdown("---")
        run_button = st.button("Run SCHON pipeline")

    # ================= Main layout: results & spectrum window =================
    # Reserve two columns: left = spectrum, right = summary & tables
    col_plot, col_info = st.columns([2, 1])

    if run_button:
        if input_path is None:
            st.error("Please provide a valid input path or upload a file.")
            return

        # Build calibration config from sidebar controls
        calib_cfg = _build_calib_cfg(
            sn_threshold=sn_threshold,
            poly_order=poly_order,
            only_CHO=only_CHO,
            min_series_length=min_series_length,
            series_ppm_tolerance=series_ppm_tolerance,
            link_ppm_tolerance=link_ppm_tolerance,
        )

        # Build per-run config object
        try:
            run_cfg = SampleRunConfig(
                input_kind=input_kind,
                sample_type=sample_type_code,  # numeric index derived from label
                calib_cfg=calib_cfg,
                formula_ppm_tolerance=formula_ppm_tolerance,
                formula_n_processes=int(formula_n_proc),
                data_root=DEFAULT_DATA_ROOT,
                peaks_output_dir=DEFAULT_PEAKS_OUTPUT_DIR,
                calibration_results_dir=DEFAULT_CALIB_OUTPUT_DIR,
                save_raw_peaks_csv=True,
                formula_element_limits=element_limits,
                ion_mode=ion_mode,
                formula_preset_name=preset_key,
            )
        except TypeError:
            # Backwards compatibility if SampleRunConfig doesn't yet accept the new fields
            run_cfg = SampleRunConfig(
                input_kind=input_kind,
                sample_type=sample_type_code,
                calib_cfg=calib_cfg,
                formula_ppm_tolerance=formula_ppm_tolerance,
                formula_n_processes=int(formula_n_proc),
                data_root=DEFAULT_DATA_ROOT,
                peaks_output_dir=DEFAULT_PEAKS_OUTPUT_DIR,
                calibration_results_dir=DEFAULT_CALIB_OUTPUT_DIR,
                save_raw_peaks_csv=True,
            )

        sample_name = input_path.stem.rstrip("/")
        st.write(f"### Running sample: `{sample_name}`")

        # Run the pipeline using SampleRun
        try:
            run = SampleRun(
                sample_name=sample_name,
                input_path=input_path,
                cfg=run_cfg,
            )

            with st.status("Loading input and extracting peaks...", expanded=True) as status:
                run.load_input()
                if run.raw_peaks_df is not None:
                    st.write(f"Raw peaks: {run.raw_peaks_df.shape[0]} rows")
                elif run.calibrated_input:
                    st.write("Input is already calibrated; using it directly.")
                status.update(label="Input loaded", state="complete")

            with st.status("Running calibration...", expanded=False) as status:
                run.run_calibration()
                status.update(label="Calibration finished", state="complete")

            with st.status("Assigning formulas on calibrated peaks...", expanded=False) as status:
                formulas_df = run.assign_formulas_on_calibrated()
                status.update(label="Formula assignment finished", state="complete")

            # Save metadata for this run (so GUI can list it later)
            meta_path = run.save_metadata()
            st.success(f"Run metadata saved to: {meta_path}")

        except Exception as e:
            st.exception(e)
            return

        # ------------------ Left: spectrum window ------------------
        with col_plot:
            st.subheader("Spectrum")

            current_df = run.current_df  # latest (with formulas if any)

            # Optional zoom controls
            mz_min = float(current_df["Calibrated m/z"].min()) if "Calibrated m/z" in current_df.columns else float(current_df["m/z"].min())
            mz_max = float(current_df["Calibrated m/z"].max()) if "Calibrated m/z" in current_df.columns else float(current_df["m/z"].max())

            mz_range = st.slider(
                "m/z range",
                min_value=float(mz_min),
                max_value=float(mz_max),
                value=(float(mz_min), float(mz_max)),
            )

            # Filter by selected m/z range
            mz_col = "Calibrated m/z" if "Calibrated m/z" in current_df.columns else "m/z"
            df_zoom = current_df[
                (current_df[mz_col] >= mz_range[0]) & (current_df[mz_col] <= mz_range[1])
            ]

            fig = _spectrum_plot(df_zoom, mz_col=mz_col)
            st.pyplot(fig, clear_figure=True)

        # ------------------ Right: summary ------------------
        with col_info:
            st.subheader("Run summary")
            st.write(f"**Input file:** {input_path}")
            st.write(f"**Input kind:** {input_kind}")
            st.write(f"**Sample type:** {sample_type_label}")
            st.write(f"**Peaks (calibrated table):** {run.calibrated_df.shape if run.calibrated_df is not None else 'N/A'}")

            if run.latest_formulas_df is not None:
                df_form = run.latest_formulas_df
                n_assigned = df_form["Formula"].notna().sum() if "Formula" in df_form.columns else 0
                st.write(f"**Formula-assigned peaks:** {n_assigned}")

                with st.expander("Show first 200 rows of formula table"):
                    st.dataframe(df_form.head(200))

            with st.expander("Calibration diagnostics", expanded=False):
                if run.calibration_debug:
                    st.json({
                        "sample_name": run.calibration_debug.get("sample_name"),
                        "n_calibrants": len(run.calibration_debug.get("calibrants", [])),
                    })
                else:
                    st.write("No calibration diagnostics available.")


if __name__ == "__main__":
    main()