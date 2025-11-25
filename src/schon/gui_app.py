from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd
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

from schon.plotting.plot import (
    plot_spectrum_with_calibrants,
    plot_raw_vs_calibrated,
    plot_spectrum_with_assigned_formulas,
)

# Defaults consistent with your Docker layout
DEFAULT_DATA_ROOT = Path("/app/data")
DEFAULT_D_SUBFOLDER = DEFAULT_DATA_ROOT / "Measurement_d-files" / "noncalibrated_d-files" / "ESI_neg_G017736-0_100_200sc_000001.d"
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
            # default to GENERIC_ESI_NEG
            default_sample_idx = sample_type_options.index("GENERIC_ESI_NEG")
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
                value=1.0,
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

    # --- Persist last successful run in session_state so UI updates (like plot view)
    # --- do not require rerunning the whole pipeline ---
    if "last_run" not in st.session_state:
        st.session_state["last_run"] = None
    if "last_input_path" not in st.session_state:
        st.session_state["last_input_path"] = None
    if "last_input_kind" not in st.session_state:
        st.session_state["last_input_kind"] = None
    if "last_sample_type_label" not in st.session_state:
        st.session_state["last_sample_type_label"] = None

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
                sample_type=sample_type_label,
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
                sample_type=sample_type_label,
                calib_cfg=calib_cfg,
                formula_ppm_tolerance=formula_ppm_tolerance,
                formula_n_processes=int(formula_n_proc),
                data_root=DEFAULT_DATA_ROOT,
                peaks_output_dir=DEFAULT_PEAKS_OUTPUT_DIR,
                calibration_results_dir=DEFAULT_CALIB_OUTPUT_DIR,
                save_raw_peaks_csv=True,
            )

        # Determine which paths to process (support bulk .d folders)
        if input_kind == "d" and input_path.is_dir():
            if input_path.suffix == ".d":
                # Single .d folder
                paths_to_process = [input_path]
            else:
                # Directory containing multiple .d folders
                paths_to_process = sorted(
                    p for p in input_path.iterdir() if p.is_dir() and p.suffix == ".d"
                )
                if not paths_to_process:
                    st.error("No .d folders found inside the selected directory.")
                    return
        else:
            # CSV and other cases: treat the input as a single item
            paths_to_process = [input_path]

        last_run: Optional[SampleRun] = None

        # Process each path
        for this_path in paths_to_process:
            sample_name = this_path.name.rstrip("/")
            st.write(f"### Running sample: `{sample_name}`")

            try:
                run = SampleRun(
                    sample_name=sample_name,
                    input_path=this_path,
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
                    _ = run.assign_formulas_on_calibrated()
                    status.update(label="Formula assignment finished", state="complete")

                # Save metadata for this run (so GUI can list it later)
                meta_path = run.save_metadata()
                st.success(f"Run metadata saved to: {meta_path}")

                # Remember the last successful run for plotting & summary
                last_run = run

            except Exception as e:
                st.error(f"Error while processing {this_path}:")
                st.exception(e)

        # After processing all paths, persist the last successful run in session state
        if last_run is not None:
            st.session_state["last_run"] = last_run
            st.session_state["last_input_path"] = input_path
            st.session_state["last_input_kind"] = input_kind
            st.session_state["last_sample_type_label"] = sample_type_label

    # ================= Rendering section (runs on every rerun) =================
    run: Optional[SampleRun] = st.session_state.get("last_run")
    if run is None:
        # Nothing to show yet (pipeline has not been run successfully)
        return

    current_df = run.current_df  # latest (with formulas if any)

    # Guard against missing or empty data
    if current_df is None or current_df.empty:
        st.warning("No data available to plot yet.")
        return

    # Decide which m/z column to use for zooming
    if "Calibrated m/z" in current_df.columns:
        mz_col = "Calibrated m/z"
    elif "m/z" in current_df.columns:
        mz_col = "m/z"
    else:
        st.warning("No m/z column found in current table.")
        return

    # Drop NaNs before computing slider bounds
    mz_series = pd.to_numeric(current_df[mz_col], errors="coerce").dropna()
    if mz_series.empty:
        st.warning("All m/z values are NaN or invalid; cannot plot spectrum.")
        return

    mz_min = float(mz_series.min())
    mz_max = float(mz_series.max())

    # Reserve two columns: left = spectrum, right = summary & controls
    col_plot, col_info = st.columns([2, 1])

    # ------------------ Right: view selector & summary ------------------
    with col_info:
        st.subheader("Run summary")
        last_input_path = st.session_state.get("last_input_path")
        last_input_kind = st.session_state.get("last_input_kind")
        last_sample_type_label = st.session_state.get("last_sample_type_label")

        st.write(f"**Input file:** {last_input_path}")
        st.write(f"**Input kind:** {last_input_kind}")
        st.write(f"**Sample type:** {last_sample_type_label}")
        st.write(
            f"**Peaks (calibrated table):** "
            f"{run.calibrated_df.shape if run.calibrated_df is not None else 'N/A'}"
        )

        # Plot view selector
        view_type = st.radio(
            "Plot view",
            options=[
                "Calibrants",
                "Raw vs calibrated",
                "Assigned formulas",
            ],
            index=0,
            help="Choose which spectrum overlay to display on the left.",
        )

        if run.latest_formulas_df is not None:
            df_form = run.latest_formulas_df
            n_assigned = (
                df_form["Formula"].notna().sum()
                if "Formula" in df_form.columns
                else 0
            )
            st.write(f"**Formula-assigned peaks:** {n_assigned}")

            with st.expander("Show first 200 rows of formula table"):
                st.dataframe(df_form.head(200))

        with st.expander("Calibration diagnostics", expanded=False):
            if run.calibration_debug:
                st.json(
                    {
                        "sample_name": run.calibration_debug.get("sample_name"),
                        "n_calibrants": len(
                            run.calibration_debug.get("calibrants", [])
                        ),
                    }
                )
            else:
                st.write("No calibration diagnostics available.")

    # ------------------ Left: spectrum window ------------------
    with col_plot:
        st.subheader("Spectrum")

        # If range collapsed, just plot without slider
        if mz_min >= mz_max:
            df_zoom = current_df.copy()
        else:
            mz_range = st.slider(
                "m/z range",
                min_value=float(mz_min),
                max_value=float(mz_max),
                value=(float(mz_min), float(mz_max)),
            )

            # Filter by selected m/z range
            df_zoom = current_df[
                (current_df[mz_col] >= mz_range[0])
                & (current_df[mz_col] <= mz_range[1])
            ]
        plots_dir = Path("/app/results/runs") / run.sample_name / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        
        # Delegate actual plotting to the plotting module
        if view_type == "Calibrants":
            fig = plot_spectrum_with_calibrants(df_zoom, output_path=plots_dir / "calibrants.png")
        elif view_type == "Raw vs calibrated":
            fig = plot_raw_vs_calibrated(df_zoom, output_path=plots_dir / "raw_vs_calibrated.png")
        else:  # "Assigned formulas"
            fig = plot_spectrum_with_assigned_formulas(df_zoom, output_path=plots_dir / "assigned_formulas.png")

        st.pyplot(fig, clear_figure=True)

    # ================= Output table preview (full width) =================
    st.markdown("---")
    st.subheader("Output table preview")
    try:
        st.dataframe(current_df.head(500), width='stretch')
    except Exception:
        st.dataframe(current_df, use_container_width=True)


if __name__ == "__main__":
    main()