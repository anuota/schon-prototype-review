# SCHON — FT-ICR-MS Processing Pipeline  
*A modular toolkit for peak extraction, calibration, and molecular formula assignment*

---

## 🌟 Overview

**SCHON** is a complete processing workflow for **FT-ICR-MS (Fourier-Transform Ion Cyclotron Resonance Mass Spectrometry)** data.  
It supports both:

- **Bruker `.d` folders** (native, no zipping required)
- **CSV peak lists** (uncalibrated or pre-calibrated)
- support for **Thermo .raw** files from Orbitrap  to be included soon 

The pipeline performs:

1. **Peak extraction** from Bruker FID / Apex `.d` folders  
2. **Internal calibration** with homologous-series recognition  
3. **Fast Molecular formula assignment** using customizable presets  
4. **Probabilistic Molecular formula assignment using KMD series** is to be included soon 
4. **Plotting** (calibrants, raw vs calibrated, assigned formulas)  
5. **Reproducible run tracking** with per-sample run history  
6. **Interactive GUI (Streamlit)**

SCHON is designed for petroleomics, environmental samples, organic geochemistry, OMICS-style FTICR datasets, and any workflow requiring robust ultra-high-resolution mass processing.

---

## ✨ Key Features

### 🔹 Peak Extraction  
- Extracts peaks from **Bruker `.d`** folders using CoreMS-based logic  
- Optional CSV input  

### 🔹 Internal Calibration  
- Supports CHO (or CHON)-restricted or full calibrant selection  
- Homologue series detection  
- Polynomial calibration (orders 1–3)  
- Per-peak calibration error (ppm)  
- Automatic tagging of peaks **used for calibration**

### 🔹 Molecular Formula Assignment  
- Fully in-memory, fast assignment engine  
- Adjustable ppm tolerance 
- Intrument mode specific presets:
    - ESI-NEG
    - ESI-POS (to be included)
    - APPI (to be included)
- Sample-specific presets:  
  - CRUDE_OIL  
  - SEDIMENTARY_ROCK  
  - COAL  
  - NATURAL_WATER  
  - GENERIC_ESI_NEG  
- Editable elemental limits in GUI  

### 🔹 Plotting  
All plots are generated via Plotly and saved to:
```
/app/results/runs/<sample_name>/<run_id>/plots/
```

 Included views:

- Spectrum with **calibrants** (red points + formula labels)  
- **Raw vs calibrated** m/z comparison  
- Spectrum with **formula-assigned peaks** (orange markers)  
- Interactive zooming and box-zoom selection  

### 🔹 Reproducible Run Management  
- Each run gets a **unique `run_id`**  
- Stored under:
```
/app/results/runs/<sample_name>/<run_id>/
```
- Includes:
  - Final CSV table  
  - run_config.json  
  - Calibration debug info  
  - Plots  
- `runs_index.json` tracks all runs per sample.

### 🔹 Streamlit GUI  
A graphical interface for:

- Selecting `.d` folders or CSV  
- Editing calibration parameters  
- Choosing ion mode and sample type  
- Adjusting formula element limits  
- Viewing interactive plots  
- Browsing run history  
- Previewing output tables  

---

## 🐳 Running SCHON in Docker

### Build (if modifying locally)

```bash
docker build -t fticr-corems .
```

### Run with GUI on localhost:8501

```bash
docker run --rm -p 8501:8501 \
  -v /Users/anya/Coding/SCHON:/app \
  fticr-corems
```
### 📁 Output Structure
```
/app/results/
    runs/
        <sample_name>/
            runs_index.json        # List of all runs
            <run_id>/
                run_config.json
                <sample>_calibrated_with_formulas.csv
                plots/
                    calibrants.png
                    raw_vs_calibrated.png
                    assigned_formulas.png
```
## 📦 Installation (development mode)
Clone:
```
git clone https://github.com/<your_name>/SCHON.git
cd SCHON
```
Run locally (no Docker):
```
pip install -r requirements.txt
streamlit run src/schon/gui_app.py
```
### 🧩 Code Structure
```
src/schon/
    gui_app.py               # Streamlit GUI
    main.py                  # CLI pipeline
    sample_run.py            # Orchestrates each run
    calibration/
        calibration.py
        calibration_config.py
    formula_assignment/
        formula_assignment.py
        formula_filters.py
        formula_presets.py
    plotting/
        plot.py
    corems_pipeline/
        generate_peaks_from_fid.py
```
## 📝 License

    SCHON uses:
        •	EUPL
        •	Streamlit (Apache 2.0)

    Streamlit Cloud hosting has restrictions (no sensitive personal data),
    but local Streamlit usage inside Docker has no restrictions.

## 🤝 Contributing

Pull requests and issues are welcome!
Areas open for contribution:
	•	ESI-POS and APPI presets
	•	Additional calibrant strategies
	•	Peak deconvolution
	•	Export formats (HDF5, Parquet, JSON-LD)
	•	Improved GUI plotting and QC workflows

⸻

## 📬 Contact

If you use SCHON for research, we’d love to hear from you.
Open an issue or reach out directly.