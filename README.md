# SCHON — FT-ICR-MS Processing Pipeline  
*A modular toolkit for peak extraction, calibration, and molecular formula assignment*

---

## 🌟 Overview

**SCHON** is a complete processing workflow for **FT-ICR-MS (Fourier-Transform Ion Cyclotron Resonance Mass Spectrometry)** data.  
It supports both:

- **Bruker `.d` folders** (native, no zipping required)
- **CSV peak lists** (uncalibrated or pre-calibrated)

The pipeline performs:

1. **Peak extraction** from Bruker FID / Apex `.d` folders  
2. **Internal calibration** with homologous-series recognition  
3. **Molecular formula assignment** using customizable presets  
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
- Supports CHO-restricted or full calibrant selection  
- Homologue series detection  
- Polynomial calibration (orders 1–3)  
- Per-peak calibration error (ppm)  
- Automatic tagging of peaks **used for calibration**

### 🔹 Molecular Formula Assignment  
- Fully in-memory, fast assignment engine  
- Adjustable ppm tolerance  
- Sample-specific presets:  
  - CRUDE_OIL  
  - SEDIMENTARY_ROCK  
  - COAL  
  - NATURAL_WATER  
  - GENERIC_ESI_NEG  
- Editable elemental limits in GUI  
- Future-ready support for ESI-POS and APPI

### 🔹 Plotting  
All plots are generated via Plotly and saved to:
