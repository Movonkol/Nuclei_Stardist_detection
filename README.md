# Bulk Microscopy Analysis (FIJI/ImageJ)

This repository contains scripts for bulk analysis of microscopy images in FIJI/ImageJ.

- **nuclei_stardist_sorted.py** – Jython macro for StarDist-based nuclei segmentation  
- **intensity.py** – Jython script for quantitative intensity measurements  
- **positiv_negativ.py** – Jython script for positive/negative classification based on intensities  

Licensed under the MIT License.

---

## 1) StarDist Batch Analysis (`nuclei_stardist_sorted.py`)

Automated segmentation and analysis of cell nuclei using StarDist and Bio-Formats.

### Features
- **Interactive setup**: model, probability/NMS thresholds, tiling, channel patterns, size filters
- **Channel detection**: e.g. DAPI via patterns like `c1-`, `dapi`
- **Size filtering**: min/max label area to remove artifacts
- **StarDist segmentation**: headless processing
- **Morphology metrics**: area, circularity per object
- **CSV export**: `nuclei_counts.csv` (summary) and `nuclei_morphology.csv` (per object)
- **Batch processing**: all supported files; multi-series aware

### Installation
1. FIJI/ImageJ (Java 8+)
2. Plugins: **Bio-Formats**, **StarDist**, **MorphoLibJ** & **IJPB-plugin**, **Jython** (bundled)
3. Place scripts in your FIJI scripts folder, then restart FIJI.

### Configuration (defaults at top of file)
```java
model = "Versatile (fluorescent nuclei)"
probability_thresh = 0.5
nms_thresh = 0.5
tiles_str = "1,1"
channel_patterns_str = "c1-,dapi"
min_label_size = 50
max_label_size = 2500
```

### Usage
1. FIJI → **Plugins > Scripting > Open...**
2. Run `nuclei_stardist_sorted.py`
3. Set parameters in dialogs
4. Select image folder and run

### Output
- Temporary annotated label windows (auto-closed)
- `nuclei_counts.csv`
- `nuclei_morphology.csv`

**CSV formats**

`nuclei_counts.csv`
| Column | Description |
| --- | --- |
| Image | Filename or series identifier |
| Nuclei_Count | Total number of detected nuclei |

`nuclei_morphology.csv`
| Column | Description |
| --- | --- |
| Image | Filename or series identifier |
| Object | Object index (per image) |
| Area | Object area (pixels) |
| Roundness | Circularity metric |

---

## 2) Intensity Measurement (`intensity.py`)

Quantifies marker-channel intensities **per series** and **aggregated per well**; optionally generates heatmaps.

### Features
- **Interactive configuration**
  - Marker channel patterns (e.g. `c3-`, `gfp`, `FITC`)
  - Thresholding: **automatic** (Otsu/Yen/Moments) with scaling factor **or** **fixed value**
  - Optional **background subtraction** (rolling-ball radius, repeats, median filter or defaults)
  - Optional **heatmap generation** (PNG) and **auto-close** of windows
- **Batch processing** for multi-series files and all supported formats
- **CSV export**
  - `intensity_measurements.csv` (per series)
  - `intensity_measurements_per_well.csv` (per well, aggregated)

### Usage
1. FIJI → **Plugins > Scripting > Run...** → select `intensity.py`
2. Dialogs:
   1) Marker channel pattern(s)  
   2) Thresholding (auto method + factor **or** fixed value)  
   3) Background subtraction (custom or default)  
   4) Heatmaps (Yes/No)  
   5) Close images after processing (Yes/No)
3. Select image folder and run

### Output
- `intensity_measurements.csv` (per series)
- `intensity_measurements_per_well.csv` (aggregated per well)
- `heatmaps/` (PNG heatmaps, if enabled)

**CSV formats**

`intensity_measurements.csv`
| Column | Description |
| --- | --- |
| Image | Filename or series identifier |
| Series | Series index/name (if applicable) |
| Well | Parsed well ID (e.g. A01) or filename prefix |
| Marker | Marker/Channel identifier (e.g. `c3-`) |
| Positive_Pixels | Pixels above threshold (count) |
| Mean_Intensity | Mean intensity within positives |
| Integrated_Intensity | Sum of intensities within positives |
| Nuclei_Count | (If available) nuclei count for normalization |
| Normalized_Intensity | Integrated_Intensity / Nuclei_Count (if provided) |

`intensity_measurements_per_well.csv`
| Column | Description |
| --- | --- |
| Well | Well ID parsed from filename/series |
| Marker | Marker/Channel identifier |
| Series_Count | Number of contributing series |
| Positive_Pixels_Total | Sum across series in well |
| Integrated_Intensity_Total | Sum across series in well |
| Nuclei_Count_Total | Sum across series in well (if available) |
| Normalized_Intensity | Integrated_Intensity_Total / Nuclei_Count_Total (if available) |

**Notes**
- Heatmaps are pseudocolor projections of marker intensity per series/well.
- Channel selection is pattern-based; prefer consistent prefixes like `c1-`, `c2-`, `c3-`.

---

## 3) Positive/Negative Classification (`positiv_negativ.py`)

Classifies nuclei as **Positive**/**Negative** based on marker intensities.

### Configure (top of script)
- **NUCLEI_CHANNEL_KEY** (e.g. `"c4-"`)
- **MARKER_CHANNEL_KEYS** (e.g. `["c3-"]`)
- **Size/shape filters**: `SIZE_MIN`, `SIZE_MAX` (pixels), `ASPECT_RATIO_MAX`
- **StarDist**: `model`, `prob`, `nms`, `tiles_str`
- **Optional background subtraction**: `ROLLING_RADIUS`, `ROLLING_REPEAT`, `MEDIAN_RADIUS`
- **Thresholding**: automatic (method + factor) or fixed (`FIXED_THRESHOLD`)

### Usage
1. FIJI → **Plugins > Scripting > Run…** → `positiv_negativ.py`
2. Select image folder
3. Answer dialogs (background subtraction, thresholding, etc.)
4. Processing steps:
   - Split channels → select nuclei (DAPI) + marker windows  
   - Segment nuclei via StarDist → filter ROIs  
   - (Optional) background subtraction + thresholding on marker  
   - Classify each nucleus → **Positive** (green) / **Negative** (red)  
   - Save:
     - `<series>_<marker>_RAW.png`
     - `<series>_<marker>_RAW_classified.png`
   - Append counts to `nuclei_counts.csv`

**CSV format (`nuclei_counts.csv`)**
| Column | Description |
| --- | --- |
| Image | Series name (e.g. `Filename_Series1_c3`) |
| Nuclei_Count | Total filtered nuclei |
| Positive_Nuclei | Positive nuclei |
| Negative_Nuclei | Negative nuclei |
## AOI-based Nuclei Pos/Neg (Fiji/Jython)

4) AOI-based Nuclei Pos/Neg (positivearea_Positivnegativ.py)
Classifies nuclei (StarDist ROIs) as positive/negative per marker **within** a thresholded AOI (“TOTAL” channel); optional background subtraction; multi-series; exports overlays and CSV.
_Based on:_ **Nuclei StarDist detection** (segments nuclei) and **Positive Area (%AREA)** (thresholds TOTAL to define AOI); combines both by evaluating ROI ∩ AOI.

Features
- Interactive configuration
- Channel patterns: nuclei/DAPI, AOI (TOTAL), marker(s) (e.g., `c4-`, `c1-`, `c3-`)
- Fixed thresholds in original bit depth (AOI + per-marker)
- Positivity rules: minimum AOI coverage (%) and positivity fraction (% pixels ≥ marker threshold in ROI∩AOI)
- Optional background subtraction (rolling-ball radius, repeats, median or defaults)
- Optional overlays (P/N labels) and contrast enhancement
- Batch processing for multi-series files; auto-close windows
- CSV export

Usage
- FIJI → Plugins > Scripting > Script Editor (Python/Jython) → open `positivearea_Positivnegativ.py` → Run
- Dialogs: nuclei channel; AOI/TOTAL pattern(s); marker pattern(s); AOI + per-marker thresholds; min AOI coverage; positivity fraction; StarDist (model, prob, NMS, tiles); background subtraction; overlays/contrast
- Select image folder and run

Output
- `nuclei_posneg_in_AOI.csv` (per image & marker summary)
- `nuclei_posneg_perROI_in_AOI.csv` (per ROI & marker details)
- `AOI_PosNeg_PNGs/`
  - `...__AOI_PosNeg__DAPI.png` (overlay on nuclei/DAPI)
  - `...__AOI_PosNeg__MARKER.png` (overlay on marker)

CSV formats

nuclei_posneg_in_AOI.csv

Column                 Description
Image                  Filename or series identifier
Series                 Series index/name (if applicable)
Marker                 Marker/channel identifier
N_in_AOI               Number of nuclei intersecting AOI
Positive               Count called positive
Negative               Count called negative
Percent_Positive       0–100 (%)

nuclei_posneg_perROI_in_AOI.csv

Column                 Description
Image                  Filename or series identifier
Series                 Series index/name
ROI_Index              Nucleus index
Marker                 Marker/channel identifier
ROI_px                 ROI pixels (area)
ROI_in_AOI_px          ROI pixels inside AOI
PosPix_in_AOI          Pixels ≥ marker threshold in ROI∩AOI
Is_Positive            Boolean (per rules)



