# Stardist Batch Analysis & Intensity Measurement

This repository contains scripts for bulk analysis of microscopy images in FIJI/ImageJ:

* **nuclei_stardist_sorted.py**: A Jython script for StarDist-based nuclei segmentation macro
* **intensity.py**: A Jython script for quantitative intensity measurements
* **positiv\_negativ.py**: A Jython script for positive/negative classification based on measured intensities

All scripts are licensed under the MIT License.

---

## StarDist Batch Analysis

Automated segmentation and analysis of cell nuclei in bulk using StarDist and Bio-Formats.

### Features

* **Automatic Parameter Configuration**: Prompts for model, probability & NMS thresholds, tiling, channel patterns, and size filters.
* **Channel Detection**: Identifies target channel(s) (e.g. DAPI) via customizable patterns.
* **Size Filtering**: Applies minimum/maximum label size to exclude artifacts.
* **StarDist Segmentation**: Runs headless StarDist for label image generation.
* **Morphology Measurement**: Measures area and circularity for each object.
* **CSV Export**: Generates `nuclei_counts.csv` (summary) and `nuclei_morphology.csv` (per-object data).
* **Batch Processing**: Processes all supported image files in a folder, including multi-series imports.

### Installation

1. Install FIJI/ImageJ (Java 8+).
2. Ensure the following plugins are installed:

   * Bio-Formats (for ND2/LIF/TIFF import)
   * StarDist plugin
   * MorphoLibJ and  IJPB-plugin (for morphological operations)
   * Jython (included with FIJI)
3. Place the scripts in your FIJI scripts directory.
4. Restart FIJI.

### Configuration

Edit the top of `nuclei_stardist_sorted.py` for hardcoded defaults:

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

1. Open FIJI > **Plugins > Scripting > Open...**
2. Select and run `nuclei_stardist_sorted.py.
3. Configure parameters in the dialog prompts.
4. Choose the folder containing your images.
5. Wait for processing; progress appears in the console.

### Output

* **Annotated label windows** (auto-closed)
* `nuclei_counts.csv` (image name, nuclei count)
* `nuclei_morphology.csv` (image name, object index, area, roundness)

### CSV Format

| Column                      | Description                     |
| --------------------------- | ------------------------------- |
| Image                       | Filename or series identifier   |
| Nuclei\_Count               | Total number of detected nuclei |
| Object (nuclei\_morphology) | Object index (per image)        |
| Area                        | Object area in pixels           |
| Roundness                   | Object circularity metric       |

---

## Intensity Measurement (`intensity.py`)

Measures integrated intensities of marker channels on a per-series and per-well basis.

### Features

* **Interactive Configuration**: Prompts for marker channel patterns, thresholding (automatic or fixed), background subtraction, heatmap generation, and window closing.
* **Background Subtraction**: Rolling-ball radius, repetitions, median filter, or defaults.
* **Thresholding**: Automatic (Otsu, Yen, Moments) with scaling factor or fixed value.
* **CSV Export**: `intensity_measurements.csv` (per-series) and `intensity_measurements_per_well.csv` (aggregated per-well).
* **Optional Heatmaps**: Pseudocolor heatmaps saved as PNGs in a `heatmaps/` subfolder.
* **Batch Processing**: Handles multi-series files and all supported formats.

### Usage

1. In FIJI, go to **Plugins > Scripting > Run...** and select `intensity.py`.
2. Follow the prompts:

   1. Marker channel pattern(s)
   2. Automatic thresholding? If yes, choose method and scaling; otherwise, enter fixed value.
   3. Background subtraction? If yes, choose custom or default settings.
   4. Generate heatmaps? Yes/No.
   5. Close images after processing? Yes/No.
3. Select the folder with your images.
4. Wait for processing; progress appears in the console.

### Output

* `intensity_measurements.csv` (per-series: positive pixels, mean intensity, integrated intensity, nuclei count, normalized intensity)
* `intensity_measurements_per_well.csv` (aggregated per-well metrics)
* Optional heatmaps in `heatmaps/` if enabled.

# Positive/Negative Classification

Configure the parameters at the top of the script before running:

- **NUCLEI_CHANNEL_KEY** (e.g., `"c4-"`): Substring to identify the nuclei/DAPI channel window  
- **MARKER_CHANNEL_KEYS** (e.g., `["c3-"]`): List of substrings for marker channel windows  
- **Size and shape filters**:  
  - `SIZE_MIN`, `SIZE_MAX` (pixel area)  
  - `ASPECT_RATIO_MAX`  
- **StarDist parameters**:  
  - `model` (e.g., "Versatile (fluorescent nuclei)")  
  - `prob` (probability threshold)  
  - `nms` (non-maximum suppression threshold)  
  - `tiles_str` (e.g., "1,1")  
- **Optional background subtraction**:  
  - `ROLLING_RADIUS`, `ROLLING_REPEAT`, `MEDIAN_RADIUS`  
- **Thresholding method**:  
  - Automatic (method and factor) or fixed value (`FIXED_THRESHOLD`)  

### Usage

1. In FIJI, go to **Plugins ▶ Scripting ▶ Run…**, then select `positiv_negativ.py`.  
2. Choose the folder containing your image files (`.tif, .tiff, .png, .jpg, .lif, .nd2`).  
3. Respond to the dialog prompts for background subtraction, thresholding, and other options.  
4. The script processes each series in the folder by:  
   1. Splitting channels and selecting nuclei (DAPI) and marker windows.  
   2. Segmenting nuclei with StarDist and filtering ROIs by size and shape.  
   3. (If enabled) Applying background subtraction and thresholding to the marker channel.  
   4. Classifying each nucleus as **Positive** (green outline) or **Negative** (red outline).  
   5. Saving two images per series and marker:  
      - `<series>_<marker>_RAW.png`  
      - `<series>_<marker>_RAW_classified.png`  
   6. Appending counts to `nuclei_counts.csv` in the selected folder.  

### Output: `nuclei_counts.csv`

| Column               | Description                                    |
|----------------------|------------------------------------------------|
| **Image**            | Series name (e.g., `Filename_Series1_c3`)      |
| **Nuclei_Count**     | Total number of filtered nuclei                |
| **Positive_Nuclei**  | Number of nuclei classified as positive        |
| **Negative_Nuclei**  | Number of nuclei classified as negative        |
