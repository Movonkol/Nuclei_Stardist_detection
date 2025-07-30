Stardist Batch Analysis Script & Intensity Measurement
This repository contains two scripts for bulk analysis of microscopy images in FIJI/ImageJ: a StarDist‑based nuclei segmentation macro and a Jython script for quantitative intensity measurements. Both scripts are licensed under the MIT License.

StarDist Batch Analysis
Automated segmentation and analysis of cell nuclei in bulk using StarDist and Bio‑Formats in FIJI/ImageJ.

Features
Automatic Parameter Configuration: Prompts for model, probability & NMS thresholds, tiling, channel patterns and size filters.

Channel Detection: Identifies target channel(s) (e.g. DAPI) via customizable search patterns.

Size Filtering: Applies minimum and maximum label size thresholds to exclude artefacts.

StarDist Segmentation: Runs headless StarDist for label image generation.

Morphology Measurement: Measures area and circularity for each detected object.

CSV Export: Generates nuclei_counts.csv and nuclei_morphology.csv with summary and per‑object data.

Batch Processing: Processes all supported image files in a selected folder, including multi‑series imports via Bio‑Formats.

Installation
Clone the repository and install FIJI/ImageJ with the required plugins:

bash
Kopieren
Bearbeiten
git clone https://github.com/yourusername/stardist-batch.git
cd stardist-batch
Install FIJI/ImageJ (Java 8+) and ensure the required plugins are installed:

For the StarDist macro: the StarDist plugin, Bio‑Formats (to import ND2/LIF/TIFF files) and MorphoLibJ (for morphological operations).

For the intensity measurement script: Bio‑Formats and Jython (included with FIJI). If you plan to perform morphological measurements or filtering in your own modifications, install MorphoLibJ as well.

Place the macro stardist_batch.ijm and the script intensity.py in your FIJI scripts directory and restart FIJI.

Configuration
Edit the first section of stardist_batch.ijm if you prefer hardcoded defaults:

groovy
Kopieren
Bearbeiten
model = "Versatile (fluorescent nuclei)"
probability_thresh = 0.5
nms_thresh = 0.5
tiles_str = "1,1"
channel_patterns_str = "c1-,dapi"
min_label_size = 50
max_label_size = 2500
Usage
Open FIJI and navigate to Plugins > Scripting > Open....

Select and run stardist_batch.ijm.

Configure parameters in the dialog prompts.

Select the folder containing your images.

Wait for processing; progress is printed to the console.

Output
Annotated label windows appear during processing (auto‑closed).

nuclei_counts.csv (image name, nuclei count)

nuclei_morphology.csv (image name, object index, area, roundness)

CSV Format
Column	Description
Image	Filename or series identifier
Nuclei_Count	Total number of detected nuclei
Object	Object index (per image)
Area	Object area in pixels
Roundness	Object circularity metric

Intensity Measurement (intensity.py)
This Jython script measures integrated intensities of marker channels on a per‑series and per‑well basis. It supports multi‑series files (ND2, TIFF, LIF, etc.), optional background subtraction and thresholding, and can generate pseudocolour heatmaps of the measured channels.

Features
Interactive Configuration: Prompts for marker channel patterns, thresholding (automatic or fixed), background subtraction (with optional custom parameters), heatmap generation and whether to close images after processing.

Background Subtraction: Configurable rolling‑ball radius, number of repetitions and median filter radius, or sensible defaults.

Thresholding: Choose between automatic thresholding (e.g. Otsu, Yen, Moments) with a scaling factor or a fixed threshold value.

CSV Export: Writes intensity_measurements.csv with per‑series measurements and intensity_measurements_per_well.csv with aggregated per‑well metrics. Normalisation relies on the nuclei_counts.csv produced by the StarDist script, so please run the segmentation macro first if you intend to compute normalised intensities.

Heatmaps (Optional): Generates pseudocolour heatmap PNGs for each marker channel, saved to a heatmaps subfolder.

Batch Processing: Processes all supported image files in a selected folder and handles Bio‑Formats multi‑series imports automatically.

Usage
Open FIJI/ImageJ and run the script via Plugins > Scripting > Run..., selecting intensity.py from this repository.

Answer the configuration prompts:

Marker channel pattern(s): Comma‑separated substrings used to identify marker channels.

Automatic thresholding? If yes, choose a threshold method and scaling factor; otherwise provide a fixed threshold value.

Subtract background? If yes, choose whether to enter custom background parameters (rolling‑ball radius, repeat count and median radius) or use defaults.

Save heatmap images? If yes, pseudocolour heatmaps will be generated for each marker channel.

Close images after processing? Select whether to close all image windows automatically when done.

Select the folder containing your images.

Wait for processing; progress and messages are printed to the console.

Output
intensity_measurements.csv (per‑series measurements: positive pixels, mean intensity, integrated intensity, nuclei count, normalised intensity)

intensity_measurements_per_well.csv (aggregated per‑well metrics)

Optional heatmaps saved in a heatmaps subfolder (if enabled)

Heatmaps are saved automatically using the Fire lookup table and do not overwrite your original images. All CSV files and heatmaps reside alongside your input data.

Positive/Negative Classification (positiv_negativ.py)
The optional positiv_negativ.py script classifies wells (or individual images) as positive or negative based on the measured intensities. It reads the aggregated results produced by intensity.py and applies a user‑defined threshold to the normalised integrated intensity.

Features
Threshold‑based classification: Uses a threshold on normalised intensity to assign each well to a positive or negative category.

CSV Input/Output: Reads intensity_measurements_per_well.csv and writes well_classification.csv with one row per well and its classification.

Simple Workflow: Prompts for a threshold value and the input folder; no additional plugins are required beyond those used for intensity.py.

Usage
Make sure you have run the StarDist segmentation and intensity measurement scripts first, so that nuclei_counts.csv and intensity_measurements_per_well.csv exist.

Run positiv_negativ.py via Plugins > Scripting > Run... in FIJI/ImageJ.

When prompted, enter the threshold for normalised intensity and select the folder containing your CSV files.

The script writes well_classification.csv to the selected folder.

Contributing
Contributions are welcome! Submit issues or pull 
