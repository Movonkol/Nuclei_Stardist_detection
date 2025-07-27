Stardist Batch Analysis Script
License: MIT IJ Macro Version

Automated segmentation and analysis of cell nuclei in bulk using StarDist and Bio-Formats in FIJI/ImageJ.

Features
Automatic Parameter Configuration: Prompts for model, probability & NMS thresholds, tiling, channel patterns, and size filters.
Channel Detection: Identifies target channel(s) (e.g., DAPI) via customizable search patterns.
Size Filtering: Applies minimum and maximum label size thresholds to exclude artifacts.
StarDist Segmentation: Runs headless StarDist for label image generation.
Morphology Measurement: Measures area and circularity for each detected object.
CSV Export: Generates `nuclei_counts.csv` and `nuclei_morphology.csv` with summary and per-object data.
Batch Processing: Processes all supported image files in a selected folder, including multi-series imports via Bio-Formats.

Installation
Clone the repository:
```
git clone https://github.com/yourusername/stardist-batch.git
cd stardist-batch
```

Install FIJI/ImageJ and required plugins:
- FIJI/ImageJ (Java 8+)
- StarDist plugin
- Bio-Formats plugin

Place the script `stardist_batch.ijm` in your FIJI scripts directory and restart FIJI.

Configuration
Edit the first section of `stardist_batch.ijm` if you prefer hardcoded defaults:
```groovy
model = "Versatile (fluorescent nuclei)"
probability_thresh = 0.5
nms_thresh = 0.5
tiles_str = "1,1"
channel_patterns_str = "c1-,dapi"
min_label_size = 50
max_label_size = 2500
```

Usage
1. Open FIJI and navigate to `Plugins > Scripting > Open...`.
2. Select and run `stardist_batch.ijm`.
3. Configure parameters in the dialog prompts.
4. Select the folder containing your images.
5. Wait for processing; progress is printed to the console.

Output
- Annotated label windows appear during processing (auto-closed).
- `nuclei_counts.csv` (image name, nuclei count)
- `nuclei_morphology.csv` (image name, object index, area, roundness)

CSV Format
| Column              | Description                        |
|---------------------|------------------------------------|
| Image               | Filename or series identifier      |
| Nuclei_Count        | Total number of detected nuclei    |
| Object              | Object index (per image)           |
| Area                | Object area in pixels              |
| Roundness           | Object circularity metric          |

Contributing
Contributions welcome! Submit issues or pull requests via GitHub.
