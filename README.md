# Nuclei‑StarDist

Batch‑process fluorescence image series in headless Fiji/ImageJ using StarDist 2D.

**Key Steps:**

* **Bio‑Formats Import:** Opens multi‑series image files (TIFF, PNG, JPG, LIF, ND2, etc.)
* **Channel Splitting & DAPI Detection:** Splits channels and identifies nuclei channel by user‑defined keywords
* **StarDist Segmentation:** Runs headless StarDist2D with configurable model, probability & NMS thresholds, and tiling
* **Size Filtering:** Applies minimum/maximum size bounds and connectivity criteria to segment labels
* **Counting & Export:** Counts nuclei per series and appends results to `nuclei_counts.csv`

## Installation

Clone or copy the script into your Fiji macros folder:

```bash
git clone https://github.com/yourusername/nuclei-stardist.git
cp nuclei-stardist/Nuclei_stardist.py ${Fiji.app}/macros/
```

Ensure you have:

* Fiji/ImageJ with CSBDeep StarDist 2D plugin
* Bio‑Formats importer
* Java 8+

## Configuration

Launch the macro to set parameters via dialog boxes:

| Parameter                    | Description                                                               | Default                        |
| ---------------------------- | ------------------------------------------------------------------------- | ------------------------------ |
| **StarDist Model**           | Pretrained model name (e.g. `Versatile (fluorescent nuclei)`)             | Versatile (fluorescent nuclei) |
| **Probability Threshold**    | Minimum probability for detection                                         | 0.5                            |
| **NMS Threshold**            | Non‑maximum suppression threshold                                         | 0.5                            |
| **Number of Tiles**          | Tiling layout (`rows,cols`) for large images                              | 1,1                            |
| **Min Cell Size**            | Lower bound for label size filtering                                      | 50                             |
| **Max Cell Size**            | Upper bound for label size filtering                                      | 2500                           |
| **Connectivity**             | Pixel connectivity for connected components (4 or 8)                      | 4                              |
| **DAPI channel identifiers** | Comma‑separated keywords to match nuclei channel titles (e.g. `c1-,dapi`) | c1-,dapi                       |

## Usage

1. In Fiji, open **Plugins ▸ New ▸ Script…**, load `Nuclei_stardist.py` and click **Run**.
2. Configure settings when prompted.
3. Select the folder containing your images.
4. The script will process each series and produce `nuclei_counts.csv` in the same folder.

## Output

* **nuclei\_counts.csv**: CSV table with columns:

  * `Image`: series identifier (e.g. `sample.tif_Series1`)
  * `Nuclei_Count`: detected nuclei count

## License

MIT © Your Name
