# Nuclei_Stardist_detection – Fiji/Jython Guide

This repository provides a **Fiji (ImageJ) + Jython** workflow for **nuclei detection with StarDist** on microscopy data. It is designed to work with **RAW microscope files** (e.g., **`.lif`**, **`.nd2`**) through **Bio-Formats**.

---

## 1) What it does

- Detects nuclei in a **DAPI** (nuclear) channel using **StarDist 2D**.
- Generates **per-nucleus ROIs** (ROI Manager), **overlay previews**, and **CSV tables** (counts & measurements).
- Works across **multi-series / multi-channel** datasets read by **Bio-Formats**.

---

## 2) Requirements

Install **Fiji** and ensure these plugins/components are present:

- **Bio-Formats** (bundled with Fiji) — *required* to open **`.lif` / `.nd2` / `.tif` / `.ome.tif`*.
- **StarDist 2D** — *required* for nuclei segmentation.
- **CSBDeep** — StarDist dependency (model execution).
- **TensorFlow for ImageJ** — runtime backend used by CSBDeep/StarDist (CPU build is fine).
- **IJPB-plugins** and **MorphoLibJ** — morphology/label operations used in ROI and mask postprocessing.
- **TrackMate** — if you enable **tracking** of nuclei across frames/series (optional for pure segmentation).

Enable/update via the Fiji Updater:  
**Help → Update… → Manage update sites… →** tick **CSBDeep**, **StarDist**, **TensorFlow**, **IJPB-plugins**, **MorphoLibJ**, **TrackMate** → *Close* → *Apply changes* → **restart**.

> Core ImageJ (IJ1) components used and bundled in Fiji: **ROI Manager**, **ResultsTable/Measurements**, **ParticleAnalyzer**, overlays/LUTs, etc.

---

## 3) Supported input formats

- **Direct via Bio-Formats:** **`.lif`** (Leica), **`.nd2`** (Nikon), **`.tif/.tiff`**, **OME‑TIFF**.  
- If a file won’t open, **update Bio‑Formats** in Fiji (*Help → Update…*).

> Proprietary formats are read through Bio‑Formats’ series/reader. If metadata is unusual, exporting to **OME‑TIFF** is a robust fallback.

---

## 4) Quick start (Fiji + Jython)

1. **Open Fiji** → *Plugins → Scripting → New → Jython*.
2. Load the `Nuclei_Stardist_detection` script (paste code or open the `.py` file).
3. When prompted:
   - Select the **input folder** containing `.lif / .nd2 / .tif / .ome.tif` files.
   - Provide **channel keywords** for the **DAPI** channel (see defaults below).
4. Run. Results are written next to the data or to a specified output folder, as overlays and **CSV**.

### Default channel keywords (edit in the script if needed)

- **DAPI / nuclei channel:** `["c1-", "dapi", "blue", "nuclei"]`

These keywords match if any appears in the channel/series title. For **OME‑TIFF** with correct names, selection is usually automatic.

---

## 5) Required imports (from the script)

Typical imports used by the script (Fiji/Jython & Bio‑Formats APIs):

```python
from ij import IJ, WindowManager, ImagePlus
from ij.io import DirectoryChooser, FileSaver
from ij.plugin.frame import RoiManager
from ij.gui import Overlay, PolygonRoi, Roi
from ij.measure import ResultsTable, Measurements
from ij.plugin.filter import ParticleAnalyzer
from ij.process import ImageStatistics as IS, ImageConverter, ByteProcessor
from ij.plugin import RGBStackMerge
from java.io import File, FileWriter, BufferedWriter, PrintWriter
from java.lang import Double
from java.awt import Color
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
import time, re, jarray
```

**StarDist invocation:** via the Fiji command `de.csbdresden.stardist.StarDist2D` (headless). Ensure **CSBDeep**, **TensorFlow**, and **StarDist** update sites are enabled.

---

## 6) Outputs

- **ROI Manager** entries for detected nuclei.  
- **Overlay previews** (PNG/TIFF) with masks/labels.  
- **CSV tables** (counts/measurements) per image or batch.

> Pixel calibration is taken from metadata when present (best with **OME‑TIFF**). Verify scale bars/units before quantitative comparisons.

---

## 7) Tips & troubleshooting

- **StarDist not found** → Enable **CSBDeep** + **StarDist** (and **TensorFlow** backend) and restart Fiji.  
- **LIF/ND2 fails to open** → Update **Bio‑Formats**; if needed, export to **OME‑TIFF**.  
- **Wrong channel picked** → Edit the **DAPI keywords** in the script to match your naming.  
- **Over/under-segmentation** → Adjust StarDist model/threshold parameters (as exposed by the script) and check bit depth.  
- **Empty CSV** → Confirm at least one valid series processed and output path is writable.

---

## 8) Reproducibility

Record your environment:
- Fiji build date
- Bio‑Formats version
- StarDist/CSBDeep/TensorFlow update dates
- IJPB‑plugins / MorphoLibJ, TrackMate versions (if used)

Pinning these alongside the script makes results easier to reproduce.

---





