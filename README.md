# StarDist Batch Analysis & Intensity Measurement

This repository provides three FIJI/ImageJ scripts for high‑throughput analysis of fluorescence microscopy images:

- **StarDist Batch Analysis**: Automated nuclei segmentation and morphology measurement using StarDist.  
- **Intensity Measurement**: Quantitative per‑well and per‑series intensity analysis, with optional background subtraction, thresholding, and heatmap generation.  
- **Positive/Negative Classification** *(optional)*: Classify wells based on normalized intensity thresholds.

---

## Table of Contents

- [Features](#features)  
- [Installation](#installation)  
- [Configuration](#configuration)  
- [Usage](#usage)  
  - [StarDist Batch Analysis](#stardist-batch-analysis)  
  - [Intensity Measurement](#intensity-measurement)  
  - [Positive/Negative Classification](#positivenegative-classification)  
- [Outputs & CSV Formats](#outputs--csv-formats)  
- [Contributing](#contributing)  

---

## Features

### StarDist Batch Analysis

- **Automatic parameter prompts** for model, probability threshold, NMS threshold, tiling, channel patterns, and size filters.  
- **Automatic channel detection** (e.g., DAPI) via customizable patterns.  
- **Size filtering** to exclude artifacts.  
- **Headless StarDist segmentation** with Bio‑Formats support.  
- **Morphology measurements** (area, roundness) per object.  
- **Batch processing** of multi‑series image files (ND2, LIF, TIFF).  
- **CSV export** of nuclei counts and per‑object morphology.

### Intensity Measurement (`intensity.py`)

- **Interactive prompts** for marker channels, thresholding method (Otsu, Yen, Moments, or fixed), and background subtraction parameters.  
- **Optional background subtraction** with rolling‑ball radius and median filtering.  
- **Per‑series and per‑well intensity metrics** (positive pixels, mean, integrated, normalized intensity).  
- **Optional heatmap generation** saved as PNGs.  
- **Batch processing** with Bio‑Formats multi‑series support.

### Positive/Negative Classification (`positiv_negativ.py`)

- **Threshold‑based classification** of wells/images by normalized intensity.  
- **Simple workflow**: prompts for threshold and input folder.  
- **Generates** `well_classification.csv`.

---

## Installation

1. Clone this repository:  
   ```bash
   git clone https://github.com/yourusername/stardist-batch.git
   cd stardist-batch
