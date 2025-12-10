#!/usr/bin/env python3
"""
StarDist Nuclei Detection for VS Code/Python
Converted from ImageJ/Jython to standard Python

This script performs nuclei detection using StarDist on microscopy images.
"""

import os
import re
import csv
import argparse
from pathlib import Path
from typing import List, Tuple, Optional
import warnings

import numpy as np
from skimage import io, measure, exposure, transform
from skimage.morphology import remove_small_objects
from skimage.color import label2rgb
from scipy import ndimage
from stardist.models import StarDist2D
from csbdeep.utils import normalize
import matplotlib.pyplot as plt
from matplotlib import patches
from PIL import Image
import tifffile


class NucleiDetectionConfig:
    """Configuration for nuclei detection pipeline"""

    def __init__(self):
        # StarDist settings
        self.model_name = "2D_versatile_fluo"  # StarDist model
        self.prob_thresh = 0.5
        self.nms_thresh = 0.5
        self.n_tiles = (1, 1)

        # Channel detection patterns
        self.channel_patterns = ['c1-', 'dapi', 'blue']

        # Size filtering
        self.size_min = 50  # px^2
        self.size_max = 2500  # px^2

        # Pre-scaling
        self.scale_factor = 1.0

        # Intensity filtering
        self.sigma_above_bg = 1.5  # standard deviations above background
        self.min_mean_intensity = 0  # minimum mean intensity (0 = off)

        # QC overlay settings
        self.save_qc_overlays = True
        self.qc_sat_pct = 1.0  # contrast saturation percentage
        self.qc_stroke_width = 2
        self.qc_color = (255, 0, 0)  # RGB red


def safe_filename(filename: str, max_length: int = 180) -> str:
    """Remove invalid characters from filename"""
    cleaned = re.sub(r'[\\/:*?"<>|]+', '_', filename)
    return cleaned[:max_length]


def find_dapi_channel(image_path: Path, patterns: List[str]) -> Optional[np.ndarray]:
    """
    Load and identify DAPI channel from multi-channel image

    Args:
        image_path: Path to image file
        patterns: List of patterns to identify DAPI channel

    Returns:
        DAPI channel as 2D numpy array
    """
    # Load image
    try:
        # Try tifffile first for .tif/.tiff
        if image_path.suffix.lower() in ['.tif', '.tiff']:
            img = tifffile.imread(str(image_path))
        else:
            img = io.imread(str(image_path))
    except Exception as e:
        print(f"Error loading {image_path}: {e}")
        return None

    # Handle different image dimensions
    if img.ndim == 2:
        # Single channel grayscale
        return img
    elif img.ndim == 3:
        # Check if RGB or multi-channel
        if img.shape[2] == 3:
            # RGB image - convert to grayscale or take blue channel
            return img[:, :, 2]  # Blue channel for DAPI
        else:
            # Multi-channel - take first channel
            return img[:, :, 0]
    elif img.ndim == 4:
        # Multi-dimensional (C, Z, Y, X) or similar
        # Take first channel, first Z-slice
        return img[0, 0, :, :]

    return img


def compute_background_stats(dapi_img: np.ndarray, mask: np.ndarray) -> Tuple[float, float]:
    """
    Compute background mean and standard deviation outside nuclei

    Args:
        dapi_img: DAPI channel image
        mask: Binary mask of nuclei (True = nuclei, False = background)

    Returns:
        (bg_mean, bg_std)
    """
    background_pixels = dapi_img[~mask]
    if len(background_pixels) == 0:
        # Fallback if no background
        return np.mean(dapi_img), np.std(dapi_img)
    return np.mean(background_pixels), np.std(background_pixels)


def filter_labels_by_intensity(
    labels: np.ndarray,
    dapi_img: np.ndarray,
    config: NucleiDetectionConfig
) -> Tuple[np.ndarray, List[dict]]:
    """
    Filter detected nuclei by intensity criteria

    Args:
        labels: Label image from StarDist
        dapi_img: Original DAPI image
        config: Detection configuration

    Returns:
        (filtered_labels, properties_list)
    """
    # Get region properties
    regions = measure.regionprops(labels, intensity_image=dapi_img)

    # Create binary mask of all nuclei for background calculation
    nuclei_mask = labels > 0
    bg_mean, bg_std = compute_background_stats(dapi_img, nuclei_mask)

    # Filter regions
    filtered_labels = np.zeros_like(labels)
    properties = []
    new_label = 1

    for region in regions:
        area = region.area
        mean_intensity = region.mean_intensity

        # Size filter
        if area < config.size_min or area > config.size_max:
            continue

        # Intensity filters
        pass_abs = (config.min_mean_intensity <= 0) or (mean_intensity >= config.min_mean_intensity)
        pass_sigma = (config.sigma_above_bg <= 0) or (mean_intensity >= bg_mean + config.sigma_above_bg * bg_std)

        if pass_abs and pass_sigma:
            # Keep this nucleus
            mask = labels == region.label
            filtered_labels[mask] = new_label

            # Store properties
            properties.append({
                'label': new_label,
                'area': area,
                'mean_intensity': mean_intensity,
                'circularity': 4 * np.pi * area / (region.perimeter ** 2) if region.perimeter > 0 else 0,
                'centroid': region.centroid
            })
            new_label += 1

    return filtered_labels, properties


def create_qc_overlay(
    dapi_img: np.ndarray,
    labels: np.ndarray,
    config: NucleiDetectionConfig,
    output_path: Path
):
    """
    Create and save QC overlay image with detected nuclei contours

    Args:
        dapi_img: DAPI channel image
        labels: Filtered label image
        config: Detection configuration
        output_path: Path to save overlay image
    """
    # Normalize DAPI for display
    dapi_norm = exposure.rescale_intensity(dapi_img, out_range=(0, 255)).astype(np.uint8)

    # Apply contrast enhancement if requested
    if config.qc_sat_pct > 0:
        p_low, p_high = np.percentile(dapi_norm, [config.qc_sat_pct, 100 - config.qc_sat_pct])
        dapi_norm = exposure.rescale_intensity(dapi_norm, in_range=(p_low, p_high))

    # Create RGB image with blue DAPI
    rgb_img = np.zeros((*dapi_norm.shape, 3), dtype=np.uint8)
    rgb_img[:, :, 2] = dapi_norm  # Blue channel

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 12), dpi=150)
    ax.imshow(rgb_img)

    # Draw contours
    for region in measure.regionprops(labels):
        # Get contour
        contours = measure.find_contours(labels == region.label, 0.5)
        for contour in contours:
            ax.plot(contour[:, 1], contour[:, 0],
                   color=tuple(c/255 for c in config.qc_color),
                   linewidth=config.qc_stroke_width)

    # Add legend
    ax.text(5, 20, 'QC: detected nuclei\nContours in RED',
           color='white', fontsize=10, va='top',
           bbox=dict(boxstyle='round', facecolor='black', alpha=0.5))

    ax.axis('off')
    plt.tight_layout()

    # Save
    plt.savefig(output_path, dpi=150, bbox_inches='tight', pad_inches=0.1)
    plt.close()
    print(f"Saved QC overlay: {output_path}")


def export_counts(csv_path: Path, image_name: str, count: int):
    """Append nuclei count to CSV"""
    file_exists = csv_path.exists()

    with open(csv_path, 'a', newline='') as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(['Image', 'Total_Nuclei_Count'])
        writer.writerow([image_name, count])


def export_morphology(csv_path: Path, image_name: str, mean_area: float, mean_circularity: float):
    """Append morphology metrics to CSV"""
    file_exists = csv_path.exists()

    with open(csv_path, 'a', newline='') as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(['Image', 'Mean_Area', 'Mean_Roundness'])
        writer.writerow([image_name, f"{mean_area:.2f}", f"{mean_circularity:.4f}"])


def process_image(
    image_path: Path,
    model: StarDist2D,
    config: NucleiDetectionConfig,
    output_dir: Path,
    qc_dir: Optional[Path] = None
) -> Tuple[int, float, float]:
    """
    Process a single image for nuclei detection

    Args:
        image_path: Path to input image
        model: StarDist2D model
        config: Detection configuration
        output_dir: Output directory
        qc_dir: QC overlay directory (optional)

    Returns:
        (count, mean_area, mean_circularity)
    """
    print(f"→ Processing: {image_path.name}")

    # Load DAPI channel
    dapi_img = find_dapi_channel(image_path, config.channel_patterns)
    if dapi_img is None:
        print(f"Error: Could not load DAPI channel from {image_path.name}")
        return 0, 0.0, 0.0

    # Ensure float type for processing
    dapi_img = dapi_img.astype(np.float32)
    orig_shape = dapi_img.shape

    # Pre-scaling if needed
    if abs(config.scale_factor - 1.0) > 1e-6:
        new_shape = tuple(int(round(s * config.scale_factor)) for s in dapi_img.shape)
        dapi_scaled = transform.resize(dapi_img, new_shape, order=1, preserve_range=True, anti_aliasing=True)
    else:
        dapi_scaled = dapi_img

    # Normalize for StarDist
    dapi_norm = normalize(dapi_scaled, 0, 100)

    # Run StarDist prediction
    labels, _ = model.predict_instances(
        dapi_norm,
        prob_thresh=config.prob_thresh,
        nms_thresh=config.nms_thresh,
        n_tiles=config.n_tiles
    )

    # Scale labels back to original size if needed
    if abs(config.scale_factor - 1.0) > 1e-6:
        labels = transform.resize(labels, orig_shape, order=0, preserve_range=True, anti_aliasing=False).astype(np.int32)

    # Filter by intensity and size
    filtered_labels, properties = filter_labels_by_intensity(labels, dapi_img, config)

    # Compute statistics
    count = len(properties)
    if count > 0:
        mean_area = np.mean([p['area'] for p in properties])
        mean_circularity = np.mean([p['circularity'] for p in properties])
    else:
        mean_area = 0.0
        mean_circularity = 0.0

    print(f"→ {image_path.name}: Count={count}")

    # Create QC overlay
    if config.save_qc_overlays and qc_dir is not None and count > 0:
        qc_filename = safe_filename(f"{image_path.stem}__QC_DAPI_Blue.png")
        qc_path = qc_dir / qc_filename
        create_qc_overlay(dapi_img, filtered_labels, config, qc_path)

    return count, mean_area, mean_circularity


def main():
    """Main processing pipeline"""
    parser = argparse.ArgumentParser(
        description='StarDist nuclei detection for microscopy images'
    )
    parser.add_argument(
        'input_dir',
        type=Path,
        help='Directory containing input images'
    )
    parser.add_argument(
        '--model',
        type=str,
        default='2D_versatile_fluo',
        help='StarDist model name (default: 2D_versatile_fluo)'
    )
    parser.add_argument(
        '--prob-thresh',
        type=float,
        default=0.5,
        help='Probability threshold (default: 0.5)'
    )
    parser.add_argument(
        '--nms-thresh',
        type=float,
        default=0.5,
        help='NMS threshold (default: 0.5)'
    )
    parser.add_argument(
        '--size-min',
        type=float,
        default=50,
        help='Minimum nucleus area in px^2 (default: 50)'
    )
    parser.add_argument(
        '--size-max',
        type=float,
        default=2500,
        help='Maximum nucleus area in px^2 (default: 2500)'
    )
    parser.add_argument(
        '--sigma-above-bg',
        type=float,
        default=1.5,
        help='Minimum σ above background (default: 1.5, 0=off)'
    )
    parser.add_argument(
        '--min-intensity',
        type=float,
        default=0,
        help='Minimum mean intensity (default: 0=off)'
    )
    parser.add_argument(
        '--scale-factor',
        type=float,
        default=1.0,
        help='Pre-scaling factor (0.5=down, 2.0=up, 1.0=off, default: 1.0)'
    )
    parser.add_argument(
        '--no-qc',
        action='store_true',
        help='Disable QC overlay generation'
    )

    args = parser.parse_args()

    # Validate input directory
    if not args.input_dir.exists():
        print(f"Error: Input directory '{args.input_dir}' does not exist")
        return 1

    # Setup configuration
    config = NucleiDetectionConfig()
    config.model_name = args.model
    config.prob_thresh = args.prob_thresh
    config.nms_thresh = args.nms_thresh
    config.size_min = args.size_min
    config.size_max = args.size_max
    config.sigma_above_bg = args.sigma_above_bg
    config.min_mean_intensity = args.min_intensity
    config.scale_factor = args.scale_factor
    config.save_qc_overlays = not args.no_qc

    # Load StarDist model
    print(f"Loading StarDist model: {config.model_name}")
    try:
        model = StarDist2D.from_pretrained(config.model_name)
    except Exception as e:
        print(f"Error loading model: {e}")
        return 1

    # Setup output directories
    output_dir = args.input_dir
    csv_counts = output_dir / "nuclei_counts.csv"
    csv_morphology = output_dir / "nuclei_morphology.csv"

    qc_dir = None
    if config.save_qc_overlays:
        qc_dir = output_dir / "QC_Overlays"
        qc_dir.mkdir(exist_ok=True)

    # Find image files
    image_extensions = ['.tif', '.tiff', '.png', '.jpg', '.jpeg']
    image_files = [f for f in args.input_dir.iterdir()
                   if f.is_file() and f.suffix.lower() in image_extensions]

    if not image_files:
        print(f"No image files found in {args.input_dir}")
        return 1

    print(f"Found {len(image_files)} image(s) to process")

    # Process each image
    grand_total = 0
    for image_path in image_files:
        try:
            count, mean_area, mean_circ = process_image(
                image_path, model, config, output_dir, qc_dir
            )

            # Export results
            export_counts(csv_counts, image_path.name, count)
            export_morphology(csv_morphology, image_path.name, mean_area, mean_circ)

            grand_total += count

        except Exception as e:
            print(f"Error processing {image_path.name}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Export grand total
    export_counts(csv_counts, "TOTAL_ALL_IMAGES", grand_total)

    # Summary
    print("\n" + "="*60)
    print(f"Processing complete!")
    print(f"Total nuclei detected: {grand_total}")
    print(f"Results saved to: {csv_counts}")
    print(f"Morphology saved to: {csv_morphology}")
    if config.save_qc_overlays:
        print(f"QC overlays saved to: {qc_dir}")
    print("="*60)

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
