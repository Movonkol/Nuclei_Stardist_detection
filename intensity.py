

from ij import IJ, WindowManager
from ij.io import DirectoryChooser
from ij.io import FileSaver
from java.io import File, FileWriter, BufferedWriter, PrintWriter
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
import time


# =====================================================================
#                             SETTINGS
# =====================================================================

# Pattern(s) used to identify marker channels.  Multiple patterns can
# be comma separated.  Spaces around patterns are ignored.
marker_patterns_str = IJ.getString(
    "Marker channel pattern (e.g., 'c3-,marker')",
    "c3-,marker"
)
MARKER_CHANNEL_KEYS = [p.strip().lower() for p in marker_patterns_str.split(",") if p.strip()]

# ---------------------------------------------------------------------
# Threshold configuration
# ---------------------------------------------------------------------
# Ask whether to use automatic thresholding (Auto Threshold command) or
# a fixed value.  When automatic thresholding is chosen, the user
# will be prompted for the threshold method and scaling factor.  When
# manual thresholding is selected, only the fixed threshold value is
# requested.
_use_auto_th_str = IJ.getString(
    "Automatic thresholding? (yes/no)",
    "yes"
)
USE_AUTO_THRESHOLD = _use_auto_th_str.lower().strip() in (
    "yes", "y", "true", "t", "1"
)

if USE_AUTO_THRESHOLD:
    # Ask for the algorithm and factor to use when computing the
    # automatic threshold.  A higher factor yields a stricter
    # threshold.
    THRESHOLD_METHOD = IJ.getString(
        "Threshold method (Otsu, Yen, Moments, ...)",
        "Otsu dark"
    )
    THRESHOLD_FACTOR = IJ.getNumber(
        "Threshold factor (>1 = stricter)",
        1.0
    )
    USE_FIXED_THRESHOLD = False
    FIXED_THRESHOLD = 0.0
else:
    # When not using automatic thresholding the method and factor are
    # irrelevant.  Instead we ask for the fixed threshold value.
    THRESHOLD_METHOD = ""
    THRESHOLD_FACTOR = 1.0
    USE_FIXED_THRESHOLD = True
    FIXED_THRESHOLD = IJ.getNumber(
        "Fixed threshold value", 50.0
    )

# ---------------------------------------------------------------------
# Background subtraction behaviour
# ---------------------------------------------------------------------
# Ask whether any background subtraction should be performed at all.
_use_bg_str = IJ.getString(
    "Subtract background? (yes/no)",
    "yes"
)
APPLY_BACKGROUND = _use_bg_str.lower().strip() in (
    "yes", "y", "true", "t", "1"
)

# If background subtraction is enabled, ask whether to prompt for
# custom parameters.  When ``automatic`` is set to false the script
# still performs background subtraction but does not ask the user for
# rolling ball radius, number of repetitions or median filter radius;
# instead it uses sensible defaults defined below.
AUTO_BACKGROUND = False
if APPLY_BACKGROUND:
    _auto_bg_str = IJ.getString(
        "Automatic background subtraction parameters? (yes/no)",
        "yes"
    )
    AUTO_BACKGROUND = _auto_bg_str.lower().strip() in (
        "yes", "y", "true", "t", "1"
    )

    if AUTO_BACKGROUND:
        # When automatic parameter collection is enabled we prompt for
        # each background subtraction parameter.
        ROLLING_RADIUS = IJ.getNumber(
            "Rolling ball radius for background subtraction", 50
        )
        ROLLING_REPEAT = int(IJ.getNumber(
            "How many times to subtract background? (integer)", 1
        ))
        MEDIAN_RADIUS = IJ.getNumber(
            "Median filter radius (0 = none)",
            0
        )
    else:
        # Use default values without prompting when automatic
        # parameter collection is disabled.  You can adjust these
        # defaults here if desired.
        ROLLING_RADIUS = 50
        ROLLING_REPEAT = 1
        MEDIAN_RADIUS = 0
else:
    # When no background subtraction is requested we set all
    # background related parameters to zero.  They will not be used
    # later in the script.
    ROLLING_RADIUS = 0
    ROLLING_REPEAT = 0
    MEDIAN_RADIUS = 0

# Note: heatmap configuration (whether to save pseudocolour heatmaps) is
# now prompted after selecting the image folder.  See below.



# =====================================================================
#                  NAME CLEANING HELPER FUNCTION
# =====================================================================
def clean_name(name):
    """Return the part of a comma‑separated image name after the last comma.

    BioFormats sometimes returns names containing multiple comma
    separated values; this helper extracts the last component and
    trims any whitespace.
    """
    if name is None:
        return ""
    parts = name.split(',')
    return parts[-1].strip()


# =====================================================================
#                         EXPORT FUNCTION
# =====================================================================
def export_intensity(
    image_name, channel_name, pos_pixels, mean_int, integrated_int,
    nuclei_count, normalized_intensity, csv_path
):
    """Append a single row of measurements to a CSV file.

    Parameters
    ----------
    image_name : str
        The name of the image or series being processed.
    channel_name : str
        The name of the marker channel being measured.
    pos_pixels : int
        Number of pixels above the effective threshold.
    mean_int : float
        Mean intensity of the positive pixels.
    integrated_int : float
        Sum of the intensities of all positive pixels.
    nuclei_count : int
        Number of nuclei in the corresponding image or series.
    normalized_intensity : float
        Integrated intensity normalised by the nucleus count.
    csv_path : str
        Full path to the CSV file to which the row should be appended.
    """
    f = File(csv_path)
    first = not f.exists()
    pw = PrintWriter(BufferedWriter(FileWriter(f, True)))
    if first:
        pw.println("Image,Channel,Positive_Pixels,Mean_Intensity,Integrated_Intensity,Nuclei_Count,Normalized_Intensity")
    pw.println("%s,%s,%d,%.3f,%.3f,%d,%.6f" % (
        image_name, channel_name, pos_pixels, mean_int,
        integrated_int, nuclei_count, normalized_intensity
    ))
    pw.close()


# =====================================================================
#               SELECT FOLDER AND LOAD NUCLEI COUNTS
# =====================================================================
dc = DirectoryChooser("Choose folder with images")
folder = dc.getDirectory()
if not folder:
    IJ.error("No folder selected")
    raise SystemExit

# ---------------------------------------------------------------------
# Heatmap configuration
# ---------------------------------------------------------------------
# Ask whether to save pseudocolour heatmaps once the folder has been
# selected.  This ensures the prompt is shown even if the script
# terminates early (e.g. due to no folder selection).  Heatmaps will
# be stored in a "heatmaps" subfolder within the chosen folder.
_use_heatmap_str = IJ.getString(
    "Save heatmap images? (yes/no)",
    "no"
)
USE_HEATMAP = _use_heatmap_str.lower().strip() in (
    "yes", "y", "true", "t", "1"
)

csv_file = folder + File.separator + "intensity_measurements.csv"
well_csv_file = folder + File.separator + "intensity_measurements_per_well.csv"

counts_filename = IJ.getString(
    "Name of nuclei count CSV (leave empty to skip normalization)",
    "nuclei_counts.csv",
)
nuclei_counts_map = {}
if counts_filename:
    counts_path = folder + File.separator + counts_filename
    try:
        f_counts = File(counts_path)
        if f_counts.exists():
            with open(counts_path, 'r') as fp:
                for line in fp:
                    line = line.strip()
                    if not line or line.lower().startswith("image"):
                        continue
                    parts = line.split(',')
                    if len(parts) >= 3:
                        img_name_raw = parts[-2].strip()
                        count_str = parts[-1].strip()
                    else:
                        img_name_raw = parts[0].strip()
                        count_str = parts[1].strip()
                    clean_img = clean_name(img_name_raw)
                    try:
                        count_val = int(float(count_str))
                    except Exception:
                        count_val = 0
                    nuclei_counts_map[clean_img] = count_val
        else:
            IJ.log("Note: File '%s' not found – skipping normalization." % counts_filename)
    except Exception as err:
        IJ.log("Error loading nuclei count CSV: %s" % str(err))

# Prepare heatmap directory if required.  The directory is created
# once at the beginning rather than on every image to avoid repeated
# filesystem checks.
heatmap_dir = None
if USE_HEATMAP:
    heatmap_dir = File(folder + File.separator + "heatmaps")
    if not heatmap_dir.exists():
        heatmap_dir.mkdir()


# =====================================================================
#          MAIN LOOP: READ FILES AND ANALYZE INTENSITY
# =====================================================================
for f in File(folder).listFiles():
    if not f.isFile():
        continue
    nm = f.getName().lower()
    # Skip non-image files
    if not nm.endswith((".tif", ".tiff", ".png", ".jpg", ".lif", ".nd2")):
        continue

    total_positive = 0
    total_integrated = 0.0
    total_nuclei = 0

    opts = ImporterOptions()
    opts.setId(f.getAbsolutePath())
    opts.setOpenAllSeries(True)
    opts.setVirtual(True)
    imps = BF.openImagePlus(opts)

    for idx, imp in enumerate(imps):
        # A unique name is created for each series in a multi‑series file.  The
        # ``clean_name`` helper is used to extract a simple, consistent key
        # from this name for matching nuclei counts.
        series_name_raw = "%s_Serie%d" % (f.getName(), idx + 1)
        series_name_clean = clean_name(series_name_raw)

        imp.setTitle(series_name_raw)
        imp.show()

        IJ.run(imp, "Split Channels", "")
        time.sleep(0.5)

        marker_channels = []
        for i in range(1, WindowManager.getImageCount() + 1):
            win = WindowManager.getImage(i)
            if win is None:
                continue
            title_lower = win.getTitle().lower()
            for mkey in MARKER_CHANNEL_KEYS:
                if mkey in title_lower:
                    marker_channels.append(win)
                    break
        if len(marker_channels) == 0:
            IJ.error("No marker channels found in %s" % series_name_raw)
            WindowManager.closeAllWindows()
            continue

        for marker in marker_channels:
            channel_name = marker.getTitle().replace(" ", "_")
            marker_analysis = marker.duplicate()
            marker_analysis.show()

            IJ.run(marker_analysis, "8-bit", "")
            IJ.run(marker_analysis, "Enhance Contrast", "saturated=0.35")

            # Apply optional background subtraction.  The values for
            # ``ROLLING_RADIUS``, ``ROLLING_REPEAT`` and ``MEDIAN_RADIUS``
            # depend on the earlier user input.
            if APPLY_BACKGROUND:
                if MEDIAN_RADIUS > 0:
                    IJ.run(marker_analysis, "Median...", "radius=%d" % int(MEDIAN_RADIUS))
                for _ in range(max(1, ROLLING_REPEAT)):
                    IJ.run(marker_analysis, "Subtract Background...", "rolling=%d" % int(ROLLING_RADIUS))

            # Determine the threshold for positive pixels.  Use a fixed
            # threshold if requested, otherwise compute one based on the
            # selected method and scale it by ``THRESHOLD_FACTOR``.
            if bool(USE_FIXED_THRESHOLD):
                threshold_value = float(FIXED_THRESHOLD)
            else:
                IJ.setAutoThreshold(marker_analysis, THRESHOLD_METHOD)
                ip_tmp = marker_analysis.getProcessor()
                threshold_value = float(ip_tmp.getMinThreshold())
                # Handle NaN threshold (e.g. when all pixels are equal)
                if threshold_value != threshold_value:
                    threshold_value = 0.0
            effective_threshold = threshold_value * float(THRESHOLD_FACTOR)

            ip = marker_analysis.getProcessor()
            width, height = ip.getWidth(), ip.getHeight()
            total_intensity_ser = 0.0
            positive_count_ser = 0
            for y in range(height):
                for x in range(width):
                    v = ip.get(x, y)
                    if v > effective_threshold:
                        total_intensity_ser += v
                        positive_count_ser += 1

            mean_int = (total_intensity_ser / positive_count_ser) if positive_count_ser > 0 else 0.0
            integrated_int = total_intensity_ser

            nuclei_count = nuclei_counts_map.get(series_name_clean, 0)

            normalized_intensity = (
                integrated_int / float(nuclei_count)
                if nuclei_count > 0
                else 0.0
            )

            export_intensity(
                series_name_raw, channel_name, positive_count_ser, mean_int,
                integrated_int, nuclei_count, normalized_intensity, csv_file
            )

            total_positive += positive_count_ser
            total_integrated += integrated_int
            total_nuclei += nuclei_count

            # Mark the processed copy as unmodified so ImageJ does not
            # prompt to save changes when closing it.
            marker_analysis.changes = False

            # -----------------------------------------------------------------
            # Save heatmap image if requested before closing the original
            # marker.  A duplicate of the original marker channel is
            # converted to a "Fire" lookup table and written to disk.  The
            # name of the file encodes the series and channel.
            if USE_HEATMAP and heatmap_dir is not None:
                # Duplicate the original marker channel (still open) and
                # apply the same preprocessing steps as ``marker_analysis``.
                heatmap_imp = marker.duplicate()
                # Convert to 8‑bit
                IJ.run(heatmap_imp, "8-bit", "")
                # Enhance contrast with the same saturation level
                IJ.run(heatmap_imp, "Enhance Contrast", "saturated=0.35")
                # Apply background subtraction if enabled using the
                # previously defined parameters
                if APPLY_BACKGROUND:
                    if MEDIAN_RADIUS > 0:
                        IJ.run(heatmap_imp, "Median...", "radius=%d" % int(MEDIAN_RADIUS))
                    for _ in range(max(1, ROLLING_REPEAT)):
                        IJ.run(heatmap_imp, "Subtract Background...", "rolling=%d" % int(ROLLING_RADIUS))
                # Apply a pseudocolour LUT; "Fire" is widely used for
                # heatmaps.  Alternative LUTs can be applied here.
                IJ.run(heatmap_imp, "Fire", "")
                # Construct a filename that uniquely identifies the
                # series and channel.  Replace path separators and
                # spaces to avoid filesystem issues.
                safe_series = series_name_clean.replace(File.separator, "_").replace(" ", "_")
                safe_channel = channel_name.replace(File.separator, "_").replace(" ", "_")
                heatmap_name = "%s_%s_heatmap.png" % (safe_series, safe_channel)
                heatmap_path = heatmap_dir.getAbsolutePath() + File.separator + heatmap_name
                # Save as PNG using FileSaver to avoid any interactive
                # save dialogues.  This writes directly to the given path.
                fs = FileSaver(heatmap_imp)
                fs.saveAsPng(heatmap_path)
                # Mark the heatmap as unmodified and close the window to
                # prevent ImageJ from asking to save on exit.
                heatmap_imp.changes = False
                heatmap_imp.close()

            # Close processed and original marker images cleanly.  The
            # ``changes`` flag must be reset on each to avoid save prompts.
            marker_analysis.close()
            marker.changes = False
            marker.close()

        WindowManager.closeAllWindows()

    # After processing all series in the file, compute mean and
    # normalised intensities for the whole well.  ``total_positive`` is
    # the sum of positive pixels across all series, while
    # ``total_integrated`` is the sum of integrated intensities.  If
    # nuclei counts are available they are used to compute the per
    # nucleus intensity.
    well_mean_int = (total_integrated / total_positive) if total_positive > 0 else 0.0
    well_norm_int = (total_integrated / total_nuclei) if total_nuclei > 0 else 0.0

    export_intensity(
        f.getName(), "AllSeries", total_positive, well_mean_int,
        total_integrated, total_nuclei, well_norm_int, well_csv_file
    )
