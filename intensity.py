from ij import IJ, WindowManager
from ij.io import DirectoryChooser, FileSaver
from java.io import File, FileWriter, BufferedWriter, PrintWriter
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
import time

# =====================================================================
#                             SETTINGS
# =====================================================================

# Marker channel patterns (e.g., 'c3-,marker')
marker_patterns_str = IJ.getString(
    "Marker channel pattern (e.g., 'c3-,marker')",
    "c3-,marker"
)
MARKER_CHANNEL_KEYS = [p.strip().lower() for p in marker_patterns_str.split(",") if p.strip()]

# ---------------------------------------------------------------------
# Threshold configuration
# ---------------------------------------------------------------------
use_auto_th_str = IJ.getString(
    "Use automatic thresholding? (yes/no)",
    "yes"
)
USE_AUTO_THRESHOLD = use_auto_th_str.lower().strip() in (
    "yes", "y", "true", "t", "1"
)

if USE_AUTO_THRESHOLD:
    THRESHOLD_METHOD = IJ.getString(
        "Threshold method (Otsu, Yen, Moments, etc.)",
        "Otsu dark"
    )
    THRESHOLD_FACTOR = IJ.getNumber(
        "Threshold factor (>1 = stricter)",
        1.0
    )
    USE_FIXED_THRESHOLD = False
    FIXED_THRESHOLD = 0.0
else:
    THRESHOLD_METHOD = ""
    THRESHOLD_FACTOR = 1.0
    USE_FIXED_THRESHOLD = True
    FIXED_THRESHOLD = IJ.getNumber(
        "Fixed threshold value", 50.0
    )

# ---------------------------------------------------------------------
# Background subtraction settings
# ---------------------------------------------------------------------
apply_bg_str = IJ.getString(
    "Apply background subtraction? (yes/no)",
    "yes"
)
APPLY_BACKGROUND = apply_bg_str.lower().strip() in (
    "yes", "y", "true", "t", "1"
)
AUTO_BACKGROUND = False
if APPLY_BACKGROUND:
    auto_bg_str = IJ.getString(
        "Use custom background parameters? (yes/no)",
        "yes"
    )
    AUTO_BACKGROUND = auto_bg_str.lower().strip() in (
        "yes", "y", "true", "t", "1"
    )
    if AUTO_BACKGROUND:
        ROLLING_RADIUS = IJ.getNumber(
            "Rolling ball radius for background subtraction", 50
        )
        ROLLING_REPEAT = int(IJ.getNumber(
            "Number of background subtraction repetitions (integer)", 1
        ))
        MEDIAN_RADIUS = IJ.getNumber(
            "Median filter radius (0 = none)", 0
        )
    else:
        ROLLING_RADIUS = 50
        ROLLING_REPEAT = 1
        MEDIAN_RADIUS = 0
else:
    ROLLING_RADIUS = 0
    ROLLING_REPEAT = 0
    MEDIAN_RADIUS = 0

# =====================================================================
#          EXPORT FUNCTION
# =====================================================================
def export_intensity(
    image_name, channel_name, positive_pixels, mean_intensity,
    integrated_intensity, nuclei_count, normalized_intensity, csv_path
):
    """
    Append a row of measurement data to a CSV file.
    """
    f = File(csv_path)
    first = not f.exists()
    pw = PrintWriter(BufferedWriter(FileWriter(f, True)))
    if first:
        pw.println(
            "Image,Channel,Positive_Pixels,Mean_Intensity,Integrated_Intensity,Nuclei_Count,Normalized_Intensity"
        )
    pw.println(
        "%s,%s,%d,%.3f,%.3f,%d,%.6f" % (
            image_name, channel_name, positive_pixels,
            mean_intensity, integrated_intensity,
            nuclei_count, normalized_intensity
        )
    )
    pw.close()

# =====================================================================
#  SELECT IMAGE FOLDER AND READ NUCLEI COUNTS (BY ORDER)
# =====================================================================
dc = DirectoryChooser("Select folder with images")
folder = dc.getDirectory()
if not folder:
    IJ.error("No folder selected")
    raise SystemExit

# CSV with nuclei counts (leave blank to skip normalization)
counts_filename = IJ.getString(
    "Name of nuclei count CSV (leave empty to skip normalization)",
    "nuclei_counts.csv"
)
counts_list = []
if counts_filename:
    counts_path = folder + File.separator + counts_filename
    try:
        f_counts = File(counts_path)
        if f_counts.exists():
            with open(counts_path, 'r') as fp:
                for line in fp:
                    line = line.strip()
                    if not line or line.lower().startswith("image"): continue
                    parts = line.split(',')
                    try:
                        counts_list.append(int(float(parts[-1].strip())))
                    except:
                        counts_list.append(0)
        else:
            IJ.log(
                "Note: '%s' not found – skipping normalization." % counts_filename
            )
    except Exception as err:
        IJ.log("Error reading nuclei count CSV: %s" % str(err))

# Heatmap option
use_heatmap = IJ.getString(
    "Save heatmap images? (yes/no)",
    "no"
).lower().strip() in ("yes","y","true","t","1")
heatmap_dir = None
if use_heatmap:
    heatmap_dir = File(folder + File.separator + "heatmaps")
    if not heatmap_dir.exists():
        heatmap_dir.mkdir()

# Output CSV paths
csv_file = folder + File.separator + "intensity_measurements.csv"
well_csv_file = folder + File.separator + "intensity_measurements_per_well.csv"

# Get sorted list of image files
def is_image(file):
    name = file.getName().lower()
    return file.isFile() and name.endswith(
        (".tif", ".tiff", ".png", ".jpg", ".lif", ".nd2")
    )
files = sorted(
    [f for f in File(folder).listFiles() if is_image(f)],
    key=lambda f: f.getName().lower()
)

# =====================================================================
#           MAIN LOOP: INTENSITY ANALYSIS
# =====================================================================
for idx, f in enumerate(files):
    # assign nuclei count by file order
    nuclei_count = counts_list[idx] if idx < len(counts_list) else 0

    total_positive = 0
    total_integrated = 0.0
    total_nuclei = nuclei_count

    opts = ImporterOptions()
    opts.setId(f.getAbsolutePath())
    opts.setOpenAllSeries(True)
    opts.setVirtual(True)
    imps = BF.openImagePlus(opts)

    for sidx, imp in enumerate(imps):
        series_name = "%s_Series%d" % (f.getName(), sidx+1)
        imp.setTitle(series_name)
        imp.show()

        IJ.run(imp, "Split Channels", "")
        time.sleep(0.5)

        # find marker channels
        marker_channels = []
        for i in range(1, WindowManager.getImageCount()+1):
            win = WindowManager.getImage(i)
            if not win: continue
            title = win.getTitle().lower()
            if any(key in title for key in MARKER_CHANNEL_KEYS):
                marker_channels.append(win)

        for marker in marker_channels:
            channel = marker.duplicate()
            channel.show()
            IJ.run(channel, "8-bit", "")
            IJ.run(channel, "Enhance Contrast", "saturated=0.35")

            if APPLY_BACKGROUND:
                if MEDIAN_RADIUS > 0:
                    IJ.run(channel, "Median...", "radius=%d" % int(MEDIAN_RADIUS))
                for _ in range(max(1, ROLLING_REPEAT)):
                    IJ.run(channel, "Subtract Background...", "rolling=%d" % int(ROLLING_RADIUS))

            # threshold
            if USE_FIXED_THRESHOLD:
                thresh = float(FIXED_THRESHOLD)
            else:
                IJ.setAutoThreshold(channel, THRESHOLD_METHOD)
                ip = channel.getProcessor()
                thresh = float(ip.getMinThreshold())
                if thresh != thresh:
                    thresh = 0.0
            effective_thresh = thresh * float(THRESHOLD_FACTOR)

            # measure
            ip = channel.getProcessor()
            w, h = ip.getWidth(), ip.getHeight()
            sum_int = 0.0
            count_pos = 0
            for y in range(h):
                for x in range(w):
                    v = ip.get(x, y)
                    if v > effective_thresh:
                        sum_int += v
                        count_pos += 1

            mean_intensity = (sum_int / count_pos) if count_pos > 0 else 0.0
            integrated_intensity = sum_int
            normalized_int = (integrated_intensity / float(nuclei_count)) if nuclei_count > 0 else 0.0

            # export per channel
            export_intensity(
                series_name, marker.getTitle(), count_pos,
                mean_intensity, integrated_intensity,
                nuclei_count, normalized_int, csv_file
            )

            total_positive += count_pos
            total_integrated += integrated_intensity

            # close images cleanly
            channel.changes = False
            marker.changes = False
            channel.close()
            marker.close()

        WindowManager.closeAllWindows()

    # summary per well
    well_mean = (total_integrated / total_positive) if total_positive > 0 else 0.0
    well_normalized = (total_integrated / total_nuclei) if total_nuclei > 0 else 0.0
    export_intensity(
        f.getName(), "AllSeries", total_positive,
        well_mean, total_integrated, total_nuclei,
        well_normalized, well_csv_file
    )
