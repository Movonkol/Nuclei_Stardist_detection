from ij import IJ, WindowManager
from ij.io import DirectoryChooser, FileSaver
from java.io import File, FileWriter, BufferedWriter, PrintWriter
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
import time

# =====================================================================
#                             SETTINGS
# =====================================================================

# Comma-separated patterns used to identify the marker channel by window title
marker_patterns_str = IJ.getString(
    "Marker channel pattern (e.g. 'c3-,marker')",
    "c3-,marker"
)
MARKER_CHANNEL_KEYS = [p.strip().lower() for p in marker_patterns_str.split(",") if p.strip()]

# ---------------------------------------------------------------------
# Threshold configuration: choose between automatic (Otsu etc.) or fixed value
# ---------------------------------------------------------------------
_use_auto_th_str = IJ.getString(
    "Use automatic thresholding? (yes/no)",
    "yes"
)
USE_AUTO_THRESHOLD = _use_auto_th_str.lower().strip() in ("yes","y","true","t","1")

if USE_AUTO_THRESHOLD:
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
    THRESHOLD_METHOD = ""
    THRESHOLD_FACTOR = 1.0
    USE_FIXED_THRESHOLD = True
    FIXED_THRESHOLD = IJ.getNumber(
        "Fixed threshold value", 50.0
    )

# ---------------------------------------------------------------------
# Background subtraction: rolling ball radius and optional median pre-filter
# ---------------------------------------------------------------------
_apply_bg_str = IJ.getString(
    "Subtract background? (yes/no)",
    "yes"
)
APPLY_BACKGROUND = _apply_bg_str.lower().strip() in ("yes","y","true","t","1")
AUTO_BACKGROUND = False
if APPLY_BACKGROUND:
    _auto_bg_str = IJ.getString(
        "Use automatic parameters? (yes/no)",
        "yes"
    )
    AUTO_BACKGROUND = _auto_bg_str.lower().strip() in ("yes","y","true","t","1")
    if AUTO_BACKGROUND:
        ROLLING_RADIUS = IJ.getNumber(
            "Rolling ball radius", 50
        )
        ROLLING_REPEAT = int(IJ.getNumber(
            "Number of repetitions (int)", 1
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
#            Export function: append one measurement row to the CSV
# =====================================================================
def export_intensity(
    image_name, channel_name, pos_pixels, mean_int,
    integrated_int, nuclei_count, normalized_intensity, csv_path
):
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
#      Select input folder and read nucleus counts CSV for normalization
# =====================================================================
dc = DirectoryChooser("Select folder with images")
folder = dc.getDirectory()
if not folder:
    IJ.error("No folder selected")
    raise SystemExit

counts_filename = IJ.getString(
    "Name of nucleus-count CSV (leave empty to skip normalization)",
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
            IJ.log("Note: file '%s' not found - normalization skipped." % counts_filename)
    except Exception as err:
        IJ.log("Error reading nucleus-count CSV: %s" % str(err))

use_heatmap = IJ.getString(
    "Save heatmaps? (yes/no)",
    "no"
).lower().strip() in ("yes","y","true","t","1","ja","j")
heatmap_dir = None
if use_heatmap:
    heatmap_dir = File(folder + File.separator + "heatmaps")
    if not heatmap_dir.exists(): heatmap_dir.mkdir()

csv_file = folder + File.separator + "intensity_measurements.csv"
well_csv_file = folder + File.separator + "intensity_measurements_per_well.csv"

# Load all supported image files from the folder in alphabetical order
files = sorted(
    [f for f in File(folder).listFiles()
     if f.isFile() and f.getName().lower().endswith((".tif", ".tiff", ".png", ".jpg", ".lif", ".nd2"))],
    key=lambda f: f.getName().lower()
)

# =====================================================================
#          Main loop: measure marker intensity for each image file
# =====================================================================
for idx, f in enumerate(files):
    # Retrieve the pre-loaded nucleus count for this file, if available
    file_nuclei_count = counts_list[idx] if idx < len(counts_list) else 0

    total_positive = 0
    total_integrated = 0.0
    total_nuclei = file_nuclei_count

    opts = ImporterOptions()
    opts.setId(f.getAbsolutePath())
    opts.setOpenAllSeries(True)
    opts.setVirtual(True)
    imps = BF.openImagePlus(opts)

    for sidx, imp in enumerate(imps):
        series_name_raw = "%s_Series%d" % (f.getName(), sidx+1)
        imp.setTitle(series_name_raw)
        imp.show()

        IJ.run(imp, "Split Channels", "")
        time.sleep(0.5)

        # Find open channel windows whose titles match the configured marker patterns
        marker_channels = []
        for i in range(1, WindowManager.getImageCount()+1):
            win = WindowManager.getImage(i)
            if not win: continue
            title = win.getTitle().lower()
            if any(mkey in title for mkey in MARKER_CHANNEL_KEYS):
                marker_channels.append(win)

        for marker in marker_channels:
            chan = marker.duplicate()
            chan.show()
            

            if APPLY_BACKGROUND:
                if MEDIAN_RADIUS>0:
                    IJ.run(chan, "Median...", "radius=%d"%int(MEDIAN_RADIUS))
                for _ in range(max(1, ROLLING_REPEAT)):
                    IJ.run(chan, "Subtract Background...", "rolling=%d"%int(ROLLING_RADIUS))

            if USE_FIXED_THRESHOLD:
                th_val = float(FIXED_THRESHOLD)
            else:
                IJ.setAutoThreshold(chan, THRESHOLD_METHOD)
                ip = chan.getProcessor()
                th_val = float(ip.getMinThreshold())
                if th_val != th_val: th_val = 0.0
            eff_th = th_val * float(THRESHOLD_FACTOR)
            # --- Save heatmaps if requested ---
            if use_heatmap and heatmap_dir is not None:
                base = ("%s__%s" % (series_name_raw, marker.getTitle())).replace(" ", "_")
                try:
                    # Save false-color heatmap using the Fire LUT for visual inspection
                    hm = chan.duplicate()
                    IJ.run(hm, "8-bit", "")
                    IJ.run(hm, "Enhance Contrast", "saturated=0.35")
                    IJ.run(hm, "Fire", "")
                    IJ.run(hm, "RGB Color", "")
                    FileSaver(hm).saveAsPng(heatmap_dir.getAbsolutePath() + File.separator + base + "_heatmap.png")
                    hm.changes = False; hm.close()
                    
                    # Binary mask (> eff_th)
                    mask = chan.duplicate()
                    IJ.setThreshold(mask, float(eff_th), 1e12)
                    IJ.run(mask, "Convert to Mask", "")
                    FileSaver(mask).saveAsPng(heatmap_dir.getAbsolutePath() + File.separator + base + "_mask.png")
                    mask.changes = False; mask.close()
                except Exception as _err:
                    IJ.log("Warning: could not save heatmap: %s" % str(_err))
            # --- End heatmap save ---


            ip = chan.getProcessor()
            w,h = ip.getWidth(), ip.getHeight()
            sum_int = 0.0
            cnt = 0
            for y in range(h):
                for x in range(w):
                    v = ip.get(x,y)
                    if v>eff_th:
                        sum_int += v
                        cnt += 1

            mean_int = (sum_int/cnt) if cnt>0 else 0.0
            integ_int = sum_int
            nuclei_count = file_nuclei_count
            norm_int = (integ_int/float(nuclei_count)) if nuclei_count>0 else 0.0

            export_intensity(
                series_name_raw, marker.getTitle(), cnt, mean_int,
                integ_int, nuclei_count, norm_int, csv_file
            )

            total_positive += cnt
            total_integrated += integ_int

            chan.changes = False
            marker.changes = False
            chan.close()
            marker.close()

        WindowManager.closeAllWindows()

    # Aggregate all series into a per-well summary row
    well_mean = (total_integrated/total_positive) if total_positive>0 else 0.0
    well_norm = (total_integrated/total_nuclei) if total_nuclei>0 else 0.0
    export_intensity(
        f.getName(), "AllSeries", total_positive,
        well_mean, total_integrated, total_nuclei, well_norm,
        well_csv_file
    )
