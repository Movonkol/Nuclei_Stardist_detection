from ij import IJ, WindowManager
from ij.io import DirectoryChooser
from java.io import File, FileWriter, BufferedWriter, PrintWriter
from ij.plugin.frame import RoiManager
from ij.gui import Overlay
from java.awt import Color
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
import time

# ==========================
#         SETTINGS
# ==========================
NUCLEI_CHANNEL_KEY  = "c4-"      # Nuclei channel key (e.g., "c1-", "dapi")
MARKER_CHANNEL_KEYS = ["c3-"]    # Marker channel keys

# Ask for main settings
USE_BACKGROUND_SUBTRACTION = IJ.showMessageWithCancel("Background Subtraction", "Do you want to subtract background?")
USE_AUTO_THRESHOLD = IJ.showMessageWithCancel("Auto Threshold", "Do you want to use automatic thresholding?")

# Ask remaining parameters only if needed
if USE_AUTO_THRESHOLD:
    THRESHOLD_METHOD = IJ.getString("Threshold Method (e.g., Otsu dark, Moments, Yen)", "Otsu dark")
    THRESHOLD_FACTOR = IJ.getNumber("Threshold multiplier (>1.0 = stricter)", 1.0)
else:
    USE_FIXED_THRESHOLD = True
    FIXED_THRESHOLD = IJ.getNumber("Fixed Threshold Value", 40)
    THRESHOLD_FACTOR = 1.0  # Ensure defined for all cases

if USE_BACKGROUND_SUBTRACTION:
    ROLLING_RADIUS = int(IJ.getNumber("Rolling Ball Radius", 100))
    ROLLING_REPEAT = int(IJ.getNumber("Background Subtraction Repeats", 1))
    MEDIAN_RADIUS = int(IJ.getNumber("Median Filter Radius (0 = none)", 0))
else:
    ROLLING_RADIUS = 0
    ROLLING_REPEAT = 0
    MEDIAN_RADIUS = 0

# Nucleus size and shape
SIZE_MIN = int(IJ.getNumber("Min nucleus area (px)", 100))
SIZE_MAX = int(IJ.getNumber("Max nucleus area (px)", 1500))
ASPECT_RATIO_MAX = IJ.getNumber("Max aspect ratio (elongated objects get filtered)", 2.0)

# StarDist parameters
model     = IJ.getString("StarDist Model", "Versatile (fluorescent nuclei)")
prob      = IJ.getNumber("Probability Threshold", 0.5)
nms       = IJ.getNumber("NMS Threshold", 0.5)
tiles_str = IJ.getString("nTiles (e.g., '1,1')", "1,1")
n_tile    = int(tiles_str.split(",")[0])

def export_nuclei_classification(image_name, count, positive, negative, csv_path):
    f     = File(csv_path)
    first = not f.exists()
    pw    = PrintWriter(BufferedWriter(FileWriter(f, True)))
    if first:
        pw.println("Image,Nuclei_Count,Positive_Nuclei,Negative_Nuclei")
    pw.println("%s,%d,%d,%d" % (image_name, count, positive, negative))
    pw.close()

def close_stardist_dialogs():
    for w in WindowManager.getNonImageWindows():
        t = w.getTitle()
        if t and t.lower().startswith("stardist"):
            w.dispose()

# Select folder
dc       = DirectoryChooser("Choose folder with images")
folder   = dc.getDirectory()
if not folder:
    IJ.error("No folder selected"); exit()
csv_file = folder + File.separator + "nuclei_counts.csv"

for f in File(folder).listFiles():
    if not f.isFile(): continue
    nm = f.getName().lower()
    if not nm.endswith((".tif", ".tiff", ".png", ".jpg", ".lif", ".nd2")):
        continue

    opts = ImporterOptions()
    opts.setId(f.getAbsolutePath())
    opts.setOpenAllSeries(True)
    opts.setVirtual(True)
    imps = BF.openImagePlus(opts)

    for idx, imp in enumerate(imps):
        series_name = "%s_Serie%d" % (f.getName(), idx+1)
        imp.setTitle(series_name)
        imp.show()
        IJ.run(imp, "Split Channels", "")
        time.sleep(0.5)

        dapi = None
        marker_channels = []
        for i in range(1, WindowManager.getImageCount()+1):
            win = WindowManager.getImage(i)
            t   = win.getTitle().lower()
            if (NUCLEI_CHANNEL_KEY in t) and dapi is None:
                dapi = win
            for mkey in MARKER_CHANNEL_KEYS:
                if mkey in t:
                    marker_channels.append(win)
        if dapi is None or len(marker_channels)==0:
            IJ.error("Missing required channels in %s" % series_name)
            IJ.run("Close All"); continue

        dapi.show()
        IJ.run(dapi, "8-bit", "")
        IJ.run(dapi, "Enhance Contrast", "saturated=0.35")

        IJ.selectWindow(dapi.getTitle())
        cmd = (
            "command=[de.csbdresden.stardist.StarDist2D],"
            "args=["
              "'input':'%s',"
              "'modelChoice':'%s',"
              "'normalizeInput':'true',"
              "'percentileBottom':'0.0','percentileTop':'100.0',"
              "'probThresh':'%f','nmsThresh':'%f',"
              "'outputType':'ROI Manager',"
              "'nTiles':'%d',"
              "'excludeBoundary':'2','verbose':'false',"
              "'showCsbdeepProgress':'false','showProbAndDist':'false'"
            "],process=[false]"
        ) % (dapi.getTitle(), model, prob, nms, n_tile)
        IJ.run("Command From Macro", cmd)
        time.sleep(1)
        close_stardist_dialogs()

        rm = RoiManager.getInstance()
        if not rm: rm = RoiManager()
        rois = rm.getRoisAsArray()
        keep_idxs = []
        for idx, roi in enumerate(rois):
            dapi.setRoi(roi)
            stats = dapi.getStatistics()
            area = stats.area
            bounds = roi.getBounds()
            aspect_ratio = max(bounds.width, bounds.height) / float(min(bounds.width, bounds.height))
            if SIZE_MIN <= area <= SIZE_MAX and aspect_ratio <= ASPECT_RATIO_MAX:
                keep_idxs.append(idx)
        if rm.getCount() > 0:
            rm.runCommand("Select All")
            rm.runCommand("Delete")
        for idx in keep_idxs:
            rm.addRoi(rois[idx])
        rois = rm.getRoisAsArray()
        rm.reset()

        if not rois or len(rois)==0:
            print("No valid nuclei found in %s" % series_name)
            IJ.run("Close All"); continue

        for marker in marker_channels:
            if marker is None:
                print("Marker is None, skipping...")
                continue
            marker_name = marker.getTitle().replace(" ", "_")
            marker_original = marker.duplicate()
            IJ.saveAs(marker_original, "PNG", folder + File.separator + series_name + "_" + marker_name + "_RAW.png")

            marker_for_analysis = marker.duplicate()
            marker_for_analysis.show()
            IJ.run(marker_for_analysis, "8-bit", "")
            IJ.run(marker_for_analysis, "Enhance Contrast", "saturated=0.35")
            if USE_BACKGROUND_SUBTRACTION:
                if MEDIAN_RADIUS > 0:
                    IJ.run(marker_for_analysis, "Median...", "radius=%d" % MEDIAN_RADIUS)
                for _ in range(ROLLING_REPEAT):
                    IJ.run(marker_for_analysis, "Subtract Background...", "rolling=%d" % ROLLING_RADIUS)

            if USE_AUTO_THRESHOLD:
                IJ.setAutoThreshold(marker_for_analysis, THRESHOLD_METHOD)
                marker_ip = marker_for_analysis.getProcessor()
                threshold_value = marker_ip.getMinThreshold()
            else:
                threshold_value = FIXED_THRESHOLD

            positive = 0
            negative = 0
            overlay = Overlay()
            marker.setRoi(None)
            for roi in rois:
                marker_for_analysis.setRoi(roi)
                stats = marker_for_analysis.getStatistics()
                mean = stats.mean
                if mean > THRESHOLD_FACTOR * threshold_value:
                    roi.setStrokeColor(Color(0,255,0))
                    positive += 1
                else:
                    roi.setStrokeColor(Color(255,0,0))
                    negative += 1
                roi.setStrokeWidth(2)
                overlay.add(roi)

            count = positive + negative
            export_nuclei_classification(series_name + "_" + marker_name, count, positive, negative, csv_file)

            marker_original.setOverlay(overlay)
            marker_original.show()
            IJ.saveAs(marker_original, "PNG", folder + File.separator + series_name + "_" + marker_name + "_RAW_classified.png")

        IJ.run("Close All")
