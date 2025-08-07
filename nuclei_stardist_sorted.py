from ij import IJ, WindowManager
from ij.io import DirectoryChooser
from java.io import File, FileWriter, BufferedWriter, PrintWriter
from java.lang import Double
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from ij.measure import ResultsTable, Measurements
from ij.plugin.filter import ParticleAnalyzer
import time

# === User Settings ===
model                = IJ.getString("StarDist model name", "Versatile (fluorescent nuclei)")
prob                 = IJ.getNumber("Probability threshold", 0.5)
nms                  = IJ.getNumber("NMS threshold", 0.5)
tiles_str            = IJ.getString("nTiles (e.g. '1,1')", "1,1")
n_tiles              = int(tiles_str.split(",")[0])
channel_patterns_str = IJ.getString("Channel pattern (e.g. 'c1-,dapi')", "c1-,dapi")
channel_patterns     = [p.strip().lower() for p in channel_patterns_str.split(",")]
size_min             = IJ.getNumber("Minimum label size (in pixels)", 50)
size_max             = IJ.getNumber("Maximum label size (in pixels)", 2500)
aspect_ratio_max     = IJ.getNumber("Maximum aspect ratio", 2.0)

# Choose input folder
dc     = DirectoryChooser("Select folder containing images")
folder = dc.getDirectory()
if not folder:
    IJ.error("No folder selected")
    exit()

csv_counts     = folder + File.separator + "nuclei_counts.csv"
csv_morphology = folder + File.separator + "nuclei_morphology.csv"

# === Export functions ===
def export_counts(image_name, total_count, path):
    f = File(path)
    first = not f.exists()
    pw = PrintWriter(BufferedWriter(FileWriter(f, True)))
    if first:
        pw.println("Image,Total_Nuclei_Count")
    pw.println("%s,%d" % (image_name, total_count))
    pw.close()

def export_morphology(image_name, mean_area, mean_roundness, path):
    f = File(path)
    first = not f.exists()
    pw = PrintWriter(BufferedWriter(FileWriter(f, True)))
    if first:
        pw.println("Image,Mean_Area,Mean_Roundness")
    pw.println("%s,%.2f,%.4f" % (image_name, mean_area, mean_roundness))
    pw.close()

def close_stardist_dialogs():
    for w in WindowManager.getNonImageWindows():
        title = w.getTitle()
        if title and title.lower().startswith("stardist"):
            w.dispose()

# === Batch-processing ===
for f in File(folder).listFiles():
    if not f.isFile():
        continue
    name = f.getName().lower()
    if not name.endswith((".tif", ".tiff", ".png", ".jpg", ".lif", ".nd2")):
        continue

    total_nuclei = 0
    areas = []
    roundness_values = []

    # Open all series via Bio-Formats
    opts = ImporterOptions()
    opts.setId(f.getAbsolutePath())
    opts.setOpenAllSeries(True)
    opts.setVirtual(True)
    imps = BF.openImagePlus(opts)

    for idx, imp in enumerate(imps):
        series_title = "%s_Series%d" % (f.getName(), idx+1)
        imp.setTitle(series_title)
        imp.show()
        print("→ Processing:", series_title)

        # Split channels and select DAPI
        IJ.run(imp, "Split Channels", "")
        time.sleep(0.5)
        dapi = None
        for i in range(1, WindowManager.getImageCount()+1):
            win = WindowManager.getImage(i)
            if any(pat in win.getTitle().lower() for pat in channel_patterns):
                dapi = win
                break
        if dapi is None:
            IJ.error("Channel not found in %s" % series_title)
            IJ.run("Close All")
            continue
        dapi.show()

        # Preprocess & run StarDist
        IJ.run(dapi, "8-bit", "")
        IJ.run(dapi, "Enhance Contrast", "saturated=0.35")
        IJ.selectWindow(dapi.getTitle())
        cmd = (
            "command=[de.csbdresden.stardist.StarDist2D],"
            "args=["
              "'input':'%s','modelChoice':'%s','normalizeInput':'true',"
              "'percentileBottom':'0.0','percentileTop':'100.0',"
              "'probThresh':'%s','nmsThresh':'%s',"
              "'outputType':'Label Image','nTiles':'%d',"
              "'excludeBoundary':'2','verbose':'false',"
              "'showCsbdeepProgress':'false','showProbAndDist':'false'"
            "],process=[false]"
        ) % (dapi.getTitle(), model, str(prob), str(nms), n_tiles)
        IJ.run("Command From Macro", cmd)
        time.sleep(1)
        close_stardist_dialogs()

        # Find the resulting label image
        label = None
        for i in range(1, WindowManager.getImageCount()+1):
            win = WindowManager.getImage(i)
            if "label image" in win.getTitle().lower():
                label = win
                break
        if label is None:
            IJ.error("Label image not found for %s" % series_title)
            IJ.run("Close All")
            continue
        label.show()

        # Connected components & pixel-based size filtering
        IJ.run(label, "Connected Components Labeling", "connectivity=4")
        IJ.run("Set Scale...", "distance=1 known=1 unit=pixel global")
        # 1) Remove all objects smaller than size_min
        IJ.run(label,
               "Label Size Filtering",
               "operation=Lower_Than size=%d connectivity=4" % size_min)
        # 2) Remove all objects larger than size_max
        IJ.run(label,
               "Label Size Filtering",
               "operation=Greater_Than size=%d connectivity=4" % size_max)
        time.sleep(0.5)

        # Count remaining labels
        hist = label.getProcessor().getHistogram()
        count = sum(1 for v in hist[1:] if v > 0)
        total_nuclei += count

        # Measure morphology
        IJ.setThreshold(label, 1, Double.POSITIVE_INFINITY)
        rt = ResultsTable()
        pa = ParticleAnalyzer(ParticleAnalyzer.SHOW_NONE,
                              Measurements.AREA | Measurements.CIRCULARITY |
                              Measurements.RECT | Measurements.SHAPE_DESCRIPTORS,
                              rt, 0, Double.POSITIVE_INFINITY)
        pa.analyze(label)
        for i in range(rt.getCounter()):
            area = rt.getValue("Area", i)
            circ = rt.getValue("Circ.", i)
            w = rt.getValue("Width", i)
            h = rt.getValue("Height", i)
            ar = max(w, h) / float(min(w, h)) if min(w, h) > 0 else 9999
            if ar <= aspect_ratio_max:
                areas.append(area)
                roundness_values.append(circ)

        IJ.run("Close All")

    # Export per-file summaries
    export_counts(f.getName(), total_nuclei, csv_counts)
    if areas:
        mean_area = sum(areas) / len(areas)
        mean_roundness = sum(roundness_values) / len(roundness_values)
    else:
        mean_area = mean_roundness = 0.0
    export_morphology(f.getName(), mean_area, mean_roundness, csv_morphology)

    print("→ %s: Total nuclei=%d, Mean Area=%.2f, Mean Roundness=%.4f" %
          (f.getName(), total_nuclei, mean_area, mean_roundness))
