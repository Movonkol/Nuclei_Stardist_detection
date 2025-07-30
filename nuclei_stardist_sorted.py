from ij import IJ, WindowManager
from ij.io import DirectoryChooser
from java.io import File, FileWriter, BufferedWriter, PrintWriter
from java.lang import Double
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from ij.measure import ResultsTable, Measurements
from ij.plugin.filter import ParticleAnalyzer
import time

# === Settings ===
# StarDist parameters
model           = IJ.getString("StarDist Model", "Versatile (fluorescent nuclei)")
prob            = IJ.getNumber("Probability Threshold", 0.5)
nms             = IJ.getNumber("NMS Threshold", 0.5)
tiles_str       = IJ.getString("nTiles (e.g. '1,1')", "1,1")
n_tile          = int(tiles_str.split(",")[0])

# Channel pattern (comma-separated, regex-like)
channel_patterns_str = IJ.getString("Channel pattern (e.g. 'c1-,dapi')", "c1-,dapi")
channel_patterns     = [p.strip().lower() for p in channel_patterns_str.split(",")]

# Size and shape filter thresholds
size_min = IJ.getNumber("Min. label size", 50)
size_max = IJ.getNumber("Max. label size", 2500)
aspect_ratio_max = IJ.getNumber("Max. aspect ratio", 2.0)

# Folder and CSV output paths
dc         = DirectoryChooser("Choose folder with images")
folder     = dc.getDirectory()
if not folder:
    IJ.error("No folder selected")
    exit()
csv_file   = folder + File.separator + "nuclei_counts.csv"
morph_file = folder + File.separator + "nuclei_morphology.csv"

# === Functions ===
def export_nuclei_count(image_name, count, csv_path):
    f     = File(csv_path)
    first = not f.exists()
    pw    = PrintWriter(BufferedWriter(FileWriter(f, True)))
    if first:
        pw.println("Image,Nuclei_Count")
    pw.println("%s,%d" % (image_name, count))
    pw.close()
    print("→ %s: %d nuclei saved to %s" % (image_name, count, csv_path))

def export_morphology(image_name, obj_index, area, roundness, csv_path):
    f     = File(csv_path)
    first = not f.exists()
    pw    = PrintWriter(BufferedWriter(FileWriter(f, True)))
    if first:
        pw.println("Image,Object,Area,Roundness")
    pw.println("%s,%d,%.2f,%.4f" % (image_name, obj_index+1, area, roundness))
    pw.close()
    print("→ %s: Object %d: Area=%.2f, Roundness=%.4f saved to %s" %
          (image_name, obj_index+1, area, roundness, csv_path))

def close_stardist_dialogs():
    for w in WindowManager.getNonImageWindows():
        t = w.getTitle()
        if t and t.lower().startswith("stardist"):
            w.dispose()

# === Batch loop over all images in the folder ===
for f in File(folder).listFiles():
    if not f.isFile(): continue
    nm = f.getName().lower()
    if not nm.endswith((".tif", ".tiff", ".png", ".jpg", ".lif", ".nd2")):
        continue

    # --- Open all series with Bio-Formats ---
    opts = ImporterOptions()
    opts.setId(f.getAbsolutePath())
    opts.setOpenAllSeries(True)
    opts.setVirtual(True)
    imps = BF.openImagePlus(opts)

    for idx, imp in enumerate(imps):
        series_name = "%s_Series%d" % (f.getName(), idx+1)
        imp.setTitle(series_name)
        imp.show()
        print("→ Processing:", series_name)

        # 1) Split channels and find DAPI
        IJ.run(imp, "Split Channels", "")
        time.sleep(0.5)
        dapi = None
        for i in range(1, WindowManager.getImageCount()+1):
            win = WindowManager.getImage(i)
            t   = win.getTitle().lower()
            if any(pat in t for pat in channel_patterns):
                dapi = win
                break
        if dapi is None:
            IJ.error("Channel not found in %s" % series_name)
            IJ.run("Close All")
            continue
        dapi.show()

        # 2) Convert to 8-bit & enhance contrast
        IJ.run(dapi, "8-bit", "")
        IJ.run(dapi, "Enhance Contrast", "saturated=0.35")

        # 3) Headless StarDist
        IJ.selectWindow(dapi.getTitle())
        cmd = (
            "command=[de.csbdresden.stardist.StarDist2D],"
            "args=["
              "'input':'%s',"
              "'modelChoice':'%s',"
              "'normalizeInput':'true',"
              "'percentileBottom':'0.0','percentileTop':'100.0',"
              "'probThresh':'%f','nmsThresh':'%f',"
              "'outputType':'Label Image',"
              "'nTiles':'%d',"
              "'excludeBoundary':'2','verbose':'false',"
              "'showCsbdeepProgress':'false','showProbAndDist':'false'"
            "],process=[false]"
        ) % (dapi.getTitle(), model, prob, nms, n_tile)
        IJ.run("Command From Macro", cmd)
        time.sleep(1)
        close_stardist_dialogs()

        # 4) Find label image
        label = None
        for i in range(1, WindowManager.getImageCount()+1):
            win = WindowManager.getImage(i)
            t   = win.getTitle().lower()
            if "label" in t and "image" in t:
                label = win
                break
        if label is None:
            IJ.error("Label image not found for %s" % series_name)
            IJ.run("Close All")
            continue
        label.show()

        # 5) Connected components and size filtering
        IJ.run(label, "Connected Components Labeling", "connectivity=4")
        label = WindowManager.getCurrentImage()
        IJ.run(label, "Label Size Filtering",
               "operation=Greater_Than size=%d connectivity=4" % size_min)
        IJ.run(label, "Label Size Filtering",
               "operation=Lower_Than  size=%d connectivity=4" % size_max)
        time.sleep(0.5)

        # 6) Count and export total nuclei
        hist  = label.getProcessor().getHistogram()
        count = sum(1 for v in hist[1:] if v > 0)
        export_nuclei_count(series_name, count, csv_file)

        # 7) Morphology & aspect ratio filtering
        IJ.setThreshold(label, 1, Double.POSITIVE_INFINITY)
        rt = ResultsTable()
        pa = ParticleAnalyzer(ParticleAnalyzer.SHOW_NONE,
                              Measurements.AREA | Measurements.CIRCULARITY | Measurements.RECT | Measurements.SHAPE_DESCRIPTORS,
                              rt, 0, Double.POSITIVE_INFINITY)
        pa.analyze(label)
        filtered_count = 0
        for i in range(rt.getCounter()):
            area = rt.getValue("Area", i)
            circ = rt.getValue("Circ.", i)
            width = rt.getValue("Width", i)
            height = rt.getValue("Height", i)
            aspect_ratio = max(width, height) / float(min(width, height)) if min(width, height) > 0 else 9999
            if aspect_ratio <= aspect_ratio_max:
                export_morphology(series_name, i, area, circ, morph_file)
                filtered_count += 1

        print("→ %s: %d objects passed Aspect Ratio ≤ %.2f" % (series_name, filtered_count, aspect_ratio_max))

        # 8) Cleanup
        IJ.run("Close All")
