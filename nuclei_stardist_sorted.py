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
model             = IJ.getString("StarDist Model", "Versatile (fluorescent nuclei)")
probability_thresh = IJ.getNumber("Probability Threshold", 0.5)
nms_thresh        = IJ.getNumber("NMS Threshold", 0.5)
tiles_str         = IJ.getString("nTiles (e.g. '1,1')", "1,1")
n_tiles           = int(tiles_str.split(",")[0])

# Channel search patterns (comma-separated, case-insensitive)
channel_patterns_str = IJ.getString("Channel Search Patterns (e.g. 'c1-,dapi')", "c1-,dapi")
channel_patterns     = [p.strip().lower() for p in channel_patterns_str.split(",")]

# Size filter thresholds
min_label_size = IJ.getNumber("Minimum Label Size (px)", 50)
max_label_size = IJ.getNumber("Maximum Label Size (px)", 2500)

# Choose folder and define output CSV paths
dc         = DirectoryChooser("Select folder with images")
input_folder = dc.getDirectory()
if not input_folder:
    IJ.error("No folder selected.")
    exit()
csv_counts = input_folder + File.separator + "nuclei_counts.csv"
csv_morph  = input_folder + File.separator + "nuclei_morphology.csv"

# === Helper functions ===
def export_nuclei_count(name, count, path):
    f     = File(path)
    first = not f.exists()
    pw    = PrintWriter(BufferedWriter(FileWriter(f, True)))
    if first:
        pw.println("Image,Nuclei_Count")
    pw.println("%s,%d" % (name, count))
    pw.close()
    print("→ %s: %d nuclei in %s" % (name, count, path))


def export_morphology(name, index, area, circ, path):
    f     = File(path)
    first = not f.exists()
    pw    = PrintWriter(BufferedWriter(FileWriter(f, True)))
    if first:
        pw.println("Image,Object,Area,Roundness")
    pw.println("%s,%d,%.2f,%.4f" % (name, index+1, area, circ))
    pw.close()
    print("→ %s: Object %d: Area=%.2f, Roundness=%.4f in %s" % (name, index+1, area, circ, path))


def close_stardist_dialogs():
    for w in WindowManager.getNonImageWindows():
        title = w.getTitle()
        if title and title.lower().startswith("stardist"):
            w.dispose()

# === Batch processing ===
for file_obj in File(input_folder).listFiles():
    if not file_obj.isFile():
        continue
    fname = file_obj.getName().lower()
    if not fname.endswith((".tif", ".tiff", ".png", ".jpg", ".lif", ".nd2")):
        continue

    # Import all series via Bio-Formats
    opts = ImporterOptions()
    opts.setId(file_obj.getAbsolutePath())
    opts.setOpenAllSeries(True)
    opts.setVirtual(True)
    image_list = BF.openImagePlus(opts)

    for idx, imp in enumerate(image_list):
        title = "%s_Series%d" % (file_obj.getName(), idx+1)
        imp.setTitle(title)
        imp.show()
        print("→ Processing:", title)

        # 1) Split channels and find target channel
        IJ.run(imp, "Split Channels", "")
        time.sleep(0.5)
        target = None
        for i in range(1, WindowManager.getImageCount()+1):
            win = WindowManager.getImage(i)
            t   = win.getTitle().lower()
            if any(p in t for p in channel_patterns):
                target = win
                break
        if not target:
            IJ.error("Target channel not found in %s" % title)
            IJ.run("Close All")
            continue
        target.show()

        # 2) Convert to 8-bit and enhance contrast
        IJ.run(target, "8-bit", "")
        IJ.run(target, "Enhance Contrast", "saturated=0.35")

        # 3) Run StarDist headlessly
        IJ.selectWindow(target.getTitle())
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
        ) % (target.getTitle(), model, probability_thresh, nms_thresh, n_tiles)
        IJ.run("Command From Macro", cmd)
        time.sleep(1)
        close_stardist_dialogs()

        # 4) Locate label image
        label = None
        for i in range(1, WindowManager.getImageCount()+1):
            win = WindowManager.getImage(i)
            if "label" in win.getTitle().lower():
                label = win
                break
        if not label:
            IJ.error("Label image not found for %s" % title)
            IJ.run("Close All")
            continue
        label.show()

        # 5) Connected components and size filtering
        IJ.run(label, "Connected Components Labeling", "connectivity=4")
        label = WindowManager.getCurrentImage()
        IJ.run(label, "Label Size Filtering", "operation=Greater_Than size=%d connectivity=4" % min_label_size)
        IJ.run(label, "Label Size Filtering", "operation=Lower_Than  size=%d connectivity=4" % max_label_size)
        time.sleep(0.5)

        # 6) Count nuclei and export
        hist  = label.getProcessor().getHistogram()
        count = sum(1 for v in hist[1:] if v > 0)
        export_nuclei_count(title, count, csv_counts)

        # 7) Measure morphology and export
        IJ.setThreshold(label, 1, Double.POSITIVE_INFINITY)
        rt = ResultsTable()
        pa = ParticleAnalyzer(ParticleAnalyzer.SHOW_NONE,
                              Measurements.AREA | Measurements.CIRCULARITY,
                              rt, 0, Double.POSITIVE_INFINITY)
        pa.analyze(label)
        for i in range(rt.getCounter()):
            area = rt.getValue("Area", i)
            circ = rt.getValue("Circ.", i)
            export_morphology(title, i, area, circ, csv_morph)

        # 8) Clean up
        IJ.run("Close All")

