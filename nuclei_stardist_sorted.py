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
model           = IJ.getString("StarDist Model", "Versatile (fluorescent nuclei)")
prob            = IJ.getNumber("Probability Threshold", 0.5)
nms             = IJ.getNumber("NMS Threshold", 0.5)
tiles_str       = IJ.getString("nTiles (e.g. '1,1')", "1,1")
n_tile          = int(tiles_str.split(",")[0])
channel_patterns_str = IJ.getString("Channel pattern (e.g. 'c1-,dapi')", "c1-,dapi")
channel_patterns     = [p.strip().lower() for p in channel_patterns_str.split(",")]
size_min = IJ.getNumber("Min. label size", 50)
size_max = IJ.getNumber("Max. label size", 2500)
aspect_ratio_max = IJ.getNumber("Max. aspect ratio", 2.0)

dc         = DirectoryChooser("Choose folder with images")
folder     = dc.getDirectory()
if not folder:
    IJ.error("No folder selected"); exit()
csv_file   = folder + File.separator + "nuclei_counts.csv"
morph_file = folder + File.separator + "nuclei_morphology.csv"

# === Export-Funktionen für „gesamt“ pro Datei ===
def export_nuclei_count_summary(image_name, total_count, csv_path):
    f     = File(csv_path)
    first = not f.exists()
    pw    = PrintWriter(BufferedWriter(FileWriter(f, True)))
    if first:
        pw.println("Image,Total_Nuclei_Count")
    pw.println("%s,%d" % (image_name, total_count))
    pw.close()

def export_morphology_summary(image_name, mean_area, mean_roundness, csv_path):
    f     = File(csv_path)
    first = not f.exists()
    pw    = PrintWriter(BufferedWriter(FileWriter(f, True)))
    if first:
        pw.println("Image,Mean_Area,Mean_Roundness")
    pw.println("%s,%.2f,%.4f" % (image_name, mean_area, mean_roundness))
    pw.close()

def close_stardist_dialogs():
    for w in WindowManager.getNonImageWindows():
        t = w.getTitle()
        if t and t.lower().startswith("stardist"):
            w.dispose()

# === Batch über alle Dateien ===
for f in File(folder).listFiles():
    if not f.isFile(): continue
    nm = f.getName().lower()
    if not nm.endswith((".tif", ".tiff", ".png", ".jpg", ".lif", ".nd2")):
        continue

    # Stats für diese Datei initialisieren
    total_count    = 0
    all_areas      = []
    all_roundness = []

    # Bio-Formats: alle Serien öffnen
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

        # Split Channels & finde DAPI
        IJ.run(imp, "Split Channels", ""); time.sleep(0.5)
        dapi = None
        for i in range(1, WindowManager.getImageCount()+1):
            win = WindowManager.getImage(i)
            if any(pat in win.getTitle().lower() for pat in channel_patterns):
                dapi = win; break
        if dapi is None:
            IJ.error("Channel not found in %s" % series_name)
            IJ.run("Close All"); continue
        dapi.show()

        # Vorverarbeitung & StarDist
        IJ.run(dapi, "8-bit", "")
        IJ.run(dapi, "Enhance Contrast", "saturated=0.35")
        IJ.selectWindow(dapi.getTitle())
        cmd = (
            "command=[de.csbdresden.stardist.StarDist2D],"
            "args=["
              "'input':'%s','modelChoice':'%s','normalizeInput':'true',"
              "'percentileBottom':'0.0','percentileTop':'100.0',"
              "'probThresh':'%f','nmsThresh':'%f',"
              "'outputType':'Label Image','nTiles':'%d',"
              "'excludeBoundary':'2','verbose':'false',"
              "'showCsbdeepProgress':'false','showProbAndDist':'false'"
            "],process=[false]"
        ) % (dapi.getTitle(), model, prob, nms, n_tile)
        IJ.run("Command From Macro", cmd); time.sleep(1)
        close_stardist_dialogs()

        # Finde Label-Image
        label = None
        for i in range(1, WindowManager.getImageCount()+1):
            win = WindowManager.getImage(i)
            t   = win.getTitle().lower()
            if "label" in t and "image" in t:
                label = win; break
        if label is None:
            IJ.error("Label image not found for %s" % series_name)
            IJ.run("Close All"); continue
        label.show()

        # Connected Components & Filter Size
        IJ.run(label, "Connected Components Labeling", "connectivity=4")
        IJ.run(label, "Label Size Filtering", "operation=Greater_Than size=%d connectivity=4" % size_min)
        IJ.run(label, "Label Size Filtering", "operation=Lower_Than size=%d connectivity=4" % size_max)
        time.sleep(0.5)

        # Count aufsummieren
        hist  = label.getProcessor().getHistogram()
        count = sum(1 for v in hist[1:] if v > 0)
        total_count += count

        # Morphologie messen
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
            w    = rt.getValue("Width", i)
            h    = rt.getValue("Height", i)
            ar   = max(w,h)/float(min(w,h)) if min(w,h)>0 else 9999
            if ar <= aspect_ratio_max:
                all_areas.append(area)
                all_roundness.append(circ)

        IJ.run("Close All")

    # ==== Nach allen Serien: Export „gesamt“ pro Datei ====
    export_nuclei_count_summary(f.getName(), total_count, csv_file)
    if all_areas:
        mean_area      = sum(all_areas) / len(all_areas)
        mean_roundness = sum(all_roundness) / len(all_roundness)
    else:
        mean_area = mean_roundness = 0.0
    export_morphology_summary(f.getName(), mean_area, mean_roundness, morph_file)

    print("→ %s: TOTAL nuclei=%d, Mean Area=%.2f, Mean Roundness=%.4f"
          % (f.getName(), total_count, mean_area, mean_roundness))
