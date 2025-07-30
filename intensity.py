from ij import IJ, WindowManager
from ij.io import DirectoryChooser
from java.io import File, FileWriter, BufferedWriter, PrintWriter
from java.lang import Double
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from ij.measure import ResultsTable, Measurements
from ij.plugin.filter import ParticleAnalyzer
import time

# === Einstellungen ===
# StarDist‑Parameter
model           = IJ.getString("StarDist‑Modell", "Versatile (fluorescent nuclei)")
prob            = IJ.getNumber("Probability Threshold", 0.5)
nms             = IJ.getNumber("NMS Threshold", 0.5)
tiles_str       = IJ.getString("nTiles (z.B. '1,1')", "1,1")
n_tile          = int(tiles_str.split(",")[0])

# Kanal‑Suchmuster (Komma‑separiert, regex-ähnlich)
channel_patterns_str = IJ.getString("Kanal‑Suchmuster (z.B. 'c1-,dapi')", "c1-,dapi")
channel_patterns     = [p.strip().lower() for p in channel_patterns_str.split(",")]

# Size‑Filter Schwellenwerte
size_min = IJ.getNumber("Min. Label‑Größe", 50)
size_max = IJ.getNumber("Max. Label‑Größe", 2500)

# Ordner‑ und CSV‑Pfade
dc         = DirectoryChooser("Wähle Ordner mit Bildern")
folder     = dc.getDirectory()
if not folder:
    IJ.error("Kein Ordner ausgewählt")
    exit()
csv_file   = folder + File.separator + "nuclei_counts.csv"
morph_file = folder + File.separator + "nuclei_morphology.csv"

# === Funktionen ===
def export_nuclei_count(image_name, count, csv_path):
    f     = File(csv_path)
    first = not f.exists()
    pw    = PrintWriter(BufferedWriter(FileWriter(f, True)))
    if first:
        pw.println("Image,Nuclei_Count")
    pw.println("%s,%d" % (image_name, count))
    pw.close()
    print("→ %s: %d Kerne in %s" % (image_name, count, csv_path))

def export_morphology(image_name, obj_index, area, roundness, csv_path):
    f     = File(csv_path)
    first = not f.exists()
    pw    = PrintWriter(BufferedWriter(FileWriter(f, True)))
    if first:
        pw.println("Image,Object,Area,Roundness")
    pw.println("%s,%d,%.2f,%.4f" % (image_name, obj_index+1, area, roundness))
    pw.close()
    print("→ %s: Objekt %d: Fläche=%.2f, Rundheit=%.4f in %s" %
          (image_name, obj_index+1, area, roundness, csv_path))

def close_stardist_dialogs():
    for w in WindowManager.getNonImageWindows():
        t = w.getTitle()
        if t and t.lower().startswith("stardist"):
            w.dispose()

# === Batch‑Loop über alle Bilder im Ordner ===
for f in File(folder).listFiles():
    if not f.isFile(): continue
    nm = f.getName().lower()
    if not nm.endswith((".tif", ".tiff", ".png", ".jpg", ".lif", ".nd2")):
        continue

    # --- Bio-Formats: alle Serien öffnen ---
    opts = ImporterOptions()
    opts.setId(f.getAbsolutePath())
    opts.setOpenAllSeries(True)
    opts.setVirtual(True)
    imps = BF.openImagePlus(opts)

    for idx, imp in enumerate(imps):
        series_name = "%s_Serie%d" % (f.getName(), idx+1)
        imp.setTitle(series_name)
        imp.show()
        print("→ Verarbeite:", series_name)

        # 1) Split Channels & DAPI finden
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
            IJ.error("Kanal nicht gefunden in %s" % series_name)
            IJ.run("Close All")
            continue
        dapi.show()

        # 2) 8‑Bit & Kontrast
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

        # 4) Label Image finden
        label = None
        for i in range(1, WindowManager.getImageCount()+1):
            win = WindowManager.getImage(i)
            t   = win.getTitle().lower()
            if "label" in t and "image" in t:
                label = win
                break
        if label is None:
            IJ.error("Label‑Bild nicht gefunden für %s" % series_name)
            IJ.run("Close All")
            continue
        label.show()

        # 5) Komponenten verbinden & size filtering
        IJ.run(label, "Connected Components Labeling", "connectivity=4")
        label = WindowManager.getCurrentImage()
        IJ.run(label, "Label Size Filtering",
               "operation=Greater_Than size=%d connectivity=4" % size_min)
        IJ.run(label, "Label Size Filtering",
               "operation=Lower_Than  size=%d connectivity=4" % size_max)
        time.sleep(0.5)

        # 6) Anzahl ermitteln & exportieren
        hist  = label.getProcessor().getHistogram()
        count = sum(1 for v in hist[1:] if v > 0)
        export_nuclei_count(series_name, count, csv_file)

        # 7) Morphologie messen & exportieren
        IJ.setThreshold(label, 1, Double.POSITIVE_INFINITY)
        rt = ResultsTable()
        pa = ParticleAnalyzer(ParticleAnalyzer.SHOW_NONE,
                              Measurements.AREA | Measurements.CIRCULARITY,
                              rt, 0, Double.POSITIVE_INFINITY)
        pa.analyze(label)
        for i in range(rt.getCounter()):
            area = rt.getValue("Area", i)
            circ = rt.getValue("Circ.", i)
            export_morphology(series_name, i, area, circ, morph_file)

        # 8) Aufräumen
        IJ.run("Close All")
