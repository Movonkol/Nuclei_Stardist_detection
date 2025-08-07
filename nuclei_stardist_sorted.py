from ij import IJ, WindowManager
from ij.io import DirectoryChooser
from java.io import File, FileWriter, BufferedWriter, PrintWriter
from java.lang import Double
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from ij.measure import ResultsTable, Measurements
from ij.plugin.filter import ParticleAnalyzer
import time

# === Benutzereinstellungen ===
model                = IJ.getString("StarDist-Modell", "Versatile (fluorescent nuclei)")
prob                 = IJ.getNumber("Wahrscheinlichkeitsschwelle", 0.5)
nms                  = IJ.getNumber("NMS-Schwelle", 0.5)
tiles_str            = IJ.getString("nTiles (z.B. '1,1')", "1,1")
n_tiles              = int(tiles_str.split(",")[0])
channel_patterns_str = IJ.getString("Kanal-Pattern (z.B. 'c1-,dapi')", "c1-,dapi")
channel_patterns     = [p.strip().lower() for p in channel_patterns_str.split(",")]

size_min             = IJ.getNumber("Minimale Label-Größe (in Pixel²)", 50)
size_max             = IJ.getNumber("Maximale Label-Größe (in Pixel²)", 2500)
# Ordner mit Bildern wählen
dc     = DirectoryChooser("Ordner mit Bildern auswählen")
folder = dc.getDirectory()
if not folder:
    IJ.error("Kein Ordner ausgewählt")
    exit()

csv_counts     = folder + File.separator + "nuclei_counts.csv" + File.separator + "nuclei_counts.csv"
csv_morphology = folder + File.separator + "nuclei_morphology.csv"

# === Export-Funktionen ===
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

# === Batch-Verarbeitung aller Dateien ===
for f in File(folder).listFiles():
    if not f.isFile(): continue
    name = f.getName().lower()
    if not name.endswith((".tif", ".tiff", ".png", ".jpg", ".lif", ".nd2")): continue

    total_nuclei = 0
    areas = []
    roundness_values = []

    # Bio-Formats öffnen
    opts = ImporterOptions()
    opts.setId(f.getAbsolutePath()); opts.setOpenAllSeries(True); opts.setVirtual(True)
    imps = BF.openImagePlus(opts)

    for idx, imp in enumerate(imps):
        title = "%s_Series%d" % (f.getName(), idx+1)
        imp.setTitle(title); imp.show()
        print("→ Processing: %s" % title)

        # Kanal splitten & DAPI wählen
        IJ.run(imp, "Split Channels", ""); time.sleep(0.5)
        dapi = None
        for i in range(1, WindowManager.getImageCount()+1):
            win = WindowManager.getImage(i)
            if any(pat in win.getTitle().lower() for pat in channel_patterns): dapi = win; break
        if not dapi:
            IJ.error("DAPI-Kanal nicht gefunden in %s" % title); IJ.run("Close All"); continue
        IJ.selectWindow(dapi.getTitle())

        # StarDist ausführen
        cmd = ("command=[de.csbdresden.stardist.StarDist2D],"
               "args=['input':'%s','modelChoice':'%s','normalizeInput':'true'," % (dapi.getTitle(), model) +
               "'percentileBottom':'0.0','percentileTop':'100.0','probThresh':'%s','nmsThresh':'%s'," % (prob, nms) +
               "'outputType':'Label Image','nTiles':'%d','excludeBoundary':'2','verbose':'false'," % n_tiles +
               "'showCsbdeepProgress':'false','showProbAndDist':'false'],process=[false]")
        IJ.run("Command From Macro", cmd); time.sleep(1); close_stardist_dialogs()

        # Label-Image finden
        label = next((WindowManager.getImage(i) for i in range(1, WindowManager.getImageCount()+1)
                      if "label image" in WindowManager.getImage(i).getTitle().lower()), None)
        if not label:
            IJ.error("Label image nicht gefunden für %s" % title); IJ.run("Close All"); continue
        IJ.selectWindow(label.getTitle()); label.show()

                # Connected Components
        IJ.run(label, "Connected Components Labeling", "connectivity=4")
        IJ.selectWindow(label.getTitle())

        # ParticleAnalyzer direkt mit Size-Range (min,max)
        IJ.setThreshold(label, 1, Double.POSITIVE_INFINITY)
        rt = ResultsTable()
        pa = ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS,
                               Measurements.AREA | Measurements.CIRCULARITY |
                               Measurements.RECT | Measurements.SHAPE_DESCRIPTORS,
                               rt, size_min, size_max)
        pa.analyze(label)

        # Ergebnis anzeigen und zählen
        count = rt.getCounter()
        print("→ %s: Count=%d" % (label.getTitle(), count))
        IJ.selectWindow(label.getTitle())
        label.show()
        IJ.log("Displayed filtered label for visual check: %s" % label.getTitle())

        total_nuclei += count
        for i in range(count):
            areas.append(rt.getValue("Area", i))
            roundness_values.append(rt.getValue("Circ.", i))

        IJ.run("Close All")
