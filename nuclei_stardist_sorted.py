from ij import IJ, WindowManager
from ij.io import DirectoryChooser
from java.io import File, FileWriter, BufferedWriter, PrintWriter
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
import time

# === Settings ===
# StarDist parameters
model = IJ.getString("StarDist Model", "Versatile (fluorescent nuclei)")
prob = IJ.getNumber("Probability Threshold", 0.5)
nms = IJ.getNumber("NMS Threshold", 0.5)
tiles_str = IJ.getString("Number of Tiles (e.g., '1,1')", "1,1")
n_tile = int(tiles_str.split(",")[0])

# Size filter parameters
min_size = IJ.getNumber("Min Cell Size (lower bound)", 50)
max_size = IJ.getNumber("Max Cell Size (upper bound)", 2500)
connectivity = int(IJ.getNumber("Connectivity", 4))

# DAPI channel identifiers (comma-separated keywords)
dapi_keywords_str = IJ.getString("DAPI channel identifier(s) (comma-separated)", "c1-,dapi")
dapi_keywords = [kw.strip().lower() for kw in dapi_keywords_str.split(",") if kw.strip()]

# Directories and paths
dc = DirectoryChooser("Select Image Folder")
folder = dc.getDirectory()
if not folder:
    IJ.error("No folder selected")
    exit()
csv_file = folder + File.separator + "nuclei_counts.csv"

# === Functions ===

def export_nuclei_count(image_name, count, csv_path):
    f = File(csv_path)
    first = not f.exists()
    pw = PrintWriter(BufferedWriter(FileWriter(f, True)))
    if first:
        pw.println("Image,Nuclei_Count")
    pw.println("%s,%d" % (image_name, count))
    pw.close()
    print("→ %s: %d nuclei in %s" % (image_name, count, csv_path))


def close_stardist_dialogs():
    for w in WindowManager.getNonImageWindows():
        t = w.getTitle()
        if t and t.lower().startswith("stardist"):
            w.dispose()

# === Main Processing ===
for f in File(folder).listFiles():
    if not f.isFile():
        continue
    nm = f.getName().lower()
    if not nm.endswith((".tif", ".tiff", ".png", ".jpg", ".lif", ".nd2")):
        continue

    # Bio-Formats: open all series
    opts = ImporterOptions()
    opts.setId(f.getAbsolutePath())
    opts.setOpenAllSeries(True)
    opts.setVirtual(True)
    imps = BF.openImagePlus(opts)

    # iterate through each series
    for idx, imp in enumerate(imps):
        series_name = "%s_Series%d" % (f.getName(), idx+1)
        imp.setTitle(series_name)
        imp.show()
        print("→ Processing:", series_name)

        # Split channels & find DAPI channel
        IJ.run(imp, "Split Channels", "")
        time.sleep(0.5)
        dapi = None
        for i in range(1, WindowManager.getImageCount()+1):
            win = WindowManager.getImage(i)
            title_lower = win.getTitle().lower()
            if any(kw in title_lower for kw in dapi_keywords):
                dapi = win
                break
        if dapi is None:
            IJ.error("DAPI channel not found in %s" % series_name)
            IJ.run("Close All")
            continue
        dapi.show()

        # Convert to 8-bit & enhance contrast
        IJ.run(dapi, "8-bit", "")
        IJ.run(dapi, "Enhance Contrast", "saturated=0.35")

        # Run StarDist (headless) via SciJava command
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

        # Find Label Image
        label = None
        for i in range(1, WindowManager.getImageCount()+1):
            win = WindowManager.getImage(i)
            t_lower = win.getTitle().lower()
            if "label" in t_lower and "image" in t_lower:
                label = win
                break
        if label is None:
            IJ.error("Label image not found for %s" % series_name)
            IJ.run("Close All")
            continue
        label.show()

        # Connected components & size filtering
        IJ.run(label, "Connected Components Labeling", "connectivity=%d" % connectivity)
        IJ.run(label, "Label Size Filtering", "operation=Greater_Than size=%d connectivity=%d" % (min_size, connectivity))
        IJ.run(label, "Label Size Filtering", "operation=Lower_Than size=%d connectivity=%d" % (max_size, connectivity))
        time.sleep(0.5)

        # Count nuclei & export results
        hist = label.getProcessor().getHistogram()
        count = sum(1 for v in hist[1:] if v > 0)
        export_nuclei_count(series_name, count, csv_file)

        # Clean up before next series
        IJ.run("Close All")
