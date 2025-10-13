from ij import IJ, WindowManager, ImagePlus
from ij.io import DirectoryChooser, FileSaver
from ij.plugin.frame import RoiManager
from ij.gui import Overlay, TextRoi, ShapeRoi
from java.io import File, FileWriter, BufferedWriter, PrintWriter
from java.lang import Double
from java.awt import Color
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from ij.measure import ResultsTable, Measurements
from ij.plugin.filter import ParticleAnalyzer
from ij.process import ImageStatistics as IS
import time, re

# ==========================
# User settings (posneg/area-style nuclei detection)
# ==========================
model                = IJ.getString("StarDist model", "Versatile (fluorescent nuclei)")
prob                 = IJ.getNumber("Probability threshold", 0.5)
nms                  = IJ.getNumber("NMS threshold", 0.5)
tiles_str            = IJ.getString("nTiles (e.g. '1,1')", "1,1")
channel_patterns_str = IJ.getString("Channel patterns for DAPI (e.g. 'c1-,dapi,blue')", "c1-,dapi,blue")
CHANNEL_PATTERNS     = [p.strip().lower() for p in channel_patterns_str.split(",") if p.strip()]

# nucleus size filter in ORIGINAL pixel area (px^2)
size_min             = IJ.getNumber("Min label area (px^2)", 50)
size_max             = IJ.getNumber("Max label area (px^2)", 2500)

# Pre-scaling for StarDist (like posneg script)
SCALE_FACTOR         = IJ.getNumber("Pre-scaling (0.5=down, 2.0=up; 1.0=off)", 1.0)

# Intensity sanity filters (same semantics as posneg/area)
SIGMA_ABOVE_BG      = IJ.getNumber("Min sigma above background (0=off)", 1.5)
MIN_MEAN_INTENSITY  = IJ.getNumber("Min mean intensity (0=off; image units)", 0)

# QC / overlay export
SAVE_QC_OVERLAYS    = IJ.getString("Save QC overlays as PNG? (yes/no)", "yes").lower().strip() in ("yes","y","ja","j","true","1")
QC_SAT_PCT          = IJ.getNumber("QC overlay contrast saturation (%)", 1.0) if SAVE_QC_OVERLAYS else 0.0
QC_STROKE_WIDTH     = IJ.getNumber("Overlay stroke width (px)", 2)
QC_COLOR            = Color(255, 0, 0)

# ==========================
# Helpers
# ==========================
def safe_name(s):
    return re.sub(r'[\\/:*?"<>|]+', '_', s)[:180]

def list_open_images():
    imgs = []
    try:
        n = WindowManager.getImageCount()
        for i in range(1, n+1):
            im = WindowManager.getImage(i)
            if im is not None:
                imgs.append(im)
    except:
        pass
    return imgs

def pick_dapi_image(patterns):
    imgs = list_open_images()
    for im in imgs:
        title = (im.getTitle() or "").lower()
        if any(p in title for p in patterns):
            return im
    if len(imgs) == 1:
        return imgs[0]
    # fallback: brightest mean
    best = None; best_mean = -1.0
    for im in imgs:
        try:
            stats = IS.getStatistics(im.getProcessor(), IS.MEAN, im.getCalibration())
            m = float(stats.mean)
        except:
            m = -1.0
        if m > best_mean:
            best_mean = m; best = im
    return best if best is not None else (imgs[0] if imgs else None)

def close_stardist_dialogs():
    try:
        for w in WindowManager.getNonImageWindows():
            try:
                if w.getTitle() and w.getTitle().lower().startswith("stardist"):
                    w.dispose()
            except:
                pass
    except:
        pass

# ==========================
# Output setup
# ==========================
dc     = DirectoryChooser("Select folder with images")
folder = dc.getDirectory()
if not folder:
    IJ.error("No folder selected"); raise SystemExit

csv_counts      = folder + File.separator + "nuclei_counts.csv"
csv_morphology  = folder + File.separator + "nuclei_morphology.csv"

qc_dir = File(folder + File.separator + "QC_Overlays")
if SAVE_QC_OVERLAYS and not qc_dir.exists():
    qc_dir.mkdirs()

# CSV writers
from java.io import File as JFile
from java.io import PrintWriter

def export_counts(image_name, total_count, path):
    f = File(path); first = not f.exists()
    pw = PrintWriter(BufferedWriter(FileWriter(f, True)))
    try:
        if first: pw.println("Image,Total_Nuclei_Count")
        pw.println("%s,%d" % (image_name, total_count))
    finally:
        pw.close()

def export_morphology(image_name, mean_area, mean_circ, path):
    f = File(path); first = not f.exists()
    pw = PrintWriter(BufferedWriter(FileWriter(f, True)))
    try:
        if first: pw.println("Image,Mean_Area,Mean_Circularity")
        pw.println("%s,%.6f,%.6f" % (image_name, mean_area, mean_circ))
    finally:
        pw.close()

# ==========================
# Main batch
# ==========================
grand_total = 0

all_files = []
try:
    for ff in JFile(folder).listFiles():
        name = ff.getName().lower()
        if ff.isFile() and name.endswith((".tif",".tiff",".png",".jpg",".jpeg",".lif",".nd2")):
            all_files.append(ff)
except:
    pass

for f in sorted(all_files, key=lambda x: x.getName().lower()):
    opts = ImporterOptions()
    opts.setId(f.getAbsolutePath())
    opts.setOpenAllSeries(True)
    opts.setVirtual(True)
    try:
        imps = BF.openImagePlus(opts)
    except:
        IJ.log("Open failed: %s" % f.getName()); continue
    if not imps:
        IJ.log("No series in: %s" % f.getName()); continue

    for sidx, imp in enumerate(imps):
        IJ.run("Close All")
        title = ("%s_Series%d" % (f.getName(), sidx+1)) if len(imps)>1 else f.getName()
        imp.setTitle(title); imp.show()

        # Split channels to isolate DAPI (same as posneg/area pipeline)
        try: IJ.run(imp, "Split Channels", "")
        except: pass
        time.sleep(0.2)

        dapi = pick_dapi_image(CHANNEL_PATTERNS)
        if dapi is None:
            IJ.log("No DAPI channel in %s" % title)
            IJ.run("Close All"); continue

        # RGB display window for overlays (blue)
        disp = dapi.duplicate()
        try: IJ.run(disp, "Blue", "")
        except:
            try: IJ.run(disp, "Blue LUT", "")
            except: pass
        IJ.run(disp, "RGB Color", "")
        disp.setTitle(title + "_RGBview"); disp.show()

        # -------- StarDist work image (may be scaled) --------
        orig_w, orig_h = dapi.getWidth(), dapi.getHeight()
        work = dapi.duplicate(); work.setTitle(title+"_work"); work.show()

        if abs(SCALE_FACTOR - 1.0) > 1e-6:
            new_w = max(1, int(round(orig_w * SCALE_FACTOR)))
            new_h = max(1, int(round(orig_h * SCALE_FACTOR)))
            IJ.run(work, "Scale...", "x=%f y=%f width=%d height=%d interpolation=Bilinear average create" %
                   (SCALE_FACTOR, SCALE_FACTOR, new_w, new_h))
            scaled = WindowManager.getCurrentImage()
            try:
                work.changes = False; work.close()
            except:
                pass
            work = scaled
        work.show()
        time.sleep(0.2)

        # -------- StarDist (Both -> we need labels + keep ROIs like posneg/area) --------
        try:
            cmd = ("command=[de.csbdresden.stardist.StarDist2D],"
                   "args=['input':'%s','modelChoice':'%s','normalizeInput':'true',"
                   "'percentileBottom':'0.0','percentileTop':'100.0','probThresh':'%s','nmsThresh':'%s',"
                   "'outputType':'Label Image','nTiles':'%s','excludeBoundary':'2','verbose':'false',"
                   "'showCsbdeepProgress':'false','showProbAndDist':'false'],process=[false]" %
                   (work.getTitle(), model, prob, nms, tiles_str))
            IJ.run("Command From Macro", cmd)
        except Exception as e:
            IJ.log("StarDist failed in %s : %s" % (title, str(e)))
            try: work.changes=False; work.close()
            except: pass
            IJ.run("Close All"); continue
        time.sleep(0.6)
        close_stardist_dialogs()

        # Find label image
        label = None
        for i in range(1, WindowManager.getImageCount()+1):
            im = WindowManager.getImage(i)
            if im and "label image" in (im.getTitle() or "").lower():
                label = im; break
        if label is None:
            IJ.log("Label image not found for %s" % title)
            IJ.run("Close All"); continue

        # Scale label back to original geometry using nearest neighbor (mask-safe)
        if abs(SCALE_FACTOR - 1.0) > 1e-6 and (label.getWidth()!=orig_w or label.getHeight()!=orig_h):
            IJ.run(label, "Scale...", "x=%f y=%f width=%d height=%d interpolation=None create" %
                   (1.0/SCALE_FACTOR, 1.0/SCALE_FACTOR, orig_w, orig_h))
            scaled_back = WindowManager.getCurrentImage()
            try: label.changes=False; label.close()
            except: pass
            label = scaled_back

        # -------- ROI extraction (ParticleAnalyzer) --------
        rm = RoiManager.getInstance() or RoiManager()
        try: rm.reset()
        except: pass

        IJ.run(label, "Connected Components Labeling", "connectivity=4")
        IJ.setThreshold(label, 1, Double.POSITIVE_INFINITY)
        rt = ResultsTable()
        pa = ParticleAnalyzer(
            ParticleAnalyzer.ADD_TO_MANAGER | ParticleAnalyzer.SHOW_MASKS,
            Measurements.AREA | Measurements.CIRCULARITY | Measurements.RECT | Measurements.SHAPE_DESCRIPTORS,
            rt, size_min, size_max
        )
        pa.analyze(label)

        rois_array = list((RoiManager.getInstance() or RoiManager()).getRoisAsArray())

        # -------- Background stats (outside union of ROIs), like posneg/area --------
        if rois_array:
            combined = None
            for r in rois_array:
                sr = ShapeRoi(r)
                combined = sr if combined is None else combined.or(sr)
            dapi.setRoi(combined)
            IJ.run(dapi, "Make Inverse", "")
            bg_stats = dapi.getStatistics(Measurements.MEAN | Measurements.STD_DEV)
            dapi.deleteRoi()
        else:
            bg_stats = dapi.getStatistics(Measurements.MEAN | Measurements.STD_DEV)
        bg_mean, bg_sd = float(bg_stats.mean), float(bg_stats.stdDev)

        # -------- Intensity-based filtering (sigma + absolute) --------
        kept, areas_kept, circ_kept = [], [], []
        for i, roi in enumerate(rois_array):
            dapi.setRoi(roi)
            rs = dapi.getStatistics(Measurements.MEAN)
            roi_mean = float(rs.mean)
            dapi.deleteRoi()

            pass_abs   = (MIN_MEAN_INTENSITY <= 0) or (roi_mean >= MIN_MEAN_INTENSITY)
            pass_sigma = (SIGMA_ABOVE_BG   <= 0) or (roi_mean >= bg_mean + SIGMA_ABOVE_BG*bg_sd)
            if pass_abs and pass_sigma:
                kept.append(roi)
                areas_kept.append(rt.getValue("Area", i))
                circ_kept.append(rt.getValue("Circ.", i))

        rm.reset()
        for r in kept:
            rm.addRoi(r)

        # -------- Summaries --------
        count = len(kept)
        mean_area = (sum(areas_kept)/count) if count>0 else 0.0
        mean_circ = (sum(circ_kept)/count) if count>0 else 0.0
        export_counts(title, count, csv_counts)
        export_morphology(title, mean_area, mean_circ, csv_morphology)
        IJ.log("%s: nuclei=%d" % (title, count))

        # -------- QC overlay on blue DAPI --------
        if count>0:
            ov = Overlay()
            legend = TextRoi(5,5, "QC: detected nuclei\nContours in RED")
            legend.setStrokeColor(Color.white)
            ov.add(legend)
            for r in kept:
                c = r.clone(); c.setStrokeWidth(int(QC_STROKE_WIDTH)); c.setStrokeColor(QC_COLOR)
                ov.add(c)
            disp.setOverlay(ov); disp.updateAndDraw()

            if SAVE_QC_OVERLAYS:
                save_copy = disp.duplicate()
                if QC_SAT_PCT > 0:
                    IJ.run(save_copy, "Enhance Contrast", "saturated=%.3f" % float(QC_SAT_PCT))
                out_name = safe_name("%s__QC_DAPI_Blue.png" % title)
                out_path = qc_dir.getAbsolutePath() + File.separator + out_name
                FileSaver(save_copy).saveAsPng(out_path)
                try: save_copy.changes=False; save_copy.close()
                except: pass
                IJ.log("Saved QC overlay: " + out_path)

        # Cleanup per series
        try: work.changes=False; work.close()
        except: pass
        IJ.run("Close All")
        grand_total += count

# Final total line
export_counts("TOTAL_ALL_IMAGES", grand_total, csv_counts)
IJ.log("Done. QC overlays in: " + qc_dir.getAbsolutePath() if SAVE_QC_OVERLAYS else "Done. QC overlays disabled")

IJ.log("Done. QC overlays in: " + qc_dir.getAbsolutePath() if SAVE_QC_OVERLAYS else "QC overlays disabled")
