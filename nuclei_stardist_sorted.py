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
import time, re, math

# === User settings ===
model                = IJ.getString("StarDist model", "Versatile (fluorescent nuclei)")
prob                 = IJ.getNumber("Probability threshold", 0.5)
nms                  = IJ.getNumber("NMS threshold", 0.5)
tiles_str            = IJ.getString("nTiles (e.g. '1,1')", "1,1")
n_tiles              = int(tiles_str.split(",")[0])
channel_patterns_str = IJ.getString("Channel patterns for DAPI (e.g. 'c1-,dapi,blue')", "c1-,dapi,blue")
channel_patterns     = [p.strip().lower() for p in channel_patterns_str.split(",")]

size_min             = IJ.getNumber("Min label area (px^2)", 50)
size_max             = IJ.getNumber("Max label area (px^2)", 2500)

# Pre-scaling (keep; no background subtraction anymore)
SCALE_FACTOR         = IJ.getNumber("Pre-scaling (0.5=down, 2.0=up; 1.0=off)", 1.0)

# NEW: intensity filter thresholds
SIGMA_ABOVE_BG      = IJ.getNumber("Min σ above background (0=off)", 1.5)
MIN_MEAN_INTENSITY  = IJ.getNumber("Min mean intensity (0=off; image units)", 0)

# QC / overlay export
SAVE_QC_OVERLAYS     = IJ.getString("Save QC overlays as PNG? (yes/no)", "yes").lower().strip() in ("yes","y","ja","j","true","1")
QC_SAT_PCT           = IJ.getNumber("QC overlay contrast saturation (%)", 1.0) if SAVE_QC_OVERLAYS else 0.0
QC_STROKE_WIDTH      = IJ.getNumber("Overlay stroke width (px)", 2)
QC_COLOR             = Color(255, 0, 0)  # red contours

def safe_name(s):
    s = re.sub(r'[\\/:*?"<>|]+', '_', s)
    return s[:180]

def list_open_images():
    imgs = []
    for i in range(1, WindowManager.getImageCount() + 1):
        im = WindowManager.getImage(i)
        if im is not None:
            imgs.append(im)
    return imgs

def pick_dapi_image(patterns):
    imgs = list_open_images()
    for im in imgs:
        title = (im.getTitle() or "").lower()
        if any(p in title for p in patterns):
            return im
    if len(imgs) == 1:
        return imgs[0]
    best = None; best_mean = -1.0
    for im in imgs:
        try:
            stats = IS.getStatistics(im.getProcessor(), IS.MEAN, im.getCalibration())
            m = float(stats.mean)
        except:
            m = 0.0
        if m > best_mean:
            best_mean = m; best = im
    return best if best is not None else (imgs[0] if imgs else None)

def close_stardist_dialogs():
    for w in WindowManager.getNonImageWindows():
        try:
            title = w.getTitle()
            if title and title.lower().startswith("stardist"):
                w.dispose()
        except:
            pass

# Pick folder
dc     = DirectoryChooser("Select folder with images")
folder = dc.getDirectory()
if not folder:
    IJ.error("No folder selected"); exit()

# CSV targets
csv_counts     = folder + File.separator + "nuclei_counts.csv"
csv_morphology = folder + File.separator + "nuclei_morphology.csv"

# QC subfolder
qc_dir = File(folder + File.separator + "QC_Overlays")
if SAVE_QC_OVERLAYS and not qc_dir.exists():
    qc_dir.mkdirs()

def export_counts(image_name, total_count, path):
    f = File(path); first = not f.exists()
    pw = PrintWriter(BufferedWriter(FileWriter(f, True)))
    try:
        if first: pw.println("Image,Total_Nuclei_Count")
        pw.println("%s,%d" % (image_name, total_count))
    finally: pw.close()

def export_morphology(image_name, mean_area, mean_roundness, path):
    f = File(path); first = not f.exists()
    pw = PrintWriter(BufferedWriter(FileWriter(f, True)))
    try:
        if first: pw.println("Image,Mean_Area,Mean_Roundness")
        pw.println("%s,%.2f,%.4f" % (image_name, mean_area, mean_roundness))
    finally: pw.close()

grand_total = 0

# === Batch over files ===
for f in File(folder).listFiles():
    if not f.isFile(): continue
    name = f.getName().lower()
    if not name.endswith((".tif", ".tiff", ".png", ".jpg", ".jpeg", ".lif", ".nd2")): continue

    opts = ImporterOptions()
    opts.setId(f.getAbsolutePath()); opts.setOpenAllSeries(True); opts.setVirtual(True)
    imps = BF.openImagePlus(opts)

    for idx, imp in enumerate(imps):
        title = "%s_Series%d" % (f.getName(), idx + 1) if len(imps) > 1 else f.getName()
        IJ.run("Close All")
        imp.setTitle(title); imp.show()
        print("→ Processing: %s" % title)

        if imp.getType() == ImagePlus.COLOR_RGB:
            IJ.run(imp, "Split Channels", "")
        else:
            try: IJ.run(imp, "Split Channels", "")
            except: pass
        time.sleep(0.3)

        dapi = pick_dapi_image(channel_patterns)
        if not dapi:
            IJ.error("No suitable channel/image found in %s" % title)
            IJ.run("Close All"); continue
        IJ.selectWindow(dapi.getTitle())

        # --- RGB display window (blue DAPI) ---
        disp = dapi.duplicate()
        try: IJ.run(disp, "Blue", "")
        except:
            try: IJ.run(disp, "Blue LUT", "")
            except: pass
        IJ.run(disp, "RGB Color", "")
        disp.setTitle(title + "_RGBview")
        disp.show()

        # ===== Pre-processing on grayscale 'work' (no background subtraction) =====
        orig_w, orig_h = dapi.getWidth(), dapi.getHeight()
        work = dapi.duplicate(); work.setTitle(title + "_work"); work.show()

        if abs(SCALE_FACTOR - 1.0) > 1e-6:
            new_w = max(1, int(round(orig_w * SCALE_FACTOR)))
            new_h = max(1, int(round(orig_h * SCALE_FACTOR)))
            IJ.run(work, "Scale...", "x=%f y=%f width=%d height=%d interpolation=Bilinear average create" %
                   (SCALE_FACTOR, SCALE_FACTOR, new_w, new_h))
            scaled = WindowManager.getCurrentImage()
            try: work.changes=False; work.close()
            except: pass
            work = scaled

        # StarDist on 'work'
        cmd = ("command=[de.csbdresden.stardist.StarDist2D],"
               "args=['input':'%s','modelChoice':'%s','normalizeInput':'true',"
               "'percentileBottom':'0.0','percentileTop':'100.0','probThresh':'%s','nmsThresh':'%s',"
               "'outputType':'Label Image','nTiles':'%d','excludeBoundary':'2','verbose':'false',"
               "'showCsbdeepProgress':'false','showProbAndDist':'false'],process=[false]" %
               (work.getTitle(), model, prob, nms, n_tiles))
        IJ.run("Command From Macro", cmd)
        time.sleep(1)

        for _ in range(3):
            close_stardist_dialogs(); time.sleep(0.1)

        # Find label image
        label = None
        for i in range(1, WindowManager.getImageCount() + 1):
            im = WindowManager.getImage(i)
            if im and "label image" in (im.getTitle() or "").lower():
                label = im; break
        if not label:
            IJ.error("Label image not found for %s" % title)
            IJ.run("Close All"); continue

        # Scale back to original if needed
        if abs(SCALE_FACTOR - 1.0) > 1e-6 and (label.getWidth() != orig_w or label.getHeight() != orig_h):
            IJ.selectWindow(label.getTitle())
            IJ.run(label, "Scale...", "x=%f y=%f width=%d height=%d interpolation=None create" %
                   (1.0/SCALE_FACTOR, 1.0/SCALE_FACTOR, orig_w, orig_h))
            scaled_back = WindowManager.getCurrentImage()
            try: label.changes=False; label.close()
            except: pass
            label = scaled_back

        IJ.selectWindow(label.getTitle()); label.show()

        # ROI extraction
        rm = RoiManager.getInstance() or RoiManager()
        rm.reset()
        IJ.run(label, "Connected Components Labeling", "connectivity=4")
        IJ.setThreshold(label, 1, Double.POSITIVE_INFINITY)
        rt = ResultsTable()
        pa = ParticleAnalyzer(
            ParticleAnalyzer.ADD_TO_MANAGER | ParticleAnalyzer.SHOW_MASKS,
            Measurements.AREA | Measurements.CIRCULARITY | Measurements.RECT | Measurements.SHAPE_DESCRIPTORS,
            rt, size_min, size_max
        )
        pa.analyze(label)

        # ===== NEW: compute background (outside nuclei) and filter ROIs =====
        rois_array = rm.getRoisAsArray()
        # Build union of ROIs for background selection
        combined = None
        for r in rois_array:
            sr = ShapeRoi(r)
            combined = sr if combined is None else combined.or(sr)

        if combined is not None:
            dapi.setRoi(combined)
            IJ.run(dapi, "Make Inverse", "")
            bg_stats = dapi.getStatistics(Measurements.MEAN | Measurements.STD_DEV)
            dapi.deleteRoi()
        else:
            bg_stats = dapi.getStatistics(Measurements.MEAN | Measurements.STD_DEV)

        bg_mean, bg_sd = bg_stats.mean, bg_stats.stdDev

        kept_rois, areas_kept, round_kept = [], [], []
        for i, roi in enumerate(rois_array):
            dapi.setRoi(roi)
            rs = dapi.getStatistics(Measurements.MEAN)
            roi_mean = rs.mean
            dapi.deleteRoi()

            pass_abs   = (MIN_MEAN_INTENSITY <= 0) or (roi_mean >= MIN_MEAN_INTENSITY)
            pass_sigma = (SIGMA_ABOVE_BG <= 0) or (roi_mean >= bg_mean + SIGMA_ABOVE_BG * bg_sd)

            if pass_abs and pass_sigma:
                kept_rois.append(roi)
                areas_kept.append(rt.getValue("Area", i))
                round_kept.append(rt.getValue("Circ.", i))

        # Replace ROIs with filtered set
        rm.reset()
        for r in kept_rois:
            rm.addRoi(r)

        # Recompute summary on kept ROIs
        count = len(kept_rois)
        if count > 0:
            mean_area  = sum(areas_kept) / count
            mean_round = sum(round_kept) / count
        else:
            mean_area = 0.0; mean_round = 0.0

        # CSV export
        export_counts(title, count, csv_counts)
        export_morphology(title, mean_area, mean_round, csv_morphology)
        print("→ %s: Count=%d" % (label.getTitle(), count))

        # --- Overlay on RGB view & optional save ---
        if rm.getCount() > 0:
            ov = Overlay()
            legend = TextRoi(5, 5, "QC: detected nuclei\nContours in RED")
            legend.setStrokeColor(Color.white)
            ov.add(legend)
            for roi in rm.getRoisAsArray():
                c = roi.clone()
                c.setStrokeWidth(int(QC_STROKE_WIDTH))
                c.setStrokeColor(QC_COLOR)
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

        # Cleanup
        try: work.changes=False; work.close()
        except: pass
        IJ.run("Close All")
        grand_total += count

# Append grand total
export_counts("TOTAL_ALL_IMAGES", grand_total, csv_counts)
IJ.log("Done. QC overlays in: " + qc_dir.getAbsolutePath() if SAVE_QC_OVERLAYS else "QC overlays disabled")

            pass_sigma = (SIGMA_ABOVE_BG <= 0) or (roi_mean >= bg_mean + SIGMA_ABOVE_BG * bg_sd)

            if pass_abs and pass_sigma:
                kept_rois.append(roi)
                areas_kept.append(rt.getValue("Area", i))
                round_kept.append(rt.getValue("Circ.", i))

        # Replace ROIs with filtered set
        rm.reset()
        for r in kept_rois:
            rm.addRoi(r)

        # Recompute summary on kept ROIs
        count = len(kept_rois)
        if count > 0:
            mean_area  = sum(areas_kept) / count
            mean_round = sum(round_kept) / count
        else:
            mean_area = 0.0; mean_round = 0.0

        # CSV export
        export_counts(title, count, csv_counts)
        export_morphology(title, mean_area, mean_round, csv_morphology)
        print("→ %s: Count=%d" % (label.getTitle(), count))

        # --- Overlay on RGB view & optional save ---
        if rm.getCount() > 0:
            ov = Overlay()
            legend = TextRoi(5, 5, "QC: detected nuclei\nContours in RED")
            legend.setStrokeColor(Color.white)
            ov.add(legend)
            for roi in rm.getRoisAsArray():
                c = roi.clone()
                c.setStrokeWidth(int(QC_STROKE_WIDTH))
                c.setStrokeColor(QC_COLOR)
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

        # Cleanup
        try: work.changes=False; work.close()
        except: pass
        IJ.run("Close All")
        grand_total += count

# Append grand total
export_counts("TOTAL_ALL_IMAGES", grand_total, csv_counts)
IJ.log("Done. QC overlays in: " + qc_dir.getAbsolutePath() if SAVE_QC_OVERLAYS else "QC overlays disabled")
