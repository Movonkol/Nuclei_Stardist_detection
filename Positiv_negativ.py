# Fiji / Jython — Nuclei Pos/Neg (StarDist-Detektion wie im "Nuclei_stardist.py")
# - DAPI-Overlay in Blau speichern
# - Marker-Overlay in Orange speichern
# - Rest wie zuvor (Kanalmuster, Pre-Scaling, DAPI-Intensitätsfilter, Pos/Neg pro Marker)
# Autor: angepasst für Moritz (2025-08-30)

from ij import IJ, WindowManager, ImagePlus
from ij.io import DirectoryChooser, FileSaver
from ij.plugin.frame import RoiManager
from ij.gui import Overlay, TextRoi, ShapeRoi
from ij.measure import ResultsTable, Measurements
from ij.plugin.filter import ParticleAnalyzer
from ij.process import ImageStatistics as IS
from java.io import File, FileWriter, BufferedWriter, PrintWriter
from java.lang import Double
from java.awt import Color
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
import time, re

# ==========================
#        EINSTELLUNGEN
# ==========================
channel_patterns_str = IJ.getString("Channel patterns for DAPI (z.B. 'c1-,dapi,blue')", "c1-,dapi,blue")
CHANNEL_PATTERNS = [p.strip().lower() for p in channel_patterns_str.split(",") if p.strip()]

marker_keys_str = IJ.getString("Marker channel identifiers (comma-separated, z.B. 'c2-,c3-')", "c3-")
MARKER_CHANNEL_KEYS = [s.strip().lower() for s in marker_keys_str.split(',') if s.strip()][:6]

marker_fixed = []
for mkey in MARKER_CHANNEL_KEYS:
    marker_fixed.append(IJ.getNumber("Fixed threshold for %s (Original bit depth)" % mkey, 1000.0))

APPLY_BG = IJ.getString("Subtract background on marker channels? (yes/no)", "yes").lower().strip() in ("yes","y","ja","j","true","1")
if APPLY_BG:
    ROLLING_RADIUS = IJ.getNumber("Rolling Ball Radius", 50)
    ROLLING_REPEAT = int(IJ.getNumber("BG repeats (integer)", 1))
    MEDIAN_RADIUS  = IJ.getNumber("Median filter radius (0 = none)", 0)
else:
    ROLLING_RADIUS, ROLLING_REPEAT, MEDIAN_RADIUS = 0, 0, 0

SHOW_OVERLAYS_AT_END = IJ.getString("Show overlay maps at end? (yes/no)", "yes").lower().strip() in ("yes","y","ja","j","true","1")
ADD_PN_LABELS = IJ.getString("Add P/N labels to overlays? (yes/no)", "yes").lower().strip() in ("yes","y","ja","j","true","1")
ENHANCE_OVERLAY_CONTRAST = IJ.getString("Enhance contrast for saved overlays? (yes/no)", "no").lower().strip() in ("yes","y","ja","j","true","1")
DAPI_OVERLAY_SAT   = IJ.getNumber("DAPI overlay saturation (%)", 1.0)
MARKER_OVERLAY_SAT = IJ.getNumber("Marker overlay saturation (%)", 0.35)

model         = IJ.getString("StarDist model", "Versatile (fluorescent nuclei)")
prob          = IJ.getNumber("Probability threshold", 0.5)
nms           = IJ.getNumber("NMS threshold", 0.5)
tiles_str     = IJ.getString("nTiles (e.g. '1,1')", "1,1")
N_TILES       = tiles_str

SCALE_FACTOR  = IJ.getNumber("Pre-scaling (0.5=down, 2.0=up; 1.0=off)", 1.0)

SIZE_MIN = int(IJ.getNumber("Min nucleus area (px)", 50))
SIZE_MAX = int(IJ.getNumber("Max nucleus area (px)", 2500))

SIGMA_ABOVE_BG     = IJ.getNumber("Min σ above background (0=off)", 1.5)
MIN_MEAN_INTENSITY = IJ.getNumber("Min mean intensity (0=off; image units)", 0)

POS_PIXELS_FRACTION_PCT = IJ.getNumber("Min. fraction of ROI pixels ≥ marker thr (%)", 5.0)
COLOR_POS = Color(0,255,0)
COLOR_NEG = Color(255,0,0)

# ==========================
#      HILFSFUNKTIONEN
# ==========================
def close_stardist_dialogs():
    for w in WindowManager.getNonImageWindows():
        try:
            if w.getTitle() and w.getTitle().lower().startswith("stardist"):
                w.dispose()
        except:
            pass

def preprocess_marker(imp):
    dup = imp.duplicate(); dup.show()
    if APPLY_BG:
        if MEDIAN_RADIUS > 0:
            IJ.run(dup, "Median...", "radius=%d" % int(MEDIAN_RADIUS))
        for _ in range(max(1, int(ROLLING_REPEAT))):
            IJ.run(dup, "Subtract Background...", "rolling=%d" % int(ROLLING_RADIUS))
    return dup, dup.getProcessor()

def roi_pixel_iter(ip, roi):
    b = roi.getBounds()
    W = ip.getWidth(); H = ip.getHeight()
    rmask = roi.getMask()
    for yy in range(b.height):
        gy = b.y + yy
        if gy < 0 or gy >= H:
            continue
        for xx in range(b.width):
            gx = b.x + xx
            if gx < 0 or gx >= W:
                continue
            if rmask is not None and rmask.get(xx, yy) == 0:
                continue
            yield gx, gy

def classify_roi(marker_ip, roi, thr, min_pos_frac):
    total = 0
    pospix = 0
    for gx, gy in roi_pixel_iter(marker_ip, roi):
        total += 1
        if marker_ip.get(gx, gy) >= thr:
            pospix += 1
    if total == 0:
        return False, 0, 0
    return (float(pospix)/float(total) >= min_pos_frac), total, pospix

def export_counts_row(csv_path, row, header):
    f = File(csv_path); first = not f.exists()
    pw = PrintWriter(BufferedWriter(FileWriter(f, True)))
    try:
        if first: pw.println(header)
        pw.println(row)
    finally:
        pw.close()

def safe_name(s):
    s = re.sub(r'[\\/:*?"<>|]+', '_', s)
    return s[:180]

def add_label(ov, roi, text, color=Color.white):
    if not ADD_PN_LABELS:
        return
    b = roi.getBounds()
    cx = b.x + int(b.width/2); cy = b.y + int(b.height/2)
    tr = TextRoi(cx, cy, text)
    tr.setStrokeColor(color); tr.setJustification(TextRoi.CENTER)
    ov.add(tr)

def pick_dapi_from_list(img_list, patterns):
    for im in img_list:
        title = (im.getTitle() or "").lower()
        if any(p in title for p in patterns):
            return im
    best = None; best_mean = -1.0
    for im in img_list:
        try:
            stats = IS.getStatistics(im.getProcessor(), IS.MEAN, im.getCalibration())
            m = float(stats.mean)
        except:
            m = 0.0
        if m > best_mean:
            best_mean = m; best = im
    return best

# ==========================
#        HAUPTPROZESS
# ==========================
folder = DirectoryChooser("Choose folder with images").getDirectory()
if not folder:
    IJ.error("No folder selected"); raise SystemExit

png_dir = File(folder + File.separator + "PosNeg_PNGs")
if not png_dir.exists():
    png_dir.mkdirs()

csv_counts = folder + File.separator + "nuclei_posneg.csv"
csv_detail = folder + File.separator + "nuclei_posneg_perROI.csv"

from java.io import File as JFile
files = sorted([f for f in JFile(folder).listFiles()
                if f.isFile() and f.getName().lower().endswith((".tif",".tiff",".png",".jpg",".lif",".nd2"))],
               key=lambda f: f.getName().lower())

min_pos_frac = float(POS_PIXELS_FRACTION_PCT) / 100.0

for f in files:
    opts = ImporterOptions(); opts.setId(f.getAbsolutePath()); opts.setOpenAllSeries(True); opts.setVirtual(True)
    imps = BF.openImagePlus(opts)

    for sidx, imp in enumerate(imps):
        series_name = "%s_Series%d" % (f.getName(), sidx+1) if len(imps) > 1 else f.getName()
        IJ.run("Close All")
        imp.setTitle(series_name); imp.show()
        try:
            IJ.run(imp, "Split Channels", "")
        except:
            pass
        time.sleep(0.3)

        windows = [WindowManager.getImage(i) for i in range(1, WindowManager.getImageCount()+1) if WindowManager.getImage(i) is not None]
        dapi = pick_dapi_from_list(windows, CHANNEL_PATTERNS)
        markers = []
        for win in windows:
            t = (win.getTitle() or "").lower()
            for mkey in MARKER_CHANNEL_KEYS:
                if mkey in t:
                    markers.append(win); break

        if dapi is None or not markers:
            IJ.log("Missing channels in %s (DAPI and/or Markers)" % series_name)
            for w in windows:
                try: w.changes=False; w.close()
                except: pass
            try: imp.changes=False; imp.close()
            except: pass
            continue

        # ===== StarDist =====
        try:
            disp = dapi.duplicate()
            try: IJ.run(disp, "Blue", "")
            except:
                try: IJ.run(disp, "Blue LUT", "")
                except: pass
            IJ.run(disp, "RGB Color", "")
            disp.setTitle(series_name + "_RGBview"); disp.show()
        except:
            disp = None

        orig_w, orig_h = dapi.getWidth(), dapi.getHeight()
        work = dapi.duplicate(); work.setTitle(series_name + "_work"); work.show()

        if abs(SCALE_FACTOR - 1.0) > 1e-6:
            new_w = max(1, int(round(orig_w * SCALE_FACTOR)))
            new_h = max(1, int(round(orig_h * SCALE_FACTOR)))
            IJ.run(work, "Scale...", "x=%f y=%f width=%d height=%d interpolation=Bilinear average create" %
                   (SCALE_FACTOR, SCALE_FACTOR, new_w, new_h))
            scaled = WindowManager.getCurrentImage()
            try: work.changes=False; work.close()
            except: pass
            work = scaled

        cmd = ("command=[de.csbdresden.stardist.StarDist2D],"
               "args=['input':'%s','modelChoice':'%s','normalizeInput':'true',"
               "'percentileBottom':'0.0','percentileTop':'100.0','probThresh':'%s','nmsThresh':'%s',"
               "'outputType':'Label Image','nTiles':'%s','excludeBoundary':'2','verbose':'false',"
               "'showCsbdeepProgress':'false','showProbAndDist':'false'],process=[false]") % (
                    work.getTitle(), model, prob, nms, N_TILES)
        IJ.run("Command From Macro", cmd); time.sleep(1.0); close_stardist_dialogs()

        label = None
        for i in range(1, WindowManager.getImageCount()+1):
            win = WindowManager.getImage(i)
            if win and "label image" in (win.getTitle() or "").lower():
                label = win; break
        if label is None:
            IJ.log("Label image nicht gefunden für %s" % series_name)
            try: work.changes=False; work.close()
            except: pass
            for w in windows:
                try: w.changes=False; w.close()
                except: pass
            try: imp.changes=False; imp.close()
            except: pass
            continue

        if abs(SCALE_FACTOR - 1.0) > 1e-6 and (label.getWidth() != orig_w or label.getHeight() != orig_h):
            IJ.selectWindow(label.getTitle())
            IJ.run(label, "Scale...", "x=%f y=%f width=%d height=%d interpolation=None create" %
                   (1.0/SCALE_FACTOR, 1.0/SCALE_FACTOR, orig_w, orig_h))
            scaled_back = WindowManager.getCurrentImage()
            try: label.changes=False; label.close()
            except: pass
            label = scaled_back

        # ===== ROIs & Filter =====
        rm = RoiManager.getInstance() or RoiManager(); rm.reset()
        IJ.selectWindow(label.getTitle())
        IJ.run(label, "Connected Components Labeling", "connectivity=4")
        IJ.setThreshold(label, 1, Double.POSITIVE_INFINITY)
        rt = ResultsTable()
        pa = ParticleAnalyzer(
            ParticleAnalyzer.ADD_TO_MANAGER,
            Measurements.AREA | Measurements.CIRCULARITY | Measurements.RECT | Measurements.SHAPE_DESCRIPTORS,
            rt, SIZE_MIN, SIZE_MAX
        )
        pa.setHideOutputImage(True)
        pa.analyze(label)

        rois_array = rm.getRoisAsArray()

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

        kept_rois = []
        for roi in rois_array:
            dapi.setRoi(roi)
            rs = dapi.getStatistics(Measurements.MEAN)
            roi_mean = rs.mean
            dapi.deleteRoi()

        # DAPI-Filter
            pass_abs   = (MIN_MEAN_INTENSITY <= 0) or (roi_mean >= MIN_MEAN_INTENSITY)
            pass_sigma = (SIGMA_ABOVE_BG <= 0) or (roi_mean >= bg_mean + SIGMA_ABOVE_BG * bg_sd)
            if pass_abs and pass_sigma:
                kept_rois.append(roi)

        rm.reset()
        for r in kept_rois:
            rm.addRoi(r)

        valid = list(rm.getRoisAsArray())
        if not valid:
            try: work.changes=False; work.close()
            except: pass
            for w in windows:
                try: w.changes=False; w.close()
                except: pass
            try: imp.changes=False; imp.close()
            except: pass
            continue

        # ===== Marker vorbereiten =====
        marker_preps = []
        for idx, m in enumerate(markers[:len(MARKER_CHANNEL_KEYS)]):
            proc_imp, proc_ip = preprocess_marker(m)
            thr = float(marker_fixed[idx])
            marker_preps.append((m, proc_imp, proc_ip, thr))

        counts_header = "Image;Series;Marker;N_total;Positive;Negative;Percent_Positive"
        detail_header = "Image;Series;ROI_Index;Marker;ROI_px;PosPix;Is_Positive"

        # ===== Klassifikation & Overlays =====
        for midx, (m_orig, m_imp, m_ip, thr) in enumerate(marker_preps):
            pos_count = 0; neg_count = 0; total = 0
            ov = Overlay()

            legend = TextRoi(5, 5, "Marker: %s\nGreen = Positive\nRed = Negative" % m_orig.getTitle())
            legend.setStrokeColor(Color.white)
            ov.add(legend)

            for ridx, roi in enumerate(valid):
                is_pos, roi_px, pospix = classify_roi(m_ip, roi, thr, min_pos_frac)
                if roi_px == 0:
                    continue
                total += 1
                c = roi.clone(); c.setStrokeWidth(2)
                if is_pos:
                    pos_count += 1
                    c.setStrokeColor(COLOR_POS); add_label(ov, roi, "P", COLOR_POS)
                else:
                    neg_count += 1
                    c.setStrokeColor(COLOR_NEG); add_label(ov, roi, "N", COLOR_NEG)
                ov.add(c)

                detail_row = "%s;%s;%d;%s;%d;%d;%s" % (
                    f.getName(), series_name, ridx+1, m_orig.getTitle(),
                    int(roi_px), int(pospix), "TRUE" if is_pos else "FALSE"
                )
                export_counts_row(csv_detail, detail_row, detail_header)

            pct = (100.0 * pos_count/float(total)) if total>0 else 0.0
            counts_row = "%s;%s;%s;%d;%d;%d;%.6f" % (
                f.getName(), series_name, m_orig.getTitle(),
                int(total), int(pos_count), int(neg_count), pct
            )
            export_counts_row(csv_counts, counts_row, counts_header)

            # ---------- Overlays speichern ----------
            # DAPI-Overlay: BLAU
            base_dapi = dapi.duplicate()
            # Setze Blau-LUT (robust mit Fallback) und fixiere in RGB
            try: IJ.run(base_dapi, "Blue", "")
            except:
                try: IJ.run(base_dapi, "Blue LUT", "")
                except: pass
            if ENHANCE_OVERLAY_CONTRAST:
                IJ.run(base_dapi, "Enhance Contrast", "saturated=%.3f" % float(DAPI_OVERLAY_SAT))
            IJ.run(base_dapi, "RGB Color", "")
            base_dapi.setOverlay(ov)
            title_dapi = safe_name("%s__%s__PosNeg__DAPI" % (series_name, m_orig.getTitle()))
            out_path_dapi = png_dir.getAbsolutePath() + File.separator + title_dapi + ".png"
            FileSaver(base_dapi).saveAsPng(out_path_dapi)
            if SHOW_OVERLAYS_AT_END:
                base_dapi.show()
            else:
                try: base_dapi.changes=False; base_dapi.close()
                except: pass

            # Marker-Overlay: ORANGE
            base_marker = m_imp.duplicate()
            # Versuche Orange-LUTs, dann Fallbacks, danach RGB fixieren
            applied = False
            for lutName in ["Orange Hot", "mpl-magma", "mpl-inferno"]:
                try:
                    IJ.run(base_marker, lutName, "")
                    applied = True
                    break
                except:
                    pass
            if not applied:
                # breite, warme Fallbacks
                for lutName in ["Fire", "Red Hot"]:
                    try:
                        IJ.run(base_marker, lutName, "")
                        applied = True
                        break
                    except:
                        pass
            if ENHANCE_OVERLAY_CONTRAST:
                IJ.run(base_marker, "Enhance Contrast", "saturated=%.3f" % float(MARKER_OVERLAY_SAT))
            IJ.run(base_marker, "RGB Color", "")
            base_marker.setOverlay(ov)
            title_marker = safe_name("%s__%s__PosNeg__MARKER" % (series_name, m_orig.getTitle()))
            out_path_marker = png_dir.getAbsolutePath() + File.separator + title_marker + ".png"
            FileSaver(base_marker).saveAsPng(out_path_marker)
            if SHOW_OVERLAYS_AT_END:
                base_marker.show()
            else:
                try: base_marker.changes=False; base_marker.close()
                except: pass

        # Aufräumen
        try: work.changes=False; work.close()
        except: pass
        for (_, m_imp, _, _) in marker_preps:
            try: m_imp.changes=False; m_imp.close()
            except: pass
        for w in windows:
            try: w.changes=False; w.close()
            except: pass
        try: imp.changes=False; imp.close()
        except: pass

IJ.log("Fertig. CSVs geschrieben:\n - " + csv_counts + "\n - " + csv_detail + "\nPNGs -> " + png_dir.getAbsolutePath())
