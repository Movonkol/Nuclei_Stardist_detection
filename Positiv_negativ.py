

from ij import IJ, WindowManager, ImagePlus
from ij.io import DirectoryChooser, FileSaver
from ij.plugin.frame import RoiManager
from ij.gui import Overlay, TextRoi
from ij.measure import ResultsTable, Measurements
from ij.plugin.filter import ParticleAnalyzer
from java.io import File, FileWriter, BufferedWriter, PrintWriter
from java.lang import Double
from java.awt import Color
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
import time, re

# ==========================
#        EINSTELLUNGEN
# ==========================
# Kanal-IDs per Teilstring (z.B. 'c1-', 'dapi', 'asma'). Mehrfach möglich, erster Treffer zählt.
NUCLEI_CHANNEL_KEY = IJ.getString("Nuclei (DAPI) channel identifier (z.B. 'c4-', 'dapi')", "c4-").lower()

marker_keys_str = IJ.getString("Marker channel identifiers (comma-separated, z.B. 'c2-,c3-')", "c3-")
MARKER_CHANNEL_KEYS = [s.strip().lower() for s in marker_keys_str.split(',') if s.strip()][:6]

# Marker-Thresholds — immer pro Marker (Original-Bit-Tiefe)
marker_fixed = []
for mkey in MARKER_CHANNEL_KEYS:
    marker_fixed.append(IJ.getNumber(
        "Fixed threshold for %s (Original bit depth)" % mkey, 1000.0
    ))

# Optional: Hintergrundkorrektur (ohne Bit-Tiefenwechsel)
APPLY_BG = IJ.getString("Subtract background? (yes/no)", "yes").lower().strip() in ("yes","y","ja","j","true","1")
if APPLY_BG:
    ROLLING_RADIUS = IJ.getNumber("Rolling Ball Radius", 50)
    ROLLING_REPEAT = int(IJ.getNumber("BG repeats (integer)", 1))
    MEDIAN_RADIUS  = IJ.getNumber("Median filter radius (0 = none)", 0)
else:
    ROLLING_RADIUS, ROLLING_REPEAT, MEDIAN_RADIUS = 0, 0, 0

# Anzeige-/Speicheroptionen
SHOW_OVERLAYS_AT_END = IJ.getString("Show overlay maps at end? (yes/no)", "yes").lower().strip() in ("yes","y","ja","j","true","1")
ADD_PN_LABELS = IJ.getString("Add P/N labels to overlays? (yes/no)", "yes").lower().strip() in ("yes","y","ja","j","true","1")
ENHANCE_OVERLAY_CONTRAST = IJ.getString("Enhance contrast for saved overlays? (yes/no)", "no").lower().strip() in ("yes","y","ja","j","true","1")
DAPI_OVERLAY_SAT   = IJ.getNumber("DAPI overlay saturation (%)", 1.0)
MARKER_OVERLAY_SAT = IJ.getNumber("Marker overlay saturation (%)", 0.35)

# StarDist (Analyse) — Normalisierung enger für „Sättigung“ in der ANALYSE
model   = IJ.getString("StarDist Model", "Versatile (fluorescent nuclei)")
prob    = IJ.getNumber("Probability Threshold", 0.4)
nms     = IJ.getNumber("NMS Threshold", 0.5)
DAPI_PB = IJ.getNumber("DAPI analysis lower percentile", 1.0)
DAPI_PT = IJ.getNumber("DAPI analysis upper percentile", 99.8)
N_TILES = IJ.getString("nTiles (e.g. '1,1')", "1,1")

# ROI-Filter (typische Kerne in 2048²@0.07 µm/px liegen ~10k–30k px)
SIZE_MIN = int(IJ.getNumber("Min nucleus area (px)", 3000))
SIZE_MAX = int(IJ.getNumber("Max nucleus area (px)", 50000))

# Positivitäts-Kriterium (ohne AOI)
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

def preprocess(imp):
    """Optional BG-Subtraktion; bleibt in Original-Bit-Tiefe (keine 8-bit-Konvertierung)."""
    dup = imp.duplicate(); dup.show()
    if APPLY_BG:
        if MEDIAN_RADIUS > 0:
            IJ.run(dup, "Median...", "radius=%d" % int(MEDIAN_RADIUS))
        for _ in range(max(1, int(ROLLING_REPEAT))):
            IJ.run(dup, "Subtract Background...", "rolling=%d" % int(ROLLING_RADIUS))
    return dup, dup.getProcessor()

def roi_pixel_iter(ip, roi):
    """Iteriere Pixelkoordinaten (global x,y) innerhalb eines ROI mit hartem Clipping an Bildgrenzen.
    Liefert (gx, gy)."""
    b = roi.getBounds()
    W = ip.getWidth(); H = ip.getHeight()
    rmask = roi.getMask()  # ByteProcessor oder None
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
    """Klassifiziere ROI basierend auf Anteil Marker-Pixel >= thr (robust, geclippt)."""
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
    f = File(csv_path)
    first = not f.exists()
    pw = PrintWriter(BufferedWriter(FileWriter(f, True)))
    if first:
        pw.println(header)
    pw.println(row)
    pw.close()

def safe_name(s):
    s = re.sub(r'[\\/:*?"<>|]+', '_', s)
    return s[:180]

def add_label(ov, roi, text, color=Color.white):
    if not ADD_PN_LABELS:
        return
    b = roi.getBounds()
    cx = b.x + int(b.width/2)
    cy = b.y + int(b.height/2)
    tr = TextRoi(cx, cy, text)
    tr.setStrokeColor(color)
    tr.setJustification(TextRoi.CENTER)
    ov.add(tr)

# ==========================
#        HAUPTPROZESS
# ==========================
folder = DirectoryChooser("Choose folder with images").getDirectory()
if not folder:
    IJ.error("No folder selected"); raise SystemExit

# Unterordner für PNGs
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
        series_name = "%s_Series%d" % (f.getName(), sidx+1)
        IJ.run("Close All")
        imp.setTitle(series_name); imp.show()
        IJ.run(imp, "Split Channels", ""); time.sleep(0.4)

        # Offene Fenster dieser Serie sammeln
        windows = [WindowManager.getImage(i) for i in range(1, WindowManager.getImageCount()+1) if WindowManager.getImage(i) is not None]

        # Kanäle zuordnen
        dapi = None
        markers = []
        for win in windows:
            t = (win.getTitle() or "").lower()
            if dapi is None and (NUCLEI_CHANNEL_KEY in t):
                dapi = win
            for mkey in MARKER_CHANNEL_KEYS:
                if mkey in t:
                    markers.append(win)

        if dapi is None or not markers:
            IJ.log("Missing channels in %s (DAPI and/or Markers)" % series_name)
            for w in windows:
                try: w.changes=False; w.close()
                except: pass
            try: imp.changes=False; imp.close()
            except: pass
            continue

        # ----- Nuklei segmentieren (StarDist -> Label Image -> CC -> Particles->ROI) -----
        dapi_proc = dapi.duplicate(); dapi_proc.setTitle("DAPI"); dapi_proc.show()
        cmd = ("command=[de.csbdresden.stardist.StarDist2D],"
               "args=['input':'%s','modelChoice':'%s','normalizeInput':'true',"
               "'percentileBottom':'%.3f','percentileTop':'%.3f',"
               "'probThresh':'%f','nmsThresh':'%f','outputType':'Label Image',"
               "'nTiles':'%s','excludeBoundary':'2','verbose':'false',"
               "'showCsbdeepProgress':'false','showProbAndDist':'false'],process=[false]") % (
                    dapi_proc.getTitle(), model, DAPI_PB, DAPI_PT, prob, nms, N_TILES)
        IJ.run("Command From Macro", cmd); time.sleep(1); close_stardist_dialogs()

        # Label-Image finden
        label = None
        for i in range(1, WindowManager.getImageCount()+1):
            win = WindowManager.getImage(i)
            if win and "label image" in (win.getTitle() or "").lower():
                label = win; break
        if label is None:
            IJ.log("Label image nicht gefunden für %s" % series_name)
            try: dapi_proc.changes=False; dapi_proc.close()
            except: pass
            for w in windows:
                try: w.changes=False; w.close()
                except: pass
            try: imp.changes=False; imp.close()
            except: pass
            continue

        # Connected Components -> ROIs in den ROI Manager (Größenfilter hier!)
        rm = RoiManager.getInstance() or RoiManager(); rm.reset()
        IJ.selectWindow(label.getTitle())
        IJ.run(label, "Connected Components Labeling", "connectivity=4")
        IJ.setThreshold(label, 1, Double.POSITIVE_INFINITY)
        rt = ResultsTable()
        pa = ParticleAnalyzer(
            ParticleAnalyzer.ADD_TO_MANAGER,
            Measurements.AREA | Measurements.CENTROID | Measurements.RECT | Measurements.SHAPE_DESCRIPTORS,
            rt, SIZE_MIN, SIZE_MAX
        )
        pa.setHideOutputImage(True)
        pa.analyze(label)
        valid = list(rm.getRoisAsArray())

        if not valid:
            try: dapi_proc.changes=False; dapi_proc.close()
            except: pass
            for w in windows:
                try: w.changes=False; w.close()
                except: pass
            try: imp.changes=False; imp.close()
            except: pass
            continue

        # ----- Marker vorbereiten (Original-Bit-Tiefe, feste Schwellen) -----
        marker_preps = []
        for idx, m in enumerate(markers[:len(MARKER_CHANNEL_KEYS)]):
            proc_imp, proc_ip = preprocess(m)  # bleibt in Original-Bit-Tiefe
            thr = float(marker_fixed[idx])     # fester Wert je Marker
            marker_preps.append((m, proc_imp, proc_ip, thr))

        # CSV-Header
        counts_header = "Image;Series;Marker;N_total;Positive;Negative;Percent_Positive"
        detail_header = "Image;Series;ROI_Index;Marker;ROI_px;PosPix;Is_Positive"

        # ----- Klassifikation & Overlays pro Marker -----
        for midx, (m_orig, m_imp, m_ip, thr) in enumerate(marker_preps):
            pos_count = 0; neg_count = 0; total = 0
            ov = Overlay()

            # kleine Legende
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
                    c.setStrokeColor(COLOR_POS)
                    add_label(ov, roi, "P", COLOR_POS)
                else:
                    neg_count += 1
                    c.setStrokeColor(COLOR_NEG)
                    add_label(ov, roi, "N", COLOR_NEG)
                ov.add(c)

                # Detailzeile
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

            # ---------- Overlays speichern (separater Ordner) ----------
            # 1) DAPI-Overlay
            base_dapi = dapi.duplicate()  # Original DAPI-Kanal (Analyse unverändert)
            if ENHANCE_OVERLAY_CONTRAST:
                IJ.run(base_dapi, "Enhance Contrast", "saturated=%.3f" % DAPI_OVERLAY_SAT)
            base_dapi.setOverlay(ov)
            title_dapi = safe_name("%s__%s__PosNeg__DAPI" % (series_name, m_orig.getTitle()))
            out_path_dapi = png_dir.getAbsolutePath() + File.separator + title_dapi + ".png"
            FileSaver(base_dapi).saveAsPng(out_path_dapi)
            if SHOW_OVERLAYS_AT_END:
                base_dapi.show()
            else:
                try: base_dapi.changes=False; base_dapi.close()
                except: pass

            # 2) Marker-Overlay
            base_marker = m_imp.duplicate()  # vorverarbeiteter Marker (BG subtrahiert, Original-Bit-Tiefe)
            if ENHANCE_OVERLAY_CONTRAST:
                IJ.run(base_marker, "Enhance Contrast", "saturated=%.3f" % MARKER_OVERLAY_SAT)
            base_marker.setOverlay(ov)
            title_marker = safe_name("%s__%s__PosNeg__MARKER" % (series_name, m_orig.getTitle()))
            out_path_marker = png_dir.getAbsolutePath() + File.separator + title_marker + ".png"
            FileSaver(base_marker).saveAsPng(out_path_marker)
            if SHOW_OVERLAYS_AT_END:
                base_marker.show()
            else:
                try: base_marker.changes=False; base_marker.close()
                except: pass

        # Aufräumen Hilfsfenster
        try: dapi_proc.changes=False; dapi_proc.close()
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
