# Fiji / Jython — Nuclei Pos/Neg nur innerhalb ausgewählter "positiver Area"
# - Feste Schwellen in Original-Bit-Tiefe (AOI & Marker), pro Marker wird ein Wert abgefragt
# - Overlays: grün=positiv, rot=negativ; optional P/N-Labels
# - PNGs werden in separatem Unterordner gespeichert
# - Für jeden Marker werden ZWEI Overlays gespeichert: auf DAPI & auf dem Markerkanal
# - Sichere Dateinamen + FileSaver
# Autor: für Moritz, 2025-08-20

from ij import IJ, WindowManager, ImagePlus
from ij.io import DirectoryChooser, FileSaver
from ij.plugin.frame import RoiManager
from ij.gui import Overlay, TextRoi
from ij.process import ByteProcessor
from java.io import File, FileWriter, BufferedWriter, PrintWriter
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

# AOI/Total-Kanal (definiert die 'positive Area' – nur dort wird klassifiziert)
aoi_keys_str = IJ.getString("AOI/Total channel identifier(s) (z.B. 'c1-,vimentin')", "c1-,vimentin")
AOI_CHANNEL_KEYS = [s.strip().lower() for s in aoi_keys_str.split(',') if s.strip()]

# AOI-Threshold in Original-Bit-Tiefe (z.B. 16-bit: 0..65535)
AOI_FIXED_THR = IJ.getNumber("AOI (Total) fixed threshold (Original bit depth)", 1000.0)

# Marker-Thresholds — immer pro Marker (Original-Bit-Tiefe)
marker_fixed = []
for mkey in MARKER_CHANNEL_KEYS:
    marker_fixed.append(IJ.getNumber(
        "Fixed threshold for %s (Original bit depth)" % mkey, float(AOI_FIXED_THR)
    ))

# Optional: Hintergrundkorrektur
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

# StarDist
model = IJ.getString("StarDist Model", "Versatile (fluorescent nuclei)")
prob  = IJ.getNumber("Probability Threshold", 0.5)
nms   = IJ.getNumber("NMS Threshold", 0.5)
n_tile = int(IJ.getString("nTiles (e.g. '1,1')", "1,1").split(',')[0])

# ROI-Filter (Nuklei)
SIZE_MIN = int(IJ.getNumber("Min nucleus area (px)", 100))
SIZE_MAX = int(IJ.getNumber("Max nucleus area (px)", 1500))
ASPECT_RATIO_MAX = IJ.getNumber("Max aspect ratio", 2.0)

# AOI-Abdeckung & Positivitäts-Kriterium
MIN_AOI_COVERAGE_PCT = IJ.getNumber("Min. ROI coverage inside AOI (%)", 50.0)
POS_PIXELS_FRACTION_PCT = IJ.getNumber("Min. fraction of AOI-ROI pixels ≥ marker thr (%)", 5.0)

# Farben
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

def make_mask_from_threshold(ip, thr):
    """Erstellt eine Set-Maske (Indices) aller Pixel >= thr."""
    w, h = ip.getWidth(), ip.getHeight()
    s = set(); k = 0
    for y in range(h):
        for x in range(w):
            if ip.get(x, y) >= thr:
                s.add(k)
            k += 1
    return s, w, h

def roi_pixel_iter(ip, roi):
    """Iteriere Pixelkoordinaten (global x,y) innerhalb eines ROI (auch nicht-rechteckig)."""
    b = roi.getBounds()
    rmask = roi.getMask()  # ByteProcessor oder None (bei Rechteck)
    w = ip.getWidth()
    for yy in range(b.height):
        gy = b.y + yy
        for xx in range(b.width):
            gx = b.x + xx
            if rmask is not None:
                if rmask.get(xx, yy) == 0:
                    continue
            yield gx, gy, gy*w + gx

def roi_area_in_mask(ip, roi, mask_set):
    """Zählt (roi_px, roi_in_aoi_px)"""
    roi_px = 0
    in_aoi = 0
    for gx, gy, k in roi_pixel_iter(ip, roi):
        roi_px += 1
        if k in mask_set:
            in_aoi += 1
    return roi_px, in_aoi

def classify_roi_in_marker(aoi_mask_set, marker_ip, roi, thr, min_pos_frac):
    """Klassifiziere ROI nur im Schnitt ROI∩AOI. Rückgabe: (is_pos, n_intersect, n_pospix)."""
    inter = 0
    pospix = 0
    for gx, gy, k in roi_pixel_iter(marker_ip, roi):
        if k in aoi_mask_set:
            inter += 1
            if marker_ip.get(gx, gy) >= thr:
                pospix += 1
    if inter == 0:
        return False, 0, 0
    return (float(pospix) / float(inter) >= min_pos_frac), inter, pospix

def export_counts_row(csv_path, row, header):
    f = File(csv_path)
    first = not f.exists()
    pw = PrintWriter(BufferedWriter(FileWriter(f, True)))
    if first:
        pw.println(header)
    pw.println(row)
    pw.close()

def safe_name(s):
    s = re.sub(r'[\\/:*?"<>|]+', '_', s)  # Windows-ungültige Zeichen
    return s[:180]  # Pfadlängen-Schutz

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
png_dir = File(folder + File.separator + "AOI_PosNeg_PNGs")
if not png_dir.exists():
    png_dir.mkdirs()

csv_counts = folder + File.separator + "nuclei_posneg_in_AOI.csv"
csv_detail = folder + File.separator + "nuclei_posneg_perROI_in_AOI.csv"

from java.io import File as JFile
files = sorted([f for f in JFile(folder).listFiles()
                if f.isFile() and f.getName().lower().endswith((".tif",".tiff",".png",".jpg",".lif",".nd2"))],
               key=lambda f: f.getName().lower())

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
        aoi_ch = None
        for win in windows:
            t = (win.getTitle() or "").lower()
            if dapi is None and (NUCLEI_CHANNEL_KEY in t):
                dapi = win
            if aoi_ch is None and any(k in t for k in AOI_CHANNEL_KEYS):
                aoi_ch = win
            for mkey in MARKER_CHANNEL_KEYS:
                if mkey in t:
                    markers.append(win)

        if dapi is None or aoi_ch is None or not markers:
            IJ.log("Missing channels in %s (DAPI/AOI/Markers)" % series_name)
            for w in windows:
                try: w.changes=False; w.close()
                except: pass
            try: imp.changes=False; imp.close()
            except: pass
            continue

        # ----- AOI-Maske bauen -----
        aoi_img, aoi_ip = preprocess(aoi_ch)  # Original-Bit-Tiefe
        aoi_mask_set, W, H = make_mask_from_threshold(aoi_ip, float(AOI_FIXED_THR))

        # ----- Nuklei segmentieren (StarDist) -----
        dapi_proc = dapi.duplicate(); dapi_proc.setTitle("DAPI"); dapi_proc.show()
        # Keine 8-bit/Enhance nötig für StarDist
        cmd = ("command=[de.csbdresden.stardist.StarDist2D],"
               "args=['input':'%s','modelChoice':'%s','normalizeInput':'true',"
               "'percentileBottom':'0.0','percentileTop':'100.0',"
               "'probThresh':'%f','nmsThresh':'%f','outputType':'ROI Manager',"
               "'nTiles':'%d','excludeBoundary':'2','verbose':'false',"
               "'showCsbdeepProgress':'false','showProbAndDist':'false'],process=[false]") % (
                    dapi_proc.getTitle(), model, prob, nms, n_tile)
        IJ.run("Command From Macro", cmd); time.sleep(1); close_stardist_dialogs()

        rm = RoiManager.getInstance() or RoiManager()
        rois = list(rm.getRoisAsArray())

        # ROI-Qualitätsfilter
        valid = []
        for roi in rois:
            dapi_proc.setRoi(roi)
            stats = dapi_proc.getStatistics()
            area = stats.area
            b = roi.getBounds()
            ar = max(b.width, b.height) / float(max(1, min(b.width, b.height)))
            if SIZE_MIN <= area <= SIZE_MAX and ar <= ASPECT_RATIO_MAX:
                valid.append(roi)
        rm.reset()
        [rm.addRoi(r) for r in valid]

        if not valid:
            try: dapi_proc.changes=False; dapi_proc.close()
            except: pass
            try: aoi_img.changes=False; aoi_img.close()
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

        # ----- Klassifikation pro ROI (nur in AOI) -----
        min_cov = float(MIN_AOI_COVERAGE_PCT) / 100.0
        min_pos_frac = float(POS_PIXELS_FRACTION_PCT) / 100.0

        # CSV-Header
        counts_header = "Image;Series;Marker;N_in_AOI;Positive;Negative;Percent_Positive"
        detail_header = "Image;Series;ROI_Index;Marker;ROI_px;ROI_in_AOI_px;PosPix_in_AOI;Is_Positive"

        # Overlays pro Marker
        for midx, (m_orig, m_imp, m_ip, thr) in enumerate(marker_preps):
            pos_count = 0; neg_count = 0; total_in_aoi = 0
            ov = Overlay()

            # kleine Legende
            legend = TextRoi(5, 5, "Marker: %s\nGreen = Positive\nRed = Negative" % m_orig.getTitle())
            legend.setStrokeColor(Color.white)
            ov.add(legend)

            for ridx, roi in enumerate(valid):
                # AOI-Abdeckung
                roi_px, in_aoi_px = roi_area_in_mask(dapi_proc.getProcessor(), roi, aoi_mask_set)
                if roi_px == 0: continue
                coverage = float(in_aoi_px)/float(roi_px)

                if coverage < min_cov or in_aoi_px == 0:
                    continue  # liegt zu wenig in AOI

                total_in_aoi += 1

                # Positivitätsprüfung nur auf ROI∩AOI
                is_pos, inter, pospix = classify_roi_in_marker(aoi_mask_set, m_ip, roi, thr, min_pos_frac)
                c = roi.clone()
                c.setStrokeWidth(2)
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
                detail_row = "%s;%s;%d;%s;%d;%d;%d;%s" % (
                    f.getName(), series_name, ridx+1, m_orig.getTitle(),
                    int(roi_px), int(in_aoi_px), int(pospix), "TRUE" if is_pos else "FALSE"
                )
                export_counts_row(csv_detail, detail_row, detail_header)

            pct = (100.0 * pos_count/float(total_in_aoi)) if total_in_aoi>0 else 0.0
            counts_row = "%s;%s;%s;%d;%d;%d;%.6f" % (
                f.getName(), series_name, m_orig.getTitle(),
                int(total_in_aoi), int(pos_count), int(neg_count), pct
            )
            export_counts_row(csv_counts, counts_row, counts_header)

            # ---------- Overlays speichern (separater Ordner) ----------
            # 1) DAPI-Overlay
            base_dapi = dapi.duplicate()  # Original DAPI-Kanal
            if ENHANCE_OVERLAY_CONTRAST:
                IJ.run(base_dapi, "Enhance Contrast", "saturated=0.35")
            base_dapi.setOverlay(ov)
            title_dapi = safe_name("%s__%s__AOI_PosNeg__DAPI" % (series_name, m_orig.getTitle()))
            out_path_dapi = png_dir.getAbsolutePath() + File.separator + title_dapi + ".png"
            ok1 = FileSaver(base_dapi).saveAsPng(out_path_dapi)
            if SHOW_OVERLAYS_AT_END: base_dapi.show()
            else:
                try: base_dapi.changes=False; base_dapi.close()
                except: pass

            # 2) Marker-Overlay
            base_marker = m_imp.duplicate()  # vorverarbeiteter Marker (BG subtrahiert, Original-Bit-Tiefe)
            if ENHANCE_OVERLAY_CONTRAST:
                IJ.run(base_marker, "Enhance Contrast", "saturated=0.35")
            base_marker.setOverlay(ov)
            title_marker = safe_name("%s__%s__AOI_PosNeg__MARKER" % (series_name, m_orig.getTitle()))
            out_path_marker = png_dir.getAbsolutePath() + File.separator + title_marker + ".png"
            ok2 = FileSaver(base_marker).saveAsPng(out_path_marker)
            if SHOW_OVERLAYS_AT_END: base_marker.show()
            else:
                try: base_marker.changes=False; base_marker.close()
                except: pass

        # Aufräumen Hilfsfenster
        try: dapi_proc.changes=False; dapi_proc.close()
        except: pass
        try: aoi_img.changes=False; aoi_img.close()
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
