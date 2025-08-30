# Fiji / Jython — AOI & Marker Pos/Neg (2 zytoplasmatische Channels, farbige Overlays + Kombi)
# - Keine Nuclei: ROIs kommen aus AOI-Maske (Connected Components)
# - AOI (Total) -> binäre Maske -> ROIs (mit Größenfilter)
# - Marker: Pos/Neg je ROI nach Anteil Marker-Pixel ≥ MARKER_THR
# - Speichert PNGs: AOI (blau), Marker (orange), COMBO (AOI+Marker), jeweils mit P/N-Konturen/Labels
# Autor: angepasst für Moritz (2025-08-30)

from ij import IJ, WindowManager
from ij.io import DirectoryChooser, FileSaver
from ij.plugin.frame import RoiManager
from ij.gui import Overlay, TextRoi, ShapeRoi
from ij.measure import ResultsTable, Measurements
from ij.plugin.filter import ParticleAnalyzer
from java.io import File, FileWriter, BufferedWriter, PrintWriter
from java.lang import Double
from java.awt import Color
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
import time, re

# ========= Einstellungen =========
# Kanal-Erkennung (erster Treffer zählt)
aoi_keys_str    = IJ.getString("AOI/Total channel identifiers (comma, e.g. 'c1-,total')","c1-,total")
marker_keys_str = IJ.getString("Marker channel identifiers (comma, e.g. 'c2-,c3-')","c2-,marker")

AOI_KEYS    = [s.strip().lower() for s in aoi_keys_str.split(",") if s.strip()]
MARKER_KEYS = [s.strip().lower() for s in marker_keys_str.split(",") if s.strip()]

# Schwellen (Original-Bit-Tiefe)
AOI_THR     = float(IJ.getNumber("AOI fixed threshold (orig bit depth)", 1000.0))
MARKER_THR  = float(IJ.getNumber("Marker fixed threshold (orig bit depth)", 1000.0))

# Optional: BG-Subtraktion (nur AOI/Marker)
APPLY_BG = IJ.getString("Subtract background on AOI/Marker? (yes/no)","yes").lower() in ("yes","y","ja","j","true","1")
if APPLY_BG:
    ROLLING_RADIUS = int(IJ.getNumber("Rolling Ball Radius", 50))
    ROLLING_REPEAT = int(IJ.getNumber("BG repeats (integer)", 1))
    MEDIAN_RADIUS  = int(IJ.getNumber("Median filter radius (0 = none)", 0))
else:
    ROLLING_RADIUS = ROLLING_REPEAT = MEDIAN_RADIUS = 0

# ROI-Filter
SIZE_MIN = int(IJ.getNumber("Min ROI area (px^2)", 80))
SIZE_MAX = int(IJ.getNumber("Max ROI area (px^2)", 1_000_000))

# Positivitäts-Kriterium
POS_FRAC_PCT = float(IJ.getNumber("Min. fraction of ROI pixels ≥ marker thr (%)", 5.0))
POS_FRAC     = POS_FRAC_PCT / 100.0

# Overlays & Export
SHOW_AT_END  = IJ.getString("Show overlays at end? (yes/no)","yes").lower() in ("yes","y","ja","j","true","1")
ADD_LABELS   = IJ.getString("Add P/N labels? (yes/no)","yes").lower() in ("yes","y","ja","j","true","1")
SAT_AOI_PCT  = float(IJ.getNumber("AOI overlay saturation (%)", 1.0))
SAT_MRK_PCT  = float(IJ.getNumber("Marker overlay saturation (%)", 0.35))
ENHANCE      = IJ.getString("Enhance contrast before saving? (yes/no)","no").lower() in ("yes","y","ja","j","true","1")
STROKE_W     = int(IJ.getNumber("Overlay stroke width (px)", 2))

COLOR_POS = Color(0,255,0)
COLOR_NEG = Color(255,0,0)

# ========= Helpers =========
def safe_name(s):
    return re.sub(r'[\\/:*?"<>|]+', '_', s)[:180]

def choose(imgs, keys):
    for im in imgs:
        t = (im.getTitle() or "").lower()
        if any(k in t for k in keys):
            return im
    return None

def preprocess(imp):
    dup = imp.duplicate(); dup.show()
    if APPLY_BG:
        if MEDIAN_RADIUS>0: IJ.run(dup, "Median...", "radius=%d" % MEDIAN_RADIUS)
        for _ in range(max(1, ROLLING_REPEAT)):
            IJ.run(dup, "Subtract Background...", "rolling=%d" % ROLLING_RADIUS)
    return dup, dup.getProcessor()

def apply_blue_rgb(imp):
    try: IJ.run(imp, "Blue", "")
    except:
        try: IJ.run(imp, "Blue LUT", "")
        except: pass
    IJ.run(imp, "RGB Color", "")

def apply_orange_rgb(imp):
    ok=False
    for lut in ["Orange Hot","mpl-inferno","mpl-magma"]:
        try: IJ.run(imp, lut, ""); ok=True; break
        except: pass
    if not ok:
        for lut in ["Fire","Red Hot"]:
            try: IJ.run(imp, lut, ""); ok=True; break
            except: pass
    IJ.run(imp, "RGB Color", "")

def add_label(ov, roi, text, color=Color.white):
    if not ADD_LABELS: return
    b = roi.getBounds()
    tr = TextRoi(b.x + b.width//2, b.y + b.height//2, text)
    tr.setStrokeColor(color); tr.setJustification(TextRoi.CENTER)
    ov.add(tr)

# ========= Pipeline =========
folder = DirectoryChooser("Choose folder with images").getDirectory()
if not folder: IJ.error("No folder selected"); raise SystemExit

out_dir = File(folder + File.separator + "AOI_Marker_Overlays")
if not out_dir.exists(): out_dir.mkdirs()

csv_counts = folder + File.separator + "aoi_marker_posneg_counts.csv"
csv_detail = folder + File.separator + "aoi_marker_posneg_perROI.csv"

def export_csv(path, row, header):
    f = File(path); first = not f.exists()
    pw = PrintWriter(BufferedWriter(FileWriter(f, True)))
    try:
        if first: pw.println(header)
        pw.println(row)
    finally: pw.close()

from java.io import File as JFile
files = sorted([f for f in JFile(folder).listFiles()
               if f.isFile() and f.getName().lower().endswith((".tif",".tiff",".png",".jpg",".jpeg",".lif",".nd2"))],
               key=lambda f: f.getName().lower())

for f in files:
    opts = ImporterOptions(); opts.setId(f.getAbsolutePath()); opts.setOpenAllSeries(True); opts.setVirtual(True)
    imps = BF.openImagePlus(opts)

    for sidx, imp in enumerate(imps):
        series = "%s_Series%d" % (f.getName(), sidx+1) if len(imps)>1 else f.getName()
        IJ.run("Close All")
        imp.setTitle(series); imp.show()
        try: IJ.run(imp, "Split Channels", "")
        except: pass
        time.sleep(0.3)

        imgs = [WindowManager.getImage(i) for i in range(1, WindowManager.getImageCount()+1) if WindowManager.getImage(i)]

        aoi = choose(imgs, AOI_KEYS)
        marker = choose(imgs, MARKER_KEYS)

        if aoi is None or marker is None:
            IJ.log("Missing AOI/Marker in %s" % series)
            for im in imgs:
                try: im.changes=False; im.close()
                except: pass
            try: imp.changes=False; imp.close()
            except: pass
            continue

        # Preprocess AOI/Marker (optional BG)
        aoi_imp, aoi_ip   = preprocess(aoi)
        mrk_imp, mrk_ip   = preprocess(marker)

        # ===== AOI → ROIs =====
        # Binär über festen Schwellwert, dann Connected Components + Größenfilter
        # 1) AOI binär setzen
        aoi_imp.setDisplayRange(0, aoi_imp.getProcessor().getMax())
        aoi_imp.setTitle(series + "__AOIproc")
        aoi_imp.show()
        IJ.setThreshold(aoi_imp, AOI_THR, Double.POSITIVE_INFINITY)
        IJ.run(aoi_imp, "Convert to Mask", "black")
        # 2) Labeling & ROIs
        IJ.run(aoi_imp, "Connected Components Labeling", "connectivity=4")
        rm = RoiManager.getInstance() or RoiManager(); rm.reset()
        rt = ResultsTable()
        pa = ParticleAnalyzer(
            ParticleAnalyzer.ADD_TO_MANAGER,
            Measurements.AREA | Measurements.CIRCULARITY | Measurements.RECT | Measurements.SHAPE_DESCRIPTORS,
            rt, float(SIZE_MIN), float(SIZE_MAX)
        )
        pa.setHideOutputImage(True)
        pa.analyze(aoi_imp)
        rois = list(rm.getRoisAsArray())
        if not rois:
            IJ.log("No ROIs after AOI mask in %s" % series)
            # Aufräumen
            for im in [aoi_imp, mrk_imp] + imgs:
                try: im.changes=False; im.close()
                except: pass
            try: imp.changes=False; imp.close()
            except: pass
            continue

        # ===== Klassifikation pro ROI (Anteil Marker-Pixel >= MARKER_THR) =====
        def roi_iter(ip, roi):
            b = roi.getBounds(); mask = roi.getMask()
            W, H = ip.getWidth(), ip.getHeight()
            for yy in range(b.height):
                gy = b.y + yy
                if gy<0 or gy>=H: continue
                for xx in range(b.width):
                    gx = b.x + xx
                    if gx<0 or gx>=W: continue
                    if mask is not None and mask.get(xx, yy)==0: continue
                    yield gx, gy

        ov = Overlay()
        pos_cnt = neg_cnt = 0
        for idx, r in enumerate(rois, 1):
            total = pospix = 0
            for x,y in roi_iter(mrk_ip, r):
                total += 1
                if mrk_ip.get(x,y) >= MARKER_THR: pospix += 1
            if total == 0: continue
            is_pos = (float(pospix)/float(total) >= POS_FRAC)
            c = r.clone(); c.setStrokeWidth(STROKE_W)
            if is_pos:
                pos_cnt += 1; c.setStrokeColor(COLOR_POS); add_label(ov, r, "P", COLOR_POS)
            else:
                neg_cnt += 1; c.setStrokeColor(COLOR_NEG); add_label(ov, r, "N", COLOR_NEG)
            ov.add(c)

            # Detail-CSV
            export_csv(csv_detail, "%s;%s;%d;%d;%d;%s" % (
                f.getName(), series, idx, total, pospix, "TRUE" if is_pos else "FALSE"
            ), "Image;Series;ROI_Index;ROI_px;PosPix;Is_Positive")

        pct = (100.0*pos_cnt/float(pos_cnt+neg_cnt)) if (pos_cnt+neg_cnt)>0 else 0.0
        export_csv(csv_counts, "%s;%s;%d;%d;%.6f" % (
            f.getName(), series, pos_cnt, neg_cnt, pct
        ), "Image;Series;Positive;Negative;Percent_Positive")

        # ===== Farbige Overlays speichern =====
        # 1) AOI (blau)
        aoi_base = aoi_imp.duplicate()
        apply_blue_rgb(aoi_base)
        if ENHANCE: IJ.run(aoi_base, "Enhance Contrast", "saturated=%.3f" % SAT_AOI_PCT)
        aoi_base.setOverlay(ov)
        FileSaver(aoi_base).saveAsPng(out_dir.getAbsolutePath()+File.separator+safe_name(series+"__AOI_Blue.png"))

        # 2) Marker (orange)
        mrk_base = mrk_imp.duplicate()
        apply_orange_rgb(mrk_base)
        if ENHANCE: IJ.run(mrk_base, "Enhance Contrast", "saturated=%.3f" % SAT_MRK_PCT)
        mrk_base.setOverlay(ov)
        FileSaver(mrk_base).saveAsPng(out_dir.getAbsolutePath()+File.separator+safe_name(series+"__MARKER_Orange.png"))

        # 3) Kombi AOI+Marker
        aoi_rgb = aoi_imp.duplicate(); apply_blue_rgb(aoi_rgb)
        if ENHANCE: IJ.run(aoi_rgb, "Enhance Contrast", "saturated=%.3f" % SAT_AOI_PCT)
        mrk_rgb = mrk_imp.duplicate(); apply_orange_rgb(mrk_rgb)
        if ENHANCE: IJ.run(mrk_rgb, "Enhance Contrast", "saturated=%.3f" % SAT_MRK_PCT)
        aoi_rgb.setTitle("AOI_RGB_TEMP"); mrk_rgb.setTitle("MRK_RGB_TEMP")
        IJ.run("Image Calculator...", "image1=[AOI_RGB_TEMP] operation=Max image2=[MRK_RGB_TEMP] create")
        combo = WindowManager.getCurrentImage()
        combo.setTitle(series+"__COMBO")
        combo.setOverlay(ov)
        FileSaver(combo).saveAsPng(out_dir.getAbsolutePath()+File.separator+safe_name(series+"__COMBO.png"))

        # Aufräumen
        if not SHOW_AT_END:
            for im in [aoi_base, mrk_base, combo, aoi_rgb, mrk_rgb, aoi_imp, mrk_imp] + imgs:
                try: im.changes=False; im.close()
                except: pass
        else:
            aoi_base.show(); mrk_base.show(); combo.show()
            for im in [aoi_rgb, mrk_rgb, aoi_imp, mrk_imp] + imgs:
                try: im.changes=False; im.close()
                except: pass
        try: imp.changes=False; imp.close()
        except: pass

IJ.log("Fertig. CSVs:\n - %s\n - %s\nOverlays → %s" % (csv_counts, csv_detail, out_dir.getAbsolutePath()))
