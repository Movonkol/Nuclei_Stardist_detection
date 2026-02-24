# Fiji / Jython: Measure percentage area with optional auto-threshold, stack splitting (series/time/z, first/second half), CSV output, and colored mask/overlay export
from ij import IJ, WindowManager, ImagePlus
from ij.io import DirectoryChooser, FileSaver
from ij.process import ByteProcessor, ImageProcessor
from ij.gui import Overlay
from java.awt import Color
from java.io import File, FileWriter, BufferedWriter, PrintWriter
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
import time, re

# ====================== One-time user prompts at script start =====================
# Channel identifiers: comma-separated substrings matched against window titles
total_patterns_str = IJ.getString("Total/mask channel (e.g. 'c1-,vimentin')", "c1-,vimentin")
TOTAL_CHANNEL_KEYS = [p.strip().lower() for p in total_patterns_str.split(",") if p.strip()]
pos_patterns_str   = IJ.getString("Positive marker channel (e.g. 'c2-,asma')", "c2-,asma")
POS_CHANNEL_KEYS   = [p.strip().lower() for p in pos_patterns_str.split(",") if p.strip()]

# Options to split the stack by series, time, or z dimension
SPLIT_BY    = IJ.getString("Split? (none/series/time/z)", "none").lower().strip()
SPLIT_WHICH = IJ.getString("Which half? (first/second)", "first").lower().strip()

# Choose between automatic (method-based) or fixed threshold values
AUTO_THRESHOLD = IJ.getString("Use auto-threshold? (yes/no)", "yes").lower().strip() in ("ja","j","yes","y","true","t","1")
if AUTO_THRESHOLD:
    T_METHOD = IJ.getString("TOTAL: auto-threshold method (e.g. 'Otsu dark')", "Otsu dark")
    T_FACTOR = IJ.getNumber("TOTAL: threshold factor (>1 = stricter)", 1.0)
    P_METHOD = IJ.getString("POS: auto-threshold method (e.g. 'Otsu dark')", "Otsu dark")
    P_FACTOR = IJ.getNumber("POS: threshold factor (>1 = stricter)", 1.0)
    T_FIXED, P_FIXED = 1000.0, 1000.0   # Fallback values (not prompted when auto mode is active)
else:
    T_FIXED = IJ.getNumber("TOTAL: fixed threshold (e.g. 16-bit: 0..65535)", 1000.0)
    P_FIXED = IJ.getNumber("POS: fixed threshold (e.g. 16-bit: 0..65535)",   1000.0)
    T_METHOD = P_METHOD = ""; T_FACTOR = P_FACTOR = 1.0

# Background correction settings: rolling ball subtraction with optional median pre-filter
apply_bg = IJ.getString("Subtract background? (yes/no)", "yes").lower().strip() in ("ja","j","yes","y","true","t","1")
if apply_bg:
    use_custom_bg = IJ.getString("Use custom background parameters? (yes/no)", "yes").lower().strip() in ("ja","j","yes","y","true","t","1")
    if use_custom_bg:
        ROLLING_RADIUS = IJ.getNumber("Rolling Ball Radius", 50)
        ROLLING_REPEAT = int(IJ.getNumber("BG repetitions (integer)", 1))
        MEDIAN_RADIUS  = IJ.getNumber("Median filter radius (0 = none)", 0)
    else:
        ROLLING_RADIUS, ROLLING_REPEAT, MEDIAN_RADIUS = 50, 1, 0
else:
    ROLLING_RADIUS, ROLLING_REPEAT, MEDIAN_RADIUS = 0, 0, 0

# Controls whether binary positivity maps and heatmap copies are saved
SAVE_BINARY_MAP = True
SAVE_HEATMAP    = IJ.getString("Also save heatmap copy? (yes/no)", "no").lower().strip() in ("ja","j","yes","y","true","t","1")

# ====================== Helper functions ======================
def export_area_row(csv_path, image_name, total_name, pos_name, total_area_px, pos_area_px, percent_pos):
    f = File(csv_path); first = not f.exists()
    pw = PrintWriter(BufferedWriter(FileWriter(f, True)))
    if first: pw.println("Image;Total_Channel;Positive_Channel;Total_Area_px;Positive_Area_px;Percent_Positive")
    percent_str = ("%.6f" % float(percent_pos)).replace(".", ",")
    pw.println("%s;%s;%s;%d;%d;%s" % (image_name, total_name, pos_name, int(total_area_px), int(pos_area_px), percent_str))
    pw.close()

def preprocess_and_threshold(imp, fixed_val, use_auto=False, method="", factor=1.0):
    ch = imp.duplicate(); ch.hide()
    if apply_bg:
        if MEDIAN_RADIUS > 0: IJ.run(ch, "Median...", "radius=%d" % int(MEDIAN_RADIUS))
        for _ in range(max(1, int(ROLLING_REPEAT))):
            IJ.run(ch, "Subtract Background...", "rolling=%d" % int(ROLLING_RADIUS))
    ip = ch.getProcessor()
    thr_eff = float(fixed_val)
    if use_auto and method.strip():
        try:
            IJ.setAutoThreshold(ch, method.strip())
            thr_auto = float(ip.getMinThreshold())
            if thr_auto == thr_auto:
                thr_eff = thr_auto * float(factor if factor else 1.0)
        except: pass
    try:
        IJ.log("Channel '%s': min=%d max=%d  (thr=%d%s)" % (
            imp.getTitle(), int(ip.getMin()), int(ip.getMax()),
            int(thr_eff), " [auto %s ×%.2f]" % (method, factor) if use_auto and method else " [fixed]"
        ))
    except: pass
    return ch, thr_eff, ip

def make_mask(ip, thr_eff):
    w, h = ip.getWidth(), ip.getHeight()
    idxset = set(); k = 0; v = int(thr_eff)
    for y in range(h):
        for x in range(w):
            if ip.get(x, y) >= v: idxset.add(k)
            k += 1
    return idxset, w, h

def find_channels_by_keys(windows, keys):
    hits = []
    for win in windows:
        title = (win.getTitle() or "").lower()
        if any(k in title for k in keys): hits.append(win)
    return hits

def is_image(file):
    n = file.getName().lower()
    return file.isFile() and n.endswith((".tif",".tiff",".png",".jpg",".lif",".nd2"))

def ensure_dir(path):
    d = File(path)
    if not d.exists(): d.mkdir()
    return d

def idxset_to_bp(idxset, w, h):
    bp = ByteProcessor(w, h)
    for k in idxset: bp.set(k % w, k // w, 255)
    return bp

def mask_to_roi(bp):
    if bp is None: return None
    imp_mask = ImagePlus("mask", bp.duplicate())
    ip = imp_mask.getProcessor()
    ip.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
    IJ.run(imp_mask, "Create Selection", "")
    roi = imp_mask.getRoi(); imp_mask.close()
    return roi

def save_masks_and_overlay(total_bp, pos_bp, out_dir, base_name):
    aoi_roi = mask_to_roi(total_bp); pos_roi = mask_to_roi(pos_bp)
    if total_bp is not None:
        w, h = total_bp.getWidth(), total_bp.getHeight()
    elif pos_bp is not None:
        w, h = pos_bp.getWidth(), pos_bp.getHeight()
    else:
        return
    def save_single(roi, fill_color, stroke_color, out_path):
        if roi is None: return
        bg = ImagePlus("bg", ByteProcessor(w, h))
        r = roi.clone() if hasattr(roi, "clone") else roi
        r.setStrokeColor(stroke_color); r.setStrokeWidth(1.5); r.setFillColor(fill_color)
        ov = Overlay(); ov.add(r); bg.setOverlay(ov)
        flat = bg.flatten(); FileSaver(flat).saveAsTiff(out_path); flat.close(); bg.close()
    out_aoi = File(out_dir, base_name + "_AOI_red.tif").getAbsolutePath()
    out_pos = File(out_dir, base_name + "_POS_orange.tif").getAbsolutePath()
    out_ovl = File(out_dir, base_name + "_overlay_compare.tif").getAbsolutePath()
    save_single(aoi_roi, Color(255,0,0,255),   Color(255,0,0),   out_aoi)
    save_single(pos_roi, Color(255,165,0,255), Color(255,165,0), out_pos)
    bg = ImagePlus("bg", ByteProcessor(w, h)); ov = Overlay()
    if aoi_roi is not None:
        tr = aoi_roi.clone() if hasattr(aoi_roi, "clone") else aoi_roi
        tr.setStrokeColor(Color(255,0,0)); tr.setStrokeWidth(1.5); tr.setFillColor(Color(255,0,0,110)); ov.add(tr)
    if pos_roi is not None:
        pr = pos_roi.clone() if hasattr(pos_roi, "clone") else pos_roi
        pr.setStrokeColor(Color(255,165,0)); pr.setStrokeWidth(1.5); pr.setFillColor(Color(255,165,0,230)); ov.add(pr)
    bg.setOverlay(ov); flat = bg.flatten(); FileSaver(flat).saveAsTiff(out_ovl); flat.close(); bg.close()

def save_pos_map(total_mask, pos_ip, p_thr, w, h, out_dir, base_name):
    bp = ByteProcessor(w, h); k = 0; pos_px = 0; v = int(p_thr)
    for y in range(h):
        for x in range(w):
            if k in total_mask and pos_ip.get(x, y) >= v:
                bp.set(x, y, 255); pos_px += 1
            k += 1
    imp_map = ImagePlus(base_name + "_POSmap", bp)
    if SAVE_BINARY_MAP:
        FileSaver(imp_map).saveAsTiff(File(out_dir, base_name + "_POSmap_binary.tif").getAbsolutePath())
    if SAVE_HEATMAP:
        imp_heat = imp_map.duplicate(); IJ.run(imp_heat, "Fire", "")
        FileSaver(imp_heat).saveAsTiff(File(out_dir, base_name + "_POSmap_heat.tif").getAbsolutePath()); imp_heat.close()
    imp_map.close()
    return pos_px, bp

# ---- Helper to split a stack into first or second half along time or z axis ----
def duplicate_half(imp, by, which):
    by = (by or "none").lower().strip()
    which = (which or "first").lower().strip()
    if by not in ("time","z"):  # series splitting is handled outside this function
        return imp
    nC = max(1, imp.getNChannels()); nZ = max(1, imp.getNSlices()); nT = max(1, imp.getNFrames())
    if by == "time" and nT > 1:
        half = nT // 2
        rng = ("1-%d" % half) if which == "first" else ("%d-%d" % (half+1, nT))
        cmd = "title=%s__T_%s duplicate channels=1-%d slices=1-%d frames=%s" % (imp.getTitle(), which, nC, nZ, rng)
    elif by == "z" and nZ > 1:
        half = nZ // 2
        rng = ("1-%d" % half) if which == "first" else ("%d-%d" % (half+1, nZ))
        cmd = "title=%s__Z_%s duplicate channels=1-%d slices=%s frames=1-%d" % (imp.getTitle(), which, nC, rng, nT)
    else:
        return imp
    IJ.run(imp, "Duplicate...", cmd)
    part = WindowManager.getCurrentImage()
    return part if part is not None else imp

# ====================== File I/O and output path setup ======================
dc = DirectoryChooser("Select folder with images")
folder = dc.getDirectory()
if not folder:
    IJ.error("No folder selected"); raise SystemExit

files = sorted([f for f in File(folder).listFiles() if is_image(f)], key=lambda f: f.getName().lower())
csv_file      = folder + File.separator + "percent_area_measurements.csv"
well_csv_file = folder + File.separator + "percent_area_measurements_per_well.csv"
maps_dir      = ensure_dir(folder + File.separator + "maps")

# ============================ Main processing ============================
for f in files:
    well_total_area = 0; well_pos_area = 0
    well_total_name = ""; well_pos_name = ""

    opts = ImporterOptions(); opts.setId(f.getAbsolutePath()); opts.setOpenAllSeries(True); opts.setVirtual(True)
    imps = BF.openImagePlus(opts)
    if not imps: continue

    # --- Series split: select only the first or second half of all series ---
    series_indices = range(len(imps))
    if SPLIT_BY == "series" and len(imps) > 1:
        total = len(imps); cut = total // 2
        series_indices = range(0, cut) if SPLIT_WHICH == "first" else range(cut, total)

    for sidx in series_indices:
        imp = imps[sidx]
        series_name = "%s_Series%d" % (f.getName(), sidx + 1) if len(imps) > 1 else f.getName()
        imp.setTitle(series_name); imp.show()

        # --- Time/Z split: duplicate the selected sub-stack and analyse it ---
        part = duplicate_half(imp, SPLIT_BY, SPLIT_WHICH)
        if part != imp:
            try: imp.changes = False; imp.close()
            except: pass
            imp = part
            series_name = imp.getTitle()

        # --- Split channels and identify total/positive channel windows ---
        IJ.run(imp, "Split Channels", ""); time.sleep(0.3)
        windows = [WindowManager.getImage(i) for i in range(1, WindowManager.getImageCount() + 1) if WindowManager.getImage(i) is not None]
        total_chs = find_channels_by_keys(windows, TOTAL_CHANNEL_KEYS)
        pos_chs   = find_channels_by_keys(windows,   POS_CHANNEL_KEYS)
        if not total_chs or not pos_chs:
            IJ.log("Warning: total or positive channel missing in %s" % series_name)
            WindowManager.closeAllWindows(); continue
        t_win, p_win = total_chs[0], pos_chs[0]

        if well_total_name == "" and well_pos_name == "":
            well_total_name = t_win.getTitle(); well_pos_name = p_win.getTitle()

        # --- Threshold and measure the total (AOI) channel ---
        total_img, t_thr, t_ip = preprocess_and_threshold(t_win, T_FIXED, use_auto=AUTO_THRESHOLD, method=T_METHOD, factor=T_FACTOR)
        total_mask, w, h       = make_mask(t_ip, t_thr)
        total_area_px          = len(total_mask)
        total_bp               = idxset_to_bp(total_mask, w, h)

        # --- Threshold the positive marker channel and count pixels inside the total mask ---
        pos_img, p_thr, p_ip   = preprocess_and_threshold(p_win, P_FIXED, use_auto=AUTO_THRESHOLD, method=P_METHOD, factor=P_FACTOR)
        base = "%s__Total[%s]__Pos[%s]" % (series_name, t_win.getTitle(), p_win.getTitle())
        pos_in_total, pos_bp   = (0, None) if total_area_px == 0 else save_pos_map(total_mask, p_ip, p_thr, w, h, maps_dir, base)

        # --- Save colored mask images and a combined comparison overlay ---
        try: save_masks_and_overlay(total_bp, pos_bp, maps_dir, base)
        except Exception as e: IJ.log("Mask/overlay save failed in %s: %s" % (series_name, str(e)))

        # --- Write per-series row to the measurements CSV ---
        percent_pos = 0.0 if total_area_px == 0 else (100.0 * float(pos_in_total) / float(total_area_px))
        export_area_row(csv_file, series_name, t_win.getTitle(), p_win.getTitle(), total_area_px, pos_in_total, percent_pos)

        # --- Accumulate per-well totals across all series of this file ---
        well_total_area += total_area_px
        well_pos_area   += (0 if pos_in_total is None else pos_in_total)

        # --- Close temporary images for this series ---
        try: total_img.changes = False; total_img.close()
        except: pass
        try: pos_img.changes = False; pos_img.close()
        except: pass
        for win in windows:
            try: win.changes = False; win.close()
            except: pass
        try: imp.changes = False; imp.close()
        except: pass
        WindowManager.closeAllWindows()

    # --- Write per-well summary row aggregating all series of this file ---
    well_percent = 0.0 if well_total_area == 0 else (100.0 * float(well_pos_area) / float(well_total_area))
    export_area_row(well_csv_file, f.getName(), well_total_name, well_pos_name, well_total_area, well_pos_area, well_percent)

