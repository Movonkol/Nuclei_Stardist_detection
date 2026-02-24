from ij import IJ, WindowManager, ImagePlus
from ij.io import DirectoryChooser, FileSaver
from ij.plugin.frame import RoiManager
from ij.gui import Overlay, TextRoi
from ij.measure import ResultsTable, Measurements
from ij.plugin.filter import ParticleAnalyzer
from ij.process import ImageStatistics as IS, ImageConverter
from ij.plugin import RGBStackMerge
from java.io import File, FileWriter, BufferedWriter, PrintWriter
from java.lang import Double
from java.awt import Color, Font
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
import time, re

# ==========================
#          Settings
# ==========================
SPLIT_BY    = IJ.getString("Split inside script? (none/series/time/z)", "none").lower().strip()
SPLIT_WHICH = IJ.getString("Which half? (first/second)", "first").lower().strip()

model                = IJ.getString("StarDist model", "Versatile (fluorescent nuclei)")
prob                 = IJ.getNumber("Probability threshold", 0.5)
nms                  = IJ.getNumber("NMS threshold", 0.5)
tiles_str            = IJ.getString("nTiles (e.g. '1,1')", "1,1")
channel_patterns_str = IJ.getString("Channel patterns for DAPI (e.g. 'c1-,dapi,blue')", "c1-,dapi,blue")
CHANNEL_PATTERNS     = [p.strip().lower() for p in channel_patterns_str.split(",") if p.strip()]
SCALE_FACTOR         = IJ.getNumber("Pre-scaling (0.5=down, 2.0=up; 1.0=off)", 1.0)
SIGMA_ABOVE_BG       = IJ.getNumber("Min sigma above background (0=off)", 1.5)
MIN_MEAN_INTENSITY   = IJ.getNumber("Min mean intensity (0=off; image units)", 0)

# Nucleus size filter range in square pixels
size_min             = IJ.getNumber("Min label area (px^2)", 50)
size_max             = IJ.getNumber("Max label area (px^2)", 2500)

# AOI channel identifiers and thresholding parameters
aoi_keys_str   = IJ.getString("AOI/Total channel identifiers (comma, e.g. 'c1-,total')", "c1-,total")
AOI_KEYS       = [s.strip().lower() for s in aoi_keys_str.split(',') if s.strip()]
AOI_THR        = float(IJ.getNumber("AOI fixed threshold (orig bit depth)", 1000.0))
APPLY_BG_AOI   = IJ.getString("Subtract background on AOI? (yes/no)", "yes").lower().strip() in ("yes","y","ja","j","true","1")
if APPLY_BG_AOI:
    AOI_ROLLING_RADIUS = int(IJ.getNumber("AOI Rolling Ball Radius", 50))
    AOI_ROLLING_REPEAT = int(IJ.getNumber("AOI BG repeats (integer)", 1))
    AOI_MEDIAN_RADIUS  = int(IJ.getNumber("AOI Median filter radius (0 = none)", 0))
else:
    AOI_ROLLING_RADIUS = AOI_ROLLING_REPEAT = AOI_MEDIAN_RADIUS = 0
AOI_ROI_MIN_FRAC_PCT = float(IJ.getNumber("Min. overlap of nucleus with AOI (%)", 50.0))
AOI_ROI_MIN_FRAC     = AOI_ROI_MIN_FRAC_PCT / 100.0

# Marker channels and positive/negative classification thresholds
marker_keys_str = IJ.getString("Marker channel identifiers (comma-separated, e.g. 'c2-,c3-')", "c2-,marker")
MARKER_KEYS     = [s.strip().lower() for s in marker_keys_str.split(',') if s.strip()][:6]
marker_fixed    = [IJ.getNumber("Fixed threshold for %s (orig bit depth)" % mkey, 1000.0) for mkey in MARKER_KEYS]
POS_FRAC_PCT    = float(IJ.getNumber("Min. fraction of nucleus pixels >= marker thr (%)", 5.0))
POS_FRAC        = POS_FRAC_PCT / 100.0

# Overlay display and export settings
SHOW_AT_END  = IJ.getString("Show overlays at end? (yes/no)","no").lower().strip() in ("yes","y","ja","j","true","1")
ADD_LABELS   = IJ.getString("Add P/N labels? (yes/no)","yes").lower().strip() in ("yes","y","ja","j","true","1")
STROKE_W     = int(IJ.getNumber("Overlay stroke width (px)", 2))
COLOR_POS    = Color(0,255,0)
COLOR_NEG    = Color(255,0,0)
LEGEND_FONT  = Font("SansSerif", Font.PLAIN, 6)
LABEL_FONT   = Font("SansSerif", Font.BOLD, 10)

# ==========================
#      Helper functions
# ==========================
def safe_name(s): return re.sub(r'[\\/:*?"<>|]+', '_', s)[:180]

def list_open_images():
    imgs = []
    try:
        n = WindowManager.getImageCount()
        for i in range(1, n+1):
            im = WindowManager.getImage(i)
            if im is not None: imgs.append(im)
    except: pass
    return imgs

def list_series_images(series):
    key = (series or "").lower()
    return [im for im in list_open_images() if key in ((im.getTitle() or "").lower())]

def pick_dapi_image(patterns, candidates):
    for im in candidates:
        title = (im.getTitle() or "").lower()
        if any(p in title for p in patterns):
            return im
    if len(candidates)==1: return candidates[0]
    # Fallback: use the brightest image if no window title matches the DAPI patterns
    best=None; best_mean=-1.0
    for im in candidates:
        try:
            stats = IS.getStatistics(im.getProcessor(), IS.MEAN, im.getCalibration())
            m = float(stats.mean)
            if m>best_mean: best_mean=m; best=im
        except: pass
    return best if best else (candidates[0] if candidates else None)

def ensure_gray8(imp):
    try:
        if imp and imp.getType() != ImagePlus.GRAY8:
            ImageConverter(imp).convertToGray8()
    except: pass

def roi_centroid(roi):
    try:
        b = roi.getBounds(); m = roi.getMask()
        if m is None:
            return b.x + b.width/2.0, b.y + b.height/2.0
        sx = sy = n = 0
        for yy in range(b.height):
            for xx in range(b.width):
                if m.get(xx, yy) != 0:
                    sx += b.x + xx + 0.5
                    sy += b.y + yy + 0.5
                    n += 1
        if n>0:
            return float(sx)/n, float(sy)/n
        return b.x + b.width/2.0, b.y + b.height/2.0
    except:
        b = roi.getBounds(); return b.x + b.width/2.0, b.y + b.height/2.0

def add_label(ov, roi, text, color=Color.white, font=None):
    if not ADD_LABELS: return
    try:
        cx, cy = roi_centroid(roi)
        tr = TextRoi(0, 0, text)
        tr.setFont(font or LABEL_FONT); tr.setAntialiased(True)
        bb = tr.getBounds()
        tr.setLocation(int(round(cx - bb.width/2.0)), int(round(cy - bb.height/2.0)))
        tr.setStrokeColor(color); tr.setJustification(TextRoi.CENTER); tr.setPosition(1)
        ov.add(tr)
    except: pass

def add_legend(ov, title):
    try:
        lg = TextRoi(5, 5, "Marker: %s\nGreen = Positive\nRed = Negative" % str(title))
        lg.setFont(LEGEND_FONT); lg.setAntialiased(True); lg.setStrokeColor(Color.white); lg.setPosition(1)
        ov.add(lg)
    except: pass

def export_csv(path, row, header):
    f = File(path); first = not f.exists()
    pw = PrintWriter(BufferedWriter(FileWriter(f, True)))
    try:
        if first: pw.println(header)
        pw.println(row)
    finally: pw.close()

def roi_pixel_iter(ip, roi):
    if ip is None or roi is None: return
    b = roi.getBounds(); rmask = roi.getMask()
    W, H = ip.getWidth(), ip.getHeight()
    for yy in range(b.height):
        gy = b.y + yy
        if gy<0 or gy>=H: continue
        for xx in range(b.width):
            gx = b.x + xx
            if gx<0 or gx>=W: continue
            if rmask is not None and rmask.get(xx, yy)==0: continue
            yield gx, gy

def roi_area_pixels(ip, roi):
    c = 0
    for _ in roi_pixel_iter(ip, roi): c += 1
    return c

def roi_mean(ip, roi):
    s = 0; n = 0
    for gx, gy in roi_pixel_iter(ip, roi):
        s += ip.get(gx, gy); n += 1
    return (float(s)/n) if n>0 else 0.0

def circ_from_roi(roi):
    try:
        mask = roi.getMask()
        if mask is None: return Double.NaN
        imp_mask = ImagePlus("roiMask", mask)
        IJ.setThreshold(imp_mask, 1, Double.POSITIVE_INFINITY)
        rtloc = ResultsTable()
        pa = ParticleAnalyzer(ParticleAnalyzer.SHOW_NONE,
                              Measurements.CIRCULARITY,
                              rtloc, 0, Double.POSITIVE_INFINITY)
        pa.setHideOutputImage(True)
        pa.analyze(imp_mask)
        if rtloc.getCounter() > 0:
            return rtloc.getValue("Circ.", 0)
    except: pass
    return Double.NaN

def safe_duplicate(imp, suffix=""):
    try:
        if imp is None: return None
        dup = imp.duplicate()
        if suffix:
            try: dup.setTitle((imp.getTitle() or "dup")+suffix)
            except: pass
        return dup
    except: return None

def close_if_open(imp):
    try:
        if imp is not None: imp.changes=False; imp.close()
    except: pass

def scale_duplicate_to(imp, scale, suffix=""):
    if imp is None: return None
    if abs(scale - 1.0) <= 1e-6:
        return safe_duplicate(imp, suffix)
    dup = safe_duplicate(imp, suffix)
    try:
        orig_w, orig_h = dup.getWidth(), dup.getHeight()
        new_w = max(1, int(round(orig_w * scale)))
        new_h = max(1, int(round(orig_h * scale)))
        IJ.run(dup, "Scale...", "x=%f y=%f width=%d height=%d interpolation=Bilinear average create" % (scale, scale, new_w, new_h))
        scaled = WindowManager.getCurrentImage()
        if scaled is not None:
            close_if_open(dup)
            return scaled
    except:
        pass
    return dup

# ==========================
#     Pipeline: process one image series
# ==========================
def process_series(imp, series, image_name, out_dir,
                   AOI_KEYS, AOI_THR, APPLY_BG_AOI, AOI_ROLLING_RADIUS, AOI_ROLLING_REPEAT, AOI_MEDIAN_RADIUS,
                   size_min, size_max, MIN_MEAN_INTENSITY, SIGMA_ABOVE_BG,
                   POS_FRAC, marker_list, marker_fixed,
                   STROKE_W, COLOR_POS, COLOR_NEG,
                   SHOW_AT_END):

    if imp is None:
        IJ.log("Skip %s (imp=None)" % series); 
        return

    try: imp.setTitle(series); imp.show()
    except: pass
    time.sleep(0.2)

    try: IJ.run(imp, "Split Channels", "")
    except: pass
    time.sleep(0.2)

    series_wins = list_series_images(series)
    if not series_wins:
        IJ.log("No channel windows for %s" % series); 
        close_if_open(imp); 
        return

    dapi = pick_dapi_image(CHANNEL_PATTERNS, series_wins)
    if dapi is None:
        IJ.log("No DAPI in %s" % series)
        close_if_open(imp); 
        for w in series_wins: close_if_open(w)
        return

    # If no marker windows were found, fall back to all non-DAPI channels of matching size
    if not marker_list:
        for win in series_wins:
            if win == dapi: continue
            if win.getWidth()==dapi.getWidth() and win.getHeight()==dapi.getHeight():
                marker_list.append(win)
    while len(marker_fixed) < len(marker_list):
        marker_fixed.append(marker_fixed[-1] if marker_fixed else 1000.0)

    # Search for the AOI channel window using the configured AOI key patterns
    aoi = None
    if AOI_KEYS:
        for win in series_wins:
            t = (win.getTitle() or "").lower()
            if any(k in t for k in AOI_KEYS):
                aoi = win; break

    # ===== StarDist (via Command From Macro) =====
    work = safe_duplicate(dapi, "_work")
    if work is None:
        IJ.log("Work duplicate failed in %s" % series)
        for w in series_wins: close_if_open(w)
        close_if_open(imp)
        return

    # Optionally downscale the working image before running StarDist to speed up detection
    if abs(SCALE_FACTOR - 1.0) > 1e-6:
        orig_w, orig_h = work.getWidth(), work.getHeight()
        new_w = max(1, int(round(orig_w * SCALE_FACTOR)))
        new_h = max(1, int(round(orig_h * SCALE_FACTOR)))
        try:
            IJ.run(work, "Scale...", "x=%f y=%f width=%d height=%d interpolation=Bilinear average create" %
                   (SCALE_FACTOR, SCALE_FACTOR, new_w, new_h))
            scaled = WindowManager.getCurrentImage()
            if scaled is not None:
                close_if_open(work); work = scaled
        except: pass

    # Set a stable window title and ensure the window is visible so StarDist can find it by name
    try:
        work.setTitle(series + "_work")
    except: pass
    work.show()
    time.sleep(0.2)

    # Reset the ROI Manager before running StarDist to avoid leftover ROIs from previous series
    rm = RoiManager.getInstance() or RoiManager()
    try: rm.reset()
    except: pass

    # Run StarDist headless via macro command, outputting both ROIs and label image
    try:
        cmd = ("command=[de.csbdresden.stardist.StarDist2D],"
               "args=['input':'%s','modelChoice':'%s','normalizeInput':'true',"
               "'percentileBottom':'0.0','percentileTop':'100.0','probThresh':'%s','nmsThresh':'%s',"
               "'outputType':'Both','nTiles':'%s','excludeBoundary':'2','verbose':'false',"
               "'showCsbdeepProgress':'false','showProbAndDist':'false'],process=[false]" %
               (work.getTitle(), model, prob, nms, tiles_str))
        IJ.run("Command From Macro", cmd)
    except Exception as e:
        IJ.log("StarDist failed in %s : %s" % (series, str(e)))
        close_if_open(work)
        for w in series_wins: close_if_open(w)
        close_if_open(imp)
        return
    time.sleep(0.4)

    # Close any StarDist dialog windows and label images left open after detection
    try:
        for w in WindowManager.getNonImageWindows():
            try:
                if w.getTitle() and w.getTitle().lower().startswith("stardist"): w.dispose()
            except: pass
    except: pass
    try:
        for imw in list_open_images():
            tt = (imw.getTitle() or "").lower()
            if "label image" in tt: close_if_open(imw)
    except: pass

    rois_array = list((RoiManager.getInstance() or RoiManager()).getRoisAsArray())
    if not rois_array:
        IJ.log("No ROIs from StarDist in %s" % series)
        close_if_open(work)
        for w in series_wins: close_if_open(w)
        close_if_open(imp)
        IJ.run("Collect Garbage", "")
        return

    # ===== Compute background intensity outside nuclei and filter ROIs by size/intensity =====
    dapi_ip = work.getProcessor()
    bgstats = IS.getStatistics(dapi_ip, IS.MEAN | IS.STD_DEV, work.getCalibration())
    bg_mean, bg_sd = float(bgstats.mean), float(bgstats.stdDev)

    area_scale = (SCALE_FACTOR * SCALE_FACTOR) if abs(SCALE_FACTOR - 1.0) > 1e-6 else 1.0
    size_min_eff = float(size_min) * area_scale
    size_max_eff = float(size_max) * area_scale

    kept = []
    for roi in rois_array:
        area_px = roi_area_pixels(dapi_ip, roi)
        if area_px < size_min_eff or area_px > size_max_eff: continue
        mval = roi_mean(dapi_ip, roi)
        if (MIN_MEAN_INTENSITY>0 and mval < MIN_MEAN_INTENSITY): continue
        if (SIGMA_ABOVE_BG>0 and mval < bg_mean + SIGMA_ABOVE_BG*bg_sd): continue
        cval = circ_from_roi(roi)
        kept.append((roi, {'Circ': cval}))

    rm.reset()
    for r,_ in kept:
        try: rm.addRoi(r)
        except: pass
    if not kept:
        IJ.log("No nuclei after filtering in %s" % series)
        close_if_open(work)
        for w in series_wins: close_if_open(w)
        close_if_open(imp)
        IJ.run("Collect Garbage", "")
        return

    # ===== Build binary AOI mask at StarDist working scale for overlap testing =====
    aoi_mask_ip = None
    if aoi is not None:
        aoi_proc = scale_duplicate_to(aoi, SCALE_FACTOR, "__AOIproc")
        if aoi_proc is not None:
            aoi_proc.show(); ensure_gray8(aoi_proc)
            try:
                if APPLY_BG_AOI:
                    if AOI_MEDIAN_RADIUS>0: IJ.run(aoi_proc, "Median...", "radius=%d" % AOI_MEDIAN_RADIUS)
                    for _ in range(max(1, AOI_ROLLING_REPEAT)):
                        IJ.run(aoi_proc, "Subtract Background...", "rolling=%d" % AOI_ROLLING_RADIUS)
                IJ.setThreshold(aoi_proc, AOI_THR, Double.POSITIVE_INFINITY)
                IJ.run(aoi_proc, "Convert to Mask", "black")
                aoi_mask_ip = aoi_proc.getProcessor()
            except: pass

    # ===== Define CSV column headers for all output files =====
    COUNTS_HDR = "Image;Series;Marker;N_total;Positive;Negative;Percent_Positive;Pos_in_AOI;Neg_in_AOI;Percent_Positive_in_AOI"
    DETAIL_HDR = "Image;Series;ROI_Index;Marker;ROI_px;PosPix;Is_Positive;In_AOI"
    MORPH_PERROI_HDR = "Image;Series;ROI_Index;Marker;Circ;Is_Positive;In_AOI"
    MORPH_SUM_HDR    = "Image;Series;N_ROIs;Mean_Circ"
    MORPH_SUM_MARKER_HDR = "Image;Series;Marker;N_ROIs;Mean_Circ;N_Positive;Mean_Circ_Positive;N_In_AOI;Mean_Circ_In_AOI;N_Positive_In_AOI;Mean_Circ_Positive_In_AOI"

    # Global morphology summary across all nuclei, independent of marker classification
    n_circ = 0; sum_circ_all = 0.0
    for _, m in kept:
        if not Double.isNaN(m['Circ']):
            n_circ += 1; sum_circ_all += m['Circ']
    mean_circ = (sum_circ_all / n_circ) if n_circ>0 else 0.0
    export_csv(csv_morph_summary, "%s;%s;%d;%.6f" % (image_name, series, n_circ, mean_circ), MORPH_SUM_HDR)

    # ===== Iterate over each marker channel for positive/negative classification =====
    if not marker_list:
        IJ.log("No marker channels in %s" % series)

    for midx, mwin in enumerate(marker_list[:len(marker_fixed)]):
        thr = float(marker_fixed[midx])
        # Duplicate and rescale the marker image to match the StarDist working resolution
        m_imp = scale_duplicate_to(mwin, SCALE_FACTOR, "__MARKER_%d_GRAY" % (midx+1)); ensure_gray8(m_imp)
        if m_imp is None:
            IJ.log("Marker duplicate failed in %s" % series)
            continue
        m_imp.show()
        m_ip  = m_imp.getProcessor()

        ov = Overlay()
        add_legend(ov, (mwin.getTitle() or "marker"))

        pos_cnt = neg_cnt = 0
        pos_aoi = neg_aoi = 0

        n_all_m = n_pos = n_in = n_posin = 0
        sum_all_m = sum_pos = sum_in = sum_posin = 0.0

        for ridx, (roi, m) in enumerate(kept, 1):
            total = pospix = 0
            for gx, gy in roi_pixel_iter(m_ip, roi):
                total += 1
                if m_ip.get(gx, gy) >= thr: pospix += 1
            if total == 0: 
                continue
            is_pos = (float(pospix)/float(total) >= POS_FRAC)

            in_aoi = False
            if aoi_mask_ip is not None:
                inside = allpix = 0
                for gx, gy in roi_pixel_iter(aoi_mask_ip, roi):
                    allpix += 1
                    if aoi_mask_ip.get(gx, gy) > 0: inside += 1
                in_aoi = (allpix>0 and float(inside)/float(allpix) >= AOI_ROI_MIN_FRAC)

            if is_pos: pos_cnt += 1
            else:      neg_cnt += 1
            if in_aoi:
                if is_pos: pos_aoi += 1
                else:      neg_aoi += 1

            cval = m['Circ']
            if not Double.isNaN(cval):
                n_all_m += 1; sum_all_m += cval
                if is_pos: n_pos += 1; sum_pos += cval
                if in_aoi:
                    n_in += 1; sum_in += cval
                    if is_pos: n_posin += 1; sum_posin += cval

            # Draw nucleus outline on the overlay only if the nucleus overlaps the AOI
            if in_aoi:
                c = roi.clone(); c.setStrokeWidth(STROKE_W)
                if is_pos:
                    c.setStrokeColor(COLOR_POS); add_label(ov, roi, "P", COLOR_POS)
                else:
                    c.setStrokeColor(COLOR_NEG); add_label(ov, roi, "N", COLOR_NEG)
                ov.add(c)

            # Write per-nucleus rows to the detail and morphology CSVs
            export_csv(csv_detail, "%s;%s;%d;%s;%d;%d;%s;%s" % (
                image_name, series, ridx, (mwin.getTitle() or "marker"), total, pospix,
                ("TRUE" if is_pos else "FALSE"), ("TRUE" if in_aoi else "FALSE")
            ), DETAIL_HDR)
            circ_str = ("%.6f" % cval) if not Double.isNaN(cval) else ""
            export_csv(csv_morph_perroi, "%s;%s;%d;%s;%s;%s;%s" % (
                image_name, series, ridx, (mwin.getTitle() or "marker"), circ_str,
                ("TRUE" if is_pos else "FALSE"),
                ("TRUE" if in_aoi else "FALSE")
            ), MORPH_PERROI_HDR)

        pct_total = (100.0*pos_cnt/float(pos_cnt+neg_cnt)) if (pos_cnt+neg_cnt)>0 else 0.0
        pct_aoi   = (100.0*pos_aoi/float(pos_aoi+neg_aoi)) if (pos_aoi+neg_aoi)>0 else 0.0
        export_csv(csv_counts, "%s;%s;%s;%d;%d;%d;%.6f;%d;%d;%.6f" % (
            image_name, series, (mwin.getTitle() or "marker"), (pos_cnt+neg_cnt),
            pos_cnt, neg_cnt, pct_total, pos_aoi, neg_aoi, pct_aoi
        ), COUNTS_HDR)

        # Write per-marker morphology summary row with circularity statistics
        mean_all = (sum_all_m/n_all_m) if n_all_m>0 else 0.0
        mean_pos = (sum_pos/n_pos) if n_pos>0 else 0.0
        mean_in  = (sum_in/n_in) if n_in>0 else 0.0
        mean_posin = (sum_posin/n_posin) if n_posin>0 else 0.0
        export_csv(csv_morph_summary_marker, "%s;%s;%s;%d;%.6f;%d;%.6f;%d;%.6f;%d;%.6f" % (
            image_name, series, (mwin.getTitle() or "marker"), n_all_m, mean_all, n_pos, mean_pos, n_in, mean_in, n_posin, mean_posin
        ), MORPH_SUM_MARKER_HDR)

        # ==== PNGs ====
        flat = combo = flat_back = None
        dapi_gray = scale_duplicate_to(dapi, SCALE_FACTOR, "__DAPI_scaled")
        mrk_gray  = scale_duplicate_to(mwin, SCALE_FACTOR, "__MRK_scaled")
        ok = (dapi_gray is not None and mrk_gray is not None and
              dapi_gray.getWidth()==mrk_gray.getWidth() and dapi_gray.getHeight()==mrk_gray.getHeight())
        if ok:
            ensure_gray8(dapi_gray); ensure_gray8(mrk_gray)
            mrk_R = safe_duplicate(mrk_gray)
            mrk_G = safe_duplicate(mrk_gray)
            try: IJ.run(mrk_G, "Multiply...", "value=0.60")
            except: pass
            try:
                combo = RGBStackMerge.mergeChannels([mrk_R, mrk_G, dapi_gray], False)  # [R,G,B]
                combo.setTitle(series+"__SCALED_PosNeg")
                combo.setOverlay(ov)
                flat = combo.flatten()
                if abs(SCALE_FACTOR - 1.0) > 1e-6:
                    orig_w, orig_h = dapi.getWidth(), dapi.getHeight()
                    IJ.run(flat, "Scale...", "x=%f y=%f width=%d height=%d interpolation=Bilinear average create" % (1.0/SCALE_FACTOR, 1.0/SCALE_FACTOR, orig_w, orig_h))
                    flat_back = WindowManager.getCurrentImage()
                else:
                    flat_back = flat
                FileSaver(flat_back).saveAsPng(
                    out_dir.getAbsolutePath() + File.separator +
                    safe_name("%s__%s__ORIG_PosNeg.png" % (series, (mwin.getTitle() or "marker")))
                )
            finally:
                for imx in [flat, combo, dapi_gray, mrk_gray, mrk_R, mrk_G]:
                    if imx is not flat_back: close_if_open(imx)
                close_if_open(flat_back)
        else:
            IJ.log("Skip PNG (size mismatch) in %s" % series)

        if not SHOW_AT_END:
            close_if_open(m_imp)
            IJ.run("Collect Garbage", "")

    # Close the AOI processing window that was opened during mask generation
    try:
        for imx in list_open_images():
            if imx and ("__AOIproc" in (imx.getTitle() or "")):
                close_if_open(imx)
    except: pass

    # Release all resources for this series and run garbage collection
    close_if_open(work)
    if not SHOW_AT_END:
        for w in series_wins: close_if_open(w)
        close_if_open(imp)
    IJ.run("Collect Garbage", "")

# ==========================
#      Main process: iterate over all input files
# ==========================
folder = DirectoryChooser("Select folder with images").getDirectory()
if not folder:
    IJ.error("No folder selected"); raise SystemExit

out_dir = File(folder + File.separator + "Nuclei_AOI_PosNeg_Overlays")
if not out_dir.exists(): out_dir.mkdirs()

# CSVs
csv_counts        = folder + File.separator + "nuclei_marker_posneg_counts.csv"
csv_detail        = folder + File.separator + "nuclei_marker_posneg_perROI.csv"
csv_morph_perroi  = folder + File.separator + "nuclei_marker_morphology_perROI.csv"
csv_morph_summary = folder + File.separator + "nuclei_marker_morphology.csv"
csv_morph_summary_marker = folder + File.separator + "nuclei_marker_morphology_summary_marker.csv"

from java.io import File as JFile
all_files = []
try:
    for ff in JFile(folder).listFiles():
        name = ff.getName().lower()
        if ff.isFile() and name.endswith((".tif",".tiff",".png",".jpg",".jpeg",".lif",".nd2")):
            all_files.append(ff)
except: pass
files = sorted(all_files, key=lambda f: f.getName().lower())

for f in files:
    image_name = f.getName()
    try:
        opts = ImporterOptions()
        opts.setId(f.getAbsolutePath())
        opts.setOpenAllSeries(True)
        opts.setVirtual(True)
        imps = BF.openImagePlus(opts)
    except:
        IJ.log("Open failed: %s" % image_name)
        continue
    if not imps:
        IJ.log("No series: %s" % image_name)
        continue

    # Derive unique series names, falling back to the image window title if unavailable
    series_names = []
    seen = set()
    for i, imp in enumerate(imps):
        nm = (imp.getTitle() or "Series%d" % (i+1))
        nm = re.sub(r'[\\/:*?"<>|]+', '_', nm).strip()
        base = nm; k=1
        while nm.lower() in seen:
            k+=1; nm = "%s_%d" % (base, k)
        seen.add(nm.lower())
        series_names.append(nm)

    # Optionally process only the first or second half of the series list
    if SPLIT_BY == "series" and len(imps) > 1:
        total = len(imps)
        cut   = total // 2
        sel = range(0, cut) if SPLIT_WHICH=="first" else range(cut, total)
    else:
        sel = range(len(imps))

    for sidx in sel:
        imp = imps[sidx]
        series = series_names[sidx]

        marker_list = []

        # Split stack along time or z axis and process the selected sub-stack
        if SPLIT_BY in ("time","z"):
            nT = max(1, imp.getNFrames())
            nZ = max(1, imp.getNSlices())
            if SPLIT_BY == "time" and nT > 1:
                half = nT // 2
                rng = ("1-%d" % half) if SPLIT_WHICH=="first" else ("%d-%d" % (half+1, nT))
                ch = imp.getNChannels(); zall = imp.getNSlices()
                cmd = "title=%s__part duplicate channels=1-%d slices=1-%d frames=%s" % (series, ch, zall, rng)
                IJ.run(imp, "Duplicate...", cmd)
                part = WindowManager.getCurrentImage()
                process_series(part, series + "__T_%s" % SPLIT_WHICH, image_name, out_dir,
                               AOI_KEYS, AOI_THR, APPLY_BG_AOI, AOI_ROLLING_RADIUS, AOI_ROLLING_REPEAT, AOI_MEDIAN_RADIUS,
                               size_min, size_max, MIN_MEAN_INTENSITY, SIGMA_ABOVE_BG,
                               POS_FRAC, marker_list, marker_fixed,
                               STROKE_W, COLOR_POS, COLOR_NEG, SHOW_AT_END)
                continue
            elif SPLIT_BY == "z" and nZ > 1:
                half = nZ // 2
                rng = ("1-%d" % half) if SPLIT_WHICH=="first" else ("%d-%d" % (half+1, nZ))
                ch = imp.getNChannels(); tall = imp.getNFrames()
                cmd = "title=%s__part duplicate channels=1-%d slices=%s frames=1-%d" % (series, ch, rng, tall)
                IJ.run(imp, "Duplicate...", cmd)
                part = WindowManager.getCurrentImage()
                process_series(part, series + "__Z_%s" % SPLIT_WHICH, image_name, out_dir,
                               AOI_KEYS, AOI_THR, APPLY_BG_AOI, AOI_ROLLING_RADIUS, AOI_ROLLING_REPEAT, AOI_MEDIAN_RADIUS,
                               size_min, size_max, MIN_MEAN_INTENSITY, SIGMA_ABOVE_BG,
                               POS_FRAC, marker_list, marker_fixed,
                               STROKE_W, COLOR_POS, COLOR_NEG, SHOW_AT_END)
                continue

        # No splitting: process the full series directly
        process_series(imp, series, image_name, out_dir,
                       AOI_KEYS, AOI_THR, APPLY_BG_AOI, AOI_ROLLING_RADIUS, AOI_ROLLING_REPEAT, AOI_MEDIAN_RADIUS,
                       size_min, size_max, MIN_MEAN_INTENSITY, SIGMA_ABOVE_BG,
                       POS_FRAC, marker_list, marker_fixed,
                       STROKE_W, COLOR_POS, COLOR_NEG, SHOW_AT_END)

    IJ.run("Collect Garbage", "")

# Log final completion message with paths to all output files
IJ.log("Done. CSVs:")
IJ.log(" - " + csv_counts)
IJ.log(" - " + csv_detail)
IJ.log(" - " + csv_morph_perroi)
IJ.log(" - " + csv_morph_summary)
IJ.log(" - " + csv_morph_summary_marker)
IJ.log("Overlays -> " + out_dir.getAbsolutePath())
