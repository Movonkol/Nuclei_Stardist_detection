# Fiji/Jython — StarDist Nuclei within AOI (no marker)
# - COMBO: AOI rot (R) + DAPI blau (B) via RGBStackMerge, Nuklei-Konturen obenauf
# - DAPI-only Overlay zusätzlich
# - Robustheit: sichere LUT-Aufrufe, Index-Checks, None-Guards, GC
# - CSV (per-ROI): Area_px, Circ, In_AOI
# - CSV (Summary): N_ROIs, N_In_AOI, Mean_Area_px, Mean_Circ, Mean_Circ_In_AOI

from ij import IJ, WindowManager, ImagePlus
from ij.io import DirectoryChooser, FileSaver
from ij.plugin.frame import RoiManager
from ij.gui import Overlay, PolygonRoi, Roi
from ij.measure import ResultsTable, Measurements
from ij.plugin.filter import ParticleAnalyzer
from ij.process import ImageStatistics as IS, ImageConverter, ByteProcessor
from ij.plugin import RGBStackMerge
from java.io import File, FileWriter, BufferedWriter, PrintWriter
from java.lang import Double
from java.awt import Color
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
import time, re, jarray

# -------- Einstellungen --------
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
size_min             = IJ.getNumber("Min nucleus area (px^2)", 50)
size_max             = IJ.getNumber("Max nucleus area (px^2)", 2500)

aoi_keys_str   = IJ.getString("AOI channel identifiers (comma, e.g. 'c2-,total')", "total")
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

SAT_DAPI_PCT = float(IJ.getNumber("DAPI overlay saturation (%)", 1.0))
ENHANCE      = IJ.getString("Enhance contrast before saving? (yes/no)","no").lower().strip() in ("yes","y","ja","j","true","1")
STROKE_W     = int(IJ.getNumber("Overlay stroke width (px)", 2))
SHOW_ONLY_IN_AOI = True

# -------- Helpers --------
def safe_name(s): return re.sub(r'[\\/:*?"<>|]+', '_', s)[:180]

def list_open_images():
    imgs=[]
    for i in range(1, WindowManager.getImageCount()+1):
        im=WindowManager.getImage(i)
        if im: imgs.append(im)
    return imgs

def list_series_images(series):
    key=(series or "").lower()
    return [im for im in list_open_images() if key in ((im.getTitle() or "").lower())]

def pick_dapi_image(patterns, candidates):
    for im in candidates:
        if any(p in (im.getTitle() or "").lower() for p in patterns):
            return im
    return candidates[0] if candidates else None

def ensure_gray8(imp):
    try:
        if imp and imp.getType()!=ImagePlus.GRAY8:
            ImageConverter(imp).convertToGray8()
    except: pass

def apply_blue_rgb(imp):
    if imp is None: return
    try:
        IJ.run(imp, "Blue", "")
    except:
        try:
            IJ.run(imp, "Blue LUT", "")
        except:
            pass
    try:
        IJ.run(imp, "RGB Color", "")
    except: pass

def export_csv(path, row, header):
    f=File(path); first=not f.exists()
    pw=PrintWriter(BufferedWriter(FileWriter(f,True)))
    try:
        if first: pw.println(header)
        pw.println(row)
    finally: pw.close()

def roi_pixel_iter(ip, roi):
    if ip is None or roi is None: return
    b=roi.getBounds(); rmask=roi.getMask()
    W,H=ip.getWidth(), ip.getHeight()
    for yy in range(b.height):
        gy=b.y+yy
        if gy<0 or gy>=H: continue
        for xx in range(b.width):
            gx=b.x+xx
            if gx<0 or gx>=W: continue
            if rmask is not None and rmask.get(xx,yy)==0: continue
            yield gx,gy

def roi_area_pixels(ip, roi):
    c=0
    for _ in roi_pixel_iter(ip, roi): c+=1
    return c

def roi_mean(ip, roi):
    s=0; n=0
    for gx,gy in roi_pixel_iter(ip, roi):
        s+=ip.get(gx,gy); n+=1
    return (float(s)/n) if n>0 else 0.0

def close_if_open(imp):
    try:
        if imp is not None: imp.changes=False; imp.close()
    except: pass

# Circularity robust über ParticleAnalyzer auf ROI-Maske
def circ_from_roi(roi):
    try:
        mask = roi.getMask()
        if mask is None:
            return Double.NaN
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
    except:
        pass
    return Double.NaN

# ROI aus work-Skalierung zurück auf Originalgröße
def rescale_roi_to_original(roi, sx, sy):
    fp = roi.getFloatPolygon()
    if fp is None or fp.npoints==0:
        b = roi.getBounds()
        nx = int(round(b.x * sx)); ny = int(round(b.y * sy))
        nw = max(1, int(round(b.width * sx))); nh = max(1, int(round(b.height * sy)))
        return Roi(nx, ny, nw, nh)
    xs = jarray.array([fp.xpoints[i] * sx for i in range(fp.npoints)], 'f')
    ys = jarray.array([fp.ypoints[i] * sy for i in range(fp.npoints)], 'f')
    nr = PolygonRoi(xs, ys, fp.npoints, Roi.POLYGON)
    try: nr.setStrokeWidth(roi.getStrokeWidth())
    except: pass
    return nr

# -------- Serie verarbeiten --------
def process_series(imp, series, folder):
    if imp is None:
        IJ.log("Skip %s (imp=None)" % series); return

    imp.setTitle(series); imp.show(); time.sleep(0.2)
    try: IJ.run(imp, "Split Channels", "")
    except: pass
    time.sleep(0.2)

    series_wins=list_series_images(series)
    if not series_wins:
        IJ.log("No channels opened for %s" % series); close_if_open(imp); return

    dapi=pick_dapi_image(CHANNEL_PATTERNS, series_wins)
    if dapi is None:
        IJ.log("No DAPI in %s" % series)
        for w in series_wins: close_if_open(w); close_if_open(imp); return

    # AOI (optional)
    aoi=None
    if AOI_KEYS:
        for win in series_wins:
            if any(k in (win.getTitle() or "").lower() for k in AOI_KEYS):
                aoi=win; break

    dapi_disp=dapi.duplicate(); ensure_gray8(dapi_disp); dapi_disp.setTitle(series+"_DAPI_GRAY"); dapi_disp.show()

    # --- StarDist → ROI-Manager ---
    orig_w, orig_h = dapi.getWidth(), dapi.getHeight()
    work = dapi.duplicate(); work.setTitle(series + "_work"); work.show()
    rm = RoiManager.getInstance() or RoiManager(); rm.reset()

    work_w, work_h = orig_w, orig_h
    if abs(SCALE_FACTOR-1.0)>1e-6:
        new_w=max(1, int(round(orig_w*SCALE_FACTOR)))
        new_h=max(1, int(round(orig_h*SCALE_FACTOR)))
        IJ.run(work, "Scale...", "x=%f y=%f width=%d height=%d interpolation=Bilinear average create" %
               (SCALE_FACTOR, SCALE_FACTOR, new_w, new_h))
        scaled=WindowManager.getCurrentImage()
        if scaled is not None:
            close_if_open(work); work=scaled
        work_w, work_h = work.getWidth(), work.getHeight()

    cmd=("command=[de.csbdresden.stardist.StarDist2D],"
         "args=['input':'%s','modelChoice':'%s','normalizeInput':'true',"
         "'percentileBottom':'0.0','percentileTop':'100.0','probThresh':'%s','nmsThresh':'%s',"
         "'outputType':'ROI Manager','nTiles':'%s','excludeBoundary':'2','verbose':'false',"
         "'showCsbdeepProgress':'false','showProbAndDist':'false'],process=[false]" %
         (work.getTitle(), model, prob, nms, tiles_str))
    IJ.run("Command From Macro", cmd); time.sleep(0.4)

    rm = RoiManager.getInstance() or RoiManager()
    rois_rm = list(rm.getRoisAsArray()) if rm else []
    if not rois_rm:
        IJ.log("No ROIs from StarDist in %s" % series)
        close_if_open(work); close_if_open(dapi_disp)
        for w in series_wins: close_if_open(w); close_if_open(imp); IJ.run("Collect Garbage",""); return

    # zurückskalieren falls nötig
    if work_w != orig_w or work_h != orig_h:
        sx = float(orig_w) / float(work_w)
        sy = float(orig_h) / float(work_h)
        rois_array = [rescale_roi_to_original(r, sx, sy) for r in rois_rm]
    else:
        rois_array = rois_rm[:]

    # ROI-Manager mit Original-ROIs (nur zur Sicht)
    rm.reset()
    for r in rois_array:
        try: rm.addRoi(r)
        except: pass

    # --- Hintergrund (DAPI) ---
    dapi_ip=dapi.getProcessor()
    if dapi_ip is None:
        IJ.log("No processor for DAPI in %s" % series); return
    bgstats = IS.getStatistics(dapi_ip, IS.MEAN | IS.STD_DEV, dapi.getCalibration())
    bg_mean, bg_sd = float(bgstats.mean), float(bgstats.stdDev)

    # --- AOI: Intensität & Maske ---
    aoi_mask_ip=None; aoi_view=None
    if aoi is not None:
        aoi_proc=aoi.duplicate(); aoi_proc.setTitle(series+"__AOIproc"); aoi_proc.show(); ensure_gray8(aoi_proc)
        if APPLY_BG_AOI:
            if AOI_MEDIAN_RADIUS>0: IJ.run(aoi_proc, "Median...", "radius=%d" % AOI_MEDIAN_RADIUS)
            for _ in range(max(1, AOI_ROLLING_REPEAT)):
                IJ.run(aoi_proc, "Subtract Background...", "rolling=%d" % AOI_ROLLING_RADIUS)
        aoi_view = aoi_proc.duplicate()
        IJ.setThreshold(aoi_proc, AOI_THR, Double.POSITIVE_INFINITY)
        IJ.run(aoi_proc, "Convert to Mask", "black")
        aoi_mask_ip = aoi_proc.getProcessor()

    # --- Filter + Messungen ---
    kept=[]  # (roi, area_px, in_aoi, circ)
    for roi in rois_array:
        area_px=roi_area_pixels(dapi_ip, roi)
        if area_px<size_min or area_px>size_max: continue
        mval=roi_mean(dapi_ip, roi)
        if (MIN_MEAN_INTENSITY>0 and mval<MIN_MEAN_INTENSITY): continue
        if (SIGMA_ABOVE_BG>0 and mval<bg_mean + SIGMA_ABOVE_BG*bg_sd): continue

        in_aoi=False
        if aoi_mask_ip is not None:
            inside=allpix=0
            for gx,gy in roi_pixel_iter(aoi_mask_ip, roi):
                allpix+=1
                if aoi_mask_ip.get(gx,gy)>0: inside+=1
            in_aoi=(allpix>0 and float(inside)/float(allpix) >= AOI_ROI_MIN_FRAC)

        circ = circ_from_roi(roi)
        kept.append((roi, area_px, in_aoi, circ))

    rm.reset()
    for r,_,_,_ in kept:
        try: rm.addRoi(r)
        except: pass

    if not kept:
        IJ.log("No nuclei after filtering in %s" % series)
        close_if_open(work); close_if_open(dapi_disp)
        for w in series_wins: close_if_open(w); close_if_open(imp); IJ.run("Collect Garbage",""); return

    # --- CSVs ---
    out_dir = File(folder + File.separator + "Nuclei_Only_Overlays")
    if not out_dir.exists(): out_dir.mkdirs()
    global csv_perroi, csv_summary, image_name
    PERROI_HDR = "Image;Series;ROI_Index;Area_px;Circ;In_AOI"
    SUM_HDR    = "Image;Series;N_ROIs;N_In_AOI;Mean_Area_px;Mean_Circ;Mean_Circ_In_AOI"

    n_all=len(kept)
    n_in = sum(1 for k in kept if k[2])
    mean_area = (sum(k[1] for k in kept)/float(n_all)) if n_all>0 else 0.0

    circ_vals = []; circ_vals_in = []
    for k in kept:
        c = k[3]
        if not Double.isNaN(c):
            circ_vals.append(c)
            if k[2]: circ_vals_in.append(c)

    mean_circ    = (sum(circ_vals)    / float(len(circ_vals)))    if circ_vals    else 0.0
    mean_circ_in = (sum(circ_vals_in) / float(len(circ_vals_in))) if circ_vals_in else 0.0

    export_csv(csv_summary, "%s;%s;%d;%d;%.6f;%.6f;%.6f" % (
        image_name, series, n_all, n_in, mean_area, mean_circ, mean_circ_in
    ), SUM_HDR)

    for ridx,k in enumerate(kept,1):
        roi, area, ina, circ = k[0], k[1], k[2], k[3]
        circ_str = ("%.6f" % circ) if not Double.isNaN(circ) else ""
        export_csv(csv_perroi, "%s;%s;%d;%d;%s;%s" % (
            image_name, series, ridx, area, circ_str, ("TRUE" if ina else "FALSE")
        ), PERROI_HDR)

    # --- COMBO via RGBStackMerge: AOI rot (R) + DAPI blau (B); ROIs obenauf ---
    W, H = dapi_disp.getWidth(), dapi_disp.getHeight()
    dapi_gray = dapi_disp.duplicate(); ensure_gray8(dapi_gray)
    if ENHANCE: IJ.run(dapi_gray, "Enhance Contrast", "saturated=%.3f" % SAT_DAPI_PCT)

    if aoi_view is not None:
        aoi_gray = aoi_view.duplicate(); ensure_gray8(aoi_gray)
    else:
        aoi_gray = ImagePlus("blankR", ByteProcessor(W, H))
    blankG = ImagePlus("blankG", ByteProcessor(W, H))

    combo = RGBStackMerge.mergeChannels([aoi_gray, blankG, dapi_gray], False)  # [R,G,B]
    combo.setTitle(series + "__COMBO_AOIred_DAPIblue_NUC")

    ov = Overlay()
    for k in kept:
        roi, ina = k[0], k[2]
        if SHOW_ONLY_IN_AOI and not ina: continue
        r = roi.clone(); r.setStrokeWidth(STROKE_W); r.setStrokeColor(Color(0,255,0))
        ov.add(r)
    combo.setOverlay(ov)

    FileSaver(combo).saveAsPng(out_dir.getAbsolutePath() + File.separator +
                               safe_name("%s__COMBO_AOIred_DAPIblue_NUC.png" % series))

    # --- DAPI-only Overlay (nur Nuklei innerhalb AOI) ---
    ov_in = Overlay()
    for k in kept:
        if not k[2]: continue
        r = k[0].clone(); r.setStrokeWidth(STROKE_W); r.setStrokeColor(Color(0,255,0))
        ov_in.add(r)
    dapi_only = dapi_disp.duplicate(); ensure_gray8(dapi_only)
    if ENHANCE: IJ.run(dapi_only, "Enhance Contrast", "saturated=%.3f" % SAT_DAPI_PCT)
    apply_blue_rgb(dapi_only)
    dapi_only.setOverlay(ov_in)
    FileSaver(dapi_only).saveAsPng(out_dir.getAbsolutePath() + File.separator +
                                   safe_name("%s__DAPIonly_NUCinAOI.png" % series))

    # Aufräumen
    for imx in (dapi_only, combo, dapi_gray, aoi_gray, blankG):
        close_if_open(imx)
    try:
        for imx in list_open_images():
            if imx and ("__AOIproc" in (imx.getTitle() or "")): close_if_open(imx)
    except: pass
    close_if_open(work); close_if_open(dapi_disp)
    for w in series_wins: close_if_open(w)
    close_if_open(imp)
    IJ.run("Collect Garbage","")

# -------- Hauptprozess --------
folder = DirectoryChooser("Select folder with images").getDirectory()
if not folder:
    IJ.error("No folder selected"); raise SystemExit

csv_perroi  = folder + File.separator + "nuclei_only_perROI.csv"
csv_summary = folder + File.separator + "nuclei_only_summary.csv"

from java.io import File as JFile
files = sorted([f for f in JFile(folder).listFiles()
               if f.isFile() and f.getName().lower().endswith((".tif",".tiff",".png",".jpg",".jpeg",".lif",".nd2"))],
               key=lambda f: f.getName().lower())

for f in files:
    image_name = f.getName()
    try:
        opts = ImporterOptions()
        opts.setId(f.getAbsolutePath())
        opts.setOpenAllSeries(True)
        opts.setVirtual(True)
        imps = BF.openImagePlus(opts)
    except:
        IJ.log("Open failed: %s" % image_name); continue
    if not imps:
        IJ.log("No series: %s" % image_name); continue

    # optional: Serie/T/Z halbieren
    if SPLIT_BY == "series" and len(imps) > 1:
        total=len(imps); cut=total//2
        sel = range(0, cut) if SPLIT_WHICH=="first" else range(cut, total)
    else:
        sel = range(len(imps))

    for sidx in sel:
        if sidx<0 or sidx>=len(imps):
            IJ.log("Index %d out of range for %s (len=%d)" % (sidx, image_name, len(imps)))
            continue
        imp = imps[sidx]
        series = "%s_Series%d" % (image_name, sidx+1) if len(imps)>1 else image_name

        if SPLIT_BY in ("time","z"):
            nT=max(1, imp.getNFrames()); nZ=max(1, imp.getNSlices())
            if SPLIT_BY=="time" and nT>1:
                half=nT//2; rng=("1-%d" % half) if SPLIT_WHICH=="first" else ("%d-%d" % (half+1, nT))
                ch=imp.getNChannels(); zall=imp.getNSlices()
                IJ.run(imp, "Duplicate...", "title=%s__part duplicate channels=1-%d slices=1-%d frames=%s" % (series, ch, zall, rng))
                part=WindowManager.getCurrentImage(); process_series(part, series+"__T_"+SPLIT_WHICH, folder); continue
            elif SPLIT_BY=="z" and nZ>1:
                half=nZ//2; rng=("1-%d" % half) if SPLIT_WHICH=="first" else ("%d-%d" % (half+1, nZ))
                ch=imp.getNChannels(); tall=imp.getNFrames()
                IJ.run(imp, "Duplicate...", "title=%s__part duplicate channels=1-%d slices=%s frames=1-%d" % (series, ch, rng, tall))
                part=WindowManager.getCurrentImage(); process_series(part, series+"__Z_"+SPLIT_WHICH, folder); continue

        process_series(imp, series, folder)

IJ.log("Done. CSVs:")
IJ.log(" - " + csv_perroi)
IJ.log(" - " + csv_summary)
IJ.log("Overlays -> " + (folder + File.separator + "Nuclei_Only_Overlays"))
