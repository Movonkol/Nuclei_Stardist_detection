# StarDist Nuclei Detection - VS Code Python Version

Dieses Projekt enthält ein Python-Skript zur Nukleusdetektion mit StarDist, konvertiert von der ursprünglichen ImageJ/Jython-Version.

## Voraussetzungen

- Python 3.8 oder höher
- VS Code mit Python-Extension
- (Optional) CUDA-fähige GPU für schnellere Verarbeitung

## Installation

### 1. Virtuelle Umgebung erstellen

```bash
# Im Projektverzeichnis
python3 -m venv .venv

# Virtuelle Umgebung aktivieren
# Linux/Mac:
source .venv/bin/activate
# Windows:
.venv\Scripts\activate
```

### 2. Abhängigkeiten installieren

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

**Hinweis:** Die Installation von TensorFlow kann einige Minuten dauern.

### 3. StarDist-Modell herunterladen

Beim ersten Ausführen wird das StarDist-Modell automatisch heruntergeladen (~100 MB).

## Verwendung

### Kommandozeile

```bash
# Grundlegende Verwendung (verarbeitet Bilder im aktuellen Verzeichnis)
python nuclei_stardist_vscode.py .

# Mit benutzerdefiniertem Verzeichnis
python nuclei_stardist_vscode.py /pfad/zu/bildern

# Mit erweiterten Optionen
python nuclei_stardist_vscode.py /pfad/zu/bildern \
    --prob-thresh 0.5 \
    --nms-thresh 0.5 \
    --size-min 50 \
    --size-max 2500 \
    --sigma-above-bg 1.5 \
    --scale-factor 1.0
```

### Verfügbare Parameter

| Parameter | Beschreibung | Standard |
|-----------|-------------|----------|
| `input_dir` | Verzeichnis mit Input-Bildern | (erforderlich) |
| `--model` | StarDist-Modell | `2D_versatile_fluo` |
| `--prob-thresh` | Wahrscheinlichkeitsschwelle | 0.5 |
| `--nms-thresh` | NMS-Schwelle | 0.5 |
| `--size-min` | Min. Nukleusbereich (px²) | 50 |
| `--size-max` | Max. Nukleusbereich (px²) | 2500 |
| `--sigma-above-bg` | Min. σ über Hintergrund | 1.5 |
| `--min-intensity` | Min. mittlere Intensität | 0 (aus) |
| `--scale-factor` | Vorskalierung (0.5=klein, 2.0=groß) | 1.0 |
| `--no-qc` | QC-Overlays deaktivieren | false |

### Hilfe anzeigen

```bash
python nuclei_stardist_vscode.py --help
```

## VS Code Debugging

Das Projekt enthält vorkonfigurierte Launch-Konfigurationen:

1. **Python: Nuclei Detection** - Verarbeitet Bilder im aktuellen Verzeichnis
2. **Python: Nuclei Detection (Custom Dir)** - Fragt nach einem benutzerdefinierten Verzeichnis
3. **Python: Current File** - Führt die aktuell geöffnete Datei aus

### Debugging starten

1. Öffnen Sie `nuclei_stardist_vscode.py` in VS Code
2. Drücken Sie `F5` oder gehen Sie zu "Run and Debug"
3. Wählen Sie eine der vorkonfigurierten Konfigurationen
4. Setzen Sie Breakpoints nach Bedarf

## Ausgabe

Das Skript erstellt folgende Dateien im Input-Verzeichnis:

### CSV-Dateien

- **nuclei_counts.csv** - Anzahl der Nuklei pro Bild
  ```
  Image,Total_Nuclei_Count
  image1.tif,235
  image2.tif,189
  TOTAL_ALL_IMAGES,424
  ```

- **nuclei_morphology.csv** - Morphologische Metriken
  ```
  Image,Mean_Area,Mean_Roundness
  image1.tif,125.50,0.8234
  image2.tif,118.75,0.7891
  ```

### QC-Overlays

Wenn aktiviert (Standard), werden PNG-Bilder im Verzeichnis `QC_Overlays/` erstellt:
- DAPI-Kanal in Blau
- Erkannte Nuklei mit roten Konturen
- Hilfreich zur visuellen Validierung der Ergebnisse

## Unterschiede zur ImageJ-Version

### Hauptunterschiede

1. **Keine GUI-Interaktion** - Alle Parameter werden über Kommandozeilenargumente übergeben
2. **Batch-Verarbeitung** - Verarbeitet alle Bilder in einem Verzeichnis automatisch
3. **Bessere Fehlerbehandlung** - Fortsetzung bei Fehlern mit einzelnen Bildern
4. **Flexiblere Eingabe** - Unterstützt TIF, TIFF, PNG, JPG/JPEG

### Funktional gleichwertig

- StarDist-Nukleusdetektierung
- Größen- und Intensitätsfilterung
- Hintergrundberechnung
- QC-Overlay-Generierung
- CSV-Export

## Fehlerbehebung

### TensorFlow-Warnungen

Warnungen wie `Your CPU supports instructions that this TensorFlow binary was not compiled to use` können ignoriert werden. Sie bedeuten nur, dass TensorFlow schneller sein könnte, wenn es für Ihre CPU neu kompiliert würde.

### GPU-Unterstützung

Für GPU-Beschleunigung:

```bash
# Installieren Sie tensorflow-gpu anstelle von tensorflow
pip uninstall tensorflow
pip install tensorflow-gpu>=2.6.0,<2.16.0
```

**Hinweis:** Erfordert kompatible NVIDIA GPU und CUDA-Installation.

### Out of Memory (OOM)

Bei großen Bildern:
- Verwenden Sie `--scale-factor 0.5` für Downsampling
- Verarbeiten Sie Bilder einzeln
- Erhöhen Sie den verfügbaren RAM

### Kein DAPI-Kanal gefunden

Das Skript sucht automatisch nach DAPI-Kanälen. Wenn Ihr Bild anders benannt ist:
- Stellen Sie sicher, dass DAPI-ähnliche Muster im Dateinamen vorhanden sind
- Oder das Skript verwendet automatisch den ersten/blauen Kanal

## Unterstützte Bildformate

- TIFF/TIF (empfohlen)
- PNG
- JPEG/JPG
- Multi-Channel TIFF
- RGB-Bilder (verwendet Blau-Kanal)

## Leistung

Verarbeitungszeit hängt ab von:
- Bildgröße und Anzahl
- CPU/GPU-Leistung
- StarDist-Modell-Konfiguration

Typische Zeiten:
- 2048x2048 px Bild: ~5-30 Sekunden (CPU)
- 2048x2048 px Bild: ~1-5 Sekunden (GPU)

## Kontakt & Support

Bei Fragen oder Problemen:
1. Überprüfen Sie die Dokumentation
2. Stellen Sie sicher, dass alle Abhängigkeiten installiert sind
3. Prüfen Sie die Log-Ausgabe auf Fehlermeldungen

## Lizenz

Dieses Projekt verwendet:
- StarDist (3-Clause BSD License)
- TensorFlow (Apache 2.0 License)
- Scikit-image (BSD License)
