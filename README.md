# SWIFT: Smith–Waterman with Intelligent Filtering and Tiling

**SWIFT** is an extended and optimized version of the classic [Smith–Waterman algorithm](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) for local sequence alignment.  
It introduces **k-mer–based intelligent filtering**, **tiled alignment**, and **parallel processing** to drastically improve runtime without compromising alignment accuracy.

---

## 🚀 Features
- **Intelligent k-mer Filtering** – Quickly eliminates unlikely match regions before alignment.
- **Tiling Strategy** – Breaks candidate regions into manageable tiles for faster processing.
- **Parallel Processing** – Leverages all available CPU cores for concurrent tile alignment.
- **Merge & Refine** – Reconstructs full alignments from top-scoring tiles using standard SW.
- **Multiple Output Formats** – Export results to JSON and SAM for downstream analysis.
- **Visualization** – ASCII alignment display and `matplotlib` identity plots.
- **Benchmarking** – Built-in performance comparison with traditional SW.

---

## 📂 Project Structure
```
swift-v2.py          # Main implementation of SWIFT algorithm
test_data.txt        # Example input sequences
alignments.json      # Exported alignment results (generated)
alignments.sam       # Exported SAM file (generated)
README.md            # Documentation
```

---

## 📦 Installation

1. **Clone this repository**
```bash
git clone https://github.com/yourusername/swift-alignment.git
cd swift-alignment
```

2. **Create a virtual environment (optional but recommended)**
```bash
python -m venv venv
source venv/bin/activate   # Mac/Linux
venv\Scripts\activate      # Windows
```

3. **Install dependencies**
```bash
pip install -r requirements.txt
```

**requirements.txt**
```
numpy
matplotlib
```

---

## 🧪 Usage

### 1. Prepare Input
Create a `test_data.txt` file in this format:
```
>query
ATATGCAATAGCGCAG...
>target
ATAAGCCACCAGCTCCG...
```

### 2. Run SWIFT
```bash
python swift-v2.py
```

### 3. Example Output
```
Execution time for (swift_alignment_parallel): 3.88 seconds
[INFO] Alignments exported to alignments.json
[INFO] Exported SAM -> alignments.sam

=== Top Refined Alignment (ASCII) ===
ATATG-CAATAGCGCAG-C-GAAC...
||| | ||  ||| | | | |  ...
ATAAGCCACCAGCTCCGACGGT...
```

---

## 📊 Performance Example
| Method            | Runtime (s) | Best Score | Tiles Processed |
|-------------------|-------------|------------|------------------|
| Traditional SW    | 7.48        | 1226       | 1                |
| SWIFT (parallel)  | 3.88        | 128        | 75               |

---

## 📜 Algorithm Overview

1. **K-mer Extraction** – Generate all k-length substrings of the query for fast match lookups.
2. **Candidate Tile Identification** – Locate target regions with high k-mer hits.
3. **Parallel Tile Alignment** – Align each tile independently using SW.
4. **Merge & Refine** – Combine adjacent/overlapping tiles and re-align.
5. **Export & Visualization** – Output alignments and show match identity plots.

---

## 📄 Output Formats

- **JSON** – Stores alignment data for programmatic use.
- **SAM** – Standard format for sequence alignment tools.
- **ASCII Plot** – Terminal-friendly alignment view.
- **Matplotlib Plot** – Visual match/mismatch representation.

---

## 🧬 Example Visualization
```
ATATG-CAATAGCGCAG-C-GAAC-GAAG--AGCGAC-ATACAAATCCG
||| | ||  ||| | | | |  | || |  || | | | ||| ||| |
ATAAGCCACCAGCTCCGACGGTGCTGACGCTAGTGCCAAAACATATCGG
```

---

## ⚙️ Parameters
You can tune the algorithm in `swift-v2.py`:
```python
K = 3                 # K-mer size
TILE_SIZE = 64        # Tile length
HIT_THRESHOLD = 2     # Min. k-mer hits for a tile
WORKERS = None        # CPU cores (None = auto-detect)
```

---

## 📚 Citation
If you use SWIFT in your research, please cite:
```
Samuel PK, 2025.
SWIFT: Smith–Waterman with Intelligent Filtering and Tiling.
GitHub Repository. https://github.com/DevPrinceK/swift
```

---

## 🤝 Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss your ideas.

---

## 📜 License
This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

---
