# FragV3 — Molecular Fragment Viewer

FragV3 is a tool for mass spectrometry research. You give it a molecule, it calculates every possible fragment produced by breaking chemical bonds, and lets you search through them interactively in a 3D web viewer.

---

## What it does

1. **Takes a molecule** (as a SMILES string) as input
2. **Calculates all possible fragments** from breaking 1, 2, 3... up to N bonds simultaneously
3. **Saves the results** to the `fragment_data/` folder (organized by number of breaks)
4. **Updates a library** so you can view and search all your molecules in one place

In the browser, you can:
- Switch between any molecule you've computed from a dropdown
- Enter an observed **m/z** value, set a **charge state (z)**, and a **tolerance**
- Instantly see which fragments match, highlighted in 3D

---

## Project Structure

```
FragV3/
├── FragV3.py             # Main Python script — run this to add a new molecule
├── web/
│   ├── index.html        # The interactive web viewer
│   └── 3Dmol-min.js      # 3D rendering library (works offline)
└── fragment_data/
    ├── manifest.json     # Registry of all computed molecules
    ├── 12_breaks/        # Fragment data organized by number of bond breaks
    ├── 6_breaks/
    └── ...
```

---

## Setup

You need Python 3 with RDKit installed:

```bash
pip install rdkit
```

---

## How to use

### Step 1 — Add a molecule to your library

```bash
python FragV3.py
```

- **SMILES**: paste your molecule's SMILES string (press Enter to use the default: 5-Iodouracil)
- **Max breaks**: type `max` to explore all combinations, or a number (e.g. `6`) to limit it

Results are cached — running the same molecule again loads instantly from disk with no recalculation.

### Step 2 — View your library

Start a local web server in the FragV3 folder:

```bash
python3 -m http.server 8000
```

Then open: **http://localhost:8000/web/index.html**

### Optional — One-command shortcut (`FragV3Run`)

Add this to your `~/.zshrc` to launch everything with one command:

```bash
FragV3Run() {
  cd "/path/to/FragV3"
  python3 -m http.server 8000 &
  sleep 1
  open http://localhost:8000/web/index.html
  wait
}
```

Then just type `FragV3Run` in any terminal.

---

## Technical Notes

- **JSON format**: Fragment data is stored as compact arrays `[indices, weight, formula]` to keep file sizes small (e.g. Sucrose with 1.26M fragments → 139 MB instead of 600 MB).
- **3D coordinates**: Generated from the SMILES string directly to preserve correct bond orders (double bonds, aromatic bonds). Uses a 4-tier fallback system for complex cage molecules like C60 (fullerenes).
- **Parallel processing**: Uses all available CPU cores for fragment generation.

---

## Requirements

- Python 3.8+
- `rdkit` (`pip install rdkit`)
- A modern web browser (Chrome, Firefox, Safari)

---

## Author
Ouassim Hocine Hafiani
