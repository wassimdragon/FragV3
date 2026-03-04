# FragV3: Combinatorial Molecular Fragmentation & Mass Spectrometry Analysis

**FragV3** is a high-performance computational tool designed for advanced mass spectrometry workflows and molecular fragmentation analysis. Built to simulate bond-breaking events, it computes exact fragment masses and generates an interactive, offline-capable 3D environment for real-time peak identification.



## 🚀 Core Features

### 1. High-Performance Combinatorial Fragmentation
* **Graph-Based Bond Breaking:** Parses SMILES strings using `rdkit`, adding explicit hydrogens to construct a complete mathematical graph of atomic nodes and bond edges.
* **Disjoint Set Union (DSU):** Implements a highly optimized DSU algorithm to rapidly resolve resulting connected components (fragments) after simulating simultaneous bond breaks.

* **Parallel Processing:** Splits the immense combinatorial workload of large molecules across all available CPU cores using Python's `multiprocessing.Pool`.
* **Chemical Validation:** Automatically cross-references generated fragments against standard atomic valencies to discard physically impossible chemical states.

### 2. Smart Caching System
* **Cryptographic Hashing:** Generates an MD5 hash of the parent molecule's atomic configuration and your specific bond-break parameters.
* **Incremental Saving:** Fragments are saved incrementally to the `fragment_data/` directory. Future runs with identical parameters instantly load the pre-calculated JSON dataset, bypassing redundant CPU computation.

### 3. Interactive, Zero-Dependency 3D Viewer
* **Energy-Minimization:** Uses RDKit's ETKDGv3 algorithm and MMFF force-field optimization to generate realistic 3D coordinates for the parent molecule.
* **Self-Contained UI:** Compiles all 3D geometries and fragment data into a single, offline-capable HTML file (`3d_viewer/viewer.html`).
* **Real-Time Mass Filtering:** Enter an observed m/z from your mass spectrometer, set a charge state, and define an error tolerance. The embedded JavaScript instantly filters the database and highlights the corresponding atoms directly on the 3D model.

## 💻 Setup & Usage

FragV3 requires Python 3 and the RDKit cheminformatics library. To get started, install the required dependency, run the script, and follow the interactive prompts in your terminal:

1. **Install RDKit:**
   `pip install rdkit`
2. **Run the script:**
   `python FragV3.py`
3. **Molecule SMILES:** Enter your target sequence when prompted (defaults to `c1nc(=O)[nH]cc1I` / 5-Iodouracil).
4. **Max Bond Breaks:** Enter the maximum number of simultaneous bonds to break. Type `max` to force an unconstrained calculation of every possible combination.

Once the computation finishes, FragV3 will automatically launch `3d_viewer/viewer.html` in your default web browser so you can interactively search for your mass spectrometry peaks.

## 📝 Author
**Ouassim Hocine Hafiani, Gemini**
