# =============================================================================
# Molecular Fragment Finder
#
# This script identifies potential molecular fragments resulting from bond breaks
# within a given molecule, calculates their molecular weight, and checks if
# they match a target mass. It uses multiprocessing for parallel bond break
# analysis and includes an option to identify fragments containing "backbone"
# atoms. It also includes functionality to save and load generated fragments
# to a file to avoid repeated computation and manage memory for large molecules
# by generating and saving fragments incrementally.
#
# Author: Ouassim Hocine Hafiani (based on user's provided code)
# Date: 2025-04-22
# =============================================================================

# --- Constants ---

# Approximate atomic masses for common elements.
ATOMIC_MASSES = {
    'H': 1, 'He': 4, 'Li': 7, 'Be': 9, 'B': 11, 'C': 12, 'N': 14, 'O': 16, 'F': 19,
    'Ne': 20, 'Na': 23, 'Mg': 24, 'Al': 27, 'Si': 28, 'P': 31, 'S': 32, 'Cl': 35,
    'Ar': 40, 'Br': 80, 'I': 127

}

# Common valencies for elements. Used for basic valency checks within fragments.
ATOM_VALENCIES = {
    'H': 1, 'He': 0, 'Li': 1, 'Be': 2, 'B': 3, 'C': 4, 'N': 3, 'O': 2, 'F': 1,
    'Ne': 0, 'Na': 1, 'Mg': 2, 'Al': 3, 'Si': 4, 'P': 3, 'S': 2, 'Cl': 1,
    'Ar': 0, 'Br': 1, 'I': 1,
}

# --- User Defined Parameters ---

TARGET_SUM = 0  # Initial weight filter (now handled in browser)
ERROR_TOLERANCE = 0  # Initial tolerance (now handled in browser)

FORCED_MAX_BOND_BREAKS = "max"
CACHE_DIR = "fragment_data" # Folder for all previous results


# --- Required Libraries ---
from itertools import combinations
from collections import defaultdict
import time
from multiprocessing import Pool, cpu_count
import sys
import json
import os
import math
from rdkit import Chem
from rdkit.Chem import AllChem
import hashlib
import webbrowser

# --- Helper Functions ---

# Global variables for worker processes
GLOBAL_ATOMS = {}
GLOBAL_BONDS = []
GLOBAL_ATOM_INDICES = set()

def init_worker(atoms, bonds):
    """Initializes global state for worker processes."""
    global GLOBAL_ATOMS, GLOBAL_BONDS, GLOBAL_ATOM_INDICES
    GLOBAL_ATOMS = atoms
    GLOBAL_BONDS = bonds
    GLOBAL_ATOM_INDICES = set(atoms.keys())

def calculate_molecular_weight(fragment_indices, atoms_dict):
    """Calculates the sum of atomic masses for atoms in a fragment."""
    return sum(ATOMIC_MASSES.get(atoms_dict.get(i), 0) for i in fragment_indices)

def generate_bond_breaks(bonds_list, max_breaks):
    """Generates combinations of bonds to be broken (1 to max_breaks)."""
    num_bonds = len(bonds_list)
    
    # Handle "max" string or None
    if max_breaks == "max":
        effective_max_breaks = num_bonds
    elif max_breaks is None:
        effective_max_breaks = 0
    else:
        effective_max_breaks = max(0, int(max_breaks))

    if not bonds_list:
        return

    for k in range(1, min(effective_max_breaks, num_bonds) + 1):
        for broken_bond_indices in combinations(range(num_bonds), k):
            yield tuple(bonds_list[i] for i in broken_bond_indices)

def is_valid_fragment_valency(fragment_indices, atoms_dict, bonds_within_fragment):
    """
    Checks if the valency of each atom within the fragment is respected
    considering only the bonds *within* that fragment.
    """
    valency_count = defaultdict(int)
    fragment_set = set(fragment_indices)

    for a1, a2 in bonds_within_fragment:
        if a1 in fragment_set and a2 in fragment_set:
            valency_count[a1] += 1
            valency_count[a2] += 1

    for atom_index in fragment_indices:
        atom_symbol = atoms_dict.get(atom_index)
        if atom_symbol:
            max_val = ATOM_VALENCIES.get(atom_symbol, float('inf'))
            if valency_count[atom_index] > max_val:
                return False
    return True

def connected_components(atom_indices, bonds_list):
    """Finds connected components (fragments) using Disjoint Set Union (DSU)."""
    parent = {i: i for i in atom_indices}

    def find(i):
        if parent[i] == i:
            return i
        parent[i] = find(parent[i])
        return parent[i]

    def union(i, j):
        root_i = find(i)
        root_j = find(j)
        if root_i != root_j:
            parent[root_i] = root_j

    for a1, a2 in bonds_list:
        if a1 in parent and a2 in parent:
            union(a1, a2)

    components = defaultdict(set)
    for i in atom_indices:
        components[find(i)].add(i)

    return list(components.values())

def process_bond_break(broken_bonds_tuple):
    """Worker function: finds and validates fragments for a set of broken bonds."""
    # Access global state
    atoms = GLOBAL_ATOMS
    all_bonds = GLOBAL_BONDS
    atom_indices = GLOBAL_ATOM_INDICES

    valid_fragments = []

    broken_bonds_set = set(broken_bonds_tuple)
    remaining_bonds = [b for b in all_bonds if b not in broken_bonds_set]

    potential_fragments = connected_components(atom_indices, remaining_bonds)

    for fragment_set in potential_fragments:
        if not fragment_set:
            continue

        bonds_within_fragment = [
            (a1, a2) for a1, a2 in remaining_bonds
            if a1 in fragment_set and a2 in fragment_set
        ]

        if is_valid_fragment_valency(fragment_set, atoms, bonds_within_fragment):
            valid_fragments.append(list(sorted(fragment_set)))

    return valid_fragments

def generate_unique_fragments_incrementally(atoms_dict, bonds_list, max_breaks, save_path):
    """
    Generates unique valid fragments incrementally using multiprocessing and yields them.
    Saves fragments to file as they are generated.
    """
    start_time = time.time()

    bond_combinations_generator = generate_bond_breaks(bonds_list, max_breaks)
    
    # Calculate number of combinations mathematically to avoid realizing the list
    num_bonds = len(bonds_list)
    
    if max_breaks == "max":
        effective_max_breaks = num_bonds
    elif max_breaks is None:
        effective_max_breaks = 0
    else:
        effective_max_breaks = max(0, int(max_breaks))
        
    effective_k = min(effective_max_breaks, num_bonds)
    num_combinations = sum(math.comb(num_bonds, k) for k in range(1, effective_k + 1))
    
    print(f"Generated {num_combinations} bond break combinations (up to {max_breaks} breaks).")
    sys.stdout.flush()

    if num_combinations == 0:
        print("No bond break combinations generated.")
        return

    print(f"Starting parallel processing with {cpu_count()} workers...")
    sys.stdout.flush()

    seen_indices = set() # Track unique fragments by sorted indices
    generated_count = 0
    unique_count = 0

    processed_combinations = 0
    try:
        # Initialize worker processes with the large read-only data
        with Pool(processes=cpu_count(), initializer=init_worker, initargs=(atoms_dict, list(bonds_list))) as pool:
            # Use a fixed reasonable chunksize or estimate it. 
            # Since we have num_combinations, we can estimate.
            chunksize = max(1, num_combinations // (cpu_count() * 4))
            chunksize = min(chunksize, 1000) # Cap it to keep updates flowing
            
            results_iterator = pool.imap_unordered(process_bond_break, bond_combinations_generator, chunksize=100)

            # Write JSON with metadata including molecule hash (includes max_breaks)
            with open(save_path, 'w') as f:
                current_hash = generate_molecule_hash(atoms_dict, list(bonds_list), max_breaks)
                f.write('{\n  "metadata": {"molecule_hash": "' + current_hash + '"},\n  "fragments": [\n')
                first_fragment = True

                for sublist in results_iterator:
                    processed_combinations += 1
                    
                    # Update progress bar in terminal
                    if processed_combinations % 100 == 0 or processed_combinations == num_combinations:
                        percent = (processed_combinations / num_combinations) * 100
                        bar_length = 30
                        filled_length = int(bar_length * processed_combinations // num_combinations)
                        bar = '█' * filled_length + '-' * (bar_length - filled_length)
                        print(f"\r  Progress: |{bar}| {percent:6.2f}% ({processed_combinations}/{num_combinations})", end="", flush=True)

                    if sublist:
                        for fragment_indices in sublist:
                            generated_count += 1
                            indices_tuple = tuple(fragment_indices)

                            if indices_tuple not in seen_indices:
                                seen_indices.add(indices_tuple)
                                unique_count += 1

                                weight = calculate_molecular_weight(fragment_indices, atoms_dict)
                                formula, _ = format_fragment(fragment_indices, atoms_dict)

                                fragment_data = {
                                    'indices': fragment_indices,
                                    'weight': weight,
                                    'formula': formula
                                }

                                if not first_fragment:
                                    f.write(',\n')
                                json.dump(fragment_data, f, indent=4)
                                first_fragment = False

                                yield fragment_data # Yield the fragment data dictionary

                f.write('\n  ]\n}') # End JSON array and object
                print() # New line after progress bar
                
                print(f"Finished generating and saving {unique_count} unique fragments (from {generated_count} total generated) to {save_path}")
                sys.stdout.flush()

    except Exception as e:
        print(f"\nError during fragment generation: {e}")
        print("Consider reducing max_breaks or simplifying the molecule.")
        # The generator will stop here
        return

def format_fragment(fragment_indices, atoms_dict):
    """Formats a fragment into a chemical formula string and composition dict."""
    composition = defaultdict(int)
    for atom_index in fragment_indices:
        atom_symbol = atoms_dict.get(atom_index)
        if atom_symbol:
            composition[atom_symbol] += 1
        else:
            print(f"Warning: Atom index {atom_index} not found in atoms dictionary.")

    formula = ''.join(
        f"{element}{(count if count > 1 else '')}"
        for element, count in sorted(composition.items())
    )
    return formula, composition

def get_max_valency_from_molecule(atoms_dict, valencies_dict, bonds_list):
    """Determines the highest common valency among atoms in the molecule."""
    if not atoms_dict:
        print("Warning: Input atoms dictionary is empty.")
        return 0

    max_v = 0
    present_atom_symbols = set(atoms_dict.values())

    found_in_valencies = False
    for atom_symbol in present_atom_symbols:
        valency = valencies_dict.get(atom_symbol)
        if valency is not None:
            found_in_valencies = True
            max_v = max(max_v, valency)

    if not found_in_valencies:
        print(f"Warning: None of the atom symbols {list(present_atom_symbols)} found in ATOM_VALENCIES. Max valency is 0.")

    # Ensure at least 1 break is considered if bonds exist but no valency > 0 found
    if max_v == 0 and bonds_list and atoms_dict:
        return 1
    return max_v

def generate_molecule_hash(atoms_dict, bonds_list, max_breaks):
    """Generates a stable hash string for a specific molecule + bond break setting."""
    atoms_str = json.dumps({str(k): v for k, v in sorted(atoms_dict.items())})
    bonds_str = json.dumps(sorted([sorted(b) for b in bonds_list]))
    combined = atoms_str + "|" + bonds_str + "|breaks=" + str(max_breaks)
    return hashlib.md5(combined.encode('utf-8')).hexdigest()

def get_fragment_cache_path(formula, smiles, max_breaks):
    """Generates an efficient filename for the fragment cache."""
    # Sanitize formula and max_breaks for filename
    clean_formula = "".join([c if c.isalnum() else "_" for c in formula])
    # Use MD5 of SMILES for a unique, safe identifier
    smiles_hash = hashlib.md5(smiles.encode()).hexdigest()[:10]
    filename = f"{clean_formula}_breaks_{max_breaks}_{smiles_hash}.json"
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, filename)


def find_and_report_matching_fragments(all_fragments_data, target_sum, error_tolerance, atoms_dict, bonds_list=None):
    """
    Generates a single, fully self-contained HTML file with ALL fragment data
    embedded as JavaScript. The user can interactively filter by mass and tolerance
    directly in the browser — no server required.
    """
    parent_formula, _ = format_fragment(list(atoms_dict.keys()), atoms_dict)
    print("-" * 30)
    print(f"Molecule: {parent_formula} | Building interactive viewer...")
    sys.stdout.flush()

    # --- Generate 3D coordinates from the atoms/bonds using SMILES-aware embedding ---
    mol_block = ""
    if bonds_list:
        print("Generating 3D coordinates...")
        sys.stdout.flush()
        rw = Chem.RWMol()
        idx_map = {}
        for orig_idx in sorted(atoms_dict.keys()):
            symbol = atoms_dict[orig_idx]
            new_idx = rw.AddAtom(Chem.Atom(symbol))
            idx_map[orig_idx] = new_idx
        for a1, a2 in bonds_list:
            if a1 in idx_map and a2 in idx_map:
                i1, i2 = idx_map[a1], idx_map[a2]
                if not rw.GetBondBetweenAtoms(i1, i2):
                    rw.AddBond(i1, i2, Chem.BondType.SINGLE)
        try:
            mol_rdkit = rw.GetMol()
            Chem.SanitizeMol(mol_rdkit)
            mol_rdkit = Chem.AddHs(mol_rdkit)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42

            print("  Embedding molecule...", end="", flush=True)
            result = AllChem.EmbedMolecule(mol_rdkit, params)
            if result == -1:
                mol_rdkit = Chem.RemoveHs(mol_rdkit)
                result = AllChem.EmbedMolecule(mol_rdkit, params)
            print(" DONE")

            if result != -1:
                print("  Optimizing geometry...", end="", flush=True)
                AllChem.MMFFOptimizeMolecule(mol_rdkit, maxIters=500)
                print(" DONE")

                mol_block = Chem.MolToMolBlock(mol_rdkit)
                print("  3D coordinates successfully generated.")
        except Exception as e:
            print(f"\n  Warning: 3D generation failed ({e}). Viewer will show text only.")
        sys.stdout.flush()

    # --- Collect all fragment data (not just matching) for JS embedding ---
    all_frags_for_js = []
    match_count = 0
    total_frags = len(all_fragments_data)
    
    # Simple loop without progress bar as requested
    for frag in all_fragments_data:
        indices = frag.get('indices', [])
        weight = frag.get('weight', 0)
        formula = frag.get('formula', '')
        if indices is None or weight is None:
            continue
            
        rdkit_indices = [idx - 1 for idx in indices]
        all_frags_for_js.append({
            'weight': round(weight, 4),
            'formula': formula,
            'rdkit_indices': rdkit_indices,
        })
        
        if abs(weight - target_sum) <= error_tolerance:
            match_count += 1
            
    # Report matching fragments found
    if match_count > 0:
        print(f"Found {match_count} matching fragment(s):")
        actual_matches = [f for f in all_fragments_data if abs(f.get('weight', 0) - target_sum) <= error_tolerance]
        for m in actual_matches:
            w, f, idx = m.get('weight'), m.get('formula'), m.get('indices')
            print(f"  Weight: {w:<10.4f} Formula: {f:<10} Indices: {idx}")

    # --- Read the local 3Dmol.js file to embed inline ---
    # Try to find 3Dmol-min.js relative to this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    dmol_path = os.path.join(script_dir, "web", "3Dmol-min.js")
    dmol_js = ""
    if os.path.exists(dmol_path):
        with open(dmol_path, 'r') as f:
            dmol_js = f.read()
        print("Bundling 3Dmol.js inline (offline-capable).")
    else:
        print("Warning: web/3Dmol-min.js not found. The viewer needs internet.")

    # Escape molBlock for embedding in JS (backticks safe via replace)
    mol_block_escaped = mol_block.replace('\\', '\\\\').replace('`', '\\`')
    all_frags_json = json.dumps(all_frags_for_js)

    # --- Build fully self-contained HTML ---
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Fragment Viewer – {parent_formula}</title>
<link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;600&display=swap" rel="stylesheet">
<style>
  :root {{
    --primary: #4f46e5;
    --primary-hover: #4338ca;
    --bg-dark: #0f172a;
    --side-bg: #1e293b;
    --text-main: #f8fafc;
    --text-muted: #94a3b8;
    --border: #334155;
    --accent: #0ea5e9;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ 
    font-family: 'Outfit', sans-serif; 
    display: flex; 
    height: 100vh; 
    background: var(--bg-dark); 
    color: var(--text-main); 
    overflow: hidden;
  }}
  
  #sidebar {{ 
    width: 380px; 
    min-width: 320px; 
    display: flex; 
    flex-direction: column; 
    background: var(--side-bg); 
    border-right: 1px solid var(--border); 
    box-shadow: 10px 0 30px rgba(0,0,0,0.3); 
    z-index: 10; 
    transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
  }}
  
  .header {{
    padding: 24px;
    background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%);
    border-bottom: 1px solid var(--border);
  }}
  .header h1 {{ font-size: 20px; font-weight: 600; margin-bottom: 4px; color: var(--accent); }}
  .header .formula {{ font-size: 14px; color: var(--text-muted); letter-spacing: 0.5px; }}

  #controls {{ 
    padding: 24px; 
    background: rgba(255,255,255,0.02);
    border-bottom: 1px solid var(--border); 
  }}
  .input-group {{ margin-bottom: 16px; }}
  .row {{ display: flex; gap: 12px; margin-bottom: 12px; }}
  label {{ font-size: 12px; font-weight: 600; color: var(--text-muted); text-transform: uppercase; letter-spacing: 1px; display: block; margin-bottom: 6px; }}
  input {{ 
    width: 100%; 
    padding: 10px 14px; 
    background: #0f172a; 
    border: 1px solid var(--border); 
    border-radius: 8px; 
    font-size: 14px; 
    color: #fff;
    transition: all 0.2s;
  }}
  input:focus {{ outline: none; border-color: var(--primary); box-shadow: 0 0 0 3px rgba(79, 70, 229, 0.2); }}
  
  #applyBtn {{ 
    width: 100%; 
    padding: 12px; 
    background: var(--primary); 
    color: #fff; 
    border: none; 
    border-radius: 8px; 
    font-size: 14px; 
    font-weight: 600; 
    cursor: pointer; 
    transition: all 0.2s;
    box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06);
  }}
  #applyBtn:hover {{ background: var(--primary-hover); transform: translateY(-1px); box-shadow: 0 10px 15px -3px rgba(0,0,0,0.2); }}
  #applyBtn:active {{ transform: translateY(0); }}

  #stats {{ font-size: 13px; color: var(--accent); margin-top: 12px; font-weight: 500; }}

  #fragList {{ 
    flex: 1; 
    overflow-y: auto; 
    padding: 16px; 
    background: rgba(15, 23, 42, 0.3);
    scrollbar-width: thin;
    scrollbar-color: var(--border) transparent;
  }}
  #fragList::-webkit-scrollbar {{ width: 6px; }}
  #fragList::-webkit-scrollbar-thumb {{ background: var(--border); border-radius: 10px; }}

  .frag-btn {{ 
    display: block; 
    width: 100%; 
    margin-bottom: 10px; 
    padding: 14px 18px; 
    text-align: left; 
    cursor: pointer; 
    border: 1px solid var(--border); 
    background: rgba(255,255,255,0.03); 
    border-radius: 12px; 
    color: var(--text-main);
    font-size: 14px; 
    transition: all 0.2s cubic-bezier(0.4, 0, 0.2, 1);
    position: relative;
    overflow: hidden;
  }}
  .frag-btn:hover {{ 
    background: rgba(255,255,255,0.07); 
    border-color: var(--text-muted);
    transform: translateX(4px);
  }}
  .frag-btn.active {{ 
    background: var(--primary); 
    border-color: var(--primary); 
    box-shadow: 0 10px 15px -3px rgba(79, 70, 229, 0.4);
  }}
  .frag-btn b {{ display: block; font-weight: 600; margin-bottom: 4px; }}
  .frag-btn small {{ font-size: 12px; color: var(--text-muted); display: block; }}
  .frag-btn.active small {{ color: rgba(255,255,255,0.8); }}

  #viewer-wrap {{ flex: 1; position: relative; background: #fff; }}
  #viewer3d {{ position: absolute; inset: 0; width: 100%; height: 100%; }}
  #no-mol {{ display: flex; align-items: center; justify-content: center; height: 100%; background: #0f172a; color: var(--text-muted); font-size: 16px; }}
  
  /* Animations */
  @keyframes fadeIn {{ from {{ opacity: 0; transform: translateY(10px); }} to {{ opacity: 1; transform: translateY(0); }} }}
  .frag-btn {{ animation: fadeIn 0.4s ease forwards; }}
</style>
</head>
<body>
<div id="sidebar">
  <div class="header">
    <h1>Molecular Finder</h1>
    <div class="formula">{parent_formula}</div>
  </div>
  <div id="controls">
    <div class="row">
      <div style="flex:1.5"><label>Target m/z</label><input id="massInput" type="number" step="0.01" value="{target_sum}" placeholder="0.00"></div>
      <div style="flex:1"><label>Charge (z)</label><input id="chargeInput" type="number" step="1" value="1" min="1"></div>
    </div>
    <div class="row">
      <div style="flex:1"><label>± Tolerance</label><input id="tolInput" type="number" step="0.01" value="{error_tolerance}" placeholder="0.00"></div>
    </div>
    <button id="applyBtn" onclick="applyFilter()">Update Results</button>
    <div id="stats">Ready to search...</div>
  </div>
  <div id="fragList">
    <button class="frag-btn active" onclick="showParent(this)"><b>Parent Molecule</b><small>Structure: {parent_formula}</small></button>
  </div>
</div>
<div id="viewer-wrap">
  {'<div id="viewer3d"></div>' if mol_block else '<div id="no-mol">3D visualization unavailable for this molecule</div>'}
</div>

{'<script>' + dmol_js + '</script>' if dmol_js else '<script src="https://3Dmol.org/build/3Dmol-min.js"></script>'}
<script>
const ALL_FRAGS = {all_frags_json};
const MOL_BLOCK = `{mol_block_escaped}`;

let viewer = null;

window.addEventListener('load', function() {{
  if (!MOL_BLOCK) return;
  const el = document.getElementById('viewer3d');
  if (!el) return;
  viewer = $3Dmol.createViewer(el, {{ backgroundColor: 'white' }});
  viewer.resize();
  showParent(document.querySelector('.frag-btn'));
}});

function applyFilter() {{
  const mz = parseFloat(document.getElementById('massInput').value) || 0;
  const z  = parseInt(document.getElementById('chargeInput').value) || 1;
  const tol  = parseFloat(document.getElementById('tolInput').value) || 0;
  
  const targetMass = mz * z;
  const list = document.getElementById('fragList');
  
  const parentBtn = list.querySelector('.frag-btn');
  list.innerHTML = '';
  list.appendChild(parentBtn);

  const matches = ALL_FRAGS.filter(f => Math.abs(f.weight - targetMass) <= tol);
  document.getElementById('stats').textContent = 'Found ' + matches.length + ' matches (m = ' + targetMass.toFixed(2) + ')';
  
  matches.forEach((frag, i) => {{
    const btn = document.createElement('button');
    btn.className = 'frag-btn';
    btn.style.animationDelay = (i * 0.05) + 's';
    btn.innerHTML = '<b>Match ' + (i+1) + '</b><small>Mass: ' + frag.weight.toFixed(4) + ' | ' + frag.formula + '</small>';
    btn.onclick = function() {{ showFragment(frag, this); }};
    list.appendChild(btn);
  }});
}}

function showParent(btn) {{
  document.querySelectorAll('.frag-btn').forEach(b => b.classList.remove('active'));
  btn.classList.add('active');
  if (!viewer || !MOL_BLOCK) return;
  viewer.clear();
  viewer.addModel(MOL_BLOCK, 'mol');
  viewer.setStyle({{}}, {{stick:{{colorscheme:'Jmol'}}, sphere:{{radius:0.3}}}});
  viewer.zoomTo();
  viewer.render();
}}

function showFragment(frag, btn) {{
  document.querySelectorAll('.frag-btn').forEach(b => b.classList.remove('active'));
  btn.classList.add('active');
  if (!viewer || !MOL_BLOCK) return;
  const indices = frag.rdkit_indices || [];
  viewer.clear();
  viewer.addModel(MOL_BLOCK, 'mol');
  viewer.setStyle({{}}, {{stick:{{colorscheme:'Jmol', opacity:0.1}}, sphere:{{radius:0.2, opacity:0.1}}}});
  
  indices.forEach(function(idx) {{
    viewer.setStyle({{index: idx}}, {{stick:{{colorscheme:'Jmol', radius:0.18}}, sphere:{{colorscheme:'Jmol', radius:0.55}}}});
  }});
  viewer.zoomTo();
  viewer.render();
}}

// Apply initial filter on load
window.addEventListener('load', function() {{ applyFilter(); }});
</script>
</body>
</html>"""

    # Write self-contained HTML
    output_dir = "3d_viewer"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "viewer.html")
    with open(output_path, "w") as f:
        f.write(html)

    abs_path = os.path.abspath(output_path)
    print(f"\n---> Self-contained viewer saved: {abs_path}")
    print("     Opening in your browser now...")
    sys.stdout.flush()
    webbrowser.open(f"file://{abs_path}")

    print("-" * 30)
    if match_count == 0:
        print(f"Found 0 fragments matching {target_sum} ± {error_tolerance}.")
    else:
        print(f"Found {match_count} fragment(s) matching {target_sum} ± {error_tolerance}.")

DEFAULT_SMILES  = "C1=C(C(=O)NC(=O)N1)I"  # 5-Iodouracil (default fallback)

# -------- Main Execution --------
if __name__ == "__main__":  # Ensure multiprocessing works correctly
    start_total = time.time()

    # --- Interactive inputs ---
    print("=" * 50)
    print("  Molecular Fragment Finder")
    print("=" * 50)

    # SMILES input
    smiles_input = input(f"  Molecule SMILES [{DEFAULT_SMILES}]: ").strip()
    smiles_to_use = smiles_input if smiles_input else DEFAULT_SMILES

    # Parse SMILES with RDKit
    _mol = Chem.MolFromSmiles(smiles_to_use)
    if _mol is None:
        print(f"Error: Could not parse SMILES '{smiles_to_use}'. Please check the string.")
        sys.exit(1)
    _mol = Chem.AddHs(_mol)  # add explicit Hs so all atoms are covered in fragmentation

    # Build atoms dict and bonds list from RDKit mol (1-based indexing)
    atoms = {}
    for _a in _mol.GetAtoms():
        atoms[_a.GetIdx() + 1] = _a.GetSymbol()
    bonds = [(_b.GetBeginAtomIdx() + 1, _b.GetEndAtomIdx() + 1) for _b in _mol.GetBonds()]

    formula_str, _ = format_fragment(list(atoms.keys()), atoms)
    print(f"  Parsed molecule: {formula_str}  ({len(atoms)} atoms, {len(bonds)} bonds)")

    # Max Breaks
    breaks_input = input(f"  Max Bond Breaks ['max', 1, 2, ...] [{FORCED_MAX_BOND_BREAKS}]: ").strip()
    if breaks_input:
        if breaks_input.lower() == "max":
            FORCED_MAX_BOND_BREAKS = "max"
        else:
            try:
                FORCED_MAX_BOND_BREAKS = int(breaks_input)
            except ValueError:
                print(f"  Invalid input '{breaks_input}'. Using default {FORCED_MAX_BOND_BREAKS}.")

    print("=" * 50)
    print(f"  Running: {formula_str}")
    print(f"  Mode: Generate all possible fragments (filter in browser)")
    print(f"  Max Breaks: {FORCED_MAX_BOND_BREAKS}")
    print("=" * 50)

    # --- Validate ---
    if not isinstance(atoms, dict) or not atoms:
        print("Error: 'atoms' must be a non-empty dictionary.")
        sys.exit(1)
    if not isinstance(bonds, list):
        print("Warning: 'bonds' is not a list (or is empty).")
    if not ATOMIC_MASSES:
        print("Error: ATOMIC_MASSES dictionary is empty.")
        sys.exit(1)
    if not ATOM_VALENCIES:
        print("Warning: ATOM_VALENCIES dictionary is empty. Valency checks will be skipped.")



    all_fragments_data = [] # List to store fragment data

    # Determine max_breaks early so it's included in the cache hash
    if FORCED_MAX_BOND_BREAKS == "max":
        dynamic_max_breaks = len(bonds)  # Fully unconstrained: try every possible bond break count
        print(f"Max bond breaks set to 'max' -> using all {dynamic_max_breaks} bonds in the molecule.")
        if dynamic_max_breaks > 15:
            print(f"  WARNING: {dynamic_max_breaks} bonds is very large. Computation may take a very long time.")
    elif FORCED_MAX_BOND_BREAKS is not None:
        dynamic_max_breaks = FORCED_MAX_BOND_BREAKS
    else:
        dynamic_max_breaks = get_max_valency_from_molecule(atoms, ATOM_VALENCIES, bonds)

    current_mol_hash = generate_molecule_hash(atoms, bonds, dynamic_max_breaks)
    fragment_file = get_fragment_cache_path(formula_str, smiles_to_use, dynamic_max_breaks)

    # --- Load or Generate Fragments ---
    regenerate = False
    if os.path.exists(fragment_file):
        print(f"Loading existing results: {os.path.basename(fragment_file)}")
        sys.stdout.flush()
        try:
            with open(fragment_file, 'r') as f:
                saved_data = json.load(f)
                
            if isinstance(saved_data, dict) and 'metadata' in saved_data:
                saved_hash = saved_data['metadata'].get('molecule_hash')
                if saved_hash == current_mol_hash:
                    all_fragments_data = saved_data.get('fragments', [])
                    print(f"Cache match! Loaded {len(all_fragments_data)} fragments.")
                    sys.stdout.flush()
                else:
                    print("Cache hash mismatch. Recomputing...")
                    regenerate = True
            else:
                print("Old cache format. Recomputing...")
                regenerate = True
                
        except (IOError, json.JSONDecodeError) as e:
            print(f"Error reading cache {fragment_file}: {e}")
            regenerate = True
    else:
        regenerate = True

    if regenerate or not all_fragments_data:
        print("Calculating new fragments...")
        sys.stdout.flush()

        # max_breaks already determined above - just report it
        if FORCED_MAX_BOND_BREAKS == "max":
            print(f"Max bond breaks set to 'max' (= {dynamic_max_breaks} bonds in molecule)")
        elif FORCED_MAX_BOND_BREAKS is not None:
            print(f"Max bond breaks forced to: {dynamic_max_breaks}")
        else:
            print(f"Max bond breaks determined dynamically: {dynamic_max_breaks}")

        print("-" * 30)
        sys.stdout.flush()

        all_fragments_data = [] # Clear data if regenerating
        try:
            # Generate fragments and populate all_fragments_data list directly
            for fragment_data in generate_unique_fragments_incrementally(atoms, bonds, max_breaks=dynamic_max_breaks, save_path=fragment_file):
                 all_fragments_data.append(fragment_data)

        except Exception as e:
             print(f"An unexpected error occurred during fragment generation: {e}")
             sys.stdout.flush()
             sys.exit(1)

    # --- Filter and Report Matching Fragments ---
    find_and_report_matching_fragments(all_fragments_data, TARGET_SUM, ERROR_TOLERANCE, atoms, bonds)

    print(f"Total execution time: {time.time() - start_total:.2f}s")