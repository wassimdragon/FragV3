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

FORCED_MAX_BOND_BREAKS = "auto"  # Defaults to optimal calculation (< 2 min runtime)
CACHE_DIR = "docs/fragment_data" # Folder for all previous results

def calculate_optimal_max_breaks(num_bonds, target_max_combos=1_500_000, absolute_cap=12):
    """
    Calculates the optimal maximum bond breaks for a molecule with `num_bonds`
    to ensure computation finishes in under 2 minutes (typically < 30s).
    """
    if num_bonds <= 0:
        return 1

    total_combos = 0
    best_breaks = 1
    for k in range(1, num_bonds + 1):
        try:
            combos_at_k = math.comb(num_bonds, k)
        except AttributeError:
            combos_at_k = math.factorial(num_bonds) // (math.factorial(k) * math.factorial(num_bonds - k))
        
        if total_combos + combos_at_k > target_max_combos and k > 1:
            break
        total_combos += combos_at_k
        best_breaks = k
        if k >= absolute_cap:
            break

    return best_breaks


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
import glob


# --- Helper Functions ---

def parse_molecule_graphs(mol):
    """
    Given an RDKit Mol with explicit hydrogens (Chem.AddHs):
    Extracts all atoms/bonds, heavy atoms/bonds, and an attached hydrogens mapping.
    
    Returns:
        all_atoms: dict {idx: symbol} (1-based)
        heavy_atoms: dict {idx: symbol} (1-based)
        all_bonds: list of (idx1, idx2) (1-based)
        heavy_bonds: list of (idx1, idx2) (1-based)
        attached_h: dict {heavy_idx: [h_idx1, h_idx2, ...]}
    """
    all_atoms = {a.GetIdx() + 1: a.GetSymbol() for a in mol.GetAtoms()}
    all_bonds = [(b.GetBeginAtomIdx() + 1, b.GetEndAtomIdx() + 1) for b in mol.GetBonds()]

    heavy_atoms = {idx: sym for idx, sym in all_atoms.items() if sym != 'H'}
    
    attached_h = defaultdict(list)
    heavy_bonds = []

    for a1, a2 in all_bonds:
        sym1 = all_atoms.get(a1)
        sym2 = all_atoms.get(a2)
        if sym1 == 'H' and sym2 != 'H':
            attached_h[a2].append(a1)
        elif sym2 == 'H' and sym1 != 'H':
            attached_h[a1].append(a2)
        elif sym1 != 'H' and sym2 != 'H':
            heavy_bonds.append((a1, a2))

    for k in attached_h:
        attached_h[k].sort()

    return all_atoms, heavy_atoms, all_bonds, heavy_bonds, dict(attached_h)


def get_fragmentation_graph(all_atoms, heavy_atoms, all_bonds, heavy_bonds, attached_h):
    """
    Determines which graph to fragment on.
    If heavy bonds exist, fragments on the heavy graph and re-attaches hydrogens.
    Otherwise (e.g. methane, water, H2), falls back to fragmenting on all bonds.
    """
    if heavy_bonds:
        return heavy_atoms, heavy_bonds, attached_h
    else:
        return all_atoms, all_bonds, {}


# Global variables for worker processes
GLOBAL_ATOMS = {}
GLOBAL_BONDS = []
GLOBAL_ATOM_INDICES = set()
GLOBAL_ALL_ATOMS = {}
GLOBAL_ATTACHED_H = {}

def init_worker(atoms, bonds, all_atoms=None, attached_h=None):
    """Initializes global state for worker processes."""
    global GLOBAL_ATOMS, GLOBAL_BONDS, GLOBAL_ATOM_INDICES, GLOBAL_ALL_ATOMS, GLOBAL_ATTACHED_H
    GLOBAL_ATOMS = atoms
    GLOBAL_BONDS = bonds
    GLOBAL_ATOM_INDICES = set(atoms.keys())
    GLOBAL_ALL_ATOMS = all_atoms if all_atoms is not None else atoms
    GLOBAL_ATTACHED_H = attached_h if attached_h is not None else {}

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
    """Worker function: finds and validates pure skeletal fragments for a set of broken bonds."""
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

def get_lone_hydrogen_fragments(attached_h, all_atoms=None):
    """Returns each lone hydrogen atom as a 1-atom fragment."""
    h_atoms = set()
    if attached_h:
        for h_list in attached_h.values():
            for h in h_list:
                h_atoms.add(h)
    if all_atoms:
        for idx, sym in all_atoms.items():
            if sym == 'H':
                h_atoms.add(idx)
    return [[h] for h in sorted(h_atoms)]

def expand_fragment_hydrogens(skeletal_fragment, attached_h):
    """
    Expands a pure heavy-atom skeletal fragment into its hydrogen-decorated variants:
    1. Bare skeletal fragment (0 H)
    2. +1 H (skeletal fragment + each 1 attached hydrogen)
    3. For single heavy atoms, all hydrogen subsets (e.g., CH2, CH3, NH2)
    4. For multi-atom fragments, full hydrogen decoration (all attached H) and -1 H variants
    """
    variants = [list(skeletal_fragment)]  # Bare skeletal fragment (0 H)

    if not attached_h:
        return variants

    frag_h = []
    for a in skeletal_fragment:
        frag_h.extend(attached_h.get(a, []))

    if not frag_h:
        return variants

    # +1 H variants (add each attached hydrogen individually)
    for h in frag_h:
        variants.append(sorted(list(skeletal_fragment) + [h]))

    if len(skeletal_fragment) == 1:
        # Single heavy atom: include all combinations of attached hydrogens (e.g. CH2, CH3)
        for k in range(2, len(frag_h) + 1):
            for h_comb in combinations(frag_h, k):
                variants.append(sorted(list(skeletal_fragment) + list(h_comb)))
    else:
        # Multi-atom fragment: full attached hydrogens
        variants.append(sorted(list(skeletal_fragment) + frag_h))
        # Loss of 1 H from full decoration
        if len(frag_h) > 2:
            for h in frag_h:
                rem = [x for x in frag_h if x != h]
                variants.append(sorted(list(skeletal_fragment) + rem))

    return variants


def generate_unique_fragments_incrementally(atoms_dict, bonds_list, max_breaks, save_path, all_atoms_dict=None, attached_h=None):
    """
    Generates unique valid fragments incrementally using multiprocessing and saves them.
    Splits into multiple files if size exceeds 80MB.
    Returns a list of generated file paths.
    """
    if all_atoms_dict is None:
        all_atoms_dict = atoms_dict
    if attached_h is None:
        attached_h = {}

    MAX_PART_SIZE = 80 * 1024 * 1024 # 80MB
    parts = []
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
        return []

    print(f"Starting parallel processing with {cpu_count()} workers...")
    sys.stdout.flush()

    seen_indices = set() # Track unique fragments by sorted indices
    generated_count = 0
    unique_count = 0

    processed_combinations = 0
    current_part = 1
    current_save_path = save_path.replace(".json", f"_part{current_part}.json")
    parts.append(current_save_path)
    
    try:
        # Initialize worker processes with the large read-only data
        with Pool(processes=cpu_count(), initializer=init_worker, initargs=(atoms_dict, list(bonds_list), all_atoms_dict, attached_h)) as pool:
            chunksize = max(1, num_combinations // (cpu_count() * 4))
            chunksize = min(chunksize, 1000)
            
            results_iterator = pool.imap_unordered(process_bond_break, bond_combinations_generator, chunksize=chunksize)

            current_hash = generate_molecule_hash(atoms_dict, list(bonds_list), max_breaks)
            f = open(current_save_path, 'w')
            f.write('{\n  "metadata": {"molecule_hash": "' + current_hash + '", "part": ' + str(current_part) + '},\n  "fragments": [\n')
            first_fragment_in_part = True

            # Emit lone hydrogens first
            lone_h_frags = get_lone_hydrogen_fragments(attached_h, all_atoms_dict)
            for h_frag in lone_h_frags:
                indices_tuple = tuple(h_frag)
                if indices_tuple not in seen_indices:
                    seen_indices.add(indices_tuple)
                    unique_count += 1
                    generated_count += 1
                    weight = calculate_molecular_weight(h_frag, all_atoms_dict)
                    formula, _ = format_fragment(h_frag, all_atoms_dict)
                    fragment_data = [h_frag, weight, formula]
                    if not first_fragment_in_part:
                        f.write(',')
                    json.dump(fragment_data, f, separators=(',', ':'))
                    first_fragment_in_part = False

            for sublist in results_iterator:
                processed_combinations += 1
                
                if processed_combinations % 100 == 0 or processed_combinations == num_combinations:
                    elapsed = time.time() - start_time
                    percent = (processed_combinations / num_combinations) * 100
                    bar_length = 30
                    filled_length = int(bar_length * processed_combinations // num_combinations)
                    bar = '█' * filled_length + '-' * (bar_length - filled_length)
                    
                    # Calculate ETA
                    if processed_combinations > 0:
                        rate = processed_combinations / elapsed
                        remaining = num_combinations - processed_combinations
                        eta_seconds = remaining / rate
                        if eta_seconds > 3600:
                            eta_str = f"{int(eta_seconds // 3600)}h {int((eta_seconds % 3600) // 60)}m"
                        elif eta_seconds > 60:
                            eta_str = f"{int(eta_seconds // 60)}m {int(eta_seconds % 60)}s"
                        else:
                            eta_str = f"{int(eta_seconds)}s"
                    else:
                        eta_str = "Calculating..."
                        
                    print(f"\r  Progress: |{bar}| {percent:6.2f}% ({processed_combinations}/{num_combinations}) ETA: {eta_str}    ", end="", flush=True)

                if sublist:
                    for skeletal_frag in sublist:
                        for fragment_indices in expand_fragment_hydrogens(skeletal_frag, attached_h):
                            generated_count += 1
                            indices_tuple = tuple(fragment_indices)

                            if indices_tuple not in seen_indices:
                                seen_indices.add(indices_tuple)
                                unique_count += 1

                                # Check if we need to rotate to a new part
                                if f.tell() > MAX_PART_SIZE:
                                    f.write('\n  ]\n}')
                                    f.close()
                                    
                                    current_part += 1
                                    current_save_path = save_path.replace(".json", f"_part{current_part}.json")
                                    parts.append(current_save_path)
                                    
                                    f = open(current_save_path, 'w')
                                    f.write('{\n  "metadata": {"molecule_hash": "' + current_hash + '", "part": ' + str(current_part) + '},\n  "fragments": [\n')
                                    first_fragment_in_part = True

                                weight = calculate_molecular_weight(fragment_indices, all_atoms_dict)
                                formula, _ = format_fragment(fragment_indices, all_atoms_dict)
                                fragment_data = [fragment_indices, weight, formula]

                                if not first_fragment_in_part:
                                    f.write(',')
                                json.dump(fragment_data, f, separators=(',', ':'))
                                first_fragment_in_part = False

            f.write('\n  ]\n}')
            f.close()

    except Exception as e:
        print(f"\nError during generation: {e}")
        raise
    finally:
        print() # New line after progress bar

    print(f"Finished generating and saving {unique_count} unique fragments (from {generated_count} total generated) to {len(parts)} parts.")
    sys.stdout.flush()
    return parts

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
    combined = atoms_str + "|" + bonds_str + "|breaks=" + str(max_breaks) + "|alg=v4_h_expansion"
    return hashlib.md5(combined.encode('utf-8')).hexdigest()

def get_fragment_cache_path(formula, smiles, max_breaks):
    """Generates an efficient filename and subdirectory path for the fragment cache."""
    # Sanitize formula and max_breaks for filename
    clean_formula = "".join([c if c.isalnum() else "_" for c in formula])
    # Use MD5 of SMILES for a unique, safe identifier
    smiles_hash = hashlib.md5(smiles.encode()).hexdigest()[:10]
    filename = f"{clean_formula}_{smiles_hash}.json"
    
    # Create a subfolder based on the number of breaks
    break_folder_name = f"{max_breaks}_breaks" if str(max_breaks).isdigit() else f"{max_breaks}_breaks"
    if max_breaks == 1:
        break_folder_name = "1_break"
        
    subfolder_path = os.path.join(CACHE_DIR, break_folder_name)
    os.makedirs(subfolder_path, exist_ok=True)
    return os.path.join(subfolder_path, filename)


def update_manifest_and_launch_viewer(atoms_dict, bonds_list, formula_str, smiles_str, max_breaks, fragment_files, custom_name=""):
    """
    Generates 3D coordinates and updates the central manifest.json
    so the SPA viewer can pick up the new molecule.
    """
    print("-" * 30)
    print("Updating UI Registry...")
    sys.stdout.flush()

    # --- Generate 3D coordinates ---
    # Parse directly from SMILES to preserve correct bond orders (C=O, C=C, etc.)
    mol_block = ""
    if smiles_str:
        print("Generating 3D coordinates...")
        sys.stdout.flush()
        try:
            mol_rdkit = Chem.MolFromSmiles(smiles_str)
            if mol_rdkit is None:
                raise ValueError(f"Could not parse SMILES: {smiles_str}")
            Chem.SanitizeMol(mol_rdkit)

            params = AllChem.ETKDGv3()
            params.randomSeed = 42

            # Tier 1: Standard ETKDGv3 with explicit H
            mol_h = Chem.AddHs(mol_rdkit)
            print("  Embedding (ETKDGv3)...", end="", flush=True)
            result = AllChem.EmbedMolecule(mol_h, params)
            if result != -1:
                AllChem.MMFFOptimizeMolecule(mol_h, maxIters=500)
                mol_block = Chem.MolToMolBlock(mol_h)
                print(" DONE")
            else:
                # Tier 2: ETKDGv3 without H (helps for some large/rigid molecules)
                result = AllChem.EmbedMolecule(mol_rdkit, params)
                if result != -1:
                    AllChem.MMFFOptimizeMolecule(mol_rdkit, maxIters=500)
                    mol_block = Chem.MolToMolBlock(mol_rdkit)
                    print(" DONE (no explicit H)")
                else:
                    # Tier 3: Random coords + ETKDGv3 (for strained/cage molecules)
                    params.useRandomCoords = True
                    result = AllChem.EmbedMolecule(mol_rdkit, params)
                    if result != -1:
                        AllChem.UFFOptimizeMolecule(mol_rdkit, maxIters=2000)
                        mol_block = Chem.MolToMolBlock(mol_rdkit)
                        print(" DONE (random coords)")
                    else:
                        # Tier 4: 2D layout + random Z perturbation → UFF opt
                        # (last resort for fullerenes/cages like C60)
                        # Compute2DCoords sets z=0 (flat), so we perturb z randomly
                        # before UFF — this lets UFF find the real 3D cage geometry.
                        import random as _rng
                        AllChem.Compute2DCoords(mol_rdkit)
                        conf = mol_rdkit.GetConformer()
                        for atom_idx in range(mol_rdkit.GetNumAtoms()):
                            pos = conf.GetAtomPosition(atom_idx)
                            conf.SetAtomPosition(atom_idx, (pos.x, pos.y, _rng.uniform(-2.0, 2.0)))
                        AllChem.UFFOptimizeMolecule(mol_rdkit, maxIters=10000)
                        mol_block = Chem.MolToMolBlock(mol_rdkit)
                        print(" DONE (2D→UFF 3D fallback)")

        except Exception as e:
            print(f"\n  Warning: 3D generation failed ({e}). Viewer will show text only.")
        sys.stdout.flush()

    # --- Update manifest.json ---
    manifest_path = os.path.join(CACHE_DIR, "manifest.json")
    manifest = []
    if os.path.exists(manifest_path):
        try:
            with open(manifest_path, 'r') as f:
                manifest = json.load(f)
        except (IOError, json.JSONDecodeError):
            pass

    # Remove existing entry for this specific file or SMILES if it exists to overwrite
    manifest = [m for m in manifest if m.get("smiles") != smiles_str]

    manifest.append({
        "name": custom_name if custom_name else formula_str,
        "formula": formula_str,
        "smiles": smiles_str,
        "max_breaks": max_breaks,
        "total_bonds": len(bonds_list),
        "file": fragment_files, # This is now a list of relative paths
        "mol_block": mol_block
    })

    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)


    print("\n✅ Done! Molecule added to the library.")
    print("   Open http://localhost:8000/docs/index.html to view it.")
    sys.stdout.flush()

DEFAULT_SMILES = "C1=C(C(=O)NC(=O)N1)I"  # 5-Iodouracil (default fallback)

# -------- Main Execution --------
if __name__ == "__main__":  # Ensure multiprocessing works correctly
    start_total = time.time()

    # --- Interactive inputs ---
    print("=" * 50, flush=True)
    print("  Molecular Fragment Finder", flush=True)
    print("=" * 50, flush=True)

    # SMILES input
    sys.stdout.write(f"  Molecule SMILES [{DEFAULT_SMILES}]: ")
    sys.stdout.flush()
    smiles_input = sys.stdin.readline().strip()
    smiles_to_use = smiles_input if smiles_input else DEFAULT_SMILES

    # Custom Name input
    sys.stdout.write("  Custom Name (optional): ")
    sys.stdout.flush()
    name_input = sys.stdin.readline().strip()
    custom_name = name_input if name_input else ""

    # Parse SMILES with RDKit
    _mol = Chem.MolFromSmiles(smiles_to_use)
    if _mol is None:
        print(f"Error: Could not parse SMILES '{smiles_to_use}'. Please check the string.")
        sys.exit(1)
    _mol = Chem.AddHs(_mol)  # add explicit Hs so all atoms are covered in fragmentation

    # Build atoms dict and bonds list from RDKit mol (1-based indexing)
    all_atoms, heavy_atoms, all_bonds, heavy_bonds, attached_h = parse_molecule_graphs(_mol)
    atoms_to_break, bonds_to_break, attached_h_to_use = get_fragmentation_graph(
        all_atoms, heavy_atoms, all_bonds, heavy_bonds, attached_h
    )

    formula_str, _ = format_fragment(list(all_atoms.keys()), all_atoms)
    print(f"  Parsed molecule: {formula_str}  ({len(all_atoms)} atoms [{len(heavy_atoms)} heavy], {len(bonds_to_break)} breakable bonds [{len(all_bonds)} total])")

    # Max Breaks
    optimal_breaks = calculate_optimal_max_breaks(len(bonds_to_break))
    sys.stdout.write(f"  Max Bond Breaks ['max', 'auto', 1, 2, ...] [optimal={optimal_breaks}]: ")
    sys.stdout.flush()
    breaks_input = sys.stdin.readline().strip()
    if breaks_input:
        if breaks_input.lower() == "max":
            FORCED_MAX_BOND_BREAKS = "max"
        elif breaks_input.lower() in ("auto", "opt", "optimal"):
            FORCED_MAX_BOND_BREAKS = "auto"
        else:
            try:
                FORCED_MAX_BOND_BREAKS = int(breaks_input)
            except ValueError:
                print(f"  Invalid input '{breaks_input}'. Using optimal default {optimal_breaks}.")
                FORCED_MAX_BOND_BREAKS = "auto"
    else:
        FORCED_MAX_BOND_BREAKS = "auto"

    print("=" * 50)
    print(f"  Running: {formula_str}")
    print(f"  Mode: Generate all possible fragments (filter in browser)")
    print(f"  Max Breaks Choice: {FORCED_MAX_BOND_BREAKS}")
    print("=" * 50)

    # --- Validate ---
    if not isinstance(all_atoms, dict) or not all_atoms:
        print("Error: 'atoms' must be a non-empty dictionary.")
        sys.exit(1)
    if not isinstance(bonds_to_break, list):
        print("Warning: 'bonds' is not a list (or is empty).")
    if not ATOMIC_MASSES:
        print("Error: ATOMIC_MASSES dictionary is empty.")
        sys.exit(1)
    if not ATOM_VALENCIES:
        print("Warning: ATOM_VALENCIES dictionary is empty. Valency checks will be skipped.")



    all_fragments_data = [] # List to store fragment data

    # Determine max_breaks early so it's included in the cache hash
    if FORCED_MAX_BOND_BREAKS == "max":
        dynamic_max_breaks = len(bonds_to_break)  # Fully unconstrained: try every possible breakable bond
        print(f"Max bond breaks set to 'max' -> using ALL {dynamic_max_breaks} breakable bonds (THIS IS THE MAXIMUM).")
        if dynamic_max_breaks > 15:
            print(f"  WARNING: {dynamic_max_breaks} bonds is very large. Computation may take a very long time.")
    elif FORCED_MAX_BOND_BREAKS == "auto" or FORCED_MAX_BOND_BREAKS is None:
        dynamic_max_breaks = calculate_optimal_max_breaks(len(bonds_to_break))
        if dynamic_max_breaks < len(bonds_to_break):
            print(f"Max bond breaks set to AUTO OPTIMAL ({dynamic_max_breaks} of {len(bonds_to_break)} bonds) [THIS IS NOT THE MAX, optimized to finish in < 2 min]. Type 'max' if you want all {len(bonds_to_break)} breaks.")
        else:
            print(f"Max bond breaks set to AUTO ({dynamic_max_breaks} of {len(bonds_to_break)} bonds) [THIS IS THE MAXIMUM].")
    elif isinstance(FORCED_MAX_BOND_BREAKS, int):
        dynamic_max_breaks = FORCED_MAX_BOND_BREAKS
        if dynamic_max_breaks < len(bonds_to_break):
            print(f"Max bond breaks set to {dynamic_max_breaks} of {len(bonds_to_break)} bonds [THIS IS NOT THE MAX].")
        else:
            print(f"Max bond breaks set to {dynamic_max_breaks} of {len(bonds_to_break)} bonds [THIS IS THE MAXIMUM].")
    else:
        dynamic_max_breaks = get_max_valency_from_molecule(atoms_to_break, ATOM_VALENCIES, bonds_to_break)

    current_mol_hash = generate_molecule_hash(atoms_to_break, bonds_to_break, dynamic_max_breaks)
    fragment_file = get_fragment_cache_path(formula_str, smiles_to_use, dynamic_max_breaks)

    # --- Load or Generate Fragments ---
    regenerate = False
    fragment_files = []
    
    # Check for existing parts
    base_cache_name = os.path.basename(fragment_file).replace(".json", "")
    existing_parts = sorted(glob.glob(os.path.join(os.path.dirname(fragment_file), f"{base_cache_name}_part*.json")))
    
    if existing_parts:
        print(f"Loading existing results: {len(existing_parts)} parts found.")
        sys.stdout.flush()
        try:
            # Check the hash of the first part
            with open(existing_parts[0], 'r') as f:
                head = "".join([next(f) for x in range(5)])
            
            if '"molecule_hash"' in head and current_mol_hash in head:
                print(f"Cache match! All parts match hash.")
                fragment_files = [os.path.relpath(p, CACHE_DIR) for p in existing_parts]
                sys.stdout.flush()
            else:
                print("Cache hash mismatch. Recomputing...")
                regenerate = True
        except Exception:
            regenerate = True
    else:
        regenerate = True

    if regenerate:
        print("Calculating new fragments...")
        sys.stdout.flush()
        os.makedirs(os.path.dirname(fragment_file), exist_ok=True)
        
        # Remove old parts if recomputing
        for p in existing_parts:
            try: os.remove(p)
            except: pass

        generated_parts = generate_unique_fragments_incrementally(
            atoms_to_break,
            bonds_to_break,
            max_breaks=dynamic_max_breaks,
            save_path=fragment_file,
            all_atoms_dict=all_atoms,
            attached_h=attached_h_to_use
        )
        fragment_files = [os.path.relpath(p, CACHE_DIR) for p in generated_parts]

    # --- Update UI Registry ---
    update_manifest_and_launch_viewer(all_atoms, bonds_to_break, formula_str, smiles_to_use, dynamic_max_breaks, fragment_files, custom_name)

    print(f"Total execution time: {time.time() - start_total:.2f}s")