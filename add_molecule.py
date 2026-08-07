import sys
import os
import time
import json
import hashlib
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import FragV3

def add_molecule(smiles, custom_name, max_breaks=None):
    start_time = time.time()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")
        
    mol = Chem.AddHs(mol)
    atoms = {a.GetIdx() + 1: a.GetSymbol() for a in mol.GetAtoms()}
    bonds = [(b.GetBeginAtomIdx() + 1, b.GetEndAtomIdx() + 1) for b in mol.GetBonds()]
    formula_str, _ = FragV3.format_fragment(list(atoms.keys()), atoms)

    if max_breaks is None or max_breaks in ("auto", "opt", "optimal"):
        max_breaks = FragV3.calculate_optimal_max_breaks(len(bonds))
    elif max_breaks == "max":
        max_breaks = len(bonds)
    elif isinstance(max_breaks, str) and max_breaks.isdigit():
        max_breaks = int(max_breaks)

    is_max = (max_breaks >= len(bonds))
    status_label = "[THIS IS THE MAXIMUM]" if is_max else f"[THIS IS NOT THE MAX, optimal cap for < 2 min run]"
    print(f"Processing '{custom_name}' ({formula_str})...")
    print(f"  Atoms: {len(atoms)}, Bonds: {len(bonds)}, Max Breaks: {max_breaks} {status_label}")

    fragment_file = FragV3.get_fragment_cache_path(formula_str, smiles, max_breaks)
    os.makedirs(os.path.dirname(fragment_file), exist_ok=True)
    current_mol_hash = FragV3.generate_molecule_hash(atoms, bonds, max_breaks)

    # Single-process fragment generation for speed & stability across OS environments
    FragV3.init_worker(atoms, bonds)
    bond_combinations_generator = FragV3.generate_bond_breaks(bonds, max_breaks)

    seen_indices = set()
    unique_count = 0
    generated_count = 0

    MAX_PART_SIZE = 80 * 1024 * 1024
    parts = []
    current_part = 1
    current_save_path = fragment_file.replace(".json", f"_part{current_part}.json")
    parts.append(current_save_path)

    f = open(current_save_path, 'w')
    f.write('{\n  "metadata": {"molecule_hash": "' + current_mol_hash + '", "part": ' + str(current_part) + '},\n  "fragments": [\n')
    first_fragment_in_part = True

    for broken_bonds_tuple in bond_combinations_generator:
        sublist = FragV3.process_bond_break(broken_bonds_tuple)
        if sublist:
            for fragment_indices in sublist:
                generated_count += 1
                indices_tuple = tuple(fragment_indices)
                if indices_tuple not in seen_indices:
                    seen_indices.add(indices_tuple)
                    unique_count += 1

                    if f.tell() > MAX_PART_SIZE:
                        f.write('\n  ]\n}')
                        f.close()
                        current_part += 1
                        current_save_path = fragment_file.replace(".json", f"_part{current_part}.json")
                        parts.append(current_save_path)
                        f = open(current_save_path, 'w')
                        f.write('{\n  "metadata": {"molecule_hash": "' + current_mol_hash + '", "part": ' + str(current_part) + '},\n  "fragments": [\n')
                        first_fragment_in_part = True

                    weight = FragV3.calculate_molecular_weight(fragment_indices, atoms)
                    formula, _ = FragV3.format_fragment(fragment_indices, atoms)
                    fragment_data = [fragment_indices, weight, formula]

                    if not first_fragment_in_part:
                        f.write(',')
                    json.dump(fragment_data, f, separators=(',', ':'))
                    first_fragment_in_part = False

    f.write('\n  ]\n}')
    f.close()

    print(f"Generated {unique_count} unique fragments ({generated_count} total evaluated) in {time.time() - start_time:.2f}s.")
    fragment_files = [os.path.relpath(p, FragV3.CACHE_DIR) for p in parts]

    # Update central manifest.json
    FragV3.update_manifest_and_launch_viewer(atoms, bonds, formula_str, smiles, max_breaks, fragment_files, custom_name)
    print("FINISHED_SUCCESSFULLY")

if __name__ == '__main__':
    smiles_arg = sys.argv[1] if len(sys.argv) > 1 else "C1[C@@H]([C@H](O[C@H]1N2C=CC(=O)NC2=O)C(O)I)O"
    name_arg = sys.argv[2] if len(sys.argv) > 2 else "5'-Iododeoxyuridine"
    raw_breaks = sys.argv[3] if len(sys.argv) > 3 else "auto"
    breaks_arg = int(raw_breaks) if raw_breaks.isdigit() else raw_breaks
    add_molecule(smiles_arg, name_arg, breaks_arg)
