import numpy as np
from Bio import PDB
import string
import warnings
import copy 
warnings.filterwarnings("ignore")

# --- SETTINGS ---
INPUT_FILE = "../data/pdb3j6s.ent"
OUTPUT_FILE = "../output/denv3_final.pdb"

def get_matrices_from_header(filename):
    """
    Manually reads the PDB header to find the 60 matrices.
    Returns a list of 60 (Rotation, Translation) tuples.
    """
    matrices = {} # Store as dictionary first {matrix_id: [row1, row2, row3]}
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith("REMARK 350   BIOMT"):
                # Format: REMARK 350   BIOMT1   1  1.000000  0.000000 ...
                parts = line.split()
                # row_num is 1, 2, or 3 (from BIOMT1, BIOMT2...)
                row_label = parts[2]  # "BIOMT1"
                row_num = int(row_label[-1]) 
                matrix_id = int(parts[3]) # The ID of the transformation (1 to 60)
                
                # The floats are the last 4 numbers
                # r1, r2, r3, t
                floats = [float(x) for x in parts[4:]]
                
                if matrix_id not in matrices:
                    matrices[matrix_id] = [[], [], []]
                
                matrices[matrix_id][row_num - 1] = floats

    # Convert to standard format for Biopython
    final_list = []
    sorted_ids = sorted(matrices.keys()) # Ensure we do 1, 2, ... 60 in order
    print(f"   -> Found {len(sorted_ids)} matrices in header.")
    
    for mid in sorted_ids:
        mat_data = matrices[mid]
        # Biopython expects a 3x3 Rotation matrix and a size-3 Translation array
        rot_matrix = np.array([
            mat_data[0][:3],
            mat_data[1][:3],
            mat_data[2][:3]
        ])
        trans_vector = np.array([
            mat_data[0][3],
            mat_data[1][3],
            mat_data[2][3]
        ])
        final_list.append((rot_matrix, trans_vector))
        
    return final_list

def run():
    print(f"1. LOADING: {INPUT_FILE}...")
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("DENV3", INPUT_FILE)
    
    # Extract the original ASU (The "Seed")
    original_model = structure[0]
    asu_chains = list(original_model.get_chains())
    print(f"   -> ASU loaded. Found chains: {[c.id for c in asu_chains]}")

    # 2. GENERATING: Manual Matrix Application
    print("2. RECONSTRUCTING: Parsing header and applying 60 matrices...")
    matrices = get_matrices_from_header(INPUT_FILE)
    
    # Create a new structure to hold the full virus
    full_model = PDB.Model.Model(0)
    full_structure = PDB.Structure.Structure("DENV3_FULL")
    full_structure.add(full_model)
    
    chain_counter = 0
    
    # Generator for unique names (A-Z, a-z, 0-9)
    chars = string.ascii_uppercase + string.ascii_lowercase + string.digits
    
    # THE BIG LOOP: 60 Matrices * 3 Chains = 180 Chains
    for mat_idx, (rot, trans) in enumerate(matrices):
        
        for chain in asu_chains:
            # Deep copy the chain so we can move it without breaking the original
            new_chain = copy.deepcopy(chain)
            
            # Apply the math (Affine Transformation)
            for atom in new_chain.get_atoms():
                # atom.transform expects rotation and translation
                atom.transform(rot, trans)
            
           
            new_id_char = chars[chain_counter % len(chars)]
            # If we run out of chars, we append a number (A1, B1...) logic is safer:
            suffix = str(chain_counter // len(chars)) if chain_counter >= len(chars) else ""
            new_chain.id = new_id_char + suffix
            
            # Determine Type (A, C, or E)
            # The order in the file is A, C, E. So:
            # Loop 0 -> Chain A
            # Loop 1 -> Chain C
            # Loop 2 -> Chain E
            # We can detect this by the original chain ID
            
            if "A" in chain.id:
                type_tag = 10.0   # Pentamer
            elif "C" in chain.id:
                type_tag = 50.0   # Face
            else:
                type_tag = 90.0   # Edge

            # Inject Data into Atoms
            for atom in new_chain.get_atoms():
                x, y, z = atom.get_coord()
                dist = np.sqrt(x**2 + y**2 + z**2)
                
                atom.set_bfactor(type_tag) # For Coloring
                atom.set_occupancy(dist)   # For Heatmap
            
            # Add to the new virus
            full_model.add(new_chain)
            chain_counter += 1

    print(f"   -> Successfully generated {len(list(full_model.get_chains()))} chains.")

    # 4. SAVING
    print(f"4. SAVING: {OUTPUT_FILE}...")
    io = PDB.PDBIO()
    io.set_structure(full_structure)
    io.save(OUTPUT_FILE)
    print("DONE. The virus is built.")

if __name__ == "__main__":

    run()
