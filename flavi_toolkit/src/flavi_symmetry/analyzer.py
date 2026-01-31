import numpy as np
from Bio import PDB
import warnings
warnings.filterwarnings("ignore")

INPUT_FILE = "../output/denv3_final.pdb"

def analyze_residue(res_num):
    print(f"\n--- DATA FOR RESIDUE {res_num} ---")
    parser = PDB.PDBParser(QUIET=True)
    try:
        s = parser.get_structure("V", INPUT_FILE)
    except FileNotFoundError:
        print("ERROR: Run main_builder.py first!")
        return

    # Find one example
    for chain in s[0]:
        if res_num in chain:
            atom = chain[res_num]['CA']
            dist = atom.get_occupancy() # We stored distance here
            tag = atom.get_bfactor()    # We stored Type here
            
            # Translate Tag back to Name
            if tag == 10.0: name = "Chain A (Pentamer)"
            elif tag == 50.0: name = "Chain C (Face)"
            else: name = "Chain E (Edge)"
            
            print(f"Location: {name}")
            print(f"Distance from Center: {dist:.2f} Angstroms")
            
            if dist > 230: status = "EXPOSED (Surface Loop)"
            elif dist < 140: status = "BURIED (RNA Interface)"
            else: status = "INTERMEDIATE (Transmembrane)"
            print(f"Status: {status}")
            return

if __name__ == "__main__":
    # Standard DENV3 Epitope site
    analyze_residue(293) 
    # Transmembrane Helix
    analyze_residue(400)