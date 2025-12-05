from flavi_symmetry.biomt_io import list_biomt_ops
from pathlib import Path

if __name__ == "__main__":
    pdb_path = Path(__file__).resolve().parents[1] / "data" / "3j6s.pdb"
    list_biomt_ops(str(pdb_path))
