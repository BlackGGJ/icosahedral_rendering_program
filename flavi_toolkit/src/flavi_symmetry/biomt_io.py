import re

def list_biomt_ops(pdb_path: str) -> None:
    op_ids = set()

    with open(pdb_path, "r") as fh:
        for line in fh:
            if line.startswith("REMARK 350   BIOMT"):
                # positions follow PDB convention:
                # columns: "REMARK 350   BIOMT1 n  r11 r12 r13 t1"
                row_id = line[23]        # '1', '2', or '3'
                op_id_str = line[24:28]  # e.g. "   1"
                op_id = int(op_id_str)
                op_ids.add(op_id)

    print(f"Found {len(op_ids)} BIOMT operators: {sorted(op_ids)[:5]} ...")
