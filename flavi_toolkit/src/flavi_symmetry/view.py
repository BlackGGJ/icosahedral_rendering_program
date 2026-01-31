def make_script():
    script = """
    # LOAD
    load ../output/denv3_final.pdb, denv3
    hide all
    bg_color white
    
    # VIEW 1: The Structure (A=Red, C=Yellow, E=Orange)
    show cartoon, denv3
    select pentamers, b < 20
    select faces, b > 40 and b < 60
    select edges, b > 80
    
    color red, pentamers
    color yellow, faces
    color orange, edges
    
    # Make it look 3D
    set ray_trace_mode, 1
    zoom denv3
    """
    with open("../output/view_structure.pml", "w") as f:
        f.write(script)
    print("Created Pymol script in output folder.")

if __name__ == "__main__":
    make_script()