"""NGL Parser

This module contains functions that parse and structure data into a 
dict for use with the NGL Molecule Viewer component. 

One or multiple input data files in the PDB or .cif.gz format can be 
entered to return a dict to input as the `data` param of the component. 
"""

import glob


# Helper function to load highlights from content string
def get_highlights(string, sep, atom_indicator):
    residues_list = []
    atoms_list = []

    str_, _str = string.split(sep)
    for e in _str.split(","):
        if atom_indicator in e:
            atoms_list.append(e.replace(atom_indicator, ""))
        else:
            residues_list.append(e)

    return (str_, {"atoms": ",".join(atoms_list), "residues": ",".join(residues_list)})


# Helper function to load the data
def get_data(data_file, pdb_id, color, resetView=False):
    chain = "ALL"
    aa_range = "ALL"
    highlight_dic = {"atoms": "", "residues": ""}

    # Check if only one chain should be shown
    if "." in pdb_id:
        pdb_id, chain = pdb_id.split(".")

        highlights_sep = "@"
        atom_indicator = "a"
        # Check if only a specified amino acids range should be shown:
        if ":" in chain:
            chain, aa_range = chain.split(":")

            # Check if atoms should be highlighted
            if highlights_sep in aa_range:
                aa_range, highlight_dic = get_highlights(
                    aa_range, highlights_sep, atom_indicator
                )

        else:
            if highlights_sep in chain:
                chain, highlight_dic = get_highlights(
                    chain, highlights_sep, atom_indicator
                )

    fname = [f for f in glob.glob(data_path + pdb_id + ".*")][0]

    if "gz" in fname:
        ext = fname.split(".")[-2]
        with gzip.open(fname, "r") as f:
            content = f.read().decode("UTF-8")
    else:
        ext = fname.split(".")[-1]
        with open(fname, "r") as f:
            content = f.read()

    return {
        "filename": fname.split("/")[-1],
        "ext": ext,
        "selectedValue": data_file,
        "chain": chain,
        "aaRange": aa_range,
        "chosen": highlight_dic,
        "color": color,
        "config": {"type": "text/plain", "input": content},
        "resetView": resetView,
        "uploaded": False,
    }
