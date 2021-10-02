"""mmCIF parser

This module contains functions that can read mmCIF files and return a
JSON representation of the structural data."""

import re
import json
import os
from shutil import copy2

import parmed as pmd


def create_data(mmcif_path):
    """
    Parse the mmCIF file to generate input modelData

    @param mmcif_path
    Name of the biomolecular structure file in mmCIF format

    """

    # Create local copy of temp file
    copy2(mmcif_path, './tmp.cif')

    # Use parmed to read the bond information from temp file
    top = pmd.load_file('tmp.cif')

    # Remove the created temp file
    os.remove('tmp.cif')

    # Read mmCIF file to create atom/bond information
    with open(mmcif_path, 'r') as infile:
        # store only non-empty lines
        lines = [l.strip() for l in infile if l.strip()]

    # Initialize all variables
    var_nchains = []
    serial = []
    atm_name = []
    res_name = []
    chain = []
    res_id = []
    positions = []
    occupancy = []
    temp_factor = []
    atom_type = []
    ct = 0

    datb = {
        'atoms': [],
        'bonds': []
    }

    # Variables that store the character positions of different
    # parameters from the molecule mmCIF file
    serialpos = 1
    atm_namepos = 3
    r_namepos = 5
    chainpos = 6
    r_idpos = 16
    xpos = 10
    ypos = 11
    zpos = 12
    occupos = 13
    bfacpos = 14
    atm_typepos = 2

    for l in lines:
        line = l.split()
        if "ATOM" in line[0] or "HETATM" in line[0]:
            serial.append(int(line[serialpos]))
            atm_name.append(line[atm_namepos].strip())
            val_r_name = line[r_namepos].strip()
            res_name.append(val_r_name)
            chain_val = line[chainpos].strip()
            chain.append(chain_val)
            if chain_val not in var_nchains:
                var_nchains.append(chain_val)
            val_r_id = int(line[r_idpos])
            res_id.append(val_r_id)
            x = float(line[xpos])
            y = float(line[ypos])
            z = float(line[zpos])
            positions.append([x, y, z])
            occupancy.append(line[occupos].strip())
            temp_factor.append(line[bfacpos].strip())
            atom_type.append(line[atm_typepos].strip())
            ct += 1
    # Create list of atoms
    tmp_res = res_id[0]
    resct = 1
    for i in range(len(chain)):  # pylint: disable=consider-using-enumerate
        if tmp_res != res_id[i]:
            tmp_res = res_id[i]
            resct += 1
        datb['atoms'].append({
            "name": atm_name[i],
            "chain": chain[i],
            "positions": positions[i],
            "residue_index": resct,
            "element": atom_type[i],
            "residue_name": res_name[i] + str(res_id[i]),
            "serial": i,
        })

    # Create list of bonds using the parmed module
    for i in range(len(top.bonds)):
        bondpair = top.bonds[i].__dict__
        atom1 = re.findall(r"\[(\d+)\]", str(bondpair['atom1']))
        atom2 = re.findall(r"\[(\d+)\]", str(bondpair['atom2']))
        datb['bonds'].append({
            'atom2_index': int(atom1[0]),
            'atom1_index': int(atom2[0])
        })
    return json.dumps(datb)
