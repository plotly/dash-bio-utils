from typing import Literal, TypedDict

def create_mol3d_style(
    atoms: list[dict],
    visualization_type: Literal['stick', 'cartoon', 'sphere'] = 'stick',
    color_element: Literal['atom', 'residue', 'residue_type', 'chain'] = 'atom',
    color_scheme: dict[str, str] = None
):
    """Function to create styles input for Molecule3dViewer

    @param atoms
    A list of atoms. Each atom should be a dict with keys: 'name', 'residue_name', 'chain'
    @param visualization_type
    A type of molecule visualization graphic: 'stick' | 'cartoon' | 'sphere'.
    @param color_element
    Elements to apply color scheme to: 'atom' | 'residue' | 'residue_type' | 'chain'.
    @param color_scheme
    Color scheme used to style moleule elements.
    This should be a dict with keys being names of atoms, residues, residue types or chains,
    depending on the value of color_element argument. If no value is provided, default color schemes will be used.
    """
    
    if not visualization_type in ['stick', 'cartoon', 'sphere']:
        raise Exception("Invalid argument type: visualization_type. Should be: 'stick' | 'cartoon' | 'sphere'.")

    if not color_element in ['atom', 'residue', 'residue_type', 'chain']:                    
        raise Exception("Invalid argument type: color_element. Should be: 'atom' | 'residue' | 'residue_type' | 'chain'.")
    
    if not isinstance(atoms, list):
        raise Exception("Invalid argument type: atoms. Should be a list of dict.")
        
    if color_scheme and not isinstance(color_scheme, dict):
        raise Exception("Invalid argument type: color_scheme. Should be a dict.")
    
    default_color = '#ABABAB'
    
    if color_scheme is None:
        color_scheme = {
            'atom': ATOM_COLORS,
            'residue': RESIDUE_COLORS,
            'residue_type': RESIDUE_TYPE_COLORS,
            'chain': CHAIN_COLORS
        }[color_element]
        
    if color_element == 'residue_type':
        residue_type_colors_map = {}
        for aa_class_name, aa_class_members in AMINO_ACID_CLASSES.items():
            for aa in aa_class_members:
                residue_type_colors_map[aa] = color_scheme.get(aa_class_name, default_color)
        color_scheme = residue_type_colors_map
    
    atom_styles = []
    
    for a in atoms:
        if color_element == 'atom':
            atom_color = color_scheme.get(a['name'], default_color)
        if color_element in ['residue', 'residue_type']:
            atom_color = color_scheme.get(a['residue_name'], default_color)
        if color_element == 'chain':
            atom_color = color_scheme.get(a['chain'], default_color)
        
        atom_styles.append({
            'visualization_type': visualization_type,
            'color': atom_color
        })

    return atom_styles