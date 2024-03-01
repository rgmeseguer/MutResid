import MDAnalysis as mda
import numpy as np
from .checkFiles import check_pdb_file
from collections import OrderedDict

def parse_ter_lines(pdb_file):
    """
    Parse 'TER' lines from a PDB file.

    Parameters:
        pdb_file (str): The path to the PDB file to parse.

    Returns:
        list of str: List of lines starting with 'TER' from the PDB file.
    """
    with open(pdb_file, 'r') as f:
        return [line for line in f if line.startswith('TER')]

def load_amberpdb(pdb_file):
    """
    Load atoms from the provided PDB file into an MDAnalysis universe.
    
    MDAnalysis does not read TER lines, and AMBER requires them to define each molecule. 
    Therefore we need to set them by hand.

    Parameters:
        pdb_file (str): The path to the PDB file to load.

    Returns:
        tuple: A tuple containing:
               - MDAnalysis universe object representing the loaded atoms.
               - List of TER lines from the PDB file.
    """
    ter_lines = parse_ter_lines(pdb_file)
    universe = mda.Universe(pdb_file)
    print(f"Loaded {pdb_file}.")
    return universe, ter_lines

def select_molecules(ter_lines, universe):
    """
    Select molecules based on terminal residue information extracted from TER lines.

    Parameters:
        ter_lines (list of str): List of TER lines from a PDB file.
        universe (MDAnalysis.core.universe.Universe): The MDAnalysis universe object.

    Returns:
        dict: Dictionary containing selected atoms for each molecule, with TER lines as keys.
    """
    molecules = OrderedDict()
    for i in range(len(ter_lines)):
        # Extract the starting and ending residues of the molecule based on
        # the information of the TER lines
        start_resid = 1 if i == 0 else int(ter_lines[i-1].split()[-1]) + 1
        end_resid = int(ter_lines[i].split()[-1])
        selection = universe.select_atoms(f"resid {start_resid}-{end_resid}")
        molecules[ter_lines[i]] = selection

    return molecules

def format_amberSelection(selection):
    """
    Format the atoms in the given selection into AMBER format.

    Parameters:
        selection (mdanalysis.core.groups.AtomGroup): AtomGroup representing the atoms to format.

    Returns:
        list of str: List of strings representing the atoms formatted in AMBER format.
    """
    return [f"ATOM  {atom.id:5} {atom.name:<4} {atom.resname:>3} {atom.resid:>5}    {atom.position[0]:8.3f}{atom.position[1]:8.3f}{atom.position[2]:8.3f}{0:6.2f}{0.0:6.2f}{atom.type:>11}\n" for atom in selection]

def format_amberpdb(molecules):
    """
    Format the atoms in the given PDB into AMBER format.
    We use the information about the TER lines parsed from the original PDB to define the 
    location of the new TER lines that define the end of a molecule.

    Parameters:
        molecules (OrderedDict): Ordered dictionary where keys indicate the TER line at the end 
                                 of each molecule, and values are AtomGroup representing
                                 the atoms in each molecule.

    Returns:
        list of str: List of strings representing the atoms formatted in AMBER format.
    """
    pdb_lines = []
    for ter_line, selection in molecules.items():
        pdb_lines.extend(format_amberSelection(selection))
        pdb_lines.append(ter_line)
    return pdb_lines

def write_amberpdb(filename, pdb_lines):
    """
    Write the provided PDB lines to a file with the given filename.

    Parameters:
        filename (str): The name of the file to write to.
        pdb_lines (list of str): List of PDB lines to write to the file.
    """
    with open(filename, 'w') as s:
        s.writelines(pdb_lines)

def remove_sidechain(u, residue_index):
    """
    Remove sidechain atoms from a specific residue in the MDAnalysis universe.

    tleap can add missing atoms from a residue but will return an error if unrelated
    atoms are present so we leave only the sidechain.

    Parameters:
        u (MDAnalysis.core.universe.Universe): The MDAnalysis universe object.
        residue_index (int): The index of the residue whose sidechain atoms are to be removed.

    Returns:
        MDAnalysis.core.universe.Universe: The MDAnalysis universe object with sidechain atoms removed.
    """
    residue = u.select_atoms(f"resid {residue_index}")
    backbone_atoms_names = set(['CA','C','O','N','H','HA'])
    sidechain_atoms_indices = [atom.index for atom in residue.atoms if atom.name not in backbone_atoms_names]
    mask = np.ones(len(u.atoms), dtype=bool)
    mask[sidechain_atoms_indices] = False

    # Remove sidechain atoms
    u.atoms = u.atoms[mask]

    print(f"Removed sidechain of residue {residue_index}.")
    return u

def mutate_residue(pdb_file, residue_id, new_residue_name):
    """
    Mutate a specific residue in a PDB file to a new residue type.

    Parameters:
        pdb_file (str): The path to the input PDB file.
        residue_id (int): The ID of the residue to mutate.
        new_residue_name (str): The name of the new residue type.

    Returns:
        output_file (str): The path of the saved file.
    """
    pdb_path, pdb_name, pdb_dir = check_pdb_file(pdb_file)
    universe, terminal_resid = load_amberpdb(pdb_path)
    residue_to_mutate = universe.select_atoms(f"resid {residue_id}")

    if len(residue_to_mutate) == 0:
        print(f"Residue {residue_id} not found in the structure.")
        return
    
    #TODO: Error: residue is not part of the protein
    
    for atom in residue_to_mutate:
        atom.residue.resname = new_residue_name

    universe = remove_sidechain(universe, residue_id)
    # Select molecules based on TER lines and updated universe
    molecules = select_molecules(terminal_resid, universe)
    
    # Write the modified structure to a new PDB file
    output_file = f"{pdb_dir}/mutated_{pdb_name}.pdb"
    write_amberpdb(output_file, format_amberpdb(molecules))
    
    print(f"Structure with residue {residue_id} mutated to {new_residue_name} and sidechain atoms deleted written to {output_file}")

    return output_file

   
