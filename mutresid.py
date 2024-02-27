import MDAnalysis as mda
import numpy as np
import os
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

    return u

def check_pdb_file(pdb_file):
    """
    Check if a PDB file exists and return its absolute path, base name, and directory path.

    Parameters:
    pdb_file (str): Path to the PDB file.

    Returns:
    tuple: A tuple containing the absolute path, base name (without extension), and the directory path of the PDB file.

    Raises:
    Exception: If the PDB file does not exist.
    """
    # Check if the pdb_file exists
    if os.path.exists(pdb_file):
        # If the file exists, get its absolute path, base name (file name without extension), and the directory path
        pdb_path = os.path.abspath(pdb_file)
        pdb_name = os.path.basename(pdb_file).rsplit('.', maxsplit=1)[0]
        pdb_dir = os.path.dirname(pdb_file)
    else:
        # If the file does not exist, raise an exception
        raise Exception(f"No {pdb_file} found in {os.getcwd()}")

    # Return the absolute path, base name, and the directory path
    return pdb_path, pdb_name, pdb_dir


def check_param_files(mol_file):
    """
    Check for the existence of molecular files and their corresponding frcmod files.

    Parameters:
    mol_file (str): Path to the molecular file (.mol2).

    Returns:
    tuple: A tuple containing the absolute path of the molecular file,
           the absolute path of the corresponding frcmod file,
           the base name of the molecular file (without extension),
           and the directory path of the molecular file.

    Raises:
    Exception: If the molecular file does not exist, if it does not have a .mol2 extension,
               or if the corresponding frcmod file does not exist in the same directory.
    """
    if os.path.exists(mol_file):
        mol_path = os.path.abspath(mol_file)
        param_file = os.path.basename(mol_file)
        param_name, extension = os.path.splitext(param_file)
        param_dir = os.path.dirname(mol_file)
        
        if extension != ".mol2":
            raise Exception("The molecular file must have a .mol2 extension.")

        frcmod_file = os.path.join(param_dir, f"{param_name}.frcmod")

        # Check if the frcmod file exists
        if os.path.exists(frcmod_file):
            # Get the absolute path of the frcmod file
            frcmod_path = os.path.abspath(frcmod_file)
        else:
            raise Exception(f"No {frcmod_file} found in {os.getcwd()}, \
                             mol2 and frcmod files have to be on the same path")
    else:
        raise Exception(f"No {mol_file} found in {os.getcwd()}")

    # Return the absolute paths, base name, and directory path
    return mol_path, frcmod_path, param_name, param_dir

def mutate_residue(pdb_file, residue_id, new_residue_name):
    """
    Mutate a specific residue in a PDB file to a new residue type.

    Parameters:
        pdb_file (str): The path to the input PDB file.
        residue_id (int): The ID of the residue to mutate.
        new_residue_name (str): The name of the new residue type.

    Returns:
        None
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