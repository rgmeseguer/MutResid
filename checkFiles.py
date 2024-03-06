import os

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
        pdb_dir = os.path.dirname(os.path.abspath(pdb_file))
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
        param_dir = os.path.dirname(os.path.abspath(mol_file))
        
        if extension != ".mol2":
            raise Exception("The molecular file must have a .mol2 extension.")

        frcmod_file = os.path.join(param_dir, f"{param_name}.frcmod")

        # Check if the frcmod file exists
        if os.path.exists(frcmod_file):
            # Get the absolute path of the frcmod file
            frcmod_path = os.path.abspath(frcmod_file)
        else:
            raise Exception(f"No {frcmod_file} found in {os.getcwd()}, \
                             mol2 and frcmod files have to be on the same path and have the same name")
    else:
        raise Exception(f"No {mol_file} found in {os.getcwd()}")

    with open(mol_path,"r") as mol:
        resname = mol.readlines()[1].strip()
    # Return the absolute paths, base name, and directory path
    return mol_path, frcmod_path, param_dir, resname

def check_topncrd(top_file,traj_file):
    if os.path.exists(top_file):
        top_path = os.path.abspath(top_file)
        top_dir = os.path.dirname(os.path.abspath(top_file))
    else:
        raise Exception(f"No {top_file} found in {os.getcwd()}")
        
    if os.path.exists(os.path.exists(traj_file)):
        traj_path = os.path.abspath(traj_file)
        traj_dir = os.path.dirname(os.path.abspath(traj_file))
    else:
        raise Exception(f"No {traj_file} found in {os.getcwd()}")
                      
    return top_path, top_dir, traj_path, traj_dir

class AmberFiles():
    def __init__(self,topfile=None,crdfile=None,pdbfile=None,sus_file=None):
        self.parm = None
        self.crd  = None
        self.pdb  = None
        self.sus_mol = None
        
        if topfile!=None or crdfile!=None:
            if crdfile == None: 
                raise Exception("If Topology is defined you have to define the coordinates too.") 
            if topfile == None:
                raise Exception("If amber coordinates are defined you have to define the topology too.") 
            
            self.parm, self.parm_dir, self.crd, self.crd_dir = check_topncrd(topfile,crdfile)
            
        if pdbfile!=None:
            self.pdb, self.pdb_name, self.pdb_dir = check_pdb_file(pdbfile)
        
        if sus_file!=None:
            self.sus_mol, self.sus_frcmod, self.sus_dir, self.sus_name = check_param_files(sus_file)