import pytraj as pt
from checkFiles import check_pdb_file, check_param_files

def create_tleap_script(pdb,extra_param,extra_name,addSolvent,saveScript=False):
    """
    Generate a tleap script for creating Amber files from PDB and additional parameters.

    Parameters:
        pdb (str): The path to the input PDB file.
        extra_param (str): The path to additional parameters file.
        extra_name (str): Name of the additional parameter.
        addSolvent (bool): Flag to indicate if solvent needs to be added.
        saveScript (bool): Flag to save the tleap script.

    Returns:
        tuple: Directory path, PDB name, and tleap script.
    """
    tleap_header = f"source leaprc.protein.ff14SB\nsource leaprc.water.tip3p\n"
    tleap_solvent = ""
    tleap_xparam = ""

    if pdb == None:
         raise Exception("If you do not define the tleap script you need to provide the pdb file")
    pdb_path, pdb_name, pdb_dir = check_pdb_file(pdb)
    tleap_pdb = f'mol = loadpdb "{pdb_path}"\n'
    
    if addSolvent:
        tleap_solvent = f"solvatebox mol TIP3PBOX 12\naddions mol Cl- 0\naddions mol Na+ 0\n" 

    if extra_param!=None:
        if extra_name==None:
            raise Exception(f"you have to define the name of the sustrate")
        mol_path, frcmod_path, param_name, param_dir = check_param_files(extra_param)
        tleap_xparam = f'{extra_name} = loadmol2 "{param_dir}/{param_name}.mol2"\nsaveoff {extra_name} "{param_dir}/{param_name}.lib"\nloadamberparams "{param_dir}/{param_name}.frcmod"\ncheck {extra_name}\n'
    
    tleap_footer = f'saveamberparm mol {pdb_dir}/{pdb_name}.prmtop {pdb_dir}/{pdb_name}.inpcrd \nsavepdb mol "{pdb_dir}/{pdb_name}_amber.pdb"'
    
    # Combine all parts to generate tleap script
    tleap_script=tleap_header+tleap_xparam+tleap_pdb+tleap_solvent+tleap_footer
    
    # Write tleap script to file if requested
    if saveScript:
        with open(f"{pdb_dir}/tleap.inp","w") as savefile:
                savefile.write(tleap_script)

    return pdb_dir,pdb_name,tleap_script

def generate_amber_files(verbose=False,tleap_script=None,pdb=None,extra_param=None,extra_name=None,addSolvent=False):
    """
    Generate Amber files from input PDB and additional parameters.

    Parameters:
        verbose (bool): Flag to enable verbose output.
        tleap_script (str): String of a tleap script. If defined the remaining parameters are ignored.
        
        if telap_script is not defined pdb is required and the rest of variables are optional.
            pdb (str): Path to the input PDB file.
            extra_param (str): Path to additional parameters file.
            extra_name (str): Name of the additional parameter.
            addSolvent (bool): Flag to indicate if solvent needs to be added.
    
    Returns:
        None
    """
    # If tleap script is not provided, create one
    if tleap_script==None:
            pdb_dir,pdb_name,tleap_script=create_tleap_script(pdb,extra_param,extra_name,addSolvent)
    
    # Load the tleap script using pytraj
    pt.load_leap(tleap_script,verbose=verbose )
    
    print(f"Saved files: \n   Amber topology:{pdb_dir}/{pdb_name}.prmtop) \n   Coordinates: {pdb_dir}/{pdb_name}.inpcrd \n   PDB: {pdb_dir}/{pdb_name}_amber.pdb")
