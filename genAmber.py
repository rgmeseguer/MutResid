import pytraj as pt
from MutResid.checkFiles import check_pdb_file, check_param_files, AmberFiles

def create_tleap_script(amberfile,addSolvent,savePDB,saveScript=False):
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
    tleap_pdb = ""
    tleap_save_pdb = ""

    if amberfile.pdb == None:
         raise Exception("If you do not define the tleap script you need to provide the pdb file")
    tleap_pdb = f'mol = loadpdb "{amberfile.pdb}"\n'
    
    if addSolvent:
        tleap_solvent = f"solvatebox mol TIP3PBOX 12\naddions mol Cl- 0\naddions mol Na+ 0\n" 

    if amberfile.sus_mol != None:
        # if extra_name==None:
        #     raise Exception(f"you have to define the name of the sustrate")
        #mol_path, frcmod_path, extra_name= check_param_files(extra_param)
        #tleap_xparam = f'{extra_name} = loadmol2 {param_dir}/{param_name}.mol2 \nsaveoff {extra_name} {param_dir}/{param_name}.lib \nloadamberparams  {param_dir}/{param_name}.frcmod \ncheck {extra_name}\n'
        tleap_xparam = f'{amberfile.sus_name} = loadmol2 {amberfile.sus_mol} \nloadamberparams  {amberfile.sus_frcmod} \ncheck {amberfile.sus_name}\n'
        
    if savePDB:
        tleap_save_pdb = f"savepdb mol {amberfile.pdb_dir}/{amberfile.pdb_name}_amber.pdb"
    
    tleap_footer = f'saveamberparm mol {amberfile.pdb_dir}/{amberfile.pdb_name}.prmtop {amberfile.pdb_dir}/{amberfile.pdb_name}.inpcrd \nquit'
    
    # Combine all parts to generate tleap script
    tleap_script=tleap_header+tleap_xparam+tleap_pdb+tleap_solvent+tleap_save_pdb+tleap_footer
    
    # Write tleap script to file if requested
    if saveScript:
        with open(f"{amberfile.pdb_dir}/tleap.inp","w") as savefile:
                savefile.write(tleap_script)

    return tleap_script

def generate_amber_files(verbose=False,tleap_script=None,amberfile=None,addSolvent=False,savePDB=False,output=False):
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
            tleap_script=create_tleap_script(amberfile,addSolvent,savePDB)
    
    # Load the tleap script using pytraj
    pt.load_leap(tleap_script,verbose=verbose )
    
    inpcrd_file = f"{amberfile.pdb_dir}/{amberfile.pdb_name}.inpcrd"
    top_file    = f"{amberfile.pdb_dir}/{amberfile.pdb_name}.prmtop"
    if savePDB:
        new_amberfile = AmberFiles(topfile=top_file,crdfile=inpcrd_file,pdbfile=f"{amberfile.pdb_dir}/{amberfile.pdb_name}_amber.pdb",sus_file=amberfile.sus_mol)
    else:
        new_amberfile = AmberFiles(topfile=top_file,crdfile=inpcrd_file)

    if output:
        print(f"Saved files: \n   Amber topology:{top_file} \n   Coordinates: {inpcrd_file}")
        if savePDB:
            print(f"   PDB: {amberfile.pdb_dir}/{amberfile.pdb_name}_amber.pdb")

    return new_amberfile