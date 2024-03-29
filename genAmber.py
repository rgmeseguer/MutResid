import pytraj as pt
from MutResid.checkFiles import AmberFiles

def create_tleap_script(amberfile,addSolvent,savePDB,saveScript):
    """
    Create a tleap script based on input parameters.

    Parameters:
        amberfile (AmberFiles): AmberFiles object containing file paths.
        addSolvent (bool): Whether to solvate the system.
        savePDB (bool): Whether to save the resulting PDB file.
        saveScript (bool): Whether to save the tleap script.

    Returns:
        str: Generated tleap script.
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

def generate_amber_files(verbose=False,tleap_script=None,amberfile=None,addSolvent=False,savePDB=False,output=False,saveScript=False):
    """
    Generate AMBER files using tleap.

    Parameters:
        verbose (bool, optional): Whether to print verbose output during tleap execution.
        tleap_script (str, optional): Tleap script to use. If None, create one.
        amberfile (AmberFiles, optional): AmberFiles object containing file paths.
        addSolvent (bool, optional): Whether to solvate the system.
        savePDB (bool, optional): Whether to save the resulting PDB file.
        output (bool, optional): Whether to print the paths of saved files.
        saveScript (bool, optional): Whether to save the tleap script.

    Returns:
        AmberFiles: An AmberFiles object containing paths to generated AMBER files.
    """
    # If tleap script is not provided, create one
    if tleap_script==None:
            tleap_script=create_tleap_script(amberfile,addSolvent,savePDB,saveScript)
    
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