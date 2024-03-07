import os
import shutil
import pytraj as pt 

from MutResid.mutres import mutate_residue
from MutResid.genAmber import generate_amber_files
from MutResid.checkFiles import AmberFiles


def extract_frame(traj,frame,output):
    """
    Extracts a single frame from a trajectory and writes it to a PDB file.

    Args:
        traj (pytraj.Trajectory): Input trajectory.
        frame (int): Index of the frame to extract.
        output (str): Output PDB file path.
    """    
    pt.write_traj(output,traj,frame_indices=[frame],overwrite=True)

def calc_mutated_ene_diff(Amberfile_Wild,Amberfile_Mut,mutated_resid):
    traj_wild = pt.load(Amberfile_Wild.crd, Amberfile_Wild.parm)["!:WAT,Cl-"]
    traj_mut  = pt.load(Amberfile_Mut.crd,  Amberfile_Mut.parm )["!:WAT,Cl-"]
        
    data_wild = pt.lie(traj_wild,f':{mutated_resid}&(!@CA,C,O,N,H,HA)',options="nopbc")
    data_mut  = pt.lie(traj_mut, f':{mutated_resid}&(!@CA,C,O,N,H,HA)',options="nopbc")        
    
    elec_energy = data_wild['LIE[EELEC]'][0] - data_mut['LIE[EELEC]'][0]
    vdw_energy  = data_wild['LIE[EVDW]'][0]  - data_mut['LIE[EVDW]'][0]

    return elec_energy, vdw_energy

def test_mutation(amberfile,residue_id,new_residue_name,verboseLeap=False,clean=False,skip=1):
    """
    Tests residue mutation in a molecular dynamics trajectory.

    Args:
        trajectory_file (str): Path to the trajectory file.
        top_file (str): Path to the topology file.
        residue_id (str): ID of the residue to mutate.
        new_residue_name (str): Name of the new residue.
        sus_file (str, optional): Path to additional parameters file (if needed). Defaults to None.
        sus_name (str, optional): Name for the additional parameters (if needed). Defaults to None.
        verboseLeap (bool, optional): Whether to enable verbose output during AMBER file generation. Defaults to False.
        clean (bool, optional): Whether to clean up intermediate files. Defaults to False.

    Returns:
        tuple: Tuple containing lists of electrostatic energy (E_elec) and van der Waals energy (E_vw).
    """
    top_path = amberfile.parm
    traj_path = amberfile.crd
    wdir = amberfile.parm_dir
    
    
    traj = pt.load(traj_path,top_path)
    traj = traj["!:WAT"] 
                      
    outDirWild = f"{wdir}/str_Wild"
    outDirMut  = f"{wdir}/str_{new_residue_name}"
    
    if not os.path.exists(f'{outDirWild}'): os.mkdir(f'{outDirWild}')
    if not os.path.exists(f'{outDirMut}'):  os.mkdir(f'{outDirMut}')
    
    # Initialize energy lists
    elec_list = []
    vdw_list = []
    # Process frames from the trajectory
    for idx, frame in enumerate(traj[::skip]):
        outfile_Wild = f"{outDirWild}/frame_{idx}.pdb"
        outfile_Mut = f"{outDirMut}/mutated_frame_{idx}.pdb"
        extract_frame(traj,idx,outfile_Wild)
        mutate_residue(outfile_Wild, output_file=outfile_Mut, residue_id=residue_id, new_residue_name=new_residue_name);
        amberfile_Wild = generate_amber_files( amberfile=AmberFiles(pdbfile=outfile_Wild,sus_file=amberfile.sus_mol), addSolvent=True, verbose=verboseLeap,saveScript=True )
        amberfile_Mut  = generate_amber_files( amberfile=AmberFiles(pdbfile=outfile_Mut,sus_file=amberfile.sus_mol),  addSolvent=True, verbose=verboseLeap,saveScript=True )
        os.remove(outfile_Wild)
        os.remove(outfile_Mut)
        
        elec_energy, vdw_energy = calc_mutated_ene_diff(amberfile_Wild,amberfile_Mut,residue_id)
        elec_list.append( elec_energy )
        vdw_list.append( vdw_energy )
                      
    # Clean up intermediate files if requested
    if clean:
        shutil.rmtree(outDirWild, ignore_errors=True)
        shutil.rmtree(outDirMut, ignore_errors=True)

    
    return elec_list, vdw_list
                      
