# MutResid
This Python module harnesses the power of MDAnalysis to facilitate the mutation of individual residues within protein structures in PDB format. Additionally, it seamlessly integrates with tleap to generate the requisite Amber files for further analysis and simulation.

## Dependecies:
- numpy
- MDAnalysis
- pytraj

## Example:
```python
from mutresid import mutate_residue
from genAmber import generate_amber_files
pdb_mutated = mutate_residue(f"{yourPdb}", residue_id=2, new_residue_name="HIS")
generate_amber_files(pdb=pdb_mutated,extra_param=f"{yourParam.mol2}",extra_name=f"{resName}",addSolvent=True,verbose=False)
```
