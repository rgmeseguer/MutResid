# MutResid
A python module that uses MDanalysis and to take a protein structure in PDB and mutates one residue returning a PDB and then runs it trough tleap to return the required Amber files.

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
