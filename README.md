# MutResid
A python script that uses MDanalysis and to take a protein structure in PDB and mutates one residue returning a PDB that can be run in tleap to generate the amber required files.

## Example:
```python
from mutresid import mutate_residue
mutate_residue("protein.pdb", residue_id=264, new_residue_name="HIS")
```
