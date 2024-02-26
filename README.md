# MutResid
A python script that uses MDanalysis and pytraj to take a protein structure in PDB or amber format and mutates one residue returning a PDB and the amber required files.

## Example:
```python
from mutresid import mutate_residue
mutate_residue("protein.pdb", residue_id=264, new_residue_name="HIS")
```
