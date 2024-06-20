# **ReadMe**



## Issues to be solved

1. MCS Calculation using Rdkit

   The following code may find inappropriate maximum common substructure for molecules containing fused rings, bridged bicyclic rings, and spiro bicyclic rings.

   ```python
   mcs_SMARTS = rdFMCS.FindMCS([mol_1, mol_2], ringMatchesRingOnly=True, completeRingsOnly=True)
   ```

   For example, the MCS found for O=C(NC1(O)C(=O)C2=CC=CC=C2C1=O)C1=CC=C(C(F)(F)F)C=C1 and C12(NC(=O)c3ccc(cc3)C)CC3CC(C1)CC(C3)CC2 is Cc1ccc(C(=O)NC2CCCCCCCC2)cc1

   

