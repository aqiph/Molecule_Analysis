from rdkit import Chem

smiles1 = 'CCO'
smiles2 = 'CN1C(NC2=NC=CC=C2)=CC=C1'
smiles3 = 'CN1C(NC2=NC=CC=C2)=CC=C1'
#smiles4 = '[C]'

molecule1 = Chem.MolFromSmiles(smiles1)
molecule2 = Chem.MolFromSmiles(smiles2)
molecule3 = Chem.MolFromSmiles(smiles3)
molecule4 = Chem.MolFromSmiles(smiles4) 

#order1 = Chem.CanonicalRankAtoms(molecule1)
#order2 = Chem.CanonicalRankAtoms(molecule2)
#order3 = Chem.CanonicalRankAtoms(molecule3)

#molecule_can1 = Chem.RenumberAtoms(molecule1, order1)
#molecule_can2 = Chem.RenumberAtoms(molecule2, order2)
#molecule_can3 = Chem.RenumberAtoms(molecule3, order3)

#smiles_can1 = Chem.MolToSmiles(molecule_can1)
#smiles_can2 = Chem.MolToSmiles(molecule_can2)
#smiles_can3 = Chem.MolToSmiles(molecule_can3)

Chem.Kekulize(molecule1)
#Chem.Kekulize(molecule2)
Chem.Kekulize(molecule3)

s1 = Chem.MolToSmiles(molecule1, kekuleSmiles = True)
s2 = Chem.MolToSmiles(molecule2, kekuleSmiles = False)
s3 = Chem.MolToSmiles(molecule3, kekuleSmiles = True)

#molecule4 = Chem.AddHs(molecule4)
#s4 = Chem.MolToSmiles(molecule4)

print(smiles1, s1)
print(smiles2, s2)
print(smiles3, s3)
#print(smiles4, s4)

