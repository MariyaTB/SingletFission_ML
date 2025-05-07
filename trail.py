from rdkit import Chem

smiles = "CC(=O)N"

mol = Chem.MolFromSmiles(smiles)

for atom in mol.GetAtoms():
    print(f"Atom: {atom.GetSymbol()}")

for bond in mol.GetBonds():
    print(f"Bond: {bond.GetBeginAtom().GetSymbol()} - {bond.GetEndAtom().GetSymbol()} ({bond.GetBondTypeAsDouble()})")