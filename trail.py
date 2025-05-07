from rdkit import Chem
from rdkit.Chem import Draw  # Import Draw module for visualization
from IPython.display import display  # Import display for inline visualization

smiles = "CC(=O)N"

mol = Chem.MolFromSmiles(smiles)

for atom in mol.GetAtoms():
    print(f"Atom: {atom.GetSymbol()}")

for bond in mol.GetBonds():
    print(f"Bond: {bond.GetBeginAtom().GetSymbol()} - {bond.GetEndAtom().GetSymbol()} ({bond.GetBondTypeAsDouble()})")

# Visualize the molecule as a graph inline
img = Draw.MolToImage(mol)
display(img)