from rdkit import Chem
from rdkit.Chem import AllChem


smiles = input("Enter SMILES string: ")
mol_name = input("Enter molecule name (default: Molecule): ") or "Molecule"
charge = int(input("Enter charge (default: 0): ") or 0)
mult = int(input("Enter multiplicity (default: 1): ") or 1)
num_conformers = int(input("Enter number of conformers (default: 2): ") or 2)


mol = Chem.MolFromSmiles(smiles)
if mol is None:
    raise ValueError("Invalid SMILES string")

# Add hydrogens
mol = Chem.AddHs(mol)

# Generate conformers
ps = AllChem.ETKDG()
ps.numConformers = num_conformers
cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=ps)

# Optimize each conformer
for cid in cids:
    AllChem.MMFFOptimizeMolecule(mol, confId=cid)

# Write .gjf file
with open(f"{mol_name}_rdkit.gjf", "w") as f:
    f.write(f"%nprocshared=4\n%mem=4GB\n#p opt b3lyp/6-31g(d)\n\n{mol_name}\n\n")
    f.write(f"{charge} {mult}\n")

    for i, cid in enumerate(cids):
        conf = mol.GetConformer(cid)
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            symbol = atom.GetSymbol()
            f.write(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")
        if i < len(cids) - 1:
            f.write("\n--Link1--\n\n")
    f.write("\n")

print(f"\n✅ Generated: {mol_name}_rdkit.gjf")