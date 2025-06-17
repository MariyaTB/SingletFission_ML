from openbabel import openbabel as ob

smiles = input("Enter SMILES string: ")
mol_name = input("Enter molecule name (default: Molecule): ") or "Molecule"
charge = int(input("Enter charge (default: 0): ") or 0)
mult = int(input("Enter multiplicity (default: 1): ") or 1)
num_conformers = int(input("Enter number of conformers (default: 2): ") or 2)
chk_file = input("Enter checkpoint file path (default: None): ")

obConversion = ob.OBConversion()
obConversion.SetInFormat("smi")
mol = ob.OBMol()
obConversion.ReadString(mol, smiles)

mol.AddHydrogens()
builder = ob.OBBuilder()
builder.Build(mol)
ff = ob.OBForceField.FindForceField('mmff94')
ff.Setup(mol)
ff.ConjugateGradients(200)
ff.GetCoordinates(mol)

conf_search = ob.OBConformerSearch()
conf_search.Setup(mol, num_conformers)
conf_search.Search()
conf_search.GetConformers(mol)

filename = f"{mol_name}_openbabel.gjf"
with open(filename, 'w') as f:
    if chk_file:
        f.write(f"%chk={chk_file}\n")
    f.write("# b3lyp/6-31g(d,p) geom=connectivity\n\n")
    f.write(f"Title Card Required\n\n")
    f.write(f"{charge} {mult}\n")

    mol.SetConformer(0)
    atoms = list(ob.OBMolAtomIter(mol))
    for atom in atoms:
        pos = atom.GetVector()
        symbol = ob.GetSymbol(atom.GetAtomicNum())
        f.write(f" {symbol:<2} {pos.GetX():>18.8f} {pos.GetY():>18.8f} {pos.GetZ():>18.8f}\n")

    f.write("\n")
    for atom in ob.OBMolAtomIter(mol):
        aid = atom.GetIdx() + 1   
        bonds = []
        for bond in ob.OBAtomBondIter(atom):
            other_atom = bond.GetNbrAtomIdx(atom) + 1
            order = bond.GetBondOrder()
            bonds.append(f"{other_atom} {order}")
        if bonds:
            f.write(f" {aid} {' '.join(bonds)}\n")
        else:
            f.write(f" {aid}\n")

print(f"\n✅ Generated: {filename}")