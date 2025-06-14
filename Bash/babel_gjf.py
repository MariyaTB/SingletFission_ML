from openbabel import openbabel as ob

# Interactive user input
smiles = input("Enter SMILES string: ")
mol_name = input("Enter molecule name (default: Molecule): ") or "Molecule"
charge = int(input("Enter charge (default: 0): ") or 0)
mult = int(input("Enter multiplicity (default: 1): ") or 1)
num_conformers = int(input("Enter number of conformers (default: 2): ") or 2)

# Create molecule from SMILES using Open Babel 3.x API
obConversion = ob.OBConversion()
obConversion.SetInAndOutFormats("smi", "mol")
mol = ob.OBMol()
obConversion.ReadString(mol, smiles)

# Add hydrogens and generate 3D coordinates
mol.AddHydrogens()
builder = ob.OBBuilder()
builder.Build(mol)
ff = ob.OBForceField.FindForceField('mmff94')
ff.Setup(mol)
ff.ConjugateGradients(200)
ff.GetCoordinates(mol)

# Generate multiple conformers
ob_mol = mol  # mol is already an OBMol object
ff = ob.OBForceField.FindForceField('mmff94')
ff.Setup(ob_mol)
confgen = ob.OBConformerSearch()
confgen.Setup(ob_mol, num_conformers)
confgen.Search()
confgen.GetConformers(ob_mol)

# Write Gaussian input (.gjf)
filename = f"{mol_name}_openbabel.gjf"
with open(filename, 'w') as f:
    f.write(f"%nprocshared=4\n%mem=4GB\n#p opt b3lyp/6-31g(d)\n\n{mol_name}\n\n")
    f.write(f"{charge} {mult}\n")

    for conf_idx in range(ob_mol.NumConformers()):
        ob_mol.SetConformer(conf_idx)
        for atom in ob.OBMolAtomIter(ob_mol):
            pos = atom.GetVector()
            symbol = ob.GetSymbol(atom.GetAtomicNum())
            f.write(f"{symbol} {pos.GetX():.6f} {pos.GetY():.6f} {pos.GetZ():.6f}\n")
        if conf_idx < num_conformers - 1:
            f.write("\n--Link1--\n\n")
    f.write("\n")

print(f"\n✅ Generated: {filename}")