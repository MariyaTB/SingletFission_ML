import os
import pandas as pd
from rdkit import Chem

smi_file_path = os.path.join('Data_Augmentation', 'families', 'family5.smi')

with open(smi_file_path, 'r') as file:
    lines = file.readlines()

molecules = [Chem.MolFromSmiles(line.strip().split()[0]) for line in lines]
smiles_list = [Chem.MolToSmiles(mol) if mol else 'Invalid' for mol in molecules]

df = pd.DataFrame({'SMILES': smiles_list})

csv_file_path = os.path.join('Data_Augmentation', 'family5.csv')
df.to_csv(csv_file_path, index=False)

print(f'SMILES data saved to {csv_file_path}')