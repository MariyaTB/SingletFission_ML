import pandas as pd
from rdkit import Chem

input_file = 'duplicate_data_clean/valid_unique_smiles_02_07_2025(2).csv'
output_file = 'duplicate_data_clean/cleaned_valid_unique_smiles.csv'

df = pd.read_csv(input_file)

deleted_count = 0
kept_rows = []

for index, row in df.iterrows():
    smi = row['SMILES']
    mol = Chem.MolFromSmiles(smi)
    
    if mol is None:
        deleted_count += 1
        continue
    
    heavy_atoms = mol.GetNumHeavyAtoms()
    
    if heavy_atoms >= 15:
        kept_rows.append(row)
    else:
        deleted_count += 1

cleaned_df = pd.DataFrame(kept_rows)

cleaned_df.to_csv(output_file, index=False)

total_entries = len(df)
remaining_entries = len(cleaned_df)

print(f"Total entries in original file: {total_entries}")
print(f"Removed {deleted_count} molecules.")
print(f"Entries remaining: {remaining_entries}")