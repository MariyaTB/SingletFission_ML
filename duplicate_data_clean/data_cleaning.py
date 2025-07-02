import pandas as pd

df_a = pd.read_csv('SF_Files/valid_unique_smiles_28_02_2025.csv')
df_b = pd.read_csv('SF_Files/valid_unique_smiles_01_04_2025.csv')

smiles_col = df_a.columns[1]

duplicates = df_a[df_a[smiles_col].isin(df_b[smiles_col])]
print(f"Number of duplicates: {len(duplicates)}")

df_a_filtered = df_a[~df_a[smiles_col].isin(df_b[smiles_col])]

result = pd.concat([df_b, df_a_filtered], ignore_index=True)

result.to_csv('SF_Files/valid_unique_smiles_28_02_2025.csv', index=False)


# import pandas as pd

# df_a = pd.read_csv('duplicate_data_clean/final_filtered_smiles(2).csv')
# df_b = pd.read_csv('duplicate_data_clean/cleaned_valid_unique_smiles.csv')

# print(df_a.columns)  # Debug: see what columns are present

# # Use the correct column index or name
# smiles_col = df_a.columns[0]  # or replace 0 with the correct index, or use the column name

# duplicates = df_a[df_a[smiles_col].isin(df_b[smiles_col])]
# print(f"Number of duplicates: {len(duplicates)}")

# df_a_filtered = df_a[~df_a[smiles_col].isin(df_b[smiles_col])]

# result = pd.concat([df_a_filtered,df_b], ignore_index=True)

# result.to_csv('duplicate_data_clean/final_filtered_smiles(2).csv', index=False)