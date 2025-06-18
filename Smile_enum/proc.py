import pandas as pd
import random

df = pd.read_csv('enumerated_validated_SF_smiles_5.csv')
if 'SF' in df.columns:
    df = df.drop(columns=['SF'])
df = df.sample(frac=1, random_state=42).reset_index(drop=True)
df.to_csv('output_SMILES.csv', index=False)