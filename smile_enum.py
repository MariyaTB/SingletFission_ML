import pandas as pd
from rdkit import Chem
from rdkit import RDLogger  
import multiprocessing
from functools import partial
from tqdm import tqdm
import time
import os
import argparse  


DEFAULT_INPUT_CSV = 'SMILE.csv'
DEFAULT_OUTPUT_CSV = 'enumerated_validated_smiles.csv' 
ENUMERATION_LEVEL = 50 
NUM_WORKERS = max(1, os.cpu_count() - 2)  

RDLogger.DisableLog('rdApp.*')

def enumerate_single_smiles(input_smiles: str, num_enumerations: int) -> list[str]:
    
    enumerated = []
    try:
        mol = Chem.MolFromSmiles(input_smiles)
        if mol is not None:
            smiles_vect = Chem.MolToRandomSmilesVect(mol, num_enumerations)
            enumerated.extend(list(smiles_vect))
        else:
            print(f"Warning: Could not parse original SMILES: {input_smiles}")
            return [] 
    except Exception as e:
        print(f"Error processing original SMILES '{input_smiles}': {e}")
        return []  
    return enumerated

def run_enumeration_and_validate(input_csv: str, output_csv: str, enumeration_level: int, num_workers: int):
   
    print(f"Starting SMILES enumeration and validation...")
    print(f"Input file: {input_csv}")
    print(f"Output file: {output_csv}")
    print(f"Enumeration level per molecule: {enumeration_level}")
    print(f"Using {num_workers} worker processes.")

    start_time = time.time()

    try:
        df_input = pd.read_csv(input_csv)
        print(f"Loaded {len(df_input)} molecules from {input_csv}.")
        if 'SMILES' not in df_input.columns or 'V' not in df_input.columns:
            raise ValueError("Input CSV must contain 'SMILES' and 'V' columns.")
    except FileNotFoundError:
        print(f"Error: Input file '{input_csv}' not found.")
        return
    except Exception as e:
        print(f"Error loading input CSV: {e}")
        return

    original_smiles_list = df_input['SMILES'].tolist()
    original_velocities = df_input['V'].tolist()

    worker_func = partial(enumerate_single_smiles, num_enumerations=enumeration_level)

    results_list = [] 
    print(f"Starting parallel enumeration across {num_workers} cores...")


    with multiprocessing.Pool(processes=num_workers) as pool:
        
        results_list = list(tqdm(pool.imap(worker_func, original_smiles_list),
                                 total=len(original_smiles_list),
                                 desc="Enumerating SMILES"))

    print("Parallel enumeration finished.")

    print("Combining results and validating enumerated SMILES...")
    output_data = []
    total_generated_count = 0  
    total_valid_count = 0     
    invalid_enumerated_count = 0
    skipped_original_count = 0

    for i, enumerated_for_mol in enumerate(results_list):
        original_velocity = original_velocities[i]
        if not enumerated_for_mol: 
            skipped_original_count += 1
            print(f"Warning: No enumerated SMILES generated for original SMILES: {original_smiles_list[i]}")
            continue

        total_generated_count += len(enumerated_for_mol)

        for enum_smiles in enumerated_for_mol:

            mol_check = Chem.MolFromSmiles(enum_smiles)
            if mol_check is not None:

                output_data.append({'SMILES': enum_smiles, 'V': original_velocity})
                total_valid_count += 1
            else:

                invalid_enumerated_count += 1
                

    df_output = pd.DataFrame(output_data)
    print(f"Total enumerated SMILES generated (before validation): {total_generated_count}")
    print(f"Number of valid enumerated SMILES kept: {total_valid_count}")
    print(f"Number of invalid enumerated SMILES discarded: {invalid_enumerated_count}")
    if skipped_original_count > 0:
        print(f"Warning: Skipped {skipped_original_count} original molecules due to initial parsing/enumeration errors.")
    print(f"Final validated dataset size: {len(df_output)} rows.")



    try:
        df_output.to_csv(output_csv, index=False)
        print(f"Successfully saved validated enumerated SMILES to {output_csv}.")
    except Exception as e:
        print(f"Error saving output CSV: {e}")

    end_time = time.time()
    print(f"Total execution time: {end_time - start_time:.2f} seconds.")



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Perform SMILES enumeration and validate structures in parallel.")
    parser.add_argument("-i", "--input", type=str, default=DEFAULT_INPUT_CSV,
                        help=f"Input CSV file name (default: {DEFAULT_INPUT_CSV})")
    parser.add_argument("-o", "--output", type=str, default=DEFAULT_OUTPUT_CSV,
                        help=f"Output CSV file name (default: {DEFAULT_OUTPUT_CSV})")
    parser.add_argument("-n", "--num", type=int, default=ENUMERATION_LEVEL,
                        help=f"Number of enumerations per molecule (default: {ENUMERATION_LEVEL})")
    parser.add_argument("-w", "--workers", type=int, default=NUM_WORKERS,
                        help=f"Number of worker processes (default: {NUM_WORKERS})")

    args = parser.parse_args()

    run_enumeration_and_validate(args.input, args.output, args.num, args.workers)