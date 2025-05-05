import pandas as pd
from rdkit import Chem
from rdkit import RDLogger # To suppress RDKit warnings if needed
import multiprocessing
from functools import partial
from tqdm import tqdm
import time
import os
import argparse # For command-line arguments

# --- Configuration ---
DEFAULT_INPUT_CSV = 'SMILE.csv' # Your input file name
DEFAULT_OUTPUT_CSV = 'enumerated_validated_ST_smiles.csv' # Descriptive output name
ENUMERATION_LEVEL = 100 # Number of SMILES to generate per input molecule
NUM_WORKERS = max(1, os.cpu_count() - 2) # Use most cores, leave some for OS

# Suppress RDKit warnings (optional, can be noisy)
RDLogger.DisableLog('rdApp.*')

def enumerate_single_smiles(input_smiles: str, num_enumerations: int) -> list[str]:
    """
    Generates a specified number of random SMILES strings for a single input SMILES.

    Args:
        input_smiles: The original SMILES string.
        num_enumerations: The target number of random SMILES to generate.

    Returns:
        A list of generated random SMILES strings. Returns an empty list
        if the input SMILES is invalid or enumeration fails.
    """
    enumerated = []
    try:
        mol = Chem.MolFromSmiles(input_smiles)
        if mol is not None:
            # RDKit's MolToRandomSmilesVect handles the randomization internally
            smiles_vect = Chem.MolToRandomSmilesVect(mol, num_enumerations)
            # Ensure we only take up to num_enumerations even if it generates more somehow
            enumerated.extend(list(smiles_vect)[:num_enumerations])
        else:
            # Don't print warning here, handle in main loop based on empty list return
            return [] # Return empty list for invalid input SMILES
    except Exception as e:
        # Error during processing (less common for MolFromSmiles/MolToRandomSmilesVect)
        # print(f"Error processing original SMILES '{input_smiles}': {e}") # Can be noisy
        return [] # Return empty list on error
    return enumerated

def run_enumeration_and_validate(input_csv: str, output_csv: str, enumeration_level: int, num_workers: int):
    """
    Main function to load data (SMILES, S, T), run parallel enumeration,
    validate enumerated SMILES, and save results with original S and T values.
    """
    print(f"Starting SMILES enumeration and validation for Singlet/Triplet data...")
    print(f"Input file: {input_csv}")
    print(f"Output file: {output_csv}")
    print(f"Enumeration level per molecule: {enumeration_level}")
    print(f"Using {num_workers} worker processes.")

    start_time = time.time()

    # 1. Load Data
    try:
        df_input = pd.read_csv(input_csv)
        print(f"Loaded {len(df_input)} molecules from {input_csv}.")
        # Ensure column names are correct
        required_cols = ['SMILES', 'S', 'T']
        if not all(col in df_input.columns for col in required_cols):
            raise ValueError(f"Input CSV must contain columns: {', '.join(required_cols)}")

        # Fill potential NaN SMILES with an empty string to avoid errors later
        df_input['SMILES'] = df_input['SMILES'].fillna('')

    except FileNotFoundError:
        print(f"Error: Input file '{input_csv}' not found.")
        return
    except Exception as e:
        print(f"Error loading input CSV: {e}")
        return

    original_smiles_list = df_input['SMILES'].tolist()
    original_s_values = df_input['S'].tolist()
    original_t_values = df_input['T'].tolist()

    # 2. Parallel Enumeration Setup
    # Use partial to fix the 'num_enumerations' argument for the worker function
    worker_func = partial(enumerate_single_smiles, num_enumerations=enumeration_level)

    results_list = [] # To store lists of enumerated SMILES for each input
    print(f"Starting parallel enumeration across {num_workers} cores...")

    # Use multiprocessing Pool
    with multiprocessing.Pool(processes=num_workers) as pool:
        # Use pool.imap for potentially better memory usage with large iterables,
        # and tqdm for progress bar.
        results_list = list(tqdm(pool.imap(worker_func, original_smiles_list),
                                 total=len(original_smiles_list),
                                 desc="Enumerating SMILES"))

    print("Parallel enumeration finished.")

    # 3. Combine Results and Validate Enumerated SMILES
    print("Combining results and validating enumerated SMILES...")
    output_data = []
    total_generated_count = 0 # Count before validation
    total_valid_count = 0     # Count after validation
    invalid_enumerated_count = 0
    skipped_original_count = 0

    for i, enumerated_for_mol in enumerate(tqdm(results_list, desc="Validating SMILES")):
        original_smiles = original_smiles_list[i] # Keep track of original for reporting
        original_s = original_s_values[i]
        original_t = original_t_values[i]

        # Check if enumeration failed for the original SMILES (returned empty list)
        if not enumerated_for_mol and original_smiles: # Check original_smiles to ensure it wasn't NaN initially
             skipped_original_count += 1
             # print(f"Warning: Original SMILES '{original_smiles}' (index {i}) could not be parsed or enumerated.") # Can be noisy
             continue
        elif not original_smiles: # Handle cases where the original SMILES itself was empty/NaN
            skipped_original_count +=1
            continue


        total_generated_count += len(enumerated_for_mol)

        for enum_smiles in enumerated_for_mol:
            # --- Validation Step ---
            mol_check = Chem.MolFromSmiles(enum_smiles)
            if mol_check is not None:
                # SMILES is valid, add it to the output with original S and T
                output_data.append({'SMILES': enum_smiles, 'S': original_s, 'T': original_t})
                total_valid_count += 1
            else:
                # SMILES is invalid, skip it and count it
                invalid_enumerated_count += 1
                # Optional: print(f"  Skipping invalid enumerated SMILES: {enum_smiles}") # Can be very verbose

    df_output = pd.DataFrame(output_data)
    print(f"\n--- Summary ---")
    print(f"Processed {len(original_smiles_list)} original molecules.")
    print(f"Skipped {skipped_original_count} original molecules (invalid SMILES or enumeration error).")
    print(f"Total enumerated SMILES generated (before validation): {total_generated_count}")
    print(f"Number of valid enumerated SMILES kept: {total_valid_count}")
    print(f"Number of invalid enumerated SMILES discarded: {invalid_enumerated_count}")
    print(f"Final validated dataset size: {len(df_output)} rows.")


    # 4. Save Output
    try:
        df_output.to_csv(output_csv, index=False)
        print(f"\nSuccessfully saved validated enumerated SMILES to {output_csv}.")
    except Exception as e:
        print(f"Error saving output CSV: {e}")

    end_time = time.time()
    print(f"Total execution time: {end_time - start_time:.2f} seconds.")


# --- Main Execution Block ---
if __name__ == "__main__":
    # Setup argparse for command-line arguments
    parser = argparse.ArgumentParser(description="Perform SMILES enumeration and validate structures, carrying over S and T values.")
    parser.add_argument("-i", "--input", type=str, default=DEFAULT_INPUT_CSV,
                        help=f"Input CSV file name (must contain SMILES, S, T columns) (default: {DEFAULT_INPUT_CSV})")
    parser.add_argument("-o", "--output", type=str, default=DEFAULT_OUTPUT_CSV,
                        help=f"Output CSV file name (default: {DEFAULT_OUTPUT_CSV})")
    parser.add_argument("-n", "--num", type=int, default=ENUMERATION_LEVEL,
                        help=f"Number of enumerations per molecule (default: {ENUMERATION_LEVEL})")
    parser.add_argument("-w", "--workers", type=int, default=NUM_WORKERS,
                        help=f"Number of worker processes (default: {NUM_WORKERS})")

    args = parser.parse_args()

    run_enumeration_and_validate(args.input, args.output, args.num, args.workers)