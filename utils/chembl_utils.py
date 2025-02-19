
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
import logging

def get_chembl_smiles(limit=50000):
    molecule = new_client.molecule
    logging.info(f"Attempting to retrieve up to {limit} molecules from ChEMBL...")

    # Retrieve molecules with a preferred name
    molecules = molecule.filter(pref_name__isnull=False)[:limit]
    logging.info(f"Retrieved {len(molecules)} molecules with a preferred name.")

    # Filter molecules with valid canonical SMILES
    valid_smiles = [mol['molecule_structures']['canonical_smiles']
                    for mol in molecules
                    if mol and 'molecule_structures' in mol and mol['molecule_structures'] is not None and 'canonical_smiles' in mol['molecule_structures']]
    logging.info(f"Collected {len(valid_smiles)} valid SMILES from ChEMBL.")

    return valid_smiles