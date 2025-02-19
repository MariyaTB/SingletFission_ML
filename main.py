import torch
import logging
from utils.chembl_utils import get_chembl_smiles
from utils.dataset import SMILESDataset
from models.clm_model import CLM
from utils.train_utils import train_model, generate_smiles, evaluate_model
from utils.clm_utils import pretrain_clm, finetune_clm
import pandas as pd
from rdkit import Chem
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    logging.info("Starting main function")
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logging.info(f"Using device: {device}")

    logging.info("Fetching SMILES from ChEMBL...")
    chembl_smiles = get_chembl_smiles(limit=50000) 
    logging.info(f"Collected {len(chembl_smiles)} SMILES from ChEMBL")
    valid_chembl = []
    for smi in chembl_smiles:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                valid_chembl.append(Chem.MolToSmiles(mol))
        except Exception as e:
            logging.warning(f"Invalid SMILES string skipped: {smi}, error: {e}")
    logging.info(f"Collected {len(valid_chembl)} valid SMILES from ChEMBL")

    logging.info("Loading augmented dataset...")
    try:
        augmented_data = pd.read_csv('data/augmented.csv')
        logging.info("Augmented dataset loaded")
    except FileNotFoundError:
        logging.error("Augmented dataset not found. Please check the file path and ensure the file exists.")
        return

    augmented_smiles = augmented_data['SMILES'].tolist()
    valid_augmented = [Chem.MolToSmiles(Chem.MolFromSmiles(smi)) for smi in augmented_smiles if Chem.MolFromSmiles(smi)]
    logging.info(f"Collected {len(valid_augmented)} valid SMILES from augmented dataset")

    combined_smiles = valid_chembl + valid_augmented
    logging.info(f"Combined {len(combined_smiles)} valid SMILES from ChEMBL and augmented dataset")

    logging.info("Creating vocabulary...")
    chars = sorted(list(set(''.join(combined_smiles)))) + ['<PAD>', '<EOS>']
    logging.info(f"Vocabulary size: {len(chars)}")
    char_to_idx = {c: i for i, c in enumerate(chars)}
    idx_to_char = {i: c for i, c in enumerate(chars)}

    logging.info("Pre-training CLM on combined dataset...")
    clm_model = CLM(vocab_size=len(chars), embed_dim=128, hidden_dim=256, n_layers=3, dropout=0.2).to(device)
    pretrain_clm(clm_model, combined_smiles, char_to_idx, device, epochs=5, batch_size=128, lr=0.001)
    logging.info("CLM pre-training complete")

    logging.info("Loading smaller dataset for fine-tuning...")
    small_dataset = pd.read_csv('data/SMILE.csv')
    small_smiles = small_dataset['SMILES'].tolist()
    valid_small_smiles = [Chem.MolToSmiles(Chem.MolFromSmiles(smi)) for smi in small_smiles if Chem.MolFromSmiles(smi)]
    logging.info(f"Loaded {len(valid_small_smiles)} valid SMILES from smaller dataset")

    logging.info("Fine-tuning CLM on smaller dataset...")
    finetune_clm(clm_model, valid_small_smiles, char_to_idx, device, epochs=20, batch_size=32, lr=0.0001)
    logging.info("CLM fine-tuning complete")

    # Evaluate CLM
    logging.info("Evaluating CLM...")
    dataset = SMILESDataset(combined_smiles, char_to_idx)
    loader = torch.utils.data.DataLoader(dataset, batch_size=32, shuffle=True)
    criterion = torch.nn.CrossEntropyLoss(ignore_index=char_to_idx['<PAD>'])
    evaluate_model(clm_model, loader, criterion, device)

    # Generate new SMILES
    logging.info("Generating new SMILES...")
    new_smiles = [generate_smiles(clm_model, char_to_idx, idx_to_char) for _ in range(1000)]
    logging.info("New SMILES generated")
    pd.DataFrame(new_smiles, columns=['SMILES']).to_csv('data/generated_smiles_17_02_2025(2).csv', index=False)
    logging.info("Generated SMILES saved to 'data/generated_smiles_17_02_2025(2).csv'")

    # Filter valid and unique SMILES
    logging.info("Filtering valid and unique SMILES...")
    valid_new_smiles = []
    seen = set(valid_small_smiles)
    for smi in new_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol and smi not in seen:
            seen.add(smi)
            valid_new_smiles.append(smi)
    logging.info(f"Generated {len(valid_new_smiles)} valid and unique SMILES")


    pd.DataFrame(valid_new_smiles, columns=['SMILES']).to_csv('data/valid_unique_smiles_17_02_2025(2).csv', index=False)
    logging.info("Valid and unique SMILES saved to 'data/valid_unique_smiles_17_02_2025(2).csv'")

if __name__ == "__main__":
    main()
