# SingletFission_ML

This project aims to use machine learning techniques to generate and evaluate molecular structures for potential singlet fission materials. Singlet fission is a process where a singlet exciton splits into two triplet excitons, which can potentially double the efficiency of solar cells. 

## Project Structure

- **Data Augmentation**: This module is responsible for augmenting the dataset of SMILES strings by generating new molecules through various techniques.
- **Model Training**: This module includes scripts for pre-training and fine-tuning a Character-Level Model (CLM) on the augmented dataset.
- **Evaluation**: This module evaluates the trained model and generates new SMILES strings.

## Files

- `main.py`: The main script that orchestrates the entire workflow, including data fetching, augmentation, model training, and evaluation.
- `Data_Augmentation/Augment.py`: Contains functions for augmenting the dataset by generating new SMILES strings.
- `utils/dataset.py`: Defines the dataset class for handling SMILES strings.
- `models/clm_model.py`: Defines the Character-Level Model (CLM) architecture.
- `utils/train_utils.py`: Utility functions for training and evaluating the model.
- `utils/clm_utils.py`: Utility functions for pre-training and fine-tuning the CLM.
_ (All of this is considered as work-in-progress )

## Installation

To run this project, you need to have the following dependencies installed:

- Python 3.8+
- pandas
- cupy
- rdkit
- torch
- tqdm

You can install the required packages using pip:

```bash
pip install pandas cupy rdkit-pypi torch tqdm
