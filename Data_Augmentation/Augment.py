import pandas as pd
import cupy as cp
from rdkit import Chem
from rdkit.Chem import Recap
from rdkit.Chem import AllChem as Chem
import random

cp.cuda.set_allocator(cp.cuda.MemoryPool().malloc)

def get_rand_smi(smi, max_folds=10):
    mol = Chem.MolFromSmiles(smi)
    smi_list = []
    for _ in range(max_folds):
        try:
            s = Chem.MolToSmiles(mol, doRandom=True)
            smi_list.append(s)
        except BaseException:
            continue
    smi_list = set(smi_list)
    max_folds = min(max_folds, len(smi_list))
    smi_list = smi_list - {smi}
    smi_list = list(smi_list)[0:max_folds - 1] + [smi]
    return smi_list

def augmentation_by_smi(lst_smi, max_folds=10):
    if max_folds <= 1:
        return lst_smi
    list_of_augmented_smi = []
    for idx in range(len(lst_smi)):
        smi = lst_smi[idx]
        list_of_rand_smi = get_rand_smi(smi, max_folds)
        list_of_augmented_smi += list_of_rand_smi
    return list_of_augmented_smi

def replace_number(smi: str, n=1) -> str:
    s = ''
    for char in smi:
        new_char = str(int(char) + n) if char.isdigit() else char
        s += new_char
    return s

def combine_fragments(fragments: list, n_combs: int):
    list_of_smi = []
    singles = [frag for frag in fragments if frag.count('*') == 1]
    n_frags, n_singles = len(fragments), len(singles)
    indices = cp.random.randint(0, n_frags, n_combs).get() 
    for idx in indices:
        bone = fragments[idx]
        bone = replace_number(bone)
        for i in range(bone.count('*')):
            frag = singles[random.randint(0, n_singles - 1)]
            frag = frag.replace('*', '')
            bone = bone.replace('*', frag, 1)
        bone = bone.replace('*', '')
        try:
            mol = Chem.MolFromSmiles(bone)
            Chem.SanitizeMol(mol)
            canonical_smi = Chem.MolToSmiles(mol, canonical=True)
            list_of_smi.append(canonical_smi)
        except:
            pass
    return list_of_smi

def augmentation_by_fragment(list_of_smi: list, n: int):
    fragments = set()
    for idx, smi in enumerate(list_of_smi):
        m = Chem.MolFromSmiles(smi)
        hierarch = Recap.RecapDecompose(m)
        leaves = list(hierarch.GetLeaves().keys())
        fragments.update(leaves)
    fragments = list(fragments)
    return combine_fragments(fragments, n)

def filter_small_molecules(smiles_list, max_atoms=28):
    filtered_smiles = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol and mol.GetNumAtoms() <= max_atoms:
            filtered_smiles.append(smi)
    return filtered_smiles

core_smiles_df = pd.read_csv(r'D:\My_Computer\Meh!!!!!!!\workspace\VIT_LAB\SingletFission_ML\data\SMILE.csv')
core_smiles_list = core_smiles_df['SMILES'].tolist()

augmented_smiles = augmentation_by_fragment(core_smiles_list, 1000000)

filtered_augmented_smiles = filter_small_molecules(augmented_smiles)

unique_augmented_smiles = list(set(filtered_augmented_smiles))

augmented_df = pd.DataFrame(unique_augmented_smiles, columns=['SMILES'])
augmented_df.to_csv(r'D:\My_Computer\Meh!!!!!!!\workspace\VIT_LAB\SingletFission_ML\Data_Augmentation\augmented_family1.csv', index=False)

print(f"Generated {len(unique_augmented_smiles)} unique augmented SMILES with up to 25 atoms.")


















# filtered_augmented_smiles = filter_small_molecules(augmented_smiles)

# augmented_df = pd.DataFrame(filtered_augmented_smiles, columns=['SMILES'])
# augmented_df.to_csv(r'D:\My_Computer\Meh!!!!!!!\workspace\VIT_LAB\SingletFission_ML\Data_Augmentation\augmented_family1.csv', index=False)

# print(f"Generated {len(filtered_augmented_smiles)} augmented SMILES with up to 25 atoms.")