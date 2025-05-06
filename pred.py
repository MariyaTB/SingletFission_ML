import pandas as pd
import numpy as np
from tqdm.notebook import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem

import torch
import torch.nn.functional as F
from torch_geometric.data import Data, Dataset, DataLoader
from torch_geometric.nn import GCNConv, global_mean_pool, GATConv, GINConv
from torch.nn import Linear, Sequential, ReLU, BatchNorm1d
import torch.optim as optim
from torch.optim.lr_scheduler import ReduceLROnPlateau

from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

import matplotlib.pyplot as plt

CSV_FILE = ''
TARGET_COLUMN = ''
SMILES_COLUMN = ''
TEST_SIZE = 0.2
VALID_SIZE = 0.1
RANDOM_STATE = 42
DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {DEVICE}")

LEARNING_RATE = 0.001
BATCH_SIZE = 64
NUM_EPOCHS = 100
HIDDEN_DIM = 128
NUM_GNN_LAYERS = 3
DROPOUT_RATE = 0.1

def get_atom_features(atom):
    possible_atoms = list(range(1, 119))
    atomic_num = atom.GetAtomicNum()
    atom_one_hot = [1.0 if i == atomic_num else 0.0 for i in possible_atoms]

    features = atom_one_hot + [
        atom.GetDegree() / 10.0,
        atom.GetFormalCharge() / 5.0,
        atom.GetTotalNumHs() / 8.0,
        float(atom.GetHybridization()),
        float(atom.GetIsAromatic()),
        float(atom.IsInRing()),
        atom.GetMass() * 0.01,
        atom.GetDoubleProp('_GasteigerCharge') if atom.HasProp('_GasteigerCharge') else 0.0
    ]
    features = [0.0 if np.isnan(f) else f for f in features]
    return torch.tensor(features, dtype=torch.float)

def smiles_to_graph(smiles, y_target):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Warning: Could not parse SMILES: {smiles}")
        return None

    try:
        mol = Chem.AddHs(mol)
        AllChem.ComputeGasteigerCharges(mol)
    except Exception as e:
        print(f"Warning: Could not process molecule {smiles}: {e}")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return None

    atom_features_list = [get_atom_features(atom) for atom in mol.GetAtoms()]
    x = torch.stack(atom_features_list, dim=0)

    edge_indices = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        edge_indices.append((i, j))
        edge_indices.append((j, i))

    edge_index = torch.tensor(edge_indices, dtype=torch.long).t().contiguous()
    y = torch.tensor([[y_target]], dtype=torch.float)
    data = Data(x=x, edge_index=edge_index, y=y)
    return data

print("Loading data...")
df = pd.read_csv(CSV_FILE)
print(f"Loaded {len(df)} molecules.")

print("Converting SMILES to graphs...")
data_list = []
for index, row in tqdm(df.iterrows(), total=len(df), desc="Processing SMILES"):
    graph = smiles_to_graph(row[SMILES_COLUMN], row[TARGET_COLUMN])
    if graph is not None:
        data_list.append(graph)

print(f"Successfully converted {len(data_list)} molecules.")
if not data_list:
    raise ValueError("No valid molecules could be converted. Check SMILES strings and RDKit installation.")

num_node_features = data_list[0].num_node_features
print(f"Number of node features: {num_node_features}")

class MoleculeDataset(Dataset):
    def __init__(self, data_list):
        super().__init__()
        self.data_list = data_list

    def len(self):
        return len(self.data_list)

    def get(self, idx):
        return self.data_list[idx]

dataset = MoleculeDataset(data_list)

train_val_data, test_data = train_test_split(dataset, test_size=TEST_SIZE, random_state=RANDOM_STATE)
train_data, val_data = train_test_split(train_val_data, test_size=VALID_SIZE/(1-TEST_SIZE), random_state=RANDOM_STATE)

print(f"Train size: {len(train_data)}")
print(f"Validation size: {len(val_data)}")
print(f"Test size: {len(test_data)}")

train_loader = DataLoader(train_data, batch_size=BATCH_SIZE, shuffle=True)
val_loader = DataLoader(val_data, batch_size=BATCH_SIZE, shuffle=False)
test_loader = DataLoader(test_data, batch_size=BATCH_SIZE, shuffle=False)

class GNNModel(torch.nn.Module):
    def __init__(self, num_node_features, hidden_dim, output_dim=1, num_layers=3, dropout_rate=0.1):
        super(GNNModel, self).__init__()
        self.num_layers = num_layers
        self.dropout_rate = dropout_rate

        self.convs = torch.nn.ModuleList()
        self.batch_norms = torch.nn.ModuleList()

        self.convs.append(GCNConv(num_node_features, hidden_dim))
        self.batch_norms.append(BatchNorm1d(hidden_dim))

        for _ in range(num_layers - 1):
            self.convs.append(GCNConv(hidden_dim, hidden_dim))
            self.batch_norms.append(BatchNorm1d(hidden_dim))

        self.lin1 = Linear(hidden_dim, hidden_dim // 2)
        self.lin2 = Linear(hidden_dim // 2, output_dim)

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch

        for i in range(self.num_layers):
            x = self.convs[i](x, edge_index)
            x = self.batch_norms[i](x)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout_rate, training=self.training)

        x = global_mean_pool(x, batch)

        x = F.relu(self.lin1(x))
        x = F.dropout(x, p=self.dropout_rate, training=self.training)
        x = self.lin2(x)

        return x

def train(model, loader, optimizer, criterion, device):
    model.train()
    total_loss = 0
    for batch in loader:
        batch = batch.to(device)
        optimizer.zero_grad()
        out = model(batch)
        loss = criterion(out, batch.y.view_as(out))
        loss.backward()
        optimizer.step()
        total_loss += loss.item() * batch.num_graphs
    return total_loss / len(loader.dataset)

def evaluate(model, loader, device):
    model.eval()
    all_preds = []
    all_targets = []
    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device)
            out = model(batch)
            all_preds.append(out.cpu().numpy())
            all_targets.append(batch.y.cpu().numpy())

    all_preds = np.concatenate(all_preds).flatten()
    all_targets = np.concatenate(all_targets).flatten()

    r2 = r2_score(all_targets, all_preds)
    mae = mean_absolute_error(all_targets, all_preds)
    mse = mean_squared_error(all_targets, all_preds)
    rmse = np.sqrt(mse)

    return r2, mae, rmse, mse, all_targets, all_preds

model = GNNModel(num_node_features=num_node_features,
                   hidden_dim=HIDDEN_DIM,
                   output_dim=1,
                   num_layers=NUM_GNN_LAYERS,
                   dropout_rate=DROPOUT_RATE).to(DEVICE)

optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE, weight_decay=1e-5)
criterion = torch.nn.MSELoss()
scheduler = ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=10, verbose=True)

print("\n--- Starting Training ---")
best_val_rmse = float('inf')
patience_counter = 0
patience_limit = 25

train_losses = []
val_rmses = []

for epoch in range(1, NUM_EPOCHS + 1):
    train_loss = train(model, train_loader, optimizer, criterion, DEVICE)
    val_r2, val_mae, val_rmse, val_mse, _, _ = evaluate(model, val_loader, DEVICE)

    train_losses.append(train_loss)
    val_rmses.append(val_rmse)

    scheduler.step(val_rmse)

    print(f'Epoch: {epoch:03d}, Train Loss: {train_loss:.4f}, ',
          f'Val R2: {val_r2:.4f}, Val MAE: {val_mae:.4f}, Val RMSE: {val_rmse:.4f}')

    if val_rmse < best_val_rmse:
        best_val_rmse = val_rmse
        torch.save(model.state_dict(), 'best_model.pth')
        print(f'    New best validation RMSE: {best_val_rmse:.4f}. Model saved.')
        patience_counter = 0
    else:
        patience_counter += 1
        if patience_counter >= patience_limit:
            print(f'    Early stopping triggered after {patience_limit} epochs without improvement.')
            break

print("\n--- Training Finished ---")

print("\nLoading best model for final evaluation...")
try:
    model.load_state_dict(torch.load('best_model.pth', map_location=DEVICE))
    print("Best model loaded successfully.")
except FileNotFoundError:
    print("Warning: 'best_model.pth' not found. Evaluating with the last model state.")

print("\n--- Evaluating on Test Set ---")
test_r2, test_mae, test_rmse, test_mse, test_targets, test_preds = evaluate(model, test_loader, DEVICE)

print(f'Test R2:    {test_r2:.4f}')
print(f'Test MAE:   {test_mae:.4f}')
print(f'Test RMSE:  {test_rmse:.4f}')
print(f'Test MSE:   {test_mse:.4f}')

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(range(1, len(train_losses) + 1), train_losses, label='Train Loss')
plt.xlabel('Epoch')
plt.ylabel('MSE Loss')
plt.title('Training Loss')
plt.legend()
plt.grid(True)

plt.subplot(1, 2, 2)
plt.scatter(test_targets, test_preds, alpha=0.5, label=f'R² = {test_r2:.3f}')
plt.plot([min(test_targets), max(test_targets)], [min(test_targets), max(test_targets)], color='red', linestyle='--', label='Ideal')
plt.xlabel('True SF Value')
plt.ylabel('Predicted SF Value')
plt.title('Test Set: True vs. Predicted')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(6, 5))
plt.plot(range(1, len(val_rmses) + 1), val_rmses, label='Validation RMSE', color='orange')
plt.xlabel('Epoch')
plt.ylabel('RMSE')
plt.title('Validation RMSE')
plt.legend()
plt.grid(True)
plt.show()