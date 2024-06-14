import scvi
import scanpy as sc
from scvi.model import SCVI

# Load the dataset
adata = sc.read_h5ad('/workspaces/SIMBA-Downstream_1/data/atlas.h5ad')

# Setup the anndata for scVI
SCVI.setup_anndata(adata)

# Initialize the SCVI model
model = SCVI(adata)

# Train the SCVI model
model.train()

# Save the trained model to a file
model.save("./model_scvi.pt", save_anndata=True)