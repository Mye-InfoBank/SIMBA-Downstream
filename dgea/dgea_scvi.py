import numpy as np
import pandas as pd
import scanpy as sc
import random
import anndata as ad
import os

import torch
from scvi.model import SCANVI, SCVI

adata = sc.read_h5ad("data/atlas.h5ad")  
directory_model = "data/"
model_path = os.path.join(directory_model, "model.pt")
#weights_biases = torch.load(model_path, map_location=torch.device('cpu'))
    
SCANVI.prepare_query_anndata(adata = adata, reference_model=directory_model)
#SCANVI.setup_anndata(adata, labels_key="cell_type", batch_key="batch", layer="counts", unlabeled_category="Unknown")
scanvi_model = SCANVI.load_query_data(adata, directory_model)
print(type(scanvi_model))

groupby = "cell_type"
groups = np.array(adata.obs[groupby].unique())

idx1= groups[0]
idx2 = groups[1]

dge_change = scanvi_model.differential_expression(adata=adata, groupby=groupby, idx1=reference, idx2=alternative, mode="change")
dge_change.head()

if not isinstance(scanvi_model, SCANVI):
    raise ValueError("Loaded model is not a trained scVI model instance.")

def scanvi_dgea(adata:ad.AnnData, groupby:str):
    
    directory_model = "data/"
    model_path = os.path.join(directory_model, "model.pt")
    #weights_biases = torch.load(model_path, map_location=torch.device('cpu'))
    
    SCANVI.prepare_query_anndata(adata = adata, reference_model=directory_model)
    #SCANVI.setup_anndata(adata, labels_key="cell_type", batch_key="batch", layer="counts", unlabeled_category="Unknown")
    #scanvi_model = SCANVI.load(directory_model, adata=adata)
    #print(type(scanvi_model))
    
    
    scanvi_model = SCANVI(adata)
    
    scanvi_model.load_state_dict(weights_biases)

    if not isinstance(scanvi_model, SCANVI):
        raise ValueError("Loaded model is not a trained scVI model instance.")
    
    #scanvi_model = SCANVI.from_scvi_model(model_to_load, adata)
    
    groups = np.array(adata.obs[groupby].unique())
    
    idx1= groups[0]
    idx2 = groups[1]
    
    dge_change = scanvi_model.differential_expression(adata=adata, groupby=groupby, idx1=idx1, idx2=idx2, mode="change")
    
    return dge_change

adata = sc.read_h5ad("data/atlas.h5ad")  
#test_dge = scanvi_dgea(adata, "cell_type")
#test_dge.head()
