import numpy as np
import pandas as pd
import scanpy as sc
import random
import anndata as ad
import os

import torch
from scvi.model import SCANVI

directory_model = "data/"
model_path = os.path.join(directory_model, "model.pt")
model_to_load = torch.load(model_path, map_location=torch.device('cpu'))
#print(model_to_load) 
print(type(model_to_load)) 

#print(hasattr(model_to_load, 'genes')) 
#print(hasattr(model_to_load, 'decoder'))


def scanvi_dgea(adata:ad.AnnData, groupby:str):
    
    if not isinstance(model_to_load, SCANVI):
        raise ValueError("Loaded model is not a trained scVI model instance.")
    
    labels_key = "labels"
    unlabeled_category = "unlabeled"
    SCANVI.setup_anndata(adata, labels_key=labels_key, unlabeled_category=unlabeled_category)
    
    scanvi_model = SCANVI.from_scvi_model(model_to_load, adata)
    
    groups = np.array(adata.obs[groupby].unique())
    
    idx1= 0
    idx2 = 1
    
    dge_change = scanvi_model.differential_expression(adata=adata, groupby=groupby, idx1=idx1, idx2=idx2, mode="change")
    
    return dge_change

adata = sc.read_h5ad("data/atlas.h5ad")  
test_dge = scanvi_dgea(adata, "cell_type")
test_dge.head()
