import numpy as np
import pandas as pd
import scanpy as sc
import random
import anndata as ad
import os

import torch
from scvi.model import SCANVI, SCVI



def scanvi_dgea(adata:ad.AnnData, groupby:str):
    
    directory_model = "data/"
    model_path = os.path.join(directory_model, "model.pt")
    #weights_biases = torch.load(model_path, map_location=torch.device('cpu'))
    
    SCVI.setup_anndata(adata)
    scanvi_model = SCANVI.load(directory_model, adata=adata)
    print(type(scanvi_model))
    
    
    scanvi_model = SCANVI(adata)
    
    scanvi_model.load_state_dict(weights_biases)

    if not isinstance(scanvi_model, SCANVI):
        raise ValueError("Loaded model is not a trained scVI model instance.")
    
    #scanvi_model = SCANVI.from_scvi_model(model_to_load, adata)
    
    groups = np.array(adata.obs[groupby].unique())
    
    idx1= 0
    idx2 = 1
    
    dge_change = scanvi_model.differential_expression(adata=adata, groupby=groupby, idx1=idx1, idx2=idx2, mode="change")
    
    return dge_change

adata = sc.read_h5ad("data/atlas.h5ad")  
test_dge = scanvi_dgea(adata, "cell_type")
test_dge.head()
