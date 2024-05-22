import numpy as np
import pandas as pd
import scanpy as sc
import random
import anndata as ad
import os

import torch
from scvi.model import SCANVI, SCVI
from scvi.data import poisson_gene_selection

def scanvi_dgea(adata:ad.AnnData, groupby:str, reference:str, alternative:str):
    
    adata = sc.read_h5ad("data/atlas.h5ad")  
    directory_model = "data/"
    model_path = os.path.join(directory_model, "model.pt")
    #weights_biases = torch.load(model_path, map_location=torch.device('cpu'))
    
    #poisson_gene_selection(adata)
    #adata.var.head()
    
    SCANVI.prepare_query_anndata(adata = adata, reference_model=directory_model)
    
    scanvi_model = SCANVI.load_query_data(adata, directory_model)
    print(type(scanvi_model))
    
    groups = np.array(adata.obs[groupby].unique())

    idx1 = adata.obs[groupby] == reference
    idx2 = adata.obs[groupby] == alternative
    
    dge_change = scanvi_model.differential_expression(adata=adata, groupby=groupby, idx1=idx1, idx2=idx2, mode="change")
    dge_change["log10_pscore"] = np.log10(dge_change["proba_not_de"])
    dge_change.head()

    return dge_change

adata = sc.read_h5ad("data/atlas.h5ad")  
print(adata.obs["cell_type"].value_counts())
test_dge = scanvi_dgea(adata, "cell_type", "Epithelial", "Endothelial")
print(test_dge['proba_de'].min(), test_dge['proba_de'].max())
print(test_dge['proba_not_de'].min(), test_dge['proba_not_de'].max())
print(test_dge['-log10_pscore'].min(), test_dge['-log10_pscore'].max())
print(test_dge['lfc_mean'].min(), test_dge['lfc_mean'].max())
print(test_dge.head(5))
print(test_dge.columns)