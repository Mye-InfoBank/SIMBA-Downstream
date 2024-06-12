import numpy as np
import pandas as pd
import scanpy as sc
import random
import anndata as ad
import os

import torch
from scvi.model import SCANVI, SCVI

def scanvi_dgea(adata:ad.AnnData, groupby:str, reference:str, alternative:str):
     
    directory_model = "data/"
    model_path = os.path.join(directory_model, "model.pt")
    print(reference, alternative)
    
    SCANVI.prepare_query_anndata(adata = adata, reference_model=directory_model)
    
    scanvi_model = SCANVI.load_query_data(adata, directory_model)
    print(type(scanvi_model))
    
    groups = np.array(adata.obs[groupby].unique())

    idx1 = adata.obs[groupby] == reference
    idx2 = adata.obs[groupby] == alternative
    
    dge_change = scanvi_model.differential_expression(adata=adata, groupby=groupby, idx1=idx1, idx2=idx2, mode="change")
    
    epsilon = 1e-10
    dge_change['proba_not_de'] = np.maximum(dge_change["proba_not_de"], epsilon)
    dge_change["log10_pscore"] = np.log10(dge_change["proba_not_de"])
    dge_change["-log10_pscore"] = -np.log10(dge_change["proba_not_de"])	
    
    return dge_change

def get_normalized_counts(adata):
    print(adata.shape)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["counts"] = adata.X.copy().tocsr()
    counts = adata.layers["counts"]
    dense_matrix = counts.toarray()
    df_counts = pd.DataFrame(dense_matrix, index=adata.obs_names, columns=adata.var_names)
    return df_counts
    