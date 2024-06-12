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
    

adata = sc.read_h5ad("data/atlas.h5ad")  
#print(adata.X.shape)
#print(adata.obs["cell_type"].value_counts())
test_dge = scanvi_dgea(adata, "cell_type", "Epithelial", "Endothelial")
#test_counts = get_normalized_counts(adata)
#print(test_dge['proba_de'].min(), test_dge['proba_de'].max())
print(test_dge['proba_not_de'].min(), test_dge['proba_not_de'].max())
print(test_dge['-log10_pscore'].min(), test_dge['-log10_pscore'].max())
print(test_dge['lfc_mean'].min(), test_dge['lfc_mean'].max())
#print(test_dge.head(5))
#print(test_dge.columns)
#print(test_counts.head(5))
#genes = list(set(test_dge.index.tolist()))
#print(len(genes))
#genes_not_found = [gene for gene in genes if gene not in test_counts.columns]
#print(len(genes_not_found))
#if genes_not_found:
    #print(f"Genes not found in the DataFrame: {genes_not_found}")
    #valid_genes = [gene for gene in genes if gene in test_counts.columns]
    #print(f"Valid genes: {valid_genes}")

#filtered = test_counts.loc[:, genes]
#print(filtered.shape)
