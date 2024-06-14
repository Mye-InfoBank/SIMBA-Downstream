import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from scvi.model import SCANVI, SCVI

def scanvi_dgea(adata:ad.AnnData, groupby:str, reference:str, alternative:str, directory_model:str):
          
    if 'cell_type' in adata.obs.columns:
        model_type = SCANVI
        print('is scavi')
        
    else:
        model_type = SCVI
        print('is scanvi')
  
    model_type.prepare_query_anndata(adata = adata, reference_model=directory_model)
    
    model = model_type.load_query_data(adata, directory_model)
    
    groups = np.array(adata.obs[groupby].unique())

    idx1 = adata.obs[groupby] == reference
    idx2 = adata.obs[groupby] == alternative
    
    dge_change = model.differential_expression(adata=adata, groupby=groupby, idx1=idx1, idx2=idx2, mode="change")
    
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

if __name__ == '__main__':
    print('Running DGEA test')
    adata = sc.read_h5ad('/workspaces/SIMBA-Downstream_1/data/atlas.h5ad')
    dge_test = scanvi_dgea(adata, "cell_type", "Endothelial", "Epithelial", './data')
    print(dge_test.head())
