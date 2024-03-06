import numpy as np
import pandas as pd
import random
import anndata as ad

from rpy2.robjects import pandas2ri, Formula, r
from rpy2.robjects.packages import importr
pandas2ri.activate()
deseq = importr('DESeq2')

def pseudobulk(adata:ad.AnnData, groupby:str):
    groups = np.array(adata.obs[groupby].unique())
    groups = groups[~pd.isnull(groups)]
    reps = []
    obs_names = []
    group_dict = {}
    for group in groups:
        adata_group = adata[adata.obs[groupby]==group]
        index = adata_group.obs.index.tolist()
        random.shuffle(index)
        split_index = np.array_split(index,5)
        for n,split in enumerate(split_index):
            subset_group = adata_group[split,:]
            pb_rep = np.array(np.sum(subset_group.layers["counts"], axis=0)).ravel()
            reps.append(pb_rep)
            obs_name = f'{group}_bulk_{n}'
            group_dict[obs_name] = group
            obs_names.append(obs_name)
    pb_df = pd.DataFrame(reps).astype(np.int32)
    pb_df.columns = adata.var_names
    pb_df.index = obs_names
    return pb_df, pd.DataFrame.from_dict(group_dict, orient='index', columns=[groupby])

def get_formula(string:str):
    return Formula(string)

def get_dds(counts:pd.DataFrame, design:pd.DataFrame, formula:Formula):
    count_matrix = counts.T
    dds = deseq.DESeqDataSetFromMatrix(countData=count_matrix, 
                                        colData=design,
                                        design=formula)
    dds = deseq.DESeq(dds)
    return dds

def get_normalized_counts(dds):
    counts = r.counts(dds, normalized=True)
    counts = np.log1p(counts)
    return counts

def get_results(dds, comparison:str):
    res = deseq.results(dds, contrast=[comparison])
    res = r('function(x) data.frame(x)')(res)
    res_df = pandas2ri.rpy2py(res)
    return res_df

def get_comparisons(dds):
    comparisons = [str(comparison) for comparison in r.resultsNames(dds)]
    comparisons = [comparison for comparison in comparisons if comparison != "Intercept"]
    return comparisons