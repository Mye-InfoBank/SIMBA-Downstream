from shiny import reactive, ui, render, module
import pandas as pd
import anndata as ad
import tempfile
import numpy as np
import random
import seaborn as sns
import shinywidgets as sw
import plotly.express as px

from rpy2.robjects import pandas2ri, Formula, r
from rpy2.robjects.packages import importr
pandas2ri.activate()
deseq = importr('DESeq2')

@module.ui
def dgea_ui():
    return ui.layout_sidebar(
                     ui.sidebar(
                        ui.output_ui("contrast_selector"),
                        ui.input_action_button("run", "Run analysis"),
                        title="Select covariates"
                     ),
                     ui.card(
                         ui.card_header("Plot settings"),
                         ui.output_ui("select_comparison"),
                         ui.input_slider("alpha", "Alpha", min=0, max=0.2, step=0.001, value=0.05),
                         ui.input_slider("lfc", "Log2 fold change", min=0, max=10, step=0.1, value=1)
                     ),
                     ui.card(
                         ui.card_header("Volcano plot"),
                         sw.output_widget("plot_volcano")
                     ),
                     ui.card(
                         ui.card_header("Heatmap"),
                         ui.output_plot("plot_heatmap"),
                         ui.card_footer(
                             ui.download_button("download_dgea", "Download DESeq2 matrix"),
                             ui.download_button("download_plot", "Download plot")
                         )
                     )
                    )

@module.server
def dgea_server(input, output, session, _adata: reactive.Value[ad.AnnData]):
    _category_columns = reactive.value([])
    _numeric_columns = reactive.value([])
    _dds = reactive.value(None)
    _comparisons = reactive.value([])
    _genes_significant = reactive.value([])
    _counts = reactive.value(None)
    _significant_counts = reactive.value(None)
    _design_matrix = reactive.value(None)
    _deseq_results = reactive.value(None)
    _heatmap = reactive.value(None)

    @reactive.effect
    def update_columns():
        obs = _adata.get().obs
        obs = obs.loc[:, obs.nunique() > 1]

        _category_columns.set(obs.select_dtypes(include=["category", "object"]).columns.to_list())
        _numeric_columns.set(obs.select_dtypes(include="number").columns.to_list())
    
    @output
    @render.ui
    def contrast_selector():
        columns = _category_columns.get()

        return ui.input_select("contrast", "Contrast", choices=columns, selected=columns[0])
    
    @reactive.effect
    @reactive.event(input["run"])
    def run_deseq():
        adata = _adata.get()
        contrast = input["contrast"].get()

        df, design = pseudobulk(adata, contrast)
        design_formula = Formula(f'~ {contrast}')

        _design_matrix.set(design)

        count_matrix = df.T

        dds = deseq.DESeqDataSetFromMatrix(countData=count_matrix, 
                                        colData=design,
                                        design=design_formula)
        
        dds = deseq.DESeq(dds)
        _dds.set(dds)

        counts = r.counts(dds, normalized=True)
        counts = np.log1p(counts)

        df_counts = pd.DataFrame(counts, index=df.columns, columns=df.index)
        _counts.set(df_counts)
    
    @reactive.effect
    def update_comparisons():
        dds = _dds.get()

        if dds is None:
            return

        comparisons = [str(comparison) for comparison in r.resultsNames(dds)]
        comparisons = [comparison for comparison in comparisons if comparison != "Intercept"]
        _comparisons.set(comparisons)
    
    @output
    @render.ui
    def select_comparison():
        comparisons = _comparisons.get()

        if not comparisons:
            return ui.p("Run analysis to see comparisons")

        return ui.input_select("comparison", "Comparison", choices=comparisons, selected=comparisons[0])
    
    @reactive.effect
    def update_significant():
        dds = _dds.get()
        alpha = input["alpha"].get()
        lfc = input["lfc"].get()
        comparison = input["comparison"].get()
        contrast = input["contrast"].get()

        if dds is None:
            return None

        if not comparison.startswith(contrast):
            return None
        
        res = deseq.results(dds, contrast=[comparison])
        res = r('function(x) data.frame(x)')(res)
        res_df = pandas2ri.rpy2py(res)

        _deseq_results.set(res_df)

        res_df = res_df[res_df["padj"] < alpha]
        res_df = res_df[np.abs(res_df["log2FoldChange"]) > lfc]

        genes = res_df.index.tolist()
        _genes_significant.set(genes)

    @reactive.effect
    def update_significant_counts():
        counts_df = _counts.get()
        genes = _genes_significant.get()

        if counts_df is None or not genes:
            return None
        
        counts_df = counts_df.loc[genes]
        _significant_counts.set(counts_df)

    @output
    @render.plot
    def plot_heatmap():
        counts_df = _significant_counts.get()
        contrast = input["contrast"].get()
        comparison = input["comparison"].get()
        design_matrix = _design_matrix.get()

        if counts_df is None \
            or contrast not in design_matrix.columns:
            return None
        
        comparison = comparison[len(contrast)+1:]
        groups = comparison.split("_vs_")
        
        design_matrix = design_matrix[design_matrix[contrast].isin(groups)]
        counts_df = counts_df.loc[:, design_matrix.index]

        if counts_df.empty:
            return None

        plot = sns.clustermap(counts_df.T, cmap="viridis", figsize=(10, 10))

        _heatmap.set(plot)

        return plot
    
    @render.download(
            filename=lambda: f"deseq2_matrix_{input['comparison'].get()}.csv"
    )
    def download_dgea():
        deseq_results = _deseq_results.get()
        if deseq_results is None:
            return None
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as temp:
            deseq_results.to_csv(temp.name)
            return temp.name
        
    @render.download(
            filename=lambda: f"heatmap_{input['comparison'].get()}.png"
    )
    def download_plot():
        plot = _heatmap.get()
        if plot is None:
            return None
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as temp:
            plot.savefig(temp.name, bbox_inches="tight")
            return temp.name

    @output
    @sw.render_plotly
    def plot_volcano():
        deseq_results = _deseq_results.get()
        alpha = input["alpha"].get()
        lfc = input["lfc"].get()

        if deseq_results is None:
            return None
        
        df_plot = deseq_results.copy()
        df_plot["-log10(padj)"] = -np.log10(df_plot["padj"])
        df_plot["category"] = "Not significant"
        df_plot["gene"] = df_plot.index
        df_plot.loc[(df_plot["log2FoldChange"] < -lfc) & (df_plot["padj"] < alpha), "category"] = "Downregulated"
        df_plot.loc[(df_plot["log2FoldChange"] > lfc) & (df_plot["padj"] < alpha), "category"] = "Upregulated"

        colormap = {"Not significant": "grey", "Downregulated": "blue", "Upregulated": "red"}
        
        return px.scatter(df_plot, x="log2FoldChange", y="-log10(padj)", color="category",
                          color_discrete_map=colormap,
                          hover_name="gene",
                          labels={"log2FoldChange": "Log2 fold change",
                                  "-log10(padj)": "-log10(padj)",
                                  "category": "Category"})


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