from shiny import reactive, ui, render, module
import numpy as np
import seaborn as sns
import shinywidgets as sw
import tempfile
import plotly.express as px

@module.ui
def plot_dgea_ui():
    return [
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
    ]

@module.server
def plot_dgea_server(input, output, session,
                     _filtered_counts,
                     _contrast,
                     _reference,
                     _alternative,
                     _design_matrix,
                     _result,
                     _alpha,
                     _lfc
                     ):
    _heatmap = reactive.value(None)

    @output
    @render.plot
    def plot_heatmap():
        counts_df = _filtered_counts.get()
        contrast = _contrast.get()
        reference = _reference.get()
        alternative = _alternative.get()
        design_matrix = _design_matrix.get()

        if counts_df is None \
            or contrast not in design_matrix.columns:
            return None

        design_matrix = design_matrix[design_matrix[contrast].isin([reference, alternative])]
        counts_df = counts_df.loc[:, design_matrix.index]

        if counts_df.empty:
            return None

        plot = sns.clustermap(counts_df.T, cmap="viridis", figsize=(10, 10))

        _heatmap.set(plot)

        return plot

    @render.download(
            filename=lambda: f"deseq2_matrix_{_contrast.get()}-{_reference.get()}:{_alternative.get()}.csv"
    )
    def download_dgea():
        deseq_results = _result.get()
        if deseq_results is None:
            return None
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as temp:
            deseq_results.to_csv(temp.name)
            return temp.name
        
    @render.download(
            filename=lambda: f"heatmap_{_contrast.get()}-{_reference.get()}:{_alternative.get()}.png"
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
        deseq_results = _result.get()
        alpha = _alpha.get()
        lfc = _lfc.get()

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