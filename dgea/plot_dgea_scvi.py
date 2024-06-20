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
                ui.download_button("download_dgea", "Download DGEA matrix"),
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
                     _result,
                     _log10_p,
                     _lfc
                     ):
    _heatmap = reactive.value(None)

    @output
    @render.plot
    def plot_heatmap():
        counts_df = _filtered_counts.get()

        if counts_df is None or counts_df.empty:
            return None
        
        plot = sns.clustermap(counts_df.T, cmap="viridis", figsize=(10, 10))
        _heatmap.set(plot)

        return plot

    @render.download(
            filename=lambda: f"dgea_matrix_{_contrast.get()}-{_reference.get()}:{_alternative.get()}.csv"
    )
    def download_dgea():
        scanvi_results = _result.get()
        if scanvi_results is None:
            return None
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as temp:
            scanvi_results.to_csv(temp.name)
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
        scanvi_results = _result.get()
        log10_p = _log10_p.get()
        lfc = _lfc.get()

        if scanvi_results is None:
            return None
        
        df_plot = scanvi_results.copy()
        #df_plot["-log10_pscore"] = -np.log10(df_plot["log10_pscore"])
        min_value = 1e-100
        #df_plot["p-value"] = df_plot["padj"].apply(lambda x: x if x > min_value else f"less than {min_value}")

        df_plot["category"] = "Not significant"
        df_plot["gene"] = df_plot.index
        df_plot.loc[(df_plot["lfc_mean"] < -lfc) & (df_plot["-log10_pscore"] < log10_p), "category"] = f"High in {_alternative.get()}"
        df_plot.loc[(df_plot["lfc_mean"] > lfc) & (df_plot["-log10_pscore"] < log10_p), "category"] = f"High in {_reference.get()}"

        colormap = {"Not significant": "grey", f"High in {_alternative.get()}": "blue", f"High in {_reference.get()}": "red"}
        
        hover_data = {
            "lfc_mean": True,
            "-log10_pscore": True,
            "category": True,
        }
        
        return px.scatter(df_plot, x="lfc_mean", y="-log10_pscore", color="category",
                          color_discrete_map=colormap,
                          hover_name="gene",
                          hover_data=hover_data,
                          labels={"lfc_mean": "Log2 fold change mean",
                                  "-log10_pscore": "Negative log 10 P-value",
                                  "category": "Category"})