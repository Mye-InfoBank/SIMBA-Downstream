from shiny import reactive, ui, render, module
import anndata as ad
import pandas as pd
from dgea.deseq2 import get_results

@module.ui
def filter_dgea_ui():
    return ui.div(
        ui.output_ui("select_comparison"),
        ui.input_slider("alpha", "Alpha (significance threshold)", min=0, max=0.2, step=0.001, value=0.05),
        ui.input_slider("lfc", "Log2 fold change", min=0, max=10, step=0.1, value=1),
        ui.output_ui("open_gprofiler")
    )

@module.server
def filter_dgea_server(input, output, session, 
                       _dds, _counts, _comparisons,
                       _result, _filtered_result, _filtered_genes,
                       _filtered_counts, _comparison, _contrast,
                       _alpha, _lfc
                       ):

    @output
    @render.ui
    def select_comparison():
        comparisons = _comparisons.get()

        if not comparisons:
            return ui.p("Run analysis to see comparisons")

        return ui.input_select("comparison", "Comparison", choices=comparisons, selected=comparisons[0])

    @reactive.effect
    def update_filters():
        _comparison.set(input["comparison"].get())
        _alpha.set(input["alpha"].get())
        _lfc.set(input["lfc"].get())

    @reactive.effect
    def update_result():
        dds = _dds.get()
        comparison = input["comparison"].get()
        contrast = _contrast.get()

        if not comparison.startswith(contrast):
            return

        res_df = get_results(dds, comparison)
        _result.set(res_df)
    
    @reactive.effect
    def filter_result():
        result = _result.get()
        alpha = input["alpha"].get()
        lfc = input["lfc"].get()
        counts = _counts.get()

        if result is None:
            return None

        result = result[(result["padj"] < alpha) & (result["log2FoldChange"].abs() > lfc)]
        _filtered_result.set(result)
        genes = result.index.tolist()
        _filtered_genes.set(genes)
        _filtered_counts.set(counts.loc[genes, :])

    @render.ui
    def open_gprofiler():
        genes = _filtered_genes.get()
        if not genes:
            return None
        
        return ui.input_action_button("gprofiler", label="Significant genes - g:Profiler",
                    href=f"https://biit.cs.ut.ee/gprofiler/gost?organism=hsapiens&query={'%0A'.join(genes)}",
                    target="_blank")
