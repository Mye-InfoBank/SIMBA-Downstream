from shiny import reactive, ui, render, module
import anndata as ad
import pandas as pd
from dgea.deseq2 import get_results

@module.ui
def filter_dgea_ui():
    return ui.div(
        ui.output_ui("select_reference"),
        ui.output_ui("select_alternative"),
        ui.input_slider("alpha", "Alpha (significance threshold)", min=0, max=0.2, step=0.001, value=0.05),
        ui.input_slider("lfc", "Log2 fold change", min=0, max=10, step=0.1, value=1),
        ui.output_ui("open_gprofiler")
    )

@module.server
def filter_dgea_server(input, output, session, 
                       _dds, _counts, _uniques,
                       _result, _filtered_result, _filtered_genes,
                       _filtered_counts, _reference, _alternative, _contrast,
                       _alpha, _lfc
                       ):

    @output
    @render.ui
    def select_reference():
        uniques = _uniques.get()

        if not uniques or len(uniques) < 2:
            return ui.p("Run analysis to see options")

        return ui.input_select("reference", "Reference", choices=uniques, selected=uniques[0])
    
    @output
    @render.ui
    def select_alternative():
        uniques = _uniques.get()

        if not uniques or len(uniques) < 2:
            return ui.p("Run analysis to see options")
        
        print(uniques)

        return ui.input_select("alternative", "Alternative", choices=uniques, selected=uniques[1])


    @reactive.effect
    def update_filters():
        _reference.set(input["reference"].get())
        _alternative.set(input["alternative"].get())
        _alpha.set(input["alpha"].get())
        _lfc.set(input["lfc"].get())

    @reactive.effect
    def update_result():
        dds = _dds.get()
        reference = _reference.get()
        alternative = _alternative.get()
        contrast = _contrast.get()

        if None in (dds, reference, alternative, contrast):
            return

        res_df = get_results(dds, contrast, reference, alternative)
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
        
        return ui.input_action_button("gprofiler", label="Open g:Profiler",
                              onclick=f"window.open('https://biit.cs.ut.ee/gprofiler/gost?organism=hsapiens&query={'%0A'.join(genes)}', '_blank')")
