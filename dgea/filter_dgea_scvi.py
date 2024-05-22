from shiny import reactive, ui, render, module
import anndata as ad
import pandas as pd
from dgea.dgea_scvi import scanvi_dgea

@module.ui
def filter_dgea_ui():
    return ui.div(
        ui.output_ui("select_reference"),
        ui.output_ui("select_alternative"),
        ui.input_slider("log10_pscore", "Ropability in Reference (significance threshold)", min=0, max=1, step=0.01, value=0.8),
        ui.input_slider("lfc", "Log2 fold change", min=0, max=10, step=0.1, value=1),
        ui.output_ui("open_gprofiler")
    )

@module.server
def filter_dgea_server(input, output, session, 
                       _adata, _uniques,
                       _result, _filtered_result, _filtered_genes, 
                       _reference, _alternative, _contrast,
                       _log10_p, _lfc
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
        _log10_p.set(input["log10_pscore"].get())
        _lfc.set(input["lfc_mean"].get())

    @reactive.effect
    def update_result():
        adata = _adata.get()
        reference = _reference.get()
        alternative = _alternative.get()
        contrast = _contrast.get()

        if None in (adata, reference, alternative, contrast):
            return

        res_df = scanvi_dgea(adata, reference, alternative, contrast)
        _result.set(res_df)
    
    @reactive.effect
    def filter_result():
        result = _result.get()
        log10_p = input["log10_pscore"].get()
        lfc = input["lfc"].get()

        if result is None:
            return None

        result = result[(result["-log10_pscore"] < log10_p) & (result["lfc_mean"].abs() > lfc)]
        _filtered_result.set(result)
        genes = result.index.tolist()
        _filtered_genes.set(genes)

    @render.ui
    def open_gprofiler():
        genes = _filtered_genes.get()
        if not genes:
            return None
        
        return ui.input_action_button("gprofiler", label="Open g:Profiler",
                              onclick=f"window.open('https://biit.cs.ut.ee/gprofiler/gost?organism=hsapiens&query={'%0A'.join(genes)}', '_blank')")
