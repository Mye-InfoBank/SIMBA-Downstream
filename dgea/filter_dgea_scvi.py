from shiny import reactive, ui, render, module
import anndata as ad
import pandas as pd
from dgea.dgea_scvi import scanvi_dgea, get_normalized_counts

@module.ui
def filter_dgea_ui():
    return ui.div(
        ui.output_ui("select_values"),
        ui.output_ui("select_reference"),
        ui.output_ui("select_alternative"),
        ui.input_slider("log10_pscore", "Ropability in Reference (significance threshold)", min=0, max=20, step=0.01, value=3),
        ui.input_slider("lfc", "Log2 fold change", min=0, max=10, step=0.1, value=1),
        ui.output_ui("open_gprofiler")
    )

@module.server
def filter_dgea_server(input, output, session, 
                       _adata, _counts, _uniques,
                       _result, _filtered_result, _filtered_genes, _filtered_counts,
                       _reference, _alternative, _contrast, _sub_category, _sub_uniques, _chosen_values,
                       _log10_p, _lfc
                       ):

    @output
    @render.ui
    def select_values():
        sub_uniques = _sub_uniques.get()

        return ui.input_select("value", "Choose values:", choices=sub_uniques, selectize=True, multiple=True, selected=sub_uniques)
    
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
        _chosen_values.set(input["value"].get())
        _reference.set(input["reference"].get())
        _alternative.set(input["alternative"].get())
        _log10_p.set(input["log10_pscore"].get())
        _lfc.set(input["lfc"].get())

    '''
    @reactive.effect
    def update_result():
        adata = _adata.get()
        reference = _reference.get()
        alternative = _alternative.get()
        contrast = _contrast.get()
        sub_category = _sub_category.get()
        chosen_values = _chosen_values.get()

        if None in (reference, alternative, contrast, sub_category, chosen_values):
            return

        res_df = scanvi_dgea(adata, contrast, reference, alternative)
        res_counts = get_normalized_counts(adata)
        _result.set(res_df)
        #_counts.set(res_counts)
        adata_sub = adata.obs[sub_category].isin(chosen_values)
        print(f"adata_sub: {adata_sub}")
        print(f"adata_sub index: {adata_sub.index}")
        filtered_values_counts = res_counts.loc[adata_sub.index[adata_sub]]

        _counts.set(filtered_values_counts)
    '''
    @reactive.effect
    def filter_result():
        result = _result.get()
        log10_p = input["log10_pscore"].get()
        lfc = input["lfc"].get()
        counts = _counts.get()
        
        if result is None:
            return None

        result = result[(result["-log10_pscore"] < log10_p) & (result["lfc_mean"].abs() > lfc)]
        _filtered_result.set(result)
        genes = result.index.tolist()
        genes_not_found = [gene for gene in genes if gene not in counts.columns]
        if genes_not_found:
            print(f"Genes not found in the DataFrame: {genes_not_found}")
        _filtered_genes.set(genes)
        _filtered_counts.set(counts.loc[:, genes])

    @render.ui
    def open_gprofiler():
        genes = _filtered_genes.get()
        if not genes:
            return None
        
        return ui.input_action_button("gprofiler", label="Open g:Profiler",
                              onclick=f"window.open('https://biit.cs.ut.ee/gprofiler/gost?organism=hsapiens&query={'%0A'.join(genes)}', '_blank')")
