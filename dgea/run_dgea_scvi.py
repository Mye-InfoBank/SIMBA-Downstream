from shiny import reactive, ui, render, module
import anndata as ad
import pandas as pd

from dgea.dgea_scvi import scanvi_dgea, get_normalized_counts


@module.ui
def run_dgea_ui():
    return ui.div(ui.output_ui("contrast_selector"),
                  ui.output_ui("subset_selector"),
                ui.input_task_button("run", "Run analysis"))

@module.server
def run_dgea_server(input, output, session,
                    _adata: reactive.Value[ad.AnnData],
                    _result: reactive.Value[pd.DataFrame],
                    _counts: reactive.Value[pd.DataFrame],
                    _reference: reactive.Value[str],
                    _alternative: reactive.Value[str],
                    _uniques: reactive.Value[list],
                    _contrast: reactive.Value[str],
                    _sub_category: reactive.Value[str],
                    _chosen_values: reactive.Value[list],
                    _sub_uniques: reactive.Value[list]):
    _category_columns = reactive.value([])
    _numeric_columns = reactive.value([])

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
    
    @output
    @render.ui
    def subset_selector():
        columns = _category_columns.get()

        return ui.input_select("sub_category", "Category to subset", choices=columns, selected=columns[0])
    
    @reactive.effect
    def update_contrast():
        contrast = input["contrast"].get()
        _contrast.set(contrast)
    
    @reactive.effect    
    def update_subset():
        sub_category = input["sub_category"].get()
        _sub_category.set(sub_category)

    @reactive.effect
    @reactive.event(input["run"])
    def handle_run():
        adata = _adata.get()
        contrast = _contrast.get()
        referece = _reference.get()
        alternative = _alternative.get()
        sub_category = _sub_category.get()
        chosen_values = _chosen_values.get()
        adata_subset = adata[adata.obs[sub_category].isin(chosen_values)].copy()
        run_scanvi(adata_subset, contrast, referece, alternative)

    @reactive.effect
    def update_scanvi_result():
        dge_change, counts = run_scanvi.result()
        _result.set(dge_change)
        _counts.set(counts)
    
    @ui.bind_task_button(button_id="run")
    @reactive.extended_task
    async def run_scanvi(adata, contrast, reference, alternative):
        dge_change = scanvi_dgea(adata, contrast, reference, alternative)
        
        counts = get_normalized_counts(adata)
        return dge_change, counts

    @reactive.effect
    def update_uniques():
        adata = _adata.get()
        contrast = _contrast.get()
        if adata is None or contrast is None:
            return
        uniques = adata.obs[contrast].unique().tolist()
        _uniques.set(uniques)
        
    @reactive.effect
    def update_sub_uniques():
        adata = _adata.get()
        sub_category = _sub_category.get()
        if adata is None or sub_category is None:
            return
        sub_uniques = adata.obs[sub_category].unique().tolist()
        _sub_uniques.set(sub_uniques)
