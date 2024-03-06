from shiny import reactive, ui, render, module
import anndata as ad
import pandas as pd

from dgea.deseq2 import pseudobulk, get_formula, get_dds, get_normalized_counts, get_comparisons


@module.ui
def run_dgea_ui():
    return ui.div(ui.output_ui("contrast_selector"),
              ui.input_task_button("run", "Run analysis"))

@module.server
def run_dgea_server(input, output, session,
                    _adata: reactive.Value[ad.AnnData],
                    _dds: reactive.Value,
                    _design_matrix: reactive.Value[pd.DataFrame],
                    _counts: reactive.Value[pd.DataFrame],
                    _comparisons: reactive.Value[list],
                    _contrast: reactive.Value[str]):
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
    
    @reactive.effect
    def update_contrast():
        contrast = input["contrast"].get()
        _contrast.set(contrast)

    @reactive.effect
    @reactive.event(input["run"])
    def handle_run():
        print("Handle")
        adata = _adata.get()
        contrast = _contrast.get()
        run_deseq(adata, contrast)

    @reactive.effect
    def update_deseq_result():
        dds, df_counts, design = run_deseq.result()
        _dds.set(dds)
        _design_matrix.set(design)
        _counts.set(df_counts)
    
    @ui.bind_task_button(button_id="run")
    @reactive.extended_task
    async def run_deseq(adata, contrast):
        df, design = pseudobulk(adata, contrast)

        design_formula = get_formula(f'~ {contrast}')

        dds = get_dds(df, design, design_formula)

        counts = get_normalized_counts(dds)
        df_counts = pd.DataFrame(counts, index=df.columns, columns=df.index)

        return dds, df_counts, design

    @reactive.effect
    def update_comparisons():
        dds = _dds.get()

        if dds is None:
            return

        _comparisons.set(get_comparisons(dds))