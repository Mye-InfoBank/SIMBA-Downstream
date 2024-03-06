from shiny import reactive, ui, render, module
import matplotlib.pyplot as plt
import pandas as pd
import anndata as ad
import tempfile

@module.ui
def dgea_ui():
    return ui.layout_sidebar(
                     ui.sidebar(
                        ui.output_ui("selectors"),
                        #ui.download_button("download", "Download plot")
                        title="Select covariates"
                     ),
                     #ui.output_plot("plot_dgea")
                    )

@module.server
def dgea_server(input, output, session, _adata: reactive.Value[ad.AnnData]):
    _category_columns = reactive.value([])
    _numeric_columns = reactive.value([])

    @reactive.effect
    def update_columns():
        _category_columns.set(_adata.get().obs.select_dtypes(include=["category", "object"]).columns.to_list())
        _numeric_columns.set(_adata.get().obs.select_dtypes(include="number").columns.to_list())

    @output
    @render.ui
    def selectors():
        category_columns = _category_columns.get()
        numeric_columns = _numeric_columns.get()

        sel_category = ui.input_selectize("categorical", "Categorical", choices=category_columns, selected=category_columns[:2], multiple=True)

        if len(numeric_columns) == 0:
            return sel_category
        
        sel_numeric = ui.input_selectize("numeric", "Numeric", choices=numeric_columns, selected=numeric_columns[:2], multiple=True)

        return sel_category, sel_numeric