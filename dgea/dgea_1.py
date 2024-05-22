from shiny import reactive, ui, render, module
import anndata as ad

from dgea.run_dgea_scvi import run_dgea_ui, run_dgea_server
from dgea.filter_dgea_scvi import filter_dgea_ui, filter_dgea_server
from dgea.plot_dgea_scvi import plot_dgea_ui, plot_dgea_server

@module.ui
def dgea_ui():
    return ui.layout_sidebar(
                     ui.sidebar(
                        run_dgea_ui("run_dgea"),
                        filter_dgea_ui("filter_dgea"),
                        title="Select covariates"
                     ),
                     *plot_dgea_ui("plot_dgea")
                    )

@module.server
def dgea_server(input, output, session, _adata: reactive.Value[ad.AnnData]):
    _uniques = reactive.value([])

    _contrast = reactive.value(None)
    _reference = reactive.value(None)
    _alternative = reactive.value(None)
    _alpha = reactive.value(0.05)
    _lfc = reactive.value(1)

    _result = reactive.value(None)
    _filtered_result = reactive.value(None)
    _filtered_genes = reactive.value(None)

    run_dgea_server("run_dgea", _adata, _result, _reference, _alternative, _uniques, _contrast)
    filter_dgea_server("filter_dgea", _uniques,
                       _result, _filtered_result, _filtered_genes,
                        _reference, _alternative, _contrast, _alpha, _lfc)
    plot_dgea_server("plot_dgea", _contrast, _reference, _alternative,
                     _result, _alpha, _lfc)
