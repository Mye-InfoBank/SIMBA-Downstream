from shiny import reactive, ui, render, module
import anndata as ad

from dgea.run_dgea import run_dgea_ui, run_dgea_server
from dgea.filter_dgea import filter_dgea_ui, filter_dgea_server
from dgea.plot_dgea import plot_dgea_ui, plot_dgea_server

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
    _dds = reactive.value(None)
    _comparisons = reactive.value([])
    _counts = reactive.value(None)
    _design_matrix = reactive.value(None)

    _contrast = reactive.value(None)
    _comparison = reactive.value(None)
    _alpha = reactive.value(0.05)
    _lfc = reactive.value(1)

    _result = reactive.value(None)
    _filtered_result = reactive.value(None)
    _filtered_genes = reactive.value(None)
    _filtered_counts = reactive.value(None)

    run_dgea_server("run_dgea", _adata, _dds, _design_matrix, _counts, _comparisons, _contrast)
    filter_dgea_server("filter_dgea", _dds, _counts, _comparisons,
                       _result, _filtered_result, _filtered_genes,
                       _filtered_counts, _comparison, _contrast, _alpha, _lfc)
    plot_dgea_server("plot_dgea", _filtered_counts, _contrast, _comparison, _design_matrix,
                     _result, _alpha, _lfc)
