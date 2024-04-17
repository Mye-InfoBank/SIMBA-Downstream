from shiny import App, reactive, ui
import scanpy as sc
import os

from composition import composition_server, composition_ui
from export import export_ui, export_server
from dgea.dgea import dgea_server, dgea_ui

with open("data/input.txt") as f:
    file_list = f.readlines()
    filename = "data/" + file_list[0].strip()
    name = file_list[1].strip()


adata = sc.read_h5ad(filename)
categorical_columns = adata.obs.select_dtypes(include="category").columns.to_list()


app_ui = ui.page_navbar(
    ui.nav_panel("Composition analysis",
                 composition_ui("composition")
                 ),
    ui.nav_panel("Differential expression analysis",
                 dgea_ui("dgea")),
    ui.nav_panel("Data export",
                    export_ui("export")),
    title=f"SIMBA🦁 downstream: {name}",
    id="page"
)


def server(input, output, session):
    _dataframe = reactive.value(adata.obs)
    _adata = reactive.value(adata)
    composition_server("composition", _dataframe)
    export_server("export")
    dgea_server("dgea", _adata)

app = App(app_ui, server)