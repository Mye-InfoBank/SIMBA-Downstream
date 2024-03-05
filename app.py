from shiny import App, reactive, ui, render
import scanpy as sc

from composition import composition_server, composition_ui

with open("input.txt") as f:
    filename = f.read().strip()
    name = filename.replace(".h5ad", "")

adata = sc.read_h5ad(filename)
categorical_columns = adata.obs.select_dtypes(include="category").columns.to_list()


app_ui = ui.page_navbar(
    ui.nav_panel("Composition analysis",
                 composition_ui("composition")
                 ),
    ui.nav_panel("Differential expression analysis"),
    title=f"SIMBA🦁 downstream: {name}",
    id="page"
)


def server(input, output, session):
    _dataframe = reactive.value(adata.obs)
    composition_server("composition", _dataframe)


app = App(app_ui, server)