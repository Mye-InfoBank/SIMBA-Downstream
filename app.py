from shiny import App, reactive, ui
import scanpy as sc
import pickle
import json

from composition import composition_server, composition_ui
from export import export_ui, export_server
from dgea.dgea_scvi_define import dgea_server, dgea_ui
#from dgea.dgea import dgea_server, dgea_ui
from tree import tree_server, tree_ui

with open("data/config.json") as f:
    config = json.load(f)
    adata = sc.read_h5ad("data/" + config["adata"])
    tree = pickle.load(open("data/" + config["tree"], "rb")) if "tree" in config else None
    name = config["name"]

categorical_columns = adata.obs.select_dtypes(include="category").columns.to_list()

app_ui = ui.page_navbar(
    ui.nav_panel("Composition analysis",
                 composition_ui("composition")
                 ),
    ui.nav_panel("Differential expression analysis",
                 dgea_ui("dgea")),
    ui.nav_panel("Tree", tree_ui("tree")),
    ui.nav_panel("Data export",
                    export_ui("export")),
    title=f"SIMBA🦁 downstream: {name}",
    id="page"
)


def server(input, output, session):
    _dataframe = reactive.value(adata.obs)
    _adata = reactive.value(adata)
    _tree = reactive.value(tree)
    composition_server("composition", _dataframe)
    export_server("export")
    dgea_server("dgea", _adata)
    tree_server("tree", _tree)

app = App(app_ui, server)