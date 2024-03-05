from shiny import App, reactive, ui, render, Session
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import io
import tempfile

with open("input.txt") as f:
    filename = f.read().strip()
    name = filename.replace(".h5ad", "")

adata = sc.read_h5ad(filename)
categorical_columns = adata.obs.select_dtypes(include="category").columns.to_list()


def plot(data: pd.DataFrame, reference: str, comparison: str):
    fractions = data.groupby([reference, comparison], observed=False).size().unstack()
    fractions = fractions.div(fractions.sum(axis=1), axis=0)

    fig, ax = plt.subplots()

    fractions.plot.barh(stacked=True, title=f"Dataset composition ({comparison.capitalize()})", color=plt.cm.tab20.colors, ax=ax)
    plt.legend(bbox_to_anchor=(1, 1.03), loc="upper left")

    return fig

app_ui = ui.page_navbar(
    ui.nav_panel("Composition analysis",
                 ui.layout_sidebar(
                     ui.sidebar(
                        ui.input_select("reference", "Reference", choices=categorical_columns, selected=categorical_columns[0]),
                        ui.input_select("comparison", "Comparison", choices=categorical_columns, selected=categorical_columns[1]),
                        ui.download_button("download", "Download plot")
                     ),
                     ui.output_plot("plot_composition")
                    )
                 ),
    ui.nav_panel("Differential expression analysis"),
    title=f"SIMBA downstream: {name}",
    id="page"
)


def server(input, output, session):
    _plot = reactive.value(None)

    @reactive.effect
    def update_plot():
        _plot.set(plot(adata.obs, input["reference"].get(), input["comparison"].get()))
    
    @output
    @render.plot
    def plot_composition():
        return _plot.get()
    
    @output
    @render.download(
        filename=lambda: f"composition_{input['reference'].get()}_{input['comparison'].get()}.png"
    )
    def download():
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as temp:
            _plot.get().savefig(temp.name, bbox_inches="tight")
            return temp.name


app = App(app_ui, server)