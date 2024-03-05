from shiny import reactive, ui, render, module
import matplotlib.pyplot as plt
import pandas as pd
import tempfile

def plot(data: pd.DataFrame, reference: str, comparison: str):
    fractions = data.groupby([reference, comparison], observed=False).size().unstack()
    fractions = fractions.div(fractions.sum(axis=1), axis=0)

    fig, ax = plt.subplots()

    fractions.plot.barh(stacked=True, title=f"Dataset composition ({comparison.capitalize()})", color=plt.cm.tab20.colors, ax=ax)
    plt.legend(bbox_to_anchor=(1, 1.03), loc="upper left")

    return fig

@module.ui
def composition_ui():
    return ui.layout_sidebar(
                     ui.sidebar(
                        ui.output_ui("selectors"),
                        ui.download_button("download", "Download plot")
                     ),
                     ui.output_plot("plot_composition")
                    )

@module.server
def composition_server(input, output, session, _dataframe):
    _category_columns = reactive.value([])
    _plot = reactive.value(None)

    @reactive.effect
    def update_category_columns():
        _category_columns.set(_dataframe.get().select_dtypes(include="category").columns.to_list())

    @output
    @render.ui
    def selectors():
        category_columns = _category_columns.get()

        return ui.input_select("reference", "Reference", choices=category_columns, selected=category_columns[0]), \
                        ui.input_select("comparison", "Comparison", choices=category_columns, selected=category_columns[1])

    @reactive.effect
    def update_plot():
        _plot.set(plot(_dataframe.get(), input["reference"].get(), input["comparison"].get()))
    
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