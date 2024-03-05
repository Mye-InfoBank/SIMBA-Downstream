from shiny import App, reactive, ui, render, Session
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

with open("input.txt") as f:
    filename = f.read().strip()
    name = filename.replace(".h5ad", "")

adata = sc.read_h5ad(filename)
categorical_columns = adata.obs.select_dtypes(include="category").columns

def plot(data: pd.DataFrame, reference: str, comparison: str):
    fractions = data.groupby([reference, comparison], observed=False).size().unstack()
    fractions = fractions.div(fractions.sum(axis=1), axis=0)

    fractions.plot.barh(stacked=True, title=f"Dataset composition ({comparison.capitalize()})")
    plt.legend(bbox_to_anchor=(1, 1.03), loc="upper left")

app_ui = ui.page_navbar(
    ui.nav_panel("Composition analysis"),  
    ui.nav_panel("Differential expression analysis"),  
    title=f"SIMBA downstream: {name}",  
    id="page",  
)  


def server(input, output, session):
    pass


app = App(app_ui, server)