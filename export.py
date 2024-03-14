from shiny import reactive, ui, render, module
import os
import re

annotations_directory = "data/annotations"

def get_files(regex, directory=annotations_directory):
    return [os.path.join(directory, f) 
             for f in os.listdir(directory) 
             if regex.match(f)]

def get_gene_sets():
    regex = re.compile(r".*-gene-sets-[A-Z0-9]+.csv")
    return get_files(regex)

def get_cell_labels():
    regex = re.compile(r".*-cell-labels-[A-Z0-9]+.csv")
    return get_files(regex)


@module.ui
def export_ui():
    return ui.input_action_button("update", "Update"), \
        ui.navset_tab(
            ui.nav_panel("Gene sets",
                ui.output_ui("render_gene_sets")
            ),
            ui.nav_panel("Cell labels",
                ui.output_ui("render_cell_labels")
            )
        )


@module.server
def export_server(input, output, session):
    _gene_sets = reactive.value([])
    _cell_labels = reactive.value([])

    @reactive.effect
    @reactive.event(input["update"])
    def update_files():
        files = get_gene_sets()
        _gene_sets.set(files)

        files = get_cell_labels()
        _cell_labels.set(files)


    @render.ui
    def render_gene_sets():
        gene_sets = _gene_sets.get()

        return [ui.input_action_button(f"download_{os.path.basename(f).replace('-', '_')[:-4]}", f"Download {os.path.basename(f)}") for f in gene_sets]

    @render.ui
    def render_cell_labels():
        cell_labels = _cell_labels.get()

        return [ui.input_action_button(f"download_{os.path.basename(f).replace('-', '_')[:-4]}", f"Download {os.path.basename(f)}") for f in cell_labels]
