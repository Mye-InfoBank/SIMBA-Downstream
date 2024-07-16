from shiny import reactive, ui, render, module
import os
import re

annotations_directory = "data/annotations"


def get_files(regex, directory=annotations_directory):
    if not os.path.exists(directory):
        return []
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
        ui.layout_columns(
            ui.card(
                ui.card_header("Gene sets"),
                ui.output_ui("render_gene_sets"),
                ui.card_footer(
                    ui.download_button("download_gene_set", "Download")
                )
            ),
            ui.card(
                ui.card_header("Cell labels"),
                ui.output_ui("render_cell_labels"),
                ui.card_footer(
                    ui.download_button("download_cell_labels", "Download")
                )
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

        return ui.input_select("gene_set", "Select file", choices={f: os.path.basename(f) for f in gene_sets})

    @render.ui
    def render_cell_labels():
        cell_labels = _cell_labels.get()

        return ui.input_select("cell_labels", "Select file", choices={f: os.path.basename(f) for f in cell_labels})

    @render.download(
        filename=lambda: os.path.basename(input["gene_set"].get())
    )
    def download_gene_set():
        file_path = input["gene_set"].get()
        print(file_path)
        return file_path

    @render.download(
        filename=lambda: os.path.basename(input["cell_labels"].get())
    )
    def download_cell_labels():
        file_path = input["cell_labels"].get()
        return file_path
