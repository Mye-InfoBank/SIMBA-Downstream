from shiny import reactive, ui, render, module

@module.ui
def tree_ui():
    return ui.output_text("tree")

@module.server
def tree_server(input, output, session):
    @output
    @render.text
    def tree():
        return "This is a tree plot"