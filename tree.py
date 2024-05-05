from shiny import reactive, ui, render, module
import scHPL
import json
from typing import List

@module.ui
def tree_ui():
    return ui.output_text("tree")

def format(node: scHPL.TreeNode):
    nodes = []
    links = []

    def format_name(name: List[str]):
        return " & ".join(name)

    for current in node.walk():
        cur_name = format_name(current.name)
        nodes.append({"name": cur_name, "id": hash(cur_name)})
        for child in current.descendants:
            links.append({"source": hash(cur_name), "target": hash(format_name(child.name))})

    return {"nodes": nodes, "links": links}

@module.server
def tree_server(input, output, session, _tree: reactive.Value[scHPL.TreeNode]):
    @output
    @render.text
    def tree():
        data = format(_tree.get()[0])
        return json.dumps(data, indent=4)