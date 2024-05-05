from shiny import reactive, ui, render, module
import scHPL
from typing import List
from shiny.render.renderer import Jsonifiable, Renderer
from htmltools import HTMLDependency

tree_dep = HTMLDependency(
    "tree",
    "1.0.0",
    source={"subdir": "tree"},
    script={"src": "tree.js", "type": "module"},
    all_files=True
)

def output_tree(id):
    return ui.div(
        ui.HTML("<script src='//unpkg.com/3d-force-graph'></script>"),
        tree_dep,
        id=module.resolve_id(id),
        class_="shiny-tree-output"
    )

class render_tree(Renderer[dict]):
    def auto_output_ui(self):
        return ui.ouptut_tree(self.ouptut_name)
    
    async def transform(self, value: dict):
        return value

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

@module.ui
def tree_ui():
    return output_tree("tree")

@module.server
def tree_server(input, output, session, _tree: reactive.Value[scHPL.TreeNode]):
    @output
    @render_tree
    def tree():
        return format(_tree.get()[0])