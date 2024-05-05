if (Shiny) {
    class TreeOutputBinding extends Shiny.OutputBinding {
        find(scope) {
            return scope.find(".shiny-tree-output");
        }

        renderValue(el, payload) {
            new ForceGraph3D()(el).graphData(payload);
        }
    }

    Shiny.outputBindings.register(
        new TreeOutputBinding(),
        "shiny-tree-output"
    )
}
