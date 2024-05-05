if (Shiny) {
    class TreeOutputBinding extends Shiny.OutputBinding {
        find(scope) {
            return scope.find(".shiny-tree-output");
        }

        renderValue(el, payload) {
            const graph = new ForceGraph3D()
                .dagMode("td")
                .dagLevelDistance(200)
                .backgroundColor('#101020')
            
            graph(el).graphData(payload);
        }
    }

    Shiny.outputBindings.register(
        new TreeOutputBinding(),
        "shiny-tree-output"
    )
}
