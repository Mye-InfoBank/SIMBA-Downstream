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
                .linkColor(() => 'rgba(255,255,255,0.2)')
                .nodeOpacity(0.9)
                .linkDirectionalParticles(2)
                .linkDirectionalParticleWidth(0.8)
                .linkDirectionalParticleSpeed(0.006)

            graph(el).graphData(payload);
        }
    }

    Shiny.outputBindings.register(
        new TreeOutputBinding(),
        "shiny-tree-output"
    )
}
