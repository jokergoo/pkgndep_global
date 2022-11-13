

setwd("~/manuscript/pkgndep_global/figures")

library(cowplot)
library(pkgndep)

df = load_pkg_stat_snapshot()

pkgs = df$package[df$adjusted_heaviness_on_children >= 30]

nt = do.call(rbind, lapply(pkg, function(x) {
	upstream_dependency(x, snapshot = TRUE)
}))

nt = nt[nt[, 4] >= 20, c(1, 2, 4)]

nt = unique(nt)



library(igraph)
g = graph.edgelist(as.matrix(nt[, 1:2]))
E(g)$heaviness = nt[, 3]

V(g)$group = ifelse(V(g)$name %in% pkgs, "leaf", "others")

library(Rcy3)
cytoscapePing()

createNetworkFromIgraph(g)




library(png)
png = readPNG("hc_upstream.png")


pdf("fig_14.pdf", width = 7, height = 6)

grid.raster(png, x = unit(2, "mm"), just = "left")
grid.points(x = unit(0.1, "npc"), y = unit(0.07, "npc"), pch = 16, size = unit(5, "mm"), gp = gpar(col = "#FF9292"))
grid.text("Top packages with adjusted HC >= 30", x = unit(0.1, "npc") + unit(4, "mm"), y = unit(0.07, "npc"), just = "left", gp = gpar(fontsize = 10))
dev.off()
