

library(pkgndep)
library(ggplot2)
library(scales)
library(ggrepel)
library(igraph)
library(ggnewscale)

### global graph analysis
lt = load_all_pkg_dep()
library(igraph)
nt_all = data.frame(parents = character(0), children = character(0), heaviness = numeric(0))
for(package in names(lt)) {
    tb = parent_dependency(package, fields = c("Depends", "Imports", "LinkingTo"))[, c(1, 2, 4)]
    nt_all = rbind(nt_all, tb)
}
nt_all = unique(nt_all)


############# heavy core graph #########
nt2 = nt_all[nt_all$heaviness >= 30, ]
g1 = igraph::graph.edgelist(as.matrix(nt2[, 1:2]))
E(g1)$heaviness = nt2[, 3]

score = readRDS("~/project/development/pkgndep_analysis/pkg_stat_score.rds")

nt = data.frame(from = character(0), to = character(0), heaviness = numeric(0))
for(package in names(lt)) {
    print(package)
    v = attr(score[[package]], "values")
    tb = data.frame(from = rep(package, length(v)), to = names(v), heaviness = v)
    nt = rbind(nt, tb)
}
nt = unique(nt)
nt2 = nt[nt$heaviness >= 30, ]

g2 = igraph::graph.edgelist(as.matrix(nt2[, 1:2]))
E(g2)$heaviness = nt2[, 3]


g3 = induced_subgraph(g1, intersect(V(g1)$name, V(g2)$name))
g3 = induced_subgraph(g3, V(g3)$name[degree(g3) > 0])

### g3 is identical to g1
identical(V(g1)$name, V(g3)$name)
el1 = as_edgelist(g1)
el2 = as_edgelist(g3)
el1 = el1[order(el1[, 1], el1[, 2]), ]
el2 = el2[order(el2[, 1], el2[, 2]), ]
identical(el1, el2)

###### g1 or g3 is the core graph #######

E(g1)$betweenness = edge_betweenness(g1)
E(g1)$betweenness_cate = ifelse(E(g1)$betweenness >= 20, "high", "low")

V(g1)$out_degree = degree(g1, mode = "out")
V(g1)$label = ifelse(V(g1)$out_degree >= 30, V(g1)$name, "")

library(RCy3)
cytoscapePing()

createNetworkFromIgraph(g1)

library(png)
p1 = grid.grabExpr({
	grid.text("A) The core dependency graph", x = unit(15, "mm"), y = unit(1, "npc") - unit(4, "mm"), just = c("left", "top"), gp = gpar(fontsize = 14))
	grid.raster(readPNG("core-heavy-graph.png"), width = unit(0.85, "snpc"))
})

## distance to leaves
mm = distances(g3, to = V(g3)$name[degree(g3, mode = "out") == 0], mode = "out")
x = mm[is.finite(mm) & mm > 0]
tb1 = table(x)

p2 = ggplot(data.frame(x = x), aes(x = x)) + 
    geom_bar() +
    labs(x = "Distance to the leaf packages", y = "Count") +
    geom_text(data = data.frame(x = as.numeric(names(tb1)), y = as.vector(tb1)), mapping = aes(x = x, y = y, label = y), vjust = -0.5) +
    ggtitle("D) Distribution of distance to leaf packages") +
    scale_x_continuous(breaks = 1:7, label = 1:7)

com = components(g1)
tb = table(com$csize)

x = as.numeric(names(tb))
y = as.vector(tb)
y = y/sum(y)

fit = lm(log(y)[1:19] ~ log(x)[1:19])

p3 = ggplot(data.frame(x = log(x), y = log(y)), aes(x = x, y = y)) +
	geom_point() + geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
	scale_x_continuous(breaks = log(c(10, 100, 1000)), label = c(10, 100, 1000)) +
	scale_y_continuous(breaks = log(c(1e-3, 1e-2, 1e-1, 0.5)), label = c(1e-3, 1e-2, 0.1, 0.5)) +
	labs(x = "Graph components size (X)", y = "Pr(X = x)") +
	ggtitle("B) Distribution of graph component sizes")

l = E(g1)$betweenness >= 20
gcore = subgraph.edges(g1, which(l))

createNetworkFromIgraph(gcore)


library(png)

p4 = grid.grabExpr({
	grid.text("C) Sub-graph with edge betweenness >= 20", x = unit(15, "mm"), y = unit(1, "npc") - unit(4, "mm"), just = c("left", "top"), gp = gpar(fontsize = 14))
	grid.raster(readPNG("between-core.png"), width = unit(0.85, "snpc"))
})

library(cowplot)
pdf("fig_12.pdf", width = 12, height = 12)
p = plot_grid(p1, p3, p4, p2, nrow = 2)
print(p)
dev.off()

