

setwd("~/manuscript/pkgndep_global/figures")


library(cowplot)
library(pkgndep)
library(grid)


compute_density = function (x, w, from, to, bw = "nrd0", adjust = 1, kernel = "gaussian", n = 512) {
    nx <- length(x)
    if (is.null(w)) {
        w <- rep(1/nx, nx)
    }
    else {
        w <- w/sum(w)
    }
    if (nx < 2) {
        warn("Groups with fewer than two data points have been dropped.")
        return(ggplot2:::new_data_frame(list(x = NA_real_, density = NA_real_, 
            scaled = NA_real_, ndensity = NA_real_, count = NA_real_, 
            n = NA_integer_), n = 1))
    }
    dens <- stats::density(x, weights = w, bw = bw, adjust = adjust, 
        kernel = kernel, n = n, from = from, to = to)
    df = ggplot2:::new_data_frame(list(x = dens$x, density = dens$y, scaled = dens$y/max(dens$y, 
        na.rm = TRUE), ndensity = dens$y/max(dens$y, na.rm = TRUE), 
        count = dens$y * nx, n = nx), n = length(dens$x))
    
    q = quantile(x, 0.99)
    df = df[df[, "x"] <= q, ]
    df
}

assignInNamespace("compute_density", compute_density, ns = "ggplot2")


library(ggplot2)
library(ggnewscale)


CUTOFF = list()
CUTOFF$adjusted_max_heaviness_from_parents = c(50, 60)
CUTOFF$adjusted_total_heaviness_from_parents = c(70, 100)
CUTOFF$adjusted_heaviness_on_children = c(15, 30)
CUTOFF$adjusted_heaviness_on_indirect_downstream = c(10, 20)


df = load_pkg_stat_snapshot()
df = df[df$n_children > 0,]


df[df$adjusted_heaviness_on_children >= 30, ]

heaviness_cate = ifelse(df$adjusted_heaviness_on_children >= CUTOFF$adjusted_heaviness_on_children[2], "high", ifelse(df$adjusted_heaviness_on_children >= CUTOFF$adjusted_heaviness_on_children[1], "median", "low"))
heaviness_color = c("high" = "red")
repo = ifelse(grepl("bioconductor", df$repository), "Bioconductor", "CRAN")
df$repo = factor(repo, levels = c("CRAN", "Bioconductor"))

x = log10(df$n_children)
y = df$heaviness_on_children
max_x = max(x)
max_y = max(y)
x = x/max(x)
y = y/max(y)

df$n_neighbours = 0
for(r in c("CRAN", "Bioconductor")) {
	l = df$repo == r
	d = dist(cbind(x[l], y[l]))
	d = as.matrix(d)
	df$n_neighbours[l] = apply(d, 1, function(x) sum(x < 0.01))
}
df = df[order(df$n_neighbours), ] # put points with higher n_neighbours on top

l_high = df$adjusted_heaviness_on_children >= CUTOFF$adjusted_heaviness_on_children[2]
l_low = !(l_high)
p1 = ggplot(df[l_low, ], aes(n_children, heaviness_on_children, color = log10(n_neighbours))) +
    geom_point(size = 0.5) + scale_color_viridis_c(breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000), name = "Number of\nneighbor points") +
    scale_x_continuous(trans='log10') +
    new_scale_color() +
    geom_point(aes(n_children, heaviness_on_children, color = "high"), data = df[l_high, ]) +
    scale_color_manual(name = "Heaviness category", values = heaviness_color) +
    ggrepel::geom_text_repel(aes(n_children, heaviness_on_children, label = package), data = df[l_high, ], min.segment.length = 0, box.padding = 0.5, max.overlaps = Inf, show.legend = FALSE, size = 3, color = "red") +
    labs(x = "Number of child packages (in log scale)", y = "Heaviness on child packages") +
    ggtitle("A) Heaviness on child packages") +
    facet_wrap(vars(repo))
p1 = grid.grabExpr(print(p1), width = 11.5, height = 6)

##
df = load_pkg_stat_snapshot()
df = df[df$n_indirect_downstream > 0,]

heaviness_cate = ifelse(df$adjusted_heaviness_on_indirect_downstream >= CUTOFF$adjusted_heaviness_on_indirect_downstream[2], "high", ifelse(df$adjusted_heaviness_on_indirect_downstream >= CUTOFF$adjusted_heaviness_on_indirect_downstream[1], "median", "low"))
heaviness_color = c("high" = "red")
repo = ifelse(grepl("bioconductor", df$repository), "Bioconductor", "CRAN")
df$repo = factor(repo, levels = c("CRAN", "Bioconductor"))

x = log10(df$n_indirect_downstream)
y = df$heaviness_on_indirect_downstream
max_x = max(x)
max_y = max(y)
x = x/max(x)
y = y/max(y)

df$n_neighbours = 0
for(r in c("CRAN", "Bioconductor")) {
	l = df$repo == r
	d = dist(cbind(x[l], y[l]))
	d = as.matrix(d)
	df$n_neighbours[l] = apply(d, 1, function(x) sum(x < 0.01))
}
df = df[order(df$n_neighbours), ] # put points with higher n_neighbours on top

l_high = df$adjusted_heaviness_on_indirect_downstream >= CUTOFF$adjusted_heaviness_on_indirect_downstream[2]
l_low = !l_high
p2 = ggplot(df[l_low, ], aes(n_indirect_downstream, heaviness_on_indirect_downstream, color = log10(n_neighbours))) +
    geom_point(size = 0.5) + scale_color_viridis_c(breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000), name = "Number of\nneighbor points") +
    scale_x_continuous(trans='log10') +
    new_scale_color() +
    geom_point(aes(n_indirect_downstream, heaviness_on_indirect_downstream, color = "high"), data = df[l_high, ]) +
    scale_color_manual(name = "Heaviness category", values = heaviness_color) +
    ggrepel::geom_text_repel(aes(n_indirect_downstream, heaviness_on_indirect_downstream, label = package), data = df[l_high, ], min.segment.length = 0, box.padding = 0.5, max.overlaps = Inf, show.legend = FALSE, size = 3, color = "red") +
    labs(x = "Number of indirect downstream packages (in log scale)", y = "Heaviness on indirect downstream packages") +
    ggtitle("D) Heaviness on indirect downstream packages") +
    facet_wrap(vars(repo))
p2 = grid.grabExpr(print(p2), width = 11.5, height = 6)


correspond_between_two_rankings = function(x1, x2, name1, name2, 
    col1 = 2, col2 = 3, top_n = round(0.25*length(x1)), transparency = 0.9, 
    pt_size = unit(1, "mm"), newpage = TRUE, ratio = c(1, 1, 1)) {
    
    if(newpage) {
        grid.newpage()
    }

    if(length(x1) != length(x2)) {
        stop("Length of `x1` and `x2` should be the same.")
    }

    r1 = rank(x1, ties.method = "random")
    r2 = rank(x2, ties.method = "random")

    if(missing(name1)) name1 = deparse(substitute(x1))
    if(missing(name2)) name2 = deparse(substitute(x2))

    n = length(x1)
    text_height = grobHeight(textGrob("foo"))*2
    pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 3, widths = unit(ratio, "null")), 
        width = unit(1, "npc") - unit(2, "mm"),
        height = unit(1, "npc") - text_height - unit(1, "cm"), y = unit(1, "cm"), just = "bottom"))
    
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2, xscale = c(0, 1), yscale = c(0, n + 1)))
    l = r1 >= n - top_n | r2 >= n - top_n
    # if(sum(!l)) grid.segments(0, r1[!l], 1, r2[!l], default.units = "native", gp = gpar(col = "#EEEEEE80"))
    if(sum(l)) {
        grid.segments(0, r1[l], 1, r2[l], default.units = "native", gp = gpar(col = add_transparency("#000000", transparency)))
        # for(ind in which(l)) {
        #   grid.bezier(c(0, 1, 0, 1), c(r1[ind], r1[ind], r2[ind], r2[ind]), default.units = "native", gp = gpar(col = add_transparency("#000000", transparency)))
        # }
    }
    upViewport()

    max_x1 = max(x1)
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1, 
        xscale = c(0, max_x1), yscale = c(0, n + 1)))
    grid.segments(max_x1 - x1, r1, max_x1, r1, default.units = "native", gp = gpar(col = "#EFEFEF"))
    upViewport()

    max_x2 = max(x2)
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3, 
        xscale = c(0, max_x2), yscale = c(0, n + 1)))
    grid.segments(0, r2, x2, r2, default.units = "native", gp = gpar(col = "#EFEFEF"))
    upViewport()

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1, 
        xscale = c(0, max_x1), yscale = c(0, n + 1)))
    l = r2 >= n - top_n
    grid.points(max_x1 - x1[l], r1[l], default.units = "native", pch = 16, size = pt_size, gp = gpar(col = add_transparency(col2, 0.5)))
    grid.text(name1, x = 1, y = unit(n + 1, "native") + unit(1, "mm"), default.units = "npc", just = c("right", "bottom"))
    upViewport()

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3, 
        xscale = c(0, max_x2), yscale = c(0, n + 1)))
    l = r1 >= n - top_n
    grid.points(x2[l], r2[l], default.units = "native", pch = 16, size = pt_size, gp = gpar(col = add_transparency(col1, 0.5)))
    grid.text(name2, x = 0, y = unit(n + 1, "native") + unit(1, "mm"), default.units = "native", just = c("left", "bottom"))
    upViewport()

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2, xscale = c(0, 1), yscale = c(0, n + 1)))
    grid.segments(c(0, 1), c(n - top_n, n - top_n), c(0, 1), c(n, n), default.units = "native", gp = gpar(lwd = 4, col = c(col1, col2)))
    upViewport()

    
    upViewport()

    # add a venn diagram at the bottom
    n_intersect = length(intersect(order(x1, decreasing = TRUE)[1:top_n], order(x2, decreasing = TRUE)[1:top_n]))
    n_union = 2*top_n - n_intersect
    grid.roundrect(x = unit(0.5 - n_intersect/2/top_n*0.4, "npc"), y = unit(0.4, "cm"), width = unit(0.4, "npc"), 
        height = unit(0.4, "cm"), gp = gpar(fill = add_transparency(col2, 0.5), col = NA), just = "left")
    grid.roundrect(x = unit(0.5 + n_intersect/2/top_n*0.4, "npc"), y = unit(0.4, "cm"), width = unit(0.4, "npc"), 
        height = unit(0.4, "cm"), gp = gpar(fill = add_transparency(col1, 0.5), col = NA), just = "right")
    grid.text(qq("top @{top_n}/@{length(x1)}"), x = unit(0.5, "npc"), y = unit(0.7, "cm"), just = "bottom", gp = gpar(fontsize = 8))

}

df = load_pkg_stat_snapshot()

p3 = grid.grabExpr({
    grid.newpage()

    pushViewport(viewport(width = 0.95, height = 0.95, x = unit(-5, "mm"), just = "left"))
    correspond_between_two_rankings(x1 = df$heaviness_on_children, x2 = df$heaviness_on_downstream, 
        name1 = "HC", name2 = "HD", top_n = 500, newpage = FALSE)
    upViewport()
})

p4 = grid.grabExpr({
    grid.newpage()
    pushViewport(viewport(width = 0.95, height = 0.95))
    correspond_between_two_rankings(x1 = df$heaviness_on_children, x2 = df$heaviness_on_indirect_downstream, 
        name1 = "HC", name2 = "HID", top_n = 500, newpage = FALSE)
    upViewport()
})

p_subtitle = textGrob("C) Compare top packages with the higest HC, HC and HID values", x = unit(-5, "mm"), just = "left")

df = load_pkg_stat_snapshot()

df$child_gini = sapply(df$package, function(x) {
    gini_index(child_dependency(x)[, 4])
})

heaviness_color = c("high" = "red")

p5 = ggplot(df[df$n_children > 1, ], aes(adjusted_heaviness_on_children, child_gini)) +
    geom_point(size = 0.5) + 
    new_scale_color() +
    geom_point(aes(adjusted_heaviness_on_children, child_gini, color = "high"), data = df[df$n_children > 1 & df$adjusted_heaviness_on_children >= 30, ]) +
    scale_color_manual(name = "Heaviness category", values = heaviness_color) +
    labs(x = "Adjusted HC", y = "Gini index", col = "top packages") +
    ggtitle("B) Gini index of heaviness on child packages")
p5 = grid.grabExpr(print(p5))



p_empty = rectGrob(gp = gpar(fill = NA, col = NA))


png("figure_09_C.png", width = 1000, height = 1000, res = 200)
plot_grid(p3 ,p4)
dev.off()

library(png)
png = readPNG("figure_09_C.png")



pdf("fig_09.pdf", width = 12, height = 15)
p = plot_grid( 
    p1, 
    plot_grid(p5, plot_grid(p_subtitle, 
            rasterGrob(png, x = unit(0, "cm"), just = "left"), 
            nrow = 2, rel_heights = c(0.1, 1)), rel_widths = c(1.3, 1)),
    p2,
    nrow = 3)
print(p)
dev.off()



library(knitr)
l = df$adjusted_heaviness_on_children >= 30 | df$adjusted_heaviness_on_indirect_downstream >= 20
df[l, ] -> df2
df2$top_HC = ifelse(df2$adjusted_heaviness_on_children >= 30, "y", "")
df2$top_HID = ifelse(df2$adjusted_heaviness_on_indirect_downstream >= 20, "y", "")
df2 = df2[order(-df2$n_downstream), ]
ind = c(1, 18, 17, 25, 24, 29, 30)
kable(df2[df2$n_downstream >= 30, ind], row.names = F)



