
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

df = load_pkg_stat_snapshot()

CUTOFF = list()
# we only looked at the second cutoff
CUTOFF$adjusted_max_heaviness_from_parents = c(50, 60)


df = df[df$n_parents > 0,]

heaviness_cate = ifelse(df$adjusted_max_heaviness_from_parents >= CUTOFF$adjusted_max_heaviness_from_parents[2], "high", "low")
heaviness_color = c("high" = "red")
repo = ifelse(grepl("bioconductor", df$repository), "Bioconductor", "CRAN")
df$repo = factor(repo, levels = c("CRAN", "Bioconductor"))

x = df$n_parents
y = df$max_heaviness_from_parents
max_x = max(x)
max_y = max(y)
x = x/max(x)
y = y/max(y)

df$n_neighbours = 0
for(r in c("CRAN", "Bioconductor")) {
	l = df$repo == r
	d = dist(cbind(x[l], y[l]))
	d = as.matrix(d)
	df$n_neighbours[l] = apply(d, 1, function(x) sum(x < 0.001))
}

df = df[order(df$n_neighbours), ] # put points with higher n_neighbours on top

l_high = df$adjusted_max_heaviness_from_parents >= CUTOFF$adjusted_max_heaviness_from_parents[2]
l_low = !l_high
p3 = ggplot(df[l_low, ], aes(n_parents, max_heaviness_from_parents, color = log10(n_neighbours))) +
	geom_point(size = 0.5) + scale_color_viridis_c(breaks = c(0, 1, 2, 3), labels = c(1, 10, 100, 1000), name = "Number of\nneighbor points") +
    new_scale_color() +
    geom_point(aes(n_parents, max_heaviness_from_parents, color = "high"), data = df[l_high, ]) +
    scale_color_manual(name = "Heaviness category", values = heaviness_color) +
    ggrepel::geom_text_repel(aes(n_parents, max_heaviness_from_parents, label = package), data = df[l_high, ], min.segment.length = 0, box.padding = 0.5, max.overlaps = Inf, show.legend = FALSE, size = 3, color = "red") +
	labs(x = "Number of parent packages", y = "Max heaviness from parents") +
	ggtitle("A) Max heaviness from parents") +
	facet_wrap(vars(repo))


gini = sapply(lt, function(x) {
    v = heaviness(x)[x$which_required]
    gini_index(v)
})

gini = gini[df$package]


df2 = data.frame(heaviness = df$max_heaviness_from_parents, gini = gini, repo = df$repo, cate = ifelse(l_high, "high", "low"))
df2$repo = factor(df2$repo, levels = c("CRAN", "Bioconductor"))
df2 = df2[df2$gini > 0 & df2$gini < 1, ]

p5 = ggplot(df2, aes(x = heaviness, y = gini)) +
    geom_point(size = 0.5) + 
    new_scale_color() +
    geom_point(aes(heaviness, gini, color = "high"), data = df2[df2$cate == "high", ]) +
    scale_color_manual(name = "Heaviness category", values = heaviness_color) +
    labs(x = "Max heaviness from parents", y = "Gini index") +
    ggtitle("B) Gini index of heaviness from strong parents") +
    facet_wrap(vars(repo))
    


pdf("fig_06.pdf", width = 12, height = 10)
print(plot_grid(p3, p5, ncol = 1))
dev.off()

