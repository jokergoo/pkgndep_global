

library(pkgndep)
library(ggplot2)
library(scales)
library(ggrepel)
library(igraph)
library(ggnewscale)
library(latex2exp)
library(GetoptLong)
library(grid)

### global graph analysis
lt = load_all_pkg_dep()
library(igraph)
nt_all = data.frame(parents = character(0), children = character(0), heaviness = numeric(0))
for(package in names(lt)) {
    tb = parent_dependency(package, fields = c("Depends", "Imports", "LinkingTo"))[, c(1, 2, 4)]
    nt_all = rbind(nt_all, tb)
}
nt_all = unique(nt_all)

tb = table(nt_all$heaviness)


x = as.numeric(names(tb))
y = tb
y = y/sum(y)
r2 = NULL
for(k in seq(0.01, 1, by = 0.01)) {
    x2 = x^k
    fit = lm(log(y) ~ x2)
    fit = summary(fit)
    r2 = c(r2, fit$r.squared)
}

plot(seq(0.01, 1, by = 0.01), r2)

a = seq(0.01, 1, by = 0.01)[which.max(r2)]

plot(as.numeric(names(tb))^a, log(tb/sum(tb)))

x2 = as.numeric(names(tb))^a
fit = lm(log(tb/sum(tb)) ~ x2)
abline(a = fit$coefficients[1], b = fit$coefficients[2])


p1 = ggplot(data.frame(x = as.numeric(names(tb))^a, y = log(as.numeric(tb)/sum(tb))), aes(x = x, y = y)) +
    geom_point() + geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
    scale_y_continuous(breaks = log(c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5)), labels = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5)) +
    scale_x_continuous(breaks = c(0, 10, 100, 200)^a, labels = c(0, 10, 100, 200)) +
    labs(x = TeX(qq("Heaviness from a parent on a child package (in $h^{@{a}}$ scale)")), y = "Pr(H = h) (in log scale)") +
    ggtitle("A) Heaviness distribution")
p1 = grid.grabExpr(print(p1))

df = load_pkg_stat_snapshot()

l = df$heaviness_on_downstream >= 10
l2 = l & df$heaviness_on_downstream*df$n_downstream >= 5000
p2 = ggplot(df[l, ], aes(x = heaviness_on_downstream, y = heaviness_on_downstream*n_downstream)) +
    geom_point(size = 0.5) + scale_x_continuous(trans = "log10") +
    new_scale_color() +
    geom_point(aes(heaviness_on_downstream, heaviness_on_downstream*n_downstream, color = "high"), data = df[l2, ], show.legend = FALSE) +
    scale_color_manual(name = "Heaviness category", values = c("high" = "red")) +
    ggrepel::geom_text_repel(aes(heaviness_on_downstream, heaviness_on_downstream*n_downstream, label = paste0(package, "(", n_downstream, ")")), data = df[l2, ], min.segment.length = 0, max.overlaps=100, show.legend = FALSE, size = 3, color = "red") +
    labs(x = "Heaviness on downstream packages (in log scale)", y = "Total heaviness on downstream packages") +
    ggtitle(qq("B) @{sum(l)} packages with HD >= 10"))

### depth distribution
library(CePa)
g = graph.edgelist(as.matrix(nt_all[, 1:2]))
x = CePa::reach(g, mode = "in")
y = df$max_heaviness_from_parents
names(y) = df$package
y = y[names(x)]
boxplot(y ~ x)

l = x < 11 & x > 0
l = l & !is.na(y)
x = x[l]
y = y[l]

p3 = ggplot(data.frame(x = factor(x), y = y), aes(x = x, y = y)) +
    geom_boxplot() +
    labs(x = "Depth", y = "Max heaviness from parents") +
    ggtitle("C) Distribution of MHP at each depth")

v = tapply(y, x, max)
v2 = tapply(y, x, length)
# p3 = p3 + geom_text(data.frame(x = 1:length(v), y = as.vector(v) + 10), mapping = aes(x = x, y = y, label = v2))
p4 = ggplot(data.frame(x = 1:length(v), y = v2), aes(x = x, y = y)) +
    geom_bar(stat = "identity") +
    scale_x_continuous(breaks = 1:10) +
    labs(x = "Depth", y = "Number of packages") +
    ggtitle("D) Number of packages at each depth")

library(cowplot)
pdf("fig_11.pdf", width = 12, height = 12)
p = plot_grid(p1, p2, p3, p4, nrow = 2)
print(p)
dev.off()





