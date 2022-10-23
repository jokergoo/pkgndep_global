
library(pkgndep)
library(ggplot2)
library(cowplot)
library(latex2exp)

########
df = load_pkg_stat_snapshot()
N = nrow(df)

lta = readRDS("~/project/development/pkgndep_analysis/adjusted_heaviness_select_a.rds")
d1 = lta$children; d1$v = 1 - d1$v/N
p1 = ggplot(d1, aes(x = a, y = v)) + geom_point() + geom_line() + geom_vline(xintercept= 10, col = "red", lty =2) +
	labs(x = TeX("Value of \\textit{a}"), y = TeX("Stability of ranks of all heavinesses compared to previous \\textit{a}")) +
	ggtitle("A) Adjust heaviness on child packages")


d3 = lta$indirect_downstream; d3$v = 1 - d3$v/N
p2 = ggplot(d3, aes(x = a, y = v)) + geom_point() + geom_line() + geom_vline(xintercept= 6, col = "red", lty =2) +
	labs(x = TeX("Value of \\textit{a}"), y = TeX("Stability of ranks of all heavinesses compared to previous \\textit{a}")) +
	ggtitle("B) Adjust heaviness on indirect downstream packages")


pdf("fig_04.pdf", width = 10.5, height = 5)
print(plot_grid(p1, p2, nrow = 1))
dev.off()

