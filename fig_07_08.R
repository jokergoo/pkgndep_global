

library(ggplot2)
library(ggrepel)
library(GetoptLong)
library(pkgndep)

df = load_pkg_stat_snapshot()

repo = ifelse(grepl("bioc", df$repository), "bioc", "cran")
names(repo) = df$package

v = df$max_co_heaviness_from_parents

tb = table(v)
x = as.numeric(names(tb))
y = as.numeric(tb)

plot(x, y)


od = order(-v)
df[od[1:10], c("package", "repository", "max_co_heaviness_from_parents", "max_co_heaviness_parents_pair", "max_co_heaviness_parents_pair_type")]

l = x >= 40

plot(x[l], y[l], type = "h")


l = v >= 40
df$max_co_heaviness_parents_pair_type[df$max_co_heaviness_parents_pair_type == ""] = "No clear relation"
df$max_co_heaviness_parents_pair_type[df$max_co_heaviness_parents_pair_type == "have-common-upstream"] = "Common-upstream"
df$max_co_heaviness_parents_pair_type[df$max_co_heaviness_parents_pair_type == "parent-child"] = "Parent-child"
df$max_co_heaviness_parents_pair_type[df$max_co_heaviness_parents_pair_type == "upstream-downstream"] = "Upstream-downstream"

df$max_co_heaviness_parents_pair_type = factor(df$max_co_heaviness_parents_pair_type, levels = c("Parent-child", "Upstream-downstream", "Common-upstream", "No clear relation"))
p = ggplot(df[l, ], aes(x = max_co_heaviness_parents_pair_type)) +
  geom_bar() +
  labs(x = "Relations of the MCoHP parents", y = "Count")

pdf("fig_07.pdf", width = 6, height = 6)
print(p)
dev.off()


df2 = df[v >= 40 & grepl("bioc", df$repository), ]
pair = df2[, "max_co_heaviness_parents_pair"]

mean = tapply(df[v >= 40, "max_co_heaviness_from_parents"], df[v>=40, "max_co_heaviness_parents_pair"], mean)

tb = table(pair)
tb = tb[order(-tb)]
labels = names(tb)
labels[-(1:5)] = ""
labels[1:5] = paste0(labels[1:5], " (", round(mean[labels[1:5]]), ")")
tb1 = data.frame(x = 1:length(tb), v = as.numeric(tb), labels = labels)


p1 = ggplot(data.frame(x = 1:length(tb), v = as.numeric(tb), labels = labels), aes(x = x, y = v, label = labels)) +
  geom_linerange(aes(x = x, ymin = 0, ymax = v)) +
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf, direction = "y", hjust = -0.5) +
  labs(x = "MCoHP parents", y = "Number of child packages having the parent pair") +
  ggtitle(qq("A) @{sum(tb)} Bioconductor packages"))


pair = do.call(rbind, strsplit(pair, ","))

library(igraph)

g = graph.edgelist(pair, directed = FALSE)
d = degree(g)

mean = sapply(names(d), function(x) {
	l = pair[, 1] %in% x | pair[, 2] %in% x
	mean(df2[l, "max_co_heaviness_from_parents"])
})

nu = sapply(names(d), function(x) {
	l = pair[, 1] %in% x | pair[, 2] %in% x
	length(unique(as.vector(pair[l, ]))) - 1
})

tb = d[order(-d)]
labels = names(tb)
labels[-(1:8)] = ""
labels[1:8] = paste0(labels[1:8], " (", round(mean[labels[1:8]]), ", ", nu[labels[1:8]], " companions)")

p2 = ggplot(data.frame(x = 1:length(tb), v = as.numeric(tb), labels = labels), aes(x = x, y = v, label = labels)) +
  geom_linerange(aes(x = x, ymin = 0, ymax = v)) +
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf, direction = "y", hjust = -0.5) +
  labs(x = "Members in the MCoHP parents", y = "Number of child packages having the co-parent") +
  ggtitle(qq("B) @{sum(tb)/2} Bioconductor packages"))

#### exclude

exclude = c("AnnotationDbi", "org.Hs.eg.db", "GenomicFeatures", "AnnotationHub",
     "ExperimentHub", "org.Rn.eg.db", "BSgenome", "org.Mm.eg.db")
df2 = df[v >= 40 & grepl("bioc", df$repository), ]
df2 = df2[sapply(strsplit(df2$max_co_heaviness_parents_pair, ","), function(x) length(intersect(x, exclude)) == 0), ]
pair = df2[, "max_co_heaviness_parents_pair"]

mean = tapply(df[v >= 40, "max_co_heaviness_from_parents"], df[v>=40, "max_co_heaviness_parents_pair"], mean)

tb = table(pair)
tb = tb[order(-tb)]
labels = names(tb)
labels[-1] = ""
labels[1] = paste0(labels[1], " (", round(mean[labels[1]]), ")")
tb1 = data.frame(x = 1:length(tb), v = as.numeric(tb), labels = labels)


p3 = ggplot(data.frame(x = 1:length(tb), v = as.numeric(tb), labels = labels), aes(x = x, y = v, label = labels)) +
  geom_linerange(aes(x = x, ymin = 0, ymax = v)) +
  geom_text_repel(box.padding = 0.25, max.overlaps = Inf, direction = "y", hjust = -0.5) +
  labs(x = "MCoHP parents", y = "Number of child packages having the parent pair") +
  ggtitle(qq("C) @{sum(tb)} Bioconductor packages,\nwhere @{length(exclude)} annotation-related parent packages are removed."))


pair = do.call(rbind, strsplit(pair, ","))

library(igraph)

g = graph.edgelist(pair, directed = FALSE)
d = degree(g)

mean = sapply(names(d), function(x) {
	l = pair[, 1] %in% x | pair[, 2] %in% x
	mean(df2[l, "max_co_heaviness_from_parents"])
})

nu = sapply(names(d), function(x) {
	l = pair[, 1] %in% x | pair[, 2] %in% x
	length(unique(as.vector(pair[l, ]))) - 1
})

tb = d[order(-d)]
labels = names(tb)
labels[-(1:5)] = ""
labels[1:5] = paste0(labels[1:5], " (", round(mean[labels[1:5]]), ", ", nu[labels[1:5]], " companions)")

p4 = ggplot(data.frame(x = 1:length(tb), v = as.numeric(tb), labels = labels), aes(x = x, y = v, label = labels)) +
  geom_linerange(aes(x = x, ymin = 0, ymax = v)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8)) +
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf, direction = "y", hjust = -0.5) +
  labs(x = "Members in the MCoHP parents", y = "Number of child packages having the co-parent") +
  ggtitle(qq("D) @{sum(tb)/2} Bioconductor packages,\nwhere @{length(exclude)} annotation-related parent packages are removed."))


#### CRAN

df2 = df[v >= 40 & !grepl("bioc", df$repository), ]
pair = df2[, "max_co_heaviness_parents_pair"]

mean = tapply(df[v >= 40, "max_co_heaviness_from_parents"], df[v>=40, "max_co_heaviness_parents_pair"], mean)

tb = table(pair)
tb = tb[order(-tb)]
labels = names(tb)
labels[-(1:5)] = ""
labels[1:5] = paste0(labels[1:5], " (", round(mean[labels[1:5]]), ")")
tb1 = data.frame(x = 1:length(tb), v = as.numeric(tb), labels = labels)


p5 = ggplot(data.frame(x = 1:length(tb), v = as.numeric(tb), labels = labels), aes(x = x, y = v, label = labels)) +
  geom_linerange(aes(x = x, ymin = 0, ymax = v)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8)) +
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf, direction = "y", hjust = -0.5) +
  labs(x = "MCoHP parents", y = "Number of child packages having the parent pair") +
  ggtitle(qq("E) @{sum(tb)} CRAN packages"))


pair = do.call(rbind, strsplit(pair, ","))

library(igraph)

g = graph.edgelist(pair, directed = FALSE)
d = degree(g)

mean = sapply(names(d), function(x) {
	l = pair[, 1] %in% x | pair[, 2] %in% x
	mean(df2[l, "max_co_heaviness_from_parents"])
})

nu = sapply(names(d), function(x) {
	l = pair[, 1] %in% x | pair[, 2] %in% x
	length(unique(as.vector(pair[l, ]))) - 1
})

tb = d[order(-d)]
labels = names(tb)
labels[-(1:5)] = ""
labels[1:5] = paste0(labels[1:5], " (", round(mean[labels[1:5]]), ", ", nu[labels[1:5]], " companions)")

p6 = ggplot(data.frame(x = 1:length(tb), v = as.numeric(tb), labels = labels), aes(x = x, y = v, label = labels)) +
  geom_linerange(aes(x = x, ymin = 0, ymax = v)) +
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf, direction = "y", hjust = -0.5) +
  labs(x = "Members in the MCoHP parents", y = "Number of child packages having the co-parent") +
  ggtitle(qq("F) @{sum(tb)/2} CRAN packages"))



library(cowplot)

pdf("fig_08.pdf", width = 12, height = 16)
p = plot_grid(p1, p2, p3, p4, p5, p6, nrow = 3)
print(p)
dev.off()

