df = load_pkg_stat_snapshot()
repo = ifelse(grepl("bioconductor", df$repository), "Bioconductor", "CRAN")
df$repo = factor(repo, levels = c("CRAN", "Bioconductor"))


nm = c("n_by_strong", "n_parents", "max_heaviness_from_parents",  "max_co_heaviness_from_parents",
       "n_children", "n_children (> 0)", "heaviness_on_children", "n_indirect_downstream", "n_indirect_downstream (> 0)", "heaviness_on_indirect_downstream")

tb1 = t(sapply(nm, function(x) {
	if(x == "heaviness_on_children") {
		l = df$n_children > 0
		tapply(df[[x]][l], df$repo[l], mean)
	} else if(x == "heaviness_on_indirect_downstream") {
		l = df$n_indirect_downstream > 0
		tapply(df[[x]][l], df$repo[l], mean)
	} else if(x == "n_children (> 0)") {
		l = df$n_children > 0
		tapply(df$n_children[l], df$repo[l], mean)
	} else if(x == "n_indirect_downstream (> 0)") {
		l = df$n_indirect_downstream > 0
		tapply(df$n_indirect_downstream[l], df$repo[l], mean)
	} else {
		tapply(df[[x]], df$repo, mean)
	}
}))


rownames(tb1) = c("Number of strong dependencies",
	"Number of parents",
	"Max heaviness from parents",
	"Max co-heaviness from parents",
	"Number of children",
	"Number of children (n > 0)",
	"Heaviness on child packages",
	"Number of indirect downstream",
	"Number of indirect downstream (n > 0)",
	"Heaviness on indirect downstream packages")

tb1 = round(tb1, 1)



v = sapply(1:60, function(i) {
	l = df$n_children >= i
	x = tapply(df$heaviness_on_children[l], df$repo[l], mean)
	x[2]/x[1]
})
matplot(1:60, t(v), type = "o")


v = sapply(1:30, function(i) {
	l = df$adjusted_heaviness_on_indirect_downstream >= i
	x = tapply(df$n_indirect_downstream[l], df$repo[l], length)
	x[1]/x[2]
})
matplot(1:30, t(v), type = "o")