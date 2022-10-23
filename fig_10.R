

library(grid)
library(png)

p1 = readPNG("example_car.png")
p2 = readPNG("example_rstatix.png")
p3 = readPNG("example_tidyverse.png")


pdf("fig_10.pdf", width = 12, height = 12)

grid.newpage()
grid.raster(p1, x = 1/4, y = unit(3/4, "npc") - unit(2, "mm"), height = 1/2, default.unit = "npc")
grid.raster(p2, x = 3/4, y = 3/4, height = 1/2/(dim(p1)[1]/dim(p2)[1]), default.unit = "npc")
grid.raster(p3, x = 1/4, y = 1/4, height = 1/2/(dim(p1)[1]/dim(p3)[1]), default.unit = "npc")

grid.text("A", x = unit(4, "mm"), y = unit(1, "npc") - unit(4, "mm"), just = c("left", "top"), gp = gpar(fontsize = 14))
grid.text("B", x = unit(4, "mm") + unit(1/2, "npc"), y = unit(1, "npc") - unit(4, "mm"), just = c("left", "top"), gp = gpar(fontsize = 14))
grid.text("C", x = unit(4, "mm"), y = unit(1/2, "npc") - unit(4, "mm"), just = c("left", "top"), gp = gpar(fontsize = 14))

dev.off()
