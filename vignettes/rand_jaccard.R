
## GDS534, "ATC:mclust" method
png("~/rand_jaccard.png", width = 620, height = 620)

col_fun = colorRamp2(c(0, 1), c("white", "blue"))
ht1 = Heatmap(get_consensus(res, k = 3), name = "matrix_1", col = col_fun, show_row_names = FALSE, show_column_names = FALSE,
	show_row_dend = FALSE, show_column_dend = FALSE, column_title = "k = 3")
ht2 = Heatmap(get_consensus(res, k = 4), name = "matrix_2", col = col_fun, show_row_names = FALSE, show_column_names = FALSE,
	show_row_dend = FALSE, show_column_dend = FALSE, column_title = "k = 4")
stat = as.data.frame(get_stats(res))

grid.newpage()
par(mfrow = c(2, 2))
plot.new()
plot.new()
plot(stat$k, stat$Rand, xlab = "k", ylab = "Rand", type = "b")
plot(stat$k, stat$Jaccard, xlab = "k", ylab = "Rand", type = "b")
pushViewport(viewport(x = 0, y = 1, width = 0.5, height = 0.5, just = c("left", "top")))
draw(ht1, newpage = FALSE)
popViewport()
pushViewport(viewport(x = 0.5, y = 1, width = 0.5, height = 0.5, just = c("left", "top")))
draw(ht2, newpage = FALSE)
popViewport()
dev.off()
