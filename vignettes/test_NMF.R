library(GetoptLong)
library(cola)

library(golubEsets)
data(Golub_Merge)
m = exprs(Golub_Merge)
anno = pData(Golub_Merge)

m[m <= 1] = NA
m = log10(m)

m = adjust_matrix(m)

set.seed(123)
m2 = m[sample(nrow(m), 500), ]

NMF1 = function(mat, k, ...) {
	fit = NNLM::nnmf(A = mat, k = k, verbose = FALSE, ...)
	fit$H
}
NMF2 = function(mat, k, ...) {
	NMF::nmf.options(maxIter = 500)
	fit = NMF::nmf(mat, rank = k)
	NMF::nmf.options(maxIter = NULL)
	fit@fit@H
}

h1 = NMF1(m2, k = 3)
h2 = NMF2(m2, k = 3)


library(grid)
library(ComplexHeatmap)

png("NMF_compare.png", width = 1000, height = 300)

grid.newpage()
pushViewport(viewport(x = 0, width = 0.5, just = "left"))
ha = HeatmapAnnotation(group = anno$ALL.AML, col = list("group" = c("ALL" = "orange", "AML" = "purple")))
draw(Heatmap(h1, top_annotation=ha, show_column_names = FALSE, show_row_dend = FALSE, column_title = "H matrix, NMF package"), newpage = FALSE)
popViewport()
pushViewport(viewport(x = 0.5, width = 0.5, just = "left"))
draw(Heatmap(h2, top_annotation=ha, show_column_names = FALSE, show_row_dend = FALSE, column_title = "H matrix, NNLM package"), newpage = FALSE)
popViewport()

dev.off()
