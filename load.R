

library(ComplexHeatmap)
library(matrixStats)
library(GetoptLong)
library(circlize)
library(clue)
library(parallel)
library(RColorBrewer)
library(cluster)
library(NMF)
library(mclust)
library(skmeans)
library(cclust)
library(png)
# library(venneuler)
library(gplots)
library(Rtsne)
library(samr)
library(genefilter)
library(kohonen)

Rfiles = list.files("~/project/cola/R", full.names = TRUE)
for(rf in Rfiles) {
	source(rf)
}
