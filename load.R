

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
library(e1071)
library(skmeans)
library(cclust)
library(digest)
library(png)
# library(venneuler)
library(gplots)
library(ggplot2)
library(bioDist)
library(pryr)
library(Rtsne)

Rfiles = list.files("~/project/development/cola/R", full.names = TRUE)
for(rf in Rfiles) {
	source(rf)
}
