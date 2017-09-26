

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
library(crayon)
library(xml2)
library(Rcpp)
library(wCorr)
library(microbenchmark)

library(knitr)
library(markdown)

library(data.tree)
library(dendextend)
library(digest)

if(grepl("tbi", Sys.info()["nodename"]) & Sys.info()["user"] == "guz") {
	Rfiles = list.files("~/project/development/cola/R", full.names = TRUE)
	cpp_files = list.files("~/project/development/cola/src", pattern = "\\.cpp$", full.names = TRUE)
} else {
	Rfiles = list.files("~/project/cola/R", full.names = TRUE)
	cpp_files = list.files("~/project/cola/src", pattern = "\\.cpp$", full.names = TRUE)
}
for(rf in Rfiles) {
	source(rf)
}

for(cf in cpp_files) {
	sourceCpp(cf)
}
