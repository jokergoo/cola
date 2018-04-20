
library(microbenchmark)

print(microbenchmark(
	library(ComplexHeatmap),
	library(matrixStats),
	library(GetoptLong),
	library(circlize),
	library(clue),
	library(parallel),
	library(RColorBrewer),
	library(cluster),
	library(NNLM),
	library(mclust),
	library(skmeans),
	library(cclust),
	library(png),
#	library(venneuler),
	library(gplots),
	library(Rtsne),
	library(samr),
	library(kohonen),
	library(crayon),
	library(xml2),
	library(Rcpp),
	library(knitr),
	library(markdown),
	library(data.tree),
	library(dendextend),
	library(digest),
	library(impute),
	library(httr),
	library(brew),
	times = 1
))

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


 #                    expr         min          lq        mean      median
 # library(ComplexHeatmap)  200917.570  200917.570  200917.570  200917.570
 #    library(matrixStats)     125.257     125.257     125.257     125.257
 #     library(GetoptLong)   95686.390   95686.390   95686.390   95686.390
 #       library(circlize)  120574.197  120574.197  120574.197  120574.197
 #           library(clue)   18897.299   18897.299   18897.299   18897.299
 #       library(parallel)     229.346     229.346     229.346     229.346
 #   library(RColorBrewer)   19534.798   19534.798   19534.798   19534.798
 #        library(cluster)   12409.542   12409.542   12409.542   12409.542
 #            library(NMF) 1627012.578 1627012.578 1627012.578 1627012.578
 #         library(mclust)   75843.475   75843.475   75843.475   75843.475
 #        library(skmeans)  128632.934  128632.934  128632.934  128632.934
 #         library(cclust)   52636.440   52636.440   52636.440   52636.440
 #            library(png)   23425.893   23425.893   23425.893   23425.893
 #      library(venneuler) 1000763.414 1000763.414 1000763.414 1000763.414
 #         library(gplots)   90659.674   90659.674   90659.674   90659.674
 #          library(Rtsne)   67219.265   67219.265   67219.265   67219.265
 #           library(samr)  113343.314  113343.314  113343.314  113343.314
 #     library(genefilter) 6105952.052 6105952.052 6105952.052 6105952.052
 #        library(kohonen)   62767.785   62767.785   62767.785   62767.785
 #         library(crayon)   70653.832   70653.832   70653.832   70653.832
 #           library(xml2)  220768.711  220768.711  220768.711  220768.711
 #           library(Rcpp)   21191.630   21191.630   21191.630   21191.630
 #          library(wCorr)  113887.028  113887.028  113887.028  113887.028
 #          library(knitr)   89526.041   89526.041   89526.041   89526.041
 #       library(markdown)   64137.355   64137.355   64137.355   64137.355
 #      library(data.tree) 1435999.041 1435999.041 1435999.041 1435999.041
 #     library(dendextend) 2549086.659 2549086.659 2549086.659 2549086.659
 #         library(digest)   59219.129   59219.129   59219.129   59219.129
 #         library(impute)     373.211     373.211     373.211     373.211
 #           library(httr)   70131.480   70131.480   70131.480   70131.480
 #           library(brew)   54475.327   54475.327   54475.327   54475.327
