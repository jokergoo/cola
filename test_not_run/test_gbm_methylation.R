library(methods)
library(GetoptLong)
ncore = 1
GetoptLong(
	"ncore=i", "mc.cores"
)

library(ComplexHeatmap)
library(matrixStats)
library(circlize)
library(RColorBrewer)
library(GetoptLong)
library(GenomicRanges)

# library(cola)
source("/home/guz/project/development/cola/load.R")
register_top_value_fun(AAC = function(mat) AAC(t(mat), mc.cores = ncore))

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data("IlluminaHumanMethylation450kanno.ilmn12.hg19", 
	package = "IlluminaHumanMethylation450kanno.ilmn12.hg19")
probe = IlluminaHumanMethylation450kanno.ilmn12.hg19 # change to a short name

setwd("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/")

library(GEOquery)
if(file.exists("GSE36278_450K.RData")) {
	load("GSE36278_450K.RData")
} else {
	gset = getGEO("GSE36278")
	save(gset, file = "GSE36278_450K.RData")
}

mat = exprs(gset[[1]])
colnames(mat) = phenoData(gset[[1]])@data$title
mat = mat[rownames(getAnnotation(probe, what = "Locations")), ]

l = getAnnotation(probe, what = "Locations")$chr %in% paste0("chr", 1:22) & 
	is.na(getAnnotation(probe, what = "SNPs.137CommonSingle")$Probe_rs)
mat = mat[l, ]

mat1 = as.matrix(mat[, grep("GBM", colnames(mat))])   # tumor samples
colnames(mat1) = gsub("GBM", "dkfz", colnames(mat1))

phenotype = read.table("/icgc/dkfzlsdf/analysis/B080/guz/ComplexHeatmap_test/450K_annotation.txt", header = TRUE, sep = "\t", row.names = 1, 
	check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
phenotype = phenotype[colnames(mat1), 1:2]
colnames(phenotype) = c("dkfz_subtype", "tcga_subtype")

mat1[is.na(mat1)] = runif(sum(is.na(mat1)))
data = adjust_matrix(mat1)
register_partition_fun(
	hclust = function(mat, k, ...) {
		hc = hclust(d = dist(t(mat)), ...)
		cutree(hc, k)
	},
	kmeans = function(mat, k, ...) {
		kmeans(t(mat), centers = k, ...)
	},
	skmeans = function(mat, k, ...) {
		skmeans(x = t(mat), k = k, ...)
	},
 
	pam = function(mat, k, ...) {
		pam(t(mat), k = k, ...)
	},
	mclust = function(mat, k, ...) {
		pca = prcomp(t(mat))
		Mclust(pca$x[, 1:3], G = k, ...)$classification
	},
	som = function(mat, k, ...) {
		kr = floor(sqrt(ncol(mat)))
		somfit = som(t(mat), grid = somgrid(kr, kr, "hexagonal"), ...)
		cl = cutree(hclust(dist(somfit$codes[[1]])), k)
		group = numeric(ncol(mat))
		for(cl_unique in unique(cl)) {
			ind = which(cl == cl_unique)
			l = somfit$unit.classif %in% ind
			group[l] = cl_unique
		}
		group
	},
	scale_method = "none"
)

res = run_all_consensus_partition_methods(data, top_n = c(5000, 10000, 15000, 20000), k = 2:8, p_sampling = 0.8, 
	known_anno = phenotype, 
	known_col = list(dkfz_cluster = structure(names = c("IDH", "K27", "G34", "RTK I PDGFRA", "Mesenchymal", "RTK II Classic"), brewer.pal(6, "Set1")),
        tcga_cluster = structure(names = c("G-CIMP+", "Cluster #2", "Cluster #3"), brewer.pal(3, "Set1"))), 
	mc.cores = ncore)

saveRDS(res, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GBM_450K_subgroup.rds"))

# cmd = qq("Rscript-3.3.1 /home/guz/project/development/cola/test_not_run/test_gbm_methylation.R --ncore 4")
# cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=25G,nodes=1:ppn=4 -N gbm_methylation_subgroup' '@{cmd}'")
# system(cmd)
