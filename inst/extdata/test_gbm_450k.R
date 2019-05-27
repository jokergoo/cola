library(GetoptLong)

library(ComplexHeatmap)
library(matrixStats)
library(circlize)
library(RColorBrewer)
library(GenomicRanges)

library(cola)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data("IlluminaHumanMethylation450kanno.ilmn12.hg19", 
	package = "IlluminaHumanMethylation450kanno.ilmn12.hg19")
probe = IlluminaHumanMethylation450kanno.ilmn12.hg19 # change to a short name

setwd("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/")

library(GEOquery)
if(file.exists("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GSE36278_450K.RData")) {
	load("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GSE36278_450K.RData")
} else {
	gset = getGEO("GSE36278")
	save(gset, file = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GSE36278_450K.RData")
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

anno_col = list(dkfz_subtype = structure(names = c("IDH", "K27", "G34", "RTK I PDGFRA", "Mesenchymal", "RTK II Classic"), brewer.pal(6, "Set1")),
        tcga_subtype = structure(names = c("G-CIMP+", "Cluster #2", "Cluster #3"), brewer.pal(3, "Set1")))

rl = run_all_consensus_partition_methods(
	mat1, 
	top_n = c(5000, 10000, 15000, 20000), 
	max_k = 8,
	scale_rows = FALSE, anno = phenotype, 
	anno_col = anno_col, 
	mc.cores = 4)

saveRDS(rl, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GBM_450K_subgroup.rds"))
cola_report(rl, output_dir = qq("@{root}/project/development/cola_examples/GBM_450K/GBM_450K_subgroup_cola_report"), mc.cores = 4)


rh = hierarchical_partition(
	mat1, 
	top_n = c(5000, 10000, 15000, 20000),
	scale_rows = FALSE, 
	anno = phenotype, 
	top_method = "MAD", 
	partition_method = "kmeans",
	anno_col = anno_col,
	mc.cores = 4 
)

saveRDS(rh, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/GBM_450K_subgroup_hierarchical_partition.rds"))
cola_report(rh, output_dir = qq("@{root}/project/development/cola_examples/GBM_450K/GBM_450K_subgroup_hierarchical_partition_cola_report"), mc.cores = 4)

