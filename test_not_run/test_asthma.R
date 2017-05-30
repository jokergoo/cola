library(methods)
library(GetoptLong)
p = 0.8
ncore = 1
GetoptLong(
	"p=f", "0.8",
	"ncore=i", "mc.cores"
)

library(cola)

source("/home/guz/project/analysis/LiNA/asthma_cohort/script/asthma_head.R")

data = deseq2[gt[rownames(deseq2)] == "protein_coding", ]

res = run_all_consensus_partition_methods(data, top_n = c(2000, 4000, 6000), k = 2:6, p_sampling = p, 
	known_anno = data.frame(phenotype = SAMPLE$phenotype), mc.cores = ncore)

saveRDS(res, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/hipo35_asthma_rnaseq_subgroup_p@{p}.rds"))


################ old dataset ################3

text2 = 
"PID  asthma
024K yes
049K yes
094K yes
239K yes
525K yes
614K yes
116K yes
150K yes
137K no
158K no
463K no
348K no
278K no
228K no
282K no
509K no
431K no"
anno = read.table(textConnection(text2), colClass = c("character", "character"), row.names = 1, header = TRUE, stringsAsFactors = FALSE)
anno$phenotype = ifelse(anno$asthma == "yes", "Asthma", "Control")

expr_old = NULL

load("/icgc/dkfzlsdf/analysis/hipo/hipo_035/data_types/RNAseq/year0_gencode19_expression.RData")
mat = expression$count
colnames(mat) = gsub("^.*?(\\d\\d\\d[KM]).*$", "\\1", colnames(mat)) 
mat = mat[, colnames(mat) %in% rownames(anno), drop = FALSE]
colnames(mat) = paste0(colnames(mat), "_year0")
expr_old = mat

load("/icgc/dkfzlsdf/analysis/hipo/hipo_035/data_types/RNAseq/year1_gencode19_expression.RData")
mat = expression$count
colnames(mat) = gsub("^.*?(\\d\\d\\d[KM]).*$", "\\1", colnames(mat)) 
mat = mat[, colnames(mat) %in% rownames(anno), drop = FALSE]
colnames(mat) = paste0(colnames(mat), "_year1")
expr_old = cbind(expr_old, mat)

load("/icgc/dkfzlsdf/analysis/hipo/hipo_035/data_types/RNAseq/year3_gencode19_expression.RData")
mat = expression$count
colnames(mat) = gsub("^.*?(\\d\\d\\d[KM]).*$", "\\1", colnames(mat)) 
mat = mat[, colnames(mat) %in% rownames(anno), drop = FALSE]
colnames(mat) = paste0(colnames(mat), "_year3")
expr_old = cbind(expr_old, mat)

load("/icgc/dkfzlsdf/analysis/hipo/hipo_035/data_types/RNAseq/year4_gencode19_expression.RData")
mat = expression$count
colnames(mat) = gsub("^.*?(\\d\\d\\d[KM]).*$", "\\1", colnames(mat)) 
mat = mat[, colnames(mat) %in% rownames(anno), drop = FALSE]
colnames(mat) = paste0(colnames(mat), "_year4")
expr_old = cbind(expr_old, mat)

load("/icgc/dkfzlsdf/analysis/hipo/hipo_035/data_types/RNAseq/year5_gencode19_expression.RData")
mat = expression$count
colnames(mat) = gsub("^.*?(\\d\\d\\d[KM]).*$", "\\1", colnames(mat)) 
mat = mat[, colnames(mat) %in% rownames(anno), drop = FALSE]
colnames(mat) = paste0(colnames(mat), "_year5")
expr_old = cbind(expr_old, mat)

load("/icgc/dkfzlsdf/analysis/hipo/hipo_035/data_types/RNAseq/year6_gencode19_expression.RData")
mat = expression$count
colnames(mat) = gsub("^.*?(\\d\\d\\d[KM]).*$", "\\1", colnames(mat)) 
mat = mat[, colnames(mat) %in% rownames(anno), drop = FALSE]
colnames(mat) = paste0(colnames(mat), "_year6")
expr_old = cbind(expr_old, mat)

anno_old = data.frame(pid = c(paste0(rownames(anno), "_year0"), paste0(rownames(anno), "_year1"), paste0(rownames(anno), "_year3"),
	                        paste0(rownames(anno), "_year4"), paste0(rownames(anno), "_year5"), paste0(rownames(anno), "_year6")),
                      type = c(anno$phenotype, anno$phenotype, anno$phenotype, anno$phenotype, anno$phenotype, anno$phenotype),
                      stringsAsFactors = FALSE)
rownames(anno_old) = anno_old$pid
anno_old = anno_old[colnames(expr_old), ]
anno_old = data.frame(pid = gsub("^(\\d\\d\\dK).*$", "\\1", anno_old$pid),
	                  year = as.numeric(gsub("^.*year(\\d).*$", "\\1", anno_old$pid)),
	                  phenotype = anno_old$type,
	                  data = rep("old", nrow(anno_old)),
	                  stringsAsFactors = FALSE)
rownames(anno_old) = qq("@{anno_old$pid}_year@{anno_old$year}_@{anno_old$data}", collapse = FALSE)


default_param = list(
	"normalization.method" = "rpkm",
	"varianceStabilize" = 1,
	"normalizeTOGeneLength" = 0
)

source("/home/guz/project/development/ngspipeline2/lib_expression.R")
library(DESeq2)
expr = normalize.count(expr_old, method = "deseq2", gene_annotation, param = default_param)
	
l = anno_old$year == 0
data = expr[gt[rownames(expr)] == "protein_coding", l]
data = t(apply(data, 1, adjust_outlier))
res = run_all_consensus_partition_methods(data, top_n = c(2000, 4000, 6000), k = 2:6, p_sampling = p, 
	known_anno = data.frame(phenotype = anno_old$phenotype[l]), mc.cores = ncore)

saveRDS(res, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/hipo35_asthma_rnaseq_heidelberg_year0_subgroup_p@{p}.rds"))


	
l = anno_old$year > 1
data = expr[gt[rownames(expr)] == "protein_coding", l]
data = t(apply(data, 1, adjust_outlier))
res = run_all_consensus_partition_methods(data, top_n = c(2000, 4000, 6000), k = 2:6, p_sampling = p, 
	known_anno = data.frame(phenotype = anno_old$phenotype[l]), mc.cores = ncore)

saveRDS(res, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/hipo35_asthma_rnaseq_heidelberg_year3more_subgroup_p@{p}.rds"))
