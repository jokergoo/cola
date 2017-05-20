library(methods)
library(GetoptLong)

datatype = "cell01"
GetoptLong(
	"datatype=s", "cell01"
)

source("/home/guz/project/development/cola/load.R")

setwd("/icgc/dkfzlsdf/analysis/hipo/hipo_015/hipo15_rnaseq_cell_analysis")
default_param = list(
	"normalization.method" = "rpkm",
	"varianceStabilize" = 1,
	"normalizeTOGeneLength" = 0
)

if(grepl("cell", datatype)) {
	source("/home/guz/project/analysis/hipo15/script_for_repo/head.R")
	gene_type = sapply(gene_annotation$gtf, function(x) x$type)
	gene_name = sapply(gene_annotation$gtf, function(x) x$name)


	if(datatype == "cell01") {
		l = anno$type %in% c("cell01") & anno$phenotype == "tumor"
		data = expr[gene_type[rownames(expr)] == "protein_coding", l]
		res = run_all(data, top_n = c(2000, 4000, 6000), k = 2:6, mc.cores = 4)

		saveRDS(res, file = "hipo15_c1_subgroups.rds")
	} else if(datatype == "cell02") {
		l = anno$type %in% c("cell02") & anno$phenotype == "tumor"
		data = expr[gene_type[rownames(expr)] == "protein_coding", l]
		res = run_all(data, top_n = c(2000, 4000, 6000), k = 2:6, mc.cores = 4)

		saveRDS(res, file = "hipo15_c2_subgroups.rds")
	} else if(datatype == "cell03") {

		l = anno$type %in% c("cell03") & anno$phenotype == "tumor"
		data = expr[gene_type[rownames(expr)] == "protein_coding", l]
		res = run_all(data, top_n = c(2000, 4000, 6000), k = 2:6, mc.cores = 4)

		saveRDS(res, file = "hipo15_c3_subgroups.rds")
	}
} else if(datatype == "primary_tumor") {

	########################### primary tumor #################################

	load("/icgc/dkfzlsdf/analysis/hipo/hipo_015/data_types/RNAseq/expression_data/hipo15_rnaseq_primary_tumor_gencode19_lincRNA_expression.RData")
	load("/icgc/dkfzlsdf/analysis/B080/guz/gencode/gencode_v19_lincRNA_transcript_merged.RData")
	load("/home/guz/project/analysis/hipo15/signature_genes/moffitt_sigANDannotation.RData")
	source("/home/guz/project/development/ngspipeline2/lib_expression.R")
	library(DESeq2)

	gene_type = sapply(gene_annotation$gtf, function(x) x$type)
	expr = normalize.count(expression$count, method = "deseq2", gene_annotation, param = default_param)
	data = expr[gene_type[rownames(expr)] == "protein_coding", rownames(pro_5types_50)]
	res = run_all(data, top_n = c(2000, 4000, 6000), k = 2:6, mc.cores = 4)

	saveRDS(res, file = "hipo15_primary_tumor_subgroups.rds")
} else if(datatype == "xenograft") {
	######################## xenograft #########################

	load("/icgc/dkfzlsdf/analysis/hipo/hipo_015/data_types/RNAseq/expression_data/rnaseq_xenograft_human_v19_lincRNA_mouse_M2_expression.RData")
	load("/icgc/dkfzlsdf/analysis/B080/guz/gencode/gencode_v19_lincRNA_transcript_merged.RData")
	source("/home/guz/project/development/ngspipeline2/lib_expression.R")
	gene_type = sapply(gene_annotation$gtf, function(x) x$type)

	library(DESeq2)
	count = expression$count
	count = count[!grepl("^ENSM", rownames(count)), ]
	expr = normalize.count(count, method = "deseq2", gene_annotation, param = default_param)
	expr = expr[gene_type[rownames(expr)] == "protein_coding", ]
	res = run_all(expr, top_n = c(2000, 4000, 6000), k = 2:6, mc.cores = 4, partition_repeat = 10)

	saveRDS(res, file = "hipo15_xenograft_subgroups.rds")
}

# for(datatype in c("cell01", "cell02", "cell03", "primary_tumor", "xenograft")) {
# 	cmd = qq("Rscript-3.1.2 /home/guz/project/development/subgroup/hipo15_subgroup.R --datatype @{datatype}")
#  	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=10G,nodes=1:ppn=4 -N hipo15_subgroup_@{datatype}' '@{cmd}'")
#  	system(cmd)
# }
