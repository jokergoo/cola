 
library(methods)
library(GetoptLong)

datatype = "cell01"
GetoptLong(
	"datatype=s", "cell01"
)

# library(cola)
source("/home/guz/project/development/cola/load.R")
# register_top_value_fun(AAC = function(mat) AAC(t(mat), cor_method = "spearman", mc.cores = 4))

source("/home/guz/project/analysis/hipo15/script_for_repo/hipo15_lib.R")

setwd("/icgc/dkfzlsdf/analysis/hipo/hipo_015/hipo15_rnaseq_cell_analysis")
default_param = list(
	"normalization.method" = "rpkm",
	"varianceStabilize" = 1,
	"normalizeTOGeneLength" = 0
)

load("/icgc/dkfzlsdf/analysis/B080/guz/gencode/gencode_v19_lincRNA_transcript_merged.RData")

add_symbol_and_type = function(mat) {
	rn = rownames(mat)
	name = sapply(gene_annotation$gtf[rn], function(x) x$name)
	type = sapply(gene_annotation$gtf[rn], function(x) x$type)
	mat = data.frame(mat, name = name, type = type, stringsAsFactors = FALSE)
	return(mat)
}

library(DESeq2)
library(ComplexHeatmap)
library(circlize)

# go/pathway analysis
deseq_analysis2 = function(count, expr, df, formula, main, cutoff = 0.01) {
	qqcat("@{main}\n")
	dds = DESeqDataSetFromMatrix(countData = count, colData = df, design = formula)
	dd2 = DESeq(dds)
	res = results(dd2)
	res.table = as.data.frame(res)
	res.table = add_symbol_and_type(res.table)

	write.csv(res.table, file = qq("@{main}_deseq2_diffgene_all_genes.csv"))

	res.table = res.table[!is.na(res.table$padj), ]

	## fold change vs base mean
	max = max(abs(res.table$log2FoldChange))
	pdf(qq("@{main}_fc_vs_basemean_padj_@{cutoff}.pdf"), height = 6, width = 8)
	plot(log2(res.table$baseMean+1), res.table$log2FoldChange, pch = 16, cex = 0.5, ylim = c(-max, max),
		col = ifelse(res.table$padj < cutoff & res.table$log2FoldChange > 1, "#FF000080", ifelse(res.table$padj < cutoff & res.table$log2FoldChange < -1, "#00FF0080", "#00000080")),
		xlab = "log2(base_mean + 1)", ylab = "log2(fold change)", main = qq("@{main}, padj < @{cutoff}"))
	dev.off()

	pdf(qq("@{main}_deseq2_diffgene_padj_@{cutoff}_diff_gene_heatmap.pdf"), height = 12, width = 8)
	mat = expr[rownames(res.table[res.table$padj < cutoff, , drop = FALSE]), , drop = FALSE]
	dn = dimnames(mat)
	if(ncol(mat) > 2) {
		mat = t(apply(mat, 1, scale))
		dimnames(mat) = dn
	}
	column_dend = as.dendrogram(hclust(dist(t(mat))))
	column_dend = stats:::reorder.dendrogram(column_dend, as.numeric(factor(df$subgroup), function(w) mean(ifelse(w == 2, 0.1, w))))
	ht = Heatmap(mat, name = "scaled_deseq2", top_annotation = HeatmapAnnotation(df = df, col = list(subgroup = c("1" = "red", "2" = "blue"))), 
		col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), show_row_names = FALSE, 
		cluster_column = column_dend, column_dend_reorder = FALSE,
		column_title = qq("heatmap for differentially expressed genes (deseq2)\ncutoff: padj < @{cutoff}, all @{sum(res.table$padj < cutoff)} genes"))
	draw(ht)

	mat = expr[rownames(res.table[res.table$padj < cutoff & res.table$type == "protein_coding", ]), ]
	dn = dimnames(mat)
	if(ncol(mat) > 2) {
		mat = t(apply(mat, 1, scale))
		dimnames(mat) = dn
	}
	column_dend = as.dendrogram(hclust(dist(t(mat))))
	column_dend = stats:::reorder.dendrogram(column_dend, as.numeric(factor(df$subgroup), function(w) mean(ifelse(w == 2, 0.1, w))))
	ht = Heatmap(mat, name = "scaled_deseq2", top_annotation = HeatmapAnnotation(df = df, col = list(subgroup = c("1" = "red", "2" = "blue"))), 
		col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), show_row_names = FALSE, 
		cluster_column = column_dend, column_dend_reorder = FALSE,
		column_title = qq("heatmap for differentially expressed genes (protein_coding) (deseq2)\ncutoff: padj < @{cutoff}, all @{sum(res.table$padj < cutoff & res.table$type == 'protein_coding')} genes"))
	draw(ht)
	dev.off()
}

if(grepl("cell", datatype)) {
	source("/home/guz/project/analysis/hipo15/script_for_repo/head.R")
	deseq_analysis = deseq_analysis2
	if(datatype == "cell01") {
		res_list = readRDS("hipo15_c1_subgroups.rds")
		res = get_single_run(res_list, "AAC", "skmeans")
		count = expression$count[, colnames(res@.env$data)]
		expr = expression$deseq2[, colnames(count)]
	} else if(datatype == "cell02") {
		res_list = readRDS("hipo15_c2_subgroups.rds")
		res = get_single_run(res_list, "AAC", "skmeans")
		count = expression$count[, colnames(res@.env$data)]
		expr = expression$deseq2[, colnames(count)]
	} else if(datatype == "cell03") {
		res_list = readRDS("hipo15_c3_subgroups.rds")
		res = get_single_run(res_list, "AAC", "skmeans")
		count = expression$count[, colnames(res@.env$data)]
		expr = expression$deseq2[, colnames(count)]
	}
} else if(datatype == "primary_tumor") {

	########################### primary tumor #################################
deseq_analysis = deseq_analysis2
	load("/icgc/dkfzlsdf/analysis/hipo/hipo_015/data_types/RNAseq/expression_data/hipo15_rnaseq_primary_tumor_gencode19_lincRNA_expression.RData")
	load("/home/guz/project/analysis/hipo15/signature_genes/moffitt_sigANDannotation.RData")
	source("/home/guz/project/development/ngspipeline2/lib_expression.R")

	res_list = readRDS("hipo15_primary_tumor_subgroups.rds")
	res = get_single_run(res_list, "AAC", "skmeans")
	count = expression$count[, colnames(res@.env$data)]
	expr = normalize.count(count, method = "deseq2", gene_annotation, param = default_param)

} else if(datatype == "xenograft") {
	######################## xenograft #########################
deseq_analysis = deseq_analysis2
	load("/icgc/dkfzlsdf/analysis/hipo/hipo_015/data_types/RNAseq/expression_data/rnaseq_xenograft_human_v19_lincRNA_mouse_M2_expression.RData")
	load("/icgc/dkfzlsdf/analysis/B080/guz/gencode/gencode_v19_lincRNA_transcript_merged.RData")
	source("/home/guz/project/development/ngspipeline2/lib_expression.R")
	gene_type = sapply(gene_annotation$gtf, function(x) x$type)

	count = expression$count
	l = !grepl("^ENSM", rownames(count))
	count = count[l, ]

	res_list = readRDS("hipo15_xenograft_subgroups.rds")
	res = get_single_run(res_list, "AAC", "skmeans")
	count = count[, colnames(res@.env$data)]
	expr = normalize.count(count, method = "deseq2", gene_annotation, param = default_param)
}

cl = get_class(res, k = 2)
anno_df = data.frame(subgroup = as.character(cl$class))
deseq_analysis(count, expr, anno_df, ~ subgroup, qq("@{datatype}_2groups"), cutoff = 0.05)


# for(datatype in c("cell01", "cell02", "cell03", "primary_tumor", "xenograft")) {
# 	cmd = qq("Rscript-3.3.1 /home/guz/project/development/cola/test_not_run/hipo15_subgroups_deseq2.R --datatype @{datatype}")
#  	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=10:00:00,mem=10G -N hipo15_subgroups_@{datatype}_deseq2' '@{cmd}'")
#  	system(cmd)
# }
