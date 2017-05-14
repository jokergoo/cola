
source("/home/guz/project/development/subgroup/subgroup_by_correlation.R")

source("/home/guz/project/analysis/hipo15/script_for_repo/head.R")
setwd("/icgc/dkfzlsdf/analysis/hipo/hipo_015/hipo15_rnaseq_cell_analysis")
gene_type = sapply(gene_annotation$gtf, function(x) x$type)
gene_name = sapply(gene_annotation$gtf, function(x) x$name)

l = anno$type %in% "cell01" & anno$phenotype == "tumor"
data = expr[gene_type[rownames(expr)] == "protein_coding", l]
res = consensus_partition(data, k = 2:8, n_bg = c(2000, 4000, 6000))
res = consensus_partition(data, k = 2:8, n_bg = c(2000, 4000, 6000), method = "top_var")
res = consensus_partition(data, k = 2:8, n_bg = c(2000, 4000), method = "NMF")

for(cell_type in c("cell01", "cell02", "cell03")) {

	pdf(qq("tumor_@{cell_type}_subgroups.pdf"), width = 8, height = 8)
	l = anno$type %in% cell_type & anno$phenotype == "tumor"
	data = expr[gene_type[rownames(expr)] == "protein_coding", l]
	res = consensus_partition(data, k = 2:8, n_bg = c(2000, 4000, 6000))
	consensus_heatmap(res, k = 2)
	consensus_heatmap(res, k = 3)
	consensus_heatmap(res, k = 4)
	select_k(res)
	mat2 = get_signatures(res, k = 2)
	mat3 = get_signatures(res, k = 3)
	mat4 = get_signatures(res, k = 4)
	gplots::venn(list("k=2" = rownames(mat2), "k=3" = rownames(mat3), "k=4" = rownames(mat4)))
	title("signatures under different k")
	dev.off()

	write.csv(cbind(mat2, gene_symbol = gene_name[rownames(mat2)]), file = qq("tumor_@{cell_type}_2_subgroups_signatures.csv"))
	write.csv(cbind(mat3, gene_symbol = gene_name[rownames(mat3)]), file = qq("tumor_@{cell_type}_3_subgroups_signatures.csv"))
	write.csv(cbind(mat4, gene_symbol = gene_name[rownames(mat4)]), file = qq("tumor_@{cell_type}_4_subgroups_signatures.csv"))

	saveRDS(res, file = qq("tumor_@{cell_type}_subgroups.rds"))
}


######################################################################
library(TCGA2STAT)

rnaseq.gbm = getTCGA(disease = "GBM", data.type = "RNASeq2", type = "RPKM", clinical = TRUE)
data = rnaseq.gbm$merged.dat
rownames(data) = data[, 1]
data = t(as.matrix(data[, -(1:3)]))
data = log10(data + 1)
pdf(qq("~/tcga_gbm_rnaseq_subgroups.pdf"), width = 8, height = 8)
res = consensus_partition(data, k = 2:8, n_bg = c(2000, 4000, 6000), method = "top_var")
consensus_heatmap(res, k = 2)
consensus_heatmap(res, k = 3)
consensus_heatmap(res, k = 4)
select_k(res)
mat2 = get_signatures(res, k = 2)
mat3 = get_signatures(res, k = 3)
mat4 = get_signatures(res, k = 4)
gplots::venn(list("k=2" = rownames(mat2), "k=3" = rownames(mat3), "k=4" = rownames(mat4)))
title("signatures under different k")
dev.off()
saveRDS(res, file = qq("~/tcga_gbm_rnaseq_subgroups.rds"))

array.gbm = getTCGA(disease = "GBM", data.type = "mRNA_Array", type = "G450", clinical = TRUE)
data = array.gbm$merged.dat
data = data[!duplicated(data[, 1]), ]
rownames(data) = data[, 1]
data = t(as.matrix(data[, -(1:3)]))
data = t(apply(data, 1, function(x) {
	m = mean(x, na.rm = TRUE)
	x[is.na(x)] = m
	return(x)
}))
res = consensus_partition(data, k = 2:8, n_bg = c(2000, 4000, 6000))
res = consensus_partition(data, k = 2:8, n_bg = 16000, method = "top_var", p = 0.6, partition_repeat = 100, order_by = "mad")
res = consensus_partition(data, k = 2:8, n_bg = 16000, method = "NMF", p = 0.6, partition_repeat = 10, order_by = "mad")
res = consensus_partition(data, k = 2:8, n_bg = c(2000, 4000, 6000), n_cor_cutoff = c(20, 50, 100, 150, 200), order_by = "mad")

data = read.table("~/Broad202.txt", header = TRUE, row.names = 1)
data = as.matrix(data)
data = log10(data)
res = consensus_partition(data, k = 2:8, top_n = c(2000, 4000, 6000), p = 0.8, top_method = "sd", partition_method = "kmeans")
res = consensus_partition(data, k = 2:8, top_n = c(2000, 4000, 6000), p = 0.8, top_method = "AUC", partition_method = "kmeans")
res = consensus_partition(data, k = 2:8, top_n = c(2000, 4000, 6000), p = 0.8, top_method = "sd", partition_method = "pam")
res = consensus_partition(data, k = 2:8, top_n = c(2000, 4000, 6000), p = 0.8, top_method = "sd", 
	partition_method = "skmeans", partition_repeat = 5)

consensus_heatmap(res, k = 2)
consensus_heatmap(res, k = 3)
consensus_heatmap(res, k = 4)
select_k(res)
mat2 = get_signatures(res, k = 2)
mat3 = get_signatures(res, k = 3)
mat4 = get_signatures(res, k = 4)
gplots::venn(list("k=2" = rownames(mat2), "k=3" = rownames(mat3), "k=4" = rownames(mat4)))
title("signatures under different k")


methyl45.gbm <- getTCGA(disease="GBM", data.type="Methylation", type="450K")
l = methyl45.gbm$cpgs$Chromosome %in% as.character(1:22)
data = methyl45.gbm$dat[l, ]


############################################################
anno_text =
"id    hipo_id     subtype batch
AK015  H016-WFRL   IDH   1
AK041   H016-7EN2   IDH  1
AK066   H016-DVZSMF IDH  1
AK068   H016-1VS79M IDH  1
AK076   H016-JGR2   IDH  1
AK085   H016-GA6F   IDH  1
AK102   H016-D3VY   IDH  1
AK103   H016-XEMY   IDH  1
AK124   H016-3N2CQY IDH  1
AK199   H016-SEMSBV IDH  1
AK213   H016-C9EF5G IDH  1
AK231   H016-K82Q   IDH  1
AK005   H016-6BP868 MES  2
AK006   H016-RKB4RB MES  2
AK030   H016-AZH7   MES  1
AK055   H016-U9HSNM MES  1
AK071   H016-TELN6S MES  1
AK072   H016-STUK   MES  1
AK079   H016-N6KNX8 MES  2
AK081   H016-KG8EA4 MES  2
AK088   H016-K48U   MES  1
AK091   H016-DD22   MES  1
AK134   H016-7CCGYW MES  2
AK139   H016-V41MG6 MES  1
AK153   H016-2BH85A MES  1
AK185   H016-1SDSG7 MES  1
AK188   H016-XXJFNP MES  1
AK195   H016-77FF   MES  1
AK218   H016-82ZL4S MES  2
AK227   H016-H1M4FV MES  1
AK235   H016-7JVLAT MES  2
AK236   H016-U676   MES  1
AK256   H016-V1A93Q MES  2
AK002   H016-8S2Z4Z RTK_I  2
AK003   H016-KBJ2J5 RTK_I  1
AK043   H016-F7WG7L RTK_I  2
AK049   H016-YKZ5   RTK_I  1
AK051   H016-AYDUQX RTK_I  1
AK142   H016-6L2VZW RTK_I  1
AK149   H016-9JGR   RTK_I  1
AK156   H016-3LK6   RTK_I  1
AK165   H016-59ND   RTK_I  1
AK173   H016-9GTQ8S RTK_I  1
AK183   H016-D3H7   RTK_I  1
AK217   H016-YN4UXM RTK_I  1
AK035   H016-1AG619 RTK_II  2
AK053   H016-AZJVFM RTK_II  1
AK074   H016-4F2ZQN RTK_II  1
AK089   H016-761V   RTK_II  1
AK098   H016-BRD3LH RTK_II  1
AK099   H016-8DFN1P RTK_II  2
AK100   H016-TTQXDW RTK_II  1
AK117   H016-TCS133 RTK_II  2
AK123   H016-PNVGQB RTK_II  2
AK132   H016-L5XL   RTK_II  1
AK133   H016-9WEKX1 RTK_II  2
AK158   H016-DUHE   RTK_II  1
AK167   H016-LEZR   RTK_II  1
AK178   H016-F6J1VJ RTK_II  1
AK205   H016-K1RYMM RTK_II  2
AK216   H016-VNDF   RTK_II  1
AK226   H016-3ZCL2Y RTK_II  2
"

SAMPLE = read.table(textConnection(anno_text), header = TRUE, stringsAsFactors = FALSE)
rownames(SAMPLE) = SAMPLE$id
SAMPLE_ID = SAMPLE$id
SUBTYPE_COLOR = RColorBrewer::brewer.pal(5, "Set1")
names(SUBTYPE_COLOR) = c("IDH", "MES", "RTK_I", "RTK_II")

PROJECT_DIR = "/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/analysis_with_epik"
library(GenomicFeatures)
TXDB = loadDb(qq("@{PROJECT_DIR}/txdb/gencode19_protein_coding_txdb.sqlite"))
GENE = genes(TXDB)
GTF = qq("@{PROJECT_DIR}/txdb/gencode.v19.annotation.gtf")

####################################################
## expression data
load(qq("@{PROJECT_DIR}/expression/hipo16_rnaseq_count_rpkm.RData"))
count = count[names(GENE), SAMPLE_ID]
rpkm = rpkm[names(GENE), SAMPLE_ID]
data = log2(rpkm + 1)
res = consensus_partition(data, k = 2:8, n_bg = 16000, method = "top_var", p = 0.6, partition_repeat = 100, order_by = "mad")
res = consensus_partition(data, k = 2:8, n_bg = 16000, method = "NMF", p = 0.6, partition_repeat = 10, order_by = "mad")
res = consensus_partition(data, k = 2:8, n_bg = c(4000, 8000, 12000, 16000), n_cor_cutoff = c(20, 50, 100, 150, 200, 400, 800), order_by = "mad")

consensus_heatmap(res, k = 4, annotation = SAMPLE[, 3, drop = FALSE], 
	annotation_color = list(subtype = SUBTYPE_COLOR))
get_signatures(res, k = 2, annotation = SAMPLE[, 3, drop = FALSE], 
	annotation_color = list(subtype = SUBTYPE_COLOR))

od = order(rowSds(data), decreasing = TRUE)[1:2000]
Heatmap(data[od, ], top_annotation = HeatmapAnnotation(subtype = SAMPLE$subtype, col = list(subtype = SUBTYPE_COLOR)))


