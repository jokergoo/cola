library(GetoptLong)

tm = c("sd", "cv", "MAD", "AAC")
pm = c("hclust", "kmeans", "skmeans", "pam", "mclust", "som")
GetoptLong("pid=s", "pid",
	"tm=s{1,}", "sd|cv|MAD|AAC",
	"pm=s{1,}", "hclust|kmeans|skmeans|pam|mclust|som")

library(SummarizedExperiment)

load(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}/@{pid}_rse_gene.Rdata"))

count = assays(rse_gene)$counts
df = rowData(rse_gene)

rpkm_normalize = function(count, gene_length) {
	all_count = colSums(count)

	rpkm = matrix(0, nrow = nrow(count), ncol = ncol(count))
	rownames(rpkm) = rownames(count)
	colnames(rpkm) = colnames(count)
	for(i in seq_len(nrow(count))) {
		rpkm[i, ] = count[i, ] / gene_length[i]
	}

	for(j in seq_len(ncol(count))) {
		rpkm[, j] = rpkm[, j] / all_count[j]
	}

	rpkm = 10^9 * rpkm
	return(rpkm)
}
rpkm = rpkm_normalize(count, df$bp_length)

load("/icgc/dkfzlsdf/analysis/B080/guz/gencode/gencode_v25_transcript_merged.RData")

gt = sapply(gene_annotation$gtf, function(x) x$type)

mat = log2(rpkm[gt[rownames(rpkm)] == "protein_coding", ] + 1)

library(cola)
mat = adjust_matrix(mat)

res_list = run_all_consensus_partition_methods(mat, k = 2:8, mc.cores = 4, top_method = tm, partition_method = pm)
saveRDS(res_list, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}/@{pid}_cola_all.rds"))
res_hc = hierarchical_partition(mat)
saveRDS(res_hc, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/@{pid}/@{pid}_cola_hc.rds"))

# for(pid in dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/recount2/", pattern = "\\d{4,}")) {
# 	cmd = qq("Rscript-3.3.1 /home/guz/project/development/cola/test_not_run/test_recount.R --pid @{pid}")
# 	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=50:00:00,mem=10G,nodes=1:ppn=4 -N cola_recount_@{pid}' '@{cmd}'")
# 	system(cmd)
# }


