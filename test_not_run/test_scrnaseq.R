library(methods)
library(GetoptLong)
p = 0.8
ncore = 1
GetoptLong(
	"p=f", "0.8",
	"ncore=i", "mc.cores"
)

source("/home/guz/project/development/cola/load.R")

# load("/icgc/dkfzlsdf/analysis/cnag/cnag_MCF10CA_scRNAseq_gencode19_expression.RData")
load("/icgc/dkfzlsdf/analysis/cnag/cnag_MCF10CA_spheroids_gencode19_expression.RData")

PROJECT_DIR = "/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/analysis_with_epik"
library(GenomicFeatures)
TXDB = loadDb(qq("@{PROJECT_DIR}/txdb/gencode19_protein_coding_txdb.sqlite"))
GENE = genes(TXDB)

rpkm = as.matrix(expression$rpkm)
count = as.matrix(expression$count)

rpkm = rpkm[names(GENE), ]
count = count[names(GENE), ]
l = apply(count, 1, function(x) sum(x > 0)/length(x) > 0.5)
data = log2(rpkm[l, ] + 1)

# data = data[, !colnames(data) %in% c("1_C50", "1_C62", "1_C87", "Bulk")]
data = data[, !colnames(data) %in% c("Invasive_1", "Invasive_10", "Invasive_11", "Invasive_14", "Round_13", "Round_3", "Round_bulk", "invasive_bulk")]
cell_type = gsub("_\\d+$", "", colnames(data))

res = run_all(data, top_n = c(2000, 4000, 6000), k = 2:6, p_sampling = p, known = cell_type, mc.cores = ncore)

saveRDS(res, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/subgroup_test/scrnaseq_subgroup_p@{p}.rds"))

# for(p in c(0.2, 0.4, 0.6, 0.8)) {
# 	cmd = qq("Rscript-3.1.2 /home/guz/project/development/subgroup/test_scrnaseq.R --p @{p} --ncore 4")
# 	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=20G,nodes=1:ppn=4 -N scrnaseq_subgroup_p@{p}' '@{cmd}'")
# 	system(cmd)
# }

