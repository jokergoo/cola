library(methods)
library(GetoptLong)
ncore = 1
GetoptLong(
	"ncore=i", "mc.cores"
)

source(qq("@{dirname(get_scriptname())}/../load.R"))
# library(cola)
register_top_value_fun(AAC = function(mat) AAC(t(mat), mc.cores = ncore))

# load("/icgc/dkfzlsdf/analysis/cnag/cnag_MCF10CA_scRNAseq_gencode19_expression.RData")
load("/icgc/dkfzlsdf/analysis/hipo/hipo_043/RNAseq/hipo_043_rnaseq_gencode19_expression.RData")

library(GenomicFeatures)
TXDB = loadDb(qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/txdb/gencode19_protein_coding_txdb.sqlite"))
GENE = genes(TXDB)

rpkm = as.matrix(expression$rpkm)
count = as.matrix(expression$count)

rpkm = rpkm[names(GENE), ]
count = count[names(GENE), ]
l = apply(count, 1, function(x) sum(x > 0)/length(x) > 0.5)
data = log10(rpkm[l, ] + 1)

data = adjust_matrix(data)
res = run_all_consensus_partition_methods(data, mc.cores = ncore)

saveRDS(res, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/hipo43_rnaseq_subgroup.rds"))

	# cmd = qq("Rscript-3.3.1 /home/guz/project/development/cola/test_not_run/test_hipo43.R --ncore 4")
	# cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=20G,nodes=1:ppn=4 -N hipo43_rnaseq_subgroup' '@{cmd}'")
	# system(cmd)
