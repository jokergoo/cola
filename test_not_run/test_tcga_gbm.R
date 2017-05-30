library(methods)
library(GetoptLong)
p = 0.8
ncore = 1
GetoptLong(
	"p=f", "0.8",
	"ncore=i", "mc.cores"
)

library(cola)

#############################################
### TCGA GBM
data = read.table("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/unifiedScaled.txt", header = TRUE, row.names = 1, check.names = FALSE)
data = as.matrix(data)

subtype = read.table("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/TCGA_unified_CORE_ClaNC840.txt", sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
subtype = structure(unlist(subtype[1, -(1:2)]), names = colnames(subtype)[-(1:2)])

data = data[, names(subtype)]

data = t(apply(data, 1, adjust_outlier))
res = run_all_consensus_partition_methods(data, top_n = c(2000, 4000, 6000), k = 2:6, mc.cores = ncore,
	known_anno = data.frame(subtype = subtype))

saveRDS(res, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/TCGA_subgroup_p@{p}.rds"))

# for(p in c(0.2, 0.4, 0.6, 0.8)) {
# 	cmd = qq("Rscript-3.1.2 /home/guz/project/development/cola/test_not_run/test_tcga_gbm.R --p @{p} --ncore 4")
# 	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=20G,nodes=1:ppn=4 -N TCGA_subgroup_p@{p}' '@{cmd}'")
# 	system(cmd)
# }
