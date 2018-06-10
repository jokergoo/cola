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
data = read.table("~/cola_test/unifiedScaled.txt", header = TRUE, row.names = 1, check.names = FALSE)
data = as.matrix(data)

subtype = read.table("~/cola_test/TCGA_unified_CORE_ClaNC840.txt", sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
subtype = structure(unlist(subtype[1, -(1:2)]), names = colnames(subtype)[-(1:2)])

data = data[, names(subtype)]

data = adjust_matrix(data)
res = run_all_consensus_partition_methods(data, top_n = c(1000, 2000, 3000, 4000), max_k = 6, mc.cores = ncore,
	anno = data.frame(subtype = subtype), 
	anno_col = list(subtype = structure(seq_len(4), names = unique(subtype))))

saveRDS(res, file = qq("~/cola_test/TCGA_subgroup_p@{p}.rds"))

res = hierarchical_partition(data, top_n = c(1000, 2000, 3000, 4000), top_value_method = "AAC", partition_method = "skmeans",
	anno = data.frame(subtype = subtype), 
	anno_col = list(subtype = structure(seq_len(4), names = unique(subtype))))
saveRDS(res, file = qq("~/cola_test/TCGA_subgroup_hierarchical_partition.rds"))

# for(p in c(0.2, 0.4, 0.6,      0.8)) {
# 	cmd = qq("Rscript-3.3.1 /home/guz/project/development/cola/test_not_run/test_tcga_gbm.R --p @{p} --ncore 4")
# 	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=20G,nodes=1:ppn=4 -N TCGA_subgroup_p@{p}' '@{cmd}'")
# 	system(cmd)
# }
