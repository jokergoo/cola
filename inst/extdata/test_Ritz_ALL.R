options(showWarnCalls = TRUE, showErrorCalls = TRUE)

# root = "/home/guz"
root = "/desktop-home/guz"


library(cola)
library(GetoptLong)


library(ALL)
data(ALL)

m = exprs(ALL)
anno = pData(ALL)
anno = anno[, c("sex", "age", "BT")]

m = adjust_matrix(m)

library(preprocessCore)
cn = colnames(m)
m = normalize.quantiles(m)
colnames(m) = cn

register_NMF()

# set.seed(123)
# rl = run_all_consensus_partition_methods(
# 	m,
# 	top_n = c(1000, 2000, 3000, 4000), 
# 	mc.cores = 4, 
# 	anno = anno
# )
# saveRDS(rl, file = qq("@{root}/project/development/cola_examples/Ritz_ALL/Ritz_ALL_subgroup.rds"))
# cola_report(rl, output_dir = qq("@{root}/project/development/cola_examples/Ritz_ALL/Ritz_ALL_subgroup_cola_report"), mc.cores = 4)

set.seed(123)
rh = hierarchical_partition(
	m,
	top_n = c(1000, 2000, 3000, 4000), 
	mc.cores = 4,
	anno = anno
)
saveRDS(rh, file = qq("@{root}/project/development/cola_examples/Ritz_ALL/Ritz_ALL_subgroup_hierarchical_partition.rds"))
cola_report(rh, output_dir = qq("@{root}/project/development/cola_examples/Ritz_ALL/Ritz_ALL_subgroup_hierarchical_partition_cola_report"), mc.cores = 4)
