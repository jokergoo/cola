library(cola)

library(ALL)
data(ALL)

m = exprs(ALL)
anno = pData(ALL)
anno = anno[, c("sex", "age", "BT", "citog", "mol.biol")]

m = adjust_matrix(m)

library(preprocessCore)
cn = colnames(m)
m = normalize.quantiles(m)
colnames(m) = cn

set.seed(123)
rl = run_all_consensus_partition_methods(
	m,
	mc.cores = 4, 
	anno = anno
)
saveRDS(rl, file = "~/project/development/cola_examples/Ritz_ALL/Ritz_ALL_subgroup.rds")
cola_report(rl, output_dir = "~/project/development/cola_examples/Ritz_ALL/Ritz_ALL_subgroup_cola_report")

set.seed(123)
rh = hierarchical_partition(
	m,
	mc.cores = 4,
	anno = anno
)
saveRDS(rh, file = "~/project/development/cola_examples/Ritz_ALL/Ritz_ALL_subgroup_hierarchical_partition.rds")
cola_report(rh, output_dir = "~/project/development/cola_examples/Ritz_ALL/Ritz_ALL_subgroup_hierarchical_partition_cola_report")
