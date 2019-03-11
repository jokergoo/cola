library(cola)

library(HSMMSingleCell)

data(HSMM_expr_matrix)
data(HSMM_sample_sheet)

# `HSMM_expr_matrix` is a FPKM matrix
m = adjust_matrix(log10(HSMM_expr_matrix + 1))
anno = HSMM_sample_sheet[, c("Hours", "Media", "Pseudotime", "State")]

register_NMF()

set.seed(123)
rl = run_all_consensus_partition_methods(
	m,
	mc.cores = 4, 
	anno = anno
)

saveRDS(rl, file = "~/project/development/cola_examples/HSMM_single_cell/HSMM_single_cell_subgroup.rds")
cola_report(rl, output_dir = "~/project/development/cola_examples/HSMM_single_cell/HSMM_single_cell_subgroup_cola_report")

set.seed(123)
rh = hierarchical_partition(
	m,
	mc.cores = 4,
	anno = anno
)
saveRDS(rh, file = "~/project/development/cola_examples/HSMM_single_cell/HSMM_single_cell_subgroup_hierarchical_partition.rds")
cola_report(rh, output_dir = "~/project/development/cola_examples/HSMM_single_cell/HSMM_single_cell_subgroup_hierarchical_partition_cola_report")
