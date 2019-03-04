library(GetoptLong)
library(cola)

library(HSMMSingleCell)

data(HSMM_expr_matrix)
data(HSMM_sample_sheet)

m = adjust(log10(HSMM_expr_matrix + 1))
anno = HSMM_sample_sheet

set.seed(123)
rl = run_all_consensus_partition_methods(
	m,
	mc.cores = 4, 
	anno = anno
)

saveRDS(rl, file = qq("~/project/development/cola_examples/HSMM_single_cell_subgroup.rds"))
cola_report(rl, output_dir = "~/project/development/cola_examples/HSMM_single_cell_subgroup_cola_report")

set.seed(123)
rh = hierarchical_partition(
	m,
	top_value_method = "ATC",
	partition_method = "skmeans",
	anno = anno
)
saveRDS(rh, file = "~/project/development/cola_examples/HSMM_single_cell_subgroup_hierarchical_partition.rds")
cola_report(rh, output_dir = "~/project/development/cola_examples/HSMM_single_cell_subgroup_hierarchical_partition_cola_report")
