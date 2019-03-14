options(showWarnCalls = TRUE, showErrorCalls = TRUE)

# root = "/home/guz"
root = "/desktop-home/guz"


library(cola)
library(GetoptLong)


library(HSMMSingleCell)

data(HSMM_expr_matrix)
data(HSMM_sample_sheet)

# `HSMM_expr_matrix` is a FPKM matrix
m = adjust_matrix(log10(HSMM_expr_matrix + 1))
anno = HSMM_sample_sheet[, c("Hours", "Media", "Pseudotime", "State")]

gt = readRDS(qq("@{root}/project/development/cola_examples/HSMM_single_cell/gene_type_gencode_v17.rds"))
m = m[gt[rownames(m)] == "protein_coding", , drop = FALSE]

register_NMF()

set.seed(123)
rl = run_all_consensus_partition_methods(
	m,
	top_n = c(1000, 2000, 3000, 4000), 
	mc.cores = 4, 
	anno = anno
)

saveRDS(rl, file = qq("@{root}/project/development/cola_examples/HSMM_single_cell/HSMM_single_cell_subgroup.rds"))
cola_report(rl, output_dir = qq("@{root}/project/development/cola_examples/HSMM_single_cell/HSMM_single_cell_subgroup_cola_report"), mc.cores = 4)

set.seed(123)
rh = hierarchical_partition(
	m,
	top_n = c(1000, 2000, 3000, 4000), 
	mc.cores = 4,
	anno = anno
)
saveRDS(rh, file = qq("@{root}/project/development/cola_examples/HSMM_single_cell/HSMM_single_cell_subgroup_hierarchical_partition.rds"))
cola_report(rh, output_dir = qq("@{root}/project/development/cola_examples/HSMM_single_cell/HSMM_single_cell_subgroup_hierarchical_partition_cola_report"), mc.cores = 4)

