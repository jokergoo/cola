if(grepl("tbi", Sys.info()["nodename"])) {
	root = "/home/guz"
} else {
	root = "/desktop-home/guz"
}

library(cola)

rpkm = readRDS(qq("@{root}/project/development/cola_examples/MCF10CA_scRNAseq/MCF10CA_scRNAseq_rpkm.rds"))

m = log2(rpkm + 1)

cell_type = ifelse(grepl("round", colnames(rpkm)), "round", "aberrant")
cell_col = cell_type = c("aberrant" = "red", "round" = "blue")

m = adjust_matrix(m)

register_NMF()

set.seed(123)
rl = run_all_consensus_partition_methods(
	m, 
	top_n = c(1000, 2000, 3000),
	mc.cores = 4,
	anno = data.frame(cell_type = cell_type), 
	anno_col = list(cell_type = cell_col)
)

saveRDS(rl, file = qq("@{root}/project/development/cola_examples/MCF10CA_scRNAseq/MCF10CA_scRNAseq_subgroup.rds"))
cola_report(rl, output_dir = qq("@{root}/project/development/cola_examples/MCF10CA_scRNAseq/MCF10CA_scRNAseq_subgroup_cola_report"), mc.cores = 4)

set.seed(123)
rh = hierarchical_partition(
	m, 
	top_n = c(1000, 2000, 3000),
	mc.cores = 4,
	anno = data.frame(cell_type = cell_type), 
	anno_col = list(cell_type = cell_col)
)
saveRDS(rh, file = qq("@{root}/project/development/cola_examples/MCF10CA_scRNAseq/MCF10CA_scRNAseq_subgroup_hierarchical_partition.rds"))
cola_report(rh, output_dir = qq("@{root}/project/development/cola_examples/MCF10CA_scRNAseq/MCF10CA_scRNAseq_subgroup_hierarchical_partition_cola_report"), mc.cores = 4)
