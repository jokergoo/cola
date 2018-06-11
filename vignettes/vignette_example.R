
library(cola)

data = read.table("~/cola_test/unifiedScaled.txt", header = TRUE, row.names = 1, check.names = FALSE)
data = as.matrix(data)

subtype = read.table("~/cola_test/TCGA_unified_CORE_ClaNC840.txt", sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
subtype = structure(unlist(subtype[1, -(1:2)]), names = colnames(subtype)[-(1:2)])

data = data[, names(subtype)]

data = adjust_matrix(data)
rl = run_all_consensus_partition_methods(data, top_n = c(1000, 2000, 3000, 4000), 
	max_k = 6, mc.cores = 4, anno = data.frame(subtype = subtype), 
	anno_col = list(subtype = structure(seq_len(4), names = unique(subtype))))
saveRDS(rl, file = "~/cola_test/TCGA_subgroup.rds")

rh = hierarchical_partition(data, top_n = c(1000, 2000, 3000, 4000),
	anno = data.frame(subtype = subtype,
		consensus = get_classes(rl, k = 4)$class), 
	anno_col = list(subtype = structure(seq_len(4), names = unique(subtype)),
		consensus = structure(brewer.pal(4, "Set1"), names = 1:4)))
saveRDS(rh, file = "~/cola_test/TCGA_subgroup_hierarchical_partition.rds")


cola_report(rl, output_dir = "~/cola_test/tcga_cola_rl_report")
cola_report(rh, output_dir = "~/cola_test/tcga_cola_rh_report")
