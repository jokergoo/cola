library(cola)
library(RColorBrewer)

m = read.table("/desktop-home/guz/project/development/cola_examples/TCGA_GBM/unifiedScaled.txt", header = TRUE, row.names = 1, check.names = FALSE)
m = as.matrix(m)

subtype = read.table("/desktop-home/guz/project/development/cola_examples/TCGA_GBM/TCGA_unified_CORE_ClaNC840.txt", sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
subtype = structure(unlist(subtype[1, -(1:2)]), names = colnames(subtype)[-(1:2)])
subtype_col = structure(seq_len(4), names = unique(subtype))

m = m[, names(subtype)]
m = adjust_matrix(m)

register_NMF()

set.seed(123)
rl = run_all_consensus_partition_methods(
	m, 
	top_n = c(1000, 2000, 3000, 4000), 
	mc.cores = 4,
	anno = data.frame(subtype = subtype), 
	anno_col = list(subtype = subtype_col)
)
saveRDS(rl, file = "/desktop-home/guz/project/development/cola_examples/TCGA_GBM/TCGA_GBM_subgroup.rds")
cola_report(rl, output_dir = "/desktop-home/guz/project/development/cola_examples/TCGA_GBM/TCGA_GBM_subgroup_cola_report")

set.seed(123)
rh = hierarchical_partition(
	m, 
	top_n = c(1000, 2000, 3000, 4000),
	mc.cores = 4, 
    anno = data.frame(subtype = subtype,
        consensus = get_classes(rl, k = 4)$class), 
    anno_col = list(subtype = subtype_col,
        consensus = structure(RColorBrewer::brewer.pal(4, "Set1"), names = 1:4))
)
saveRDS(rh, file = "/desktop-home/guz/project/development/cola_examples/TCGA_GBM/TCGA_GBM_subgroup_hierarchical_partition.rds")
cola_report(rh, output_dir = "/desktop-home/guz/project/development/cola_examples/TCGA_GBM/TCGA_GBM_subgroup_hierarchical_partition_cola_report")
