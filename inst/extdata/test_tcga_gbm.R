library(GetoptLong)
library(cola)
library(RColorBrewer)

data = read.table("~/project/development/cola_examples/unifiedScaled.txt", header = TRUE, row.names = 1, check.names = FALSE)
data = as.matrix(data)

subtype = read.table("~/project/development/cola_examples/TCGA_unified_CORE_ClaNC840.txt", sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
subtype = structure(unlist(subtype[1, -(1:2)]), names = colnames(subtype)[-(1:2)])

data = data[, names(subtype)]
data = adjust_matrix(data)

set.seed(123)
rl = run_all_consensus_partition_methods(
	data, 
	top_n = c(1000, 2000, 3000, 4000), 
	mc.cores = 4,
	anno = data.frame(subtype = subtype), 
	anno_col = list(subtype = structure(seq_len(4), names = unique(subtype)))
)
saveRDS(rl, file = "~/project/development/cola_examples/TCGA_GBM_subgroup.rds")
cola_report(rl, output_dir = "~/project/development/cola_examples/TCGA_GBM_subgroup_cola_report")

set.seed(123)
rh = hierarchical_partition(data, 
	top_n = c(1000, 2000, 3000, 4000),
    top_value_method = "MAD", 
    partition_method = "kmeans",
    anno = data.frame(subtype = subtype,
        consensus = get_classes(rl, k = 4)$class), 
    anno_col = list(subtype = structure(seq_len(4), names = unique(subtype)),
        consensus = structure(RColorBrewer::brewer.pal(4, "Set1"), names = 1:4))
)
saveRDS(rh, file = "~/project/development/cola_examples/TCGA_GBM_subgroup_hierarchical_partition.rds")
cola_report(rh, output_dir = "~/project/development/cola_examples/TCGA_GBM_subgroup_hierarchical_partition_cola_report")
