library(cola)

library(golubEsets)
data(Golub_Merge)
m = exprs(Golub_Merge)
colnames(m) = paste0("sample_", colnames(m))
anno = pData(Golub_Merge)

m[m <= 1] = NA
m = log10(m)

m = adjust_matrix(m)

library(preprocessCore)
cn = colnames(m)
m = normalize.quantiles(m)
colnames(m) = cn

set.seed(123)
rl = run_all_consensus_partition_methods(
	m,
	mc.cores = 4, 
	anno = anno[, c("ALL.AML"), drop = FALSE],
	anno_col = c("ALL" = "red", "AML" = "blue")
)

saveRDS(rl, file = "~/project/development/cola_examples/Golub_leukemia/Golub_leukemia_subgroup.rds")
cola_report(rl, output_dir = "~/project/development/cola_examples/Golub_leukemia/Golub_leukemia_subgroup_cola_report")

set.seed(123)
rh = hierarchical_partition(
	m,
	top_value_method = "ATC",
	partition_method = "skmeans",
	mc.cores = 4,
	anno = anno[, c("ALL.AML"), drop = FALSE],
	anno_col = c("ALL" = "red", "AML" = "blue")
)
saveRDS(rh, file = "~/project/development/cola_examples/Golub_leukemia/Golub_leukemia_subgroup_hierarchical_partition.rds")
cola_report(rh, output_dir = "~/project/development/cola_examples/Golub_leukemia/Golub_leukemia_subgroup_hierarchical_partition_cola_report")
