if(grepl("tbi", Sys.info()["nodename"])) {
	root = "/home/guz"
} else {
	root = "/desktop-home/guz"
}

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

register_NMF()

set.seed(123)
rl = run_all_consensus_partition_methods(
	m,
	top_n = c(1000, 2000, 3000, 4000), 
	mc.cores = 4, 
	anno = anno[, c("ALL.AML"), drop = FALSE],
	anno_col = c("ALL" = "red", "AML" = "blue")
)

saveRDS(rl, file = qq("@{root}/project/development/cola_examples/Golub_leukemia/Golub_leukemia_subgroup.rds"))
cola_report(rl, output_dir = qq("@{root}/project/development/cola_examples/Golub_leukemia/Golub_leukemia_subgroup_cola_report"), mc.cores = 4)

set.seed(123)
rh = hierarchical_partition(
	m,
	top_n = c(1000, 2000, 3000, 4000), 
	top_value_method = "ATC",
	partition_method = "skmeans",
	mc.cores = 4,
	anno = anno[, c("ALL.AML"), drop = FALSE],
	anno_col = c("ALL" = "red", "AML" = "blue")
)
saveRDS(rh, file = qq("@{root}/project/development/cola_examples/Golub_leukemia/Golub_leukemia_subgroup_hierarchical_partition.rds"))
cola_report(rh, output_dir = qq("@{root}/project/development/cola_examples/Golub_leukemia/Golub_leukemia_subgroup_hierarchical_partition_cola_report"), mc.cores = 4)
