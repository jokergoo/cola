source("/home/guz/project/development/cola/load.R")

######### TCGA ###########
for(p in c(0.2, 0.4, 0.6, 0.8)) {
	res_list = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/subgroup_test/TCGA_subgroup_p@{p}.rds"))
	for(k in 2:6) {
		pdf(qq("/icgc/dkfzlsdf/analysis/B080/guz/subgroup_test/TCGA_subgroup_k@{k}_p@{p}.pdf"), width = 15, height = 10)
		for(fun in c("plot_ecdf", "consensus_heatmap", "membership_heatmap")) {
			qqcat("collecting plots for @{fun}, p@{p}, k@{k}\n")
			collect_plots(res_list, k = k, fun = get(fun))
		}
		dev.off()
	}
}

######### hipo16  ###########
for(p in c(0.2, 0.4, 0.6, 0.8)) {
	res_list = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/subgroup_test/hipo16_rnaseq_subgroup_p@{p}.rds"))
	for(k in 2:6) {
		pdf(qq("/icgc/dkfzlsdf/analysis/B080/guz/subgroup_test/hipo16_rnaseq_subgroup_k@{k}_p@{p}.pdf"), width = 15, height = 10)
		for(fun in c("plot_ecdf", "consensus_heatmap", "membership_heatmap")) {
			qqcat("collecting plots for @{fun}, p@{p}, k@{k}\n")
			collect_plots(res_list, k = k, fun = get(fun))
		}
		dev.off()
	}
}


######### scrnaseq  ###########
for(p in c(0.2, 0.4, 0.6, 0.8)) {
	res_list = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/subgroup_test/scrnaseq_subgroup_p@{p}.rds"))
	for(k in 2:6) {
		pdf(qq("/icgc/dkfzlsdf/analysis/B080/guz/subgroup_test/scrnaseq_subgroup_k@{k}_p@{p}.pdf"), width = 15, height = 10)
		for(fun in c("plot_ecdf", "consensus_heatmap", "membership_heatmap")) {
			qqcat("collecting plots for @{fun}, p@{p}, k@{k}\n")
			collect_plots(res_list, k = k, fun = get(fun))
		}
		dev.off()
	}
}
