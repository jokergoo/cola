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

load("/home/guz/project/analysis/hipo15/signature_genes/moffitt_sigANDannotation.RData")
df_cell = read.table("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/CellularComposition_Cell01Cell03")

setwd("/icgc/dkfzlsdf/analysis/hipo/hipo_015/hipo15_rnaseq_cell_analysis/")
for(type in c("c1", "c2", "c3", "primary_tumor", "xenograft")) {
	fn = qq("hipo15_@{type}_subgroups.rds")
	res_list = readRDS(fn)
	res = get_single_run(res_list, "AAC", "skmeans")

	pdf(qq("hipo15_@{type}_subgroups.pdf"))
	plot_ecdf(res)
	consensus_heatmap(res, k = 2)
	if(type == "primary_tumor") {
		col_fun = colorRamp2(c(0, 1), c("white", "purple"))
		anno = as.data.frame(pro_5types_50[colnames(res@.env$data), ])
		p = test_to_known_factors(res, k = 2, anno)
		get_signatures(res, k = 2, anno = anno,
			anno_col = list(malignant = col_fun, immune = col_fun, stromal = col_fun,
				endothelial = col_fun, normal = col_fun))
		for(an in colnames(anno)) {
			decorate_annotation(an, {
				grid.text(qq("@{round(p[1, an], 2)}"), unit(0, "npc"), just = "right")
			})
		}
	} else if(type %in% c("c1", "c3")) {
		col_fun = colorRamp2(c(0, 1), c("white", "purple"))
		anno = df_cell[colnames(res@.env$data), ]
		p = test_to_known_factors(res, k = 2, anno)
		get_signatures(res, k = 2, anno = anno,
			anno_col = list(malignant = col_fun, immune = col_fun, stromal = col_fun, normal = col_fun))
		for(an in colnames(anno)) {
			decorate_annotation(an, {
				grid.text(qq("@{round(p[1, an], 2)}"), unit(0, "npc"), just = "right")
			})
		}
	} else {
		get_signatures(res, k = 2)
	}
	dimension_reduction(res, k = 2)
	dev.off()
}

