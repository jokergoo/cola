 
if(grepl("tbi", Sys.info()["nodename"]) & Sys.info()["user"] == "guz") {
	source("~/project/development/cola/load.R")
} else {
	source("~/project/cola/load.R")
}

load("~/project/development/cola/inst/extdata/TCGA_GBM_microarray_subset.Rdata")
data = t(apply(data, 1, adjust_outlier))
res = run_all_consensus_partition_methods(data, top_n = c(500, 1000), k = 2:6,
	top_method = c("sd", "vc"), partition_method = c("kmeans", "hclust"),
	known_anno = data.frame(subtype = subtype),
	known_col = list(subtype = structure(1:4, names = unique(subtype))))

get_class(res, k = 3)
collect_plots(res, k = 3, fun = consensus_heatmap)
collect_plots(res, k = 3, fun = membership_heatmap)
collect_plots(res, k = 3, fun = get_signatures)
collect_plots(res, k = 3, fun = plot_ecdf)
collect_plots(res, k = 3, fun = dimension_reduction)
collect_classes(res, k = 3)
top_rows_overlap(res, top_n = 2000)
top_rows_heatmap(res, top_n = 2000)
get_stat(res, k = 3)
get_class(res, k = 3)
test_to_known(res, k = 3)

obj = get_single_run(res, top_method = "AAC", partition_method = "skmeans")
select_partition_number(obj)
collect_plots(obj)
collect_classes(obj)
consensus_heatmap(obj, k = 3)
membership_heatmap(obj, k = 3)
get_signatures(obj, k = 3)
dimension_reduction(obj, k = 3)

get_stat(obj)
get_class(obj, k = 3)
get_consensus(obj, k = 3)
get_membership(obj, k = 3)

test_to_known(obj, k = 3)

