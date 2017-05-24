 
if(grepl("tbi", Sys.info()["nodename"]) & Sys.info()["user"] == "guz") {
	source("~/project/development/cola/load.R")
} else {
	source("~/project/cola/load.R")
}

load("~/project/development/cola/inst/extdata/TCGA_GBM_microarray_subset.Rdata")
data = t(apply(data, 1, adjust_outlier))
res = run_all_consensus_partition_methods(data, top_n = c(500, 1000), k = 2:6,
	known_anno = data.frame(subtype = subtype))

obj = get_single_run(res, top_method = "MAD", partition_method = "pam")

