
#' Run subgroup classification in a batch
#'
#' @param data a numeric matrix where subgroups are found by columns.
#' @param top_method method which are used to extract top n rows. Allowed methods
#'        are in [ALL_TOP_VALUE_METHOD()] and can be self-added by [register_top_value_fun()].
#' @param partition_method method which are used to do partition on data columns. 
#'        Allowed methods are in [ALL_PARTITION_METHOD()] and can be self-added 
#'        by [register_partition_fun()].
#' @param mc.cores number of cores to use.
#' @param get_signatures whether to run [get_signatures()] for each partition.
#' @param ... other arguments passed to [consensus_partition()].
#'
#' @return 
#' a `consensus_partition_all_methods` class object. Following methods can be used on it: [collect_plots()],
#' [collect_classes()], [get_class()], [get_single_run()].
#' 
#' @export
#' @import GetoptLong
#' @import parallel
#' @import pryr
run_all = function(data, top_method = ALL_TOP_VALUE_METHOD(), 
	partition_method = ALL_PARTITION_METHOD(), 
	mc.cores = 1, get_signatures = TRUE, ...) {
		
	.env = new.env()
	
	if(is.data.frame(data)) data = as.matrix(data)
	if(is.null(rownames(data))) rownames(data) = seq_len(nrow(data))

	l = rowSds(data) == 0
	data = data[!l, , drop = FALSE]
	qqcat("removed @{sum(l)} rows with sd = 0\n")

	all_value_list = lapply(top_method, function(tm) {
		qqcat("calculate @{tm} score for all rows\n")
		all_value = get_top_value_fun(tm)(data)
		all_value[is.na(all_value)] = -Inf
		return(all_value)
	})
	names(all_value_list) = top_method
	.env$all_value_list = all_value_list

	.env$data = data
	res = list(list = list(), top_method = top_method, partition_method = partition_method, 
		.env = .env)
	comb = expand.grid(top_method, partition_method, stringsAsFactors = FALSE)
	comb = comb[sample(nrow(comb), nrow(comb)), ]
	res$list = mclapply(seq_len(nrow(comb)), function(i) {
		gc(verbose = FALSE)
		tm = comb[i, 1]
		pm = comb[i, 2]
		qqcat("running classification for @{tm}:@{pm}. @{i}/@{nrow(comb)}\n")
		time_used = system.time(res <- consensus_partition(top_method = tm, 
			partition_method = pm, .env = .env, ...))
		mem = pryr::mem_used()
		attr(res, "system.time") = time_used
		attr(res, "mem_used") = as.numeric(mem)
		if(get_signatures) {
			for(k in res$k) {
				try(get_signatures(res, k = k))
			}
		}
		return(res)
	}, mc.cores = mc.cores)

	names(res$list) = paste(comb[, 1], comb[, 2], sep = ":")

	class(res) = c("consensus_partition_all_methods", "list")

	return(res)
}


#' Print the consensus_partition_all_methods object
#'
#' @param x a `consensus_partition_all_methods` object
#' @param ... other arguments
#' 
#' @export
#' @import GetoptLong
print.consensus_partition_all_methods = function(x, ...) {
	qqcat("Top rows are extracted by '@{paste(x$top_method, collapse = ', ')}' methods.\n")
	qqcat("Subgroups are detected by '@{paste(x$partition_method, collapse = ', ')}' method.\n")
	qqcat("Number of partitions are tried for k = @{paste(x$list[[1]]$k, collapse = ', ')}\n")
}

#' Get object for a single combination of top method and partition method
#'
#' @param res_list object returned from [run_all()].
#' @param top_method a single string which is used in [run_all()]
#' @param partition_method a single string which is used in [run_all()]
#'
#' @return a `consensus_partition` class object.
#' @export
get_single_run = function(res_list, top_method = res_list$top_method[1], 
	partition_method = res_list$partition_method[1]) {
	nm = paste0(top_method, ":", partition_method)
	res_list$list[[nm]]
}
