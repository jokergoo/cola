
#' Run subgroup classification in a batch
#'
#' @param data 
#' @param top_method 
#' @param partition_method 
#' @param mc.cores 
#' @param get_signatures 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
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
	res = list(list = list(), top_method = top_method, partition_method = partition_method, .env = .env)
	comb = expand.grid(top_method, partition_method, stringsAsFactors = FALSE)
	comb = comb[sample(nrow(comb), nrow(comb)), ]
	res$list = mclapply(seq_len(nrow(comb)), function(i) {
		gc(verbose = FALSE)
		tm = comb[i, 1]
		pm = comb[i, 2]
		qqcat("running classification for @{tm}:@{pm}. @{i}/@{nrow(comb)}\n")
		time_used = system.time(res <- consensus_partition(top_method = tm, partition_method = pm, .env = .env, ...))
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

	class(res) = c("run_all", "list")

	return(res)
}

print.run_all = function(x, ...) {
	qqcat("Top rows are extracted by '@{paste(x$top_method, collapse = ', ')}' methods.\n")
	qqcat("Subgroups are detected by '@{paste(x$partition_method, collapse = ', ')}' method.\n")
	qqcat("Number of partitions are tried for k = @{paste(x$list[[1]]$k, collapse = ', ')}\n")
}

get_single_run = function(res, top_method = "sd", partition_method = "kmeans") {
	nm = paste0(top_method, ":", partition_method)
	res$list[[nm]]
}
