ConsensusPartitionList = setClass("ConsensusPartitionList",
	slots = list(
		list = "list",
		top_method = "character",
		partition_method = "character",
		.env = "environment"
))

#' Run subgroup classification in a batch
#'
#' @param data a numeric matrix where subgroups are found by columns.
#' @param top_method method which are used to extract top n rows. Allowed methods
#'        are in [ALL_TOP_VALUE_METHOD()] and can be self-added by [register_top_value_fun()].
#' @param partition_method method which are used to do partition on data columns. 
#'        Allowed methods are in [ALL_PARTITION_METHOD()] and can be self-added 
#'        by [register_partition_fun()].
#' @param mc.cores number of cores to use`.
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
run_all_consensus_partition_methods = function(data, top_method = ALL_TOP_VALUE_METHOD(), 
	partition_method = ALL_PARTITION_METHOD(), 
	mc.cores = 1, get_signatures = TRUE, ...) {
		
	.env = new.env()
	
	if(is.data.frame(data)) data = as.matrix(data)
	if(is.null(rownames(data))) rownames(data) = seq_len(nrow(data))

	l = rowSds(data) == 0
	data = data[!l, , drop = FALSE]
	if(sum(l)) qqcat("removed @{sum(l)} rows with sd = 0\n")

	all_value_list = lapply(top_method, function(tm) {
		qqcat("calculate @{tm} score for all rows\n")
		all_value = get_top_value_fun(tm)(data)
		all_value[is.na(all_value)] = -Inf
		return(all_value)
	})
	names(all_value_list) = top_method
	.env$all_value_list = all_value_list

	.env$data = data
	res = ConsensusPartitionList(
		list = list(), 
		top_method = top_method, 
		partition_method = partition_method, 
		.env = .env)
	comb = expand.grid(top_method, partition_method, stringsAsFactors = FALSE)
	comb = comb[sample(nrow(comb), nrow(comb)), ]
	lt = mclapply(seq_len(nrow(comb)), function(i) {
		gc(verbose = FALSE)
		tm = comb[i, 1]
		pm = comb[i, 2]
		qqcat("running classification for @{tm}:@{pm}. @{i}/@{nrow(comb)}\n")
		time_used = system.time(res <- consensus_partition(top_method = tm, 
			partition_method = pm, .env = .env, ...))
		attr(res, "system.time") = time_used
		if(get_signatures) {
			for(k in res$k) {
				try(get_signatures(res, k = k))
			}
		}
		return(res)
	}, mc.cores = mc.cores)

	cat("adjust class labels\n")
	reference_class = vector("list", length(lt[[1]]@k))
	for(i in seq_along(lt)) {
		res = lt[[i]]
		for(k in res@k)
			# relabel the class according to the class in the first object
	        ik = which(res@k == k)
	        if(is.null(reference_class[[ik]])) {
	        	reference_class[[ik]] = get_class(res, k)[, "class"]
	        } else {
	        	# following elements need to be relabeled
	        	# - res$object_list[[ik]]$classification$class
	        	# - column order of res$object_list[[ik]]$membership
	        	# - res$object_list[[ik]]$membership_each
	        	class = get_class(res, k)[, "class"]
	        	map = relabel_class(class, reference_class[[ik]])
	        	map2 = structure(names(map), names = map)
	        	res@object_list[[ik]]$class_df$class = as.numeric(map[as.character(class)])
	        	
	        	res@object_list[[ik]]$membership = res@object_list[[ik]]$membership[, as.numeric(map2[as.character(1:k)]) ]
				colnames(res@object_list[[ik]]$membership) = paste0("p", 1:k)
				
				odim = dim(res@object_list[[ik]]$membership_each)
				res@object_list[[ik]]$membership_each = as.numeric(map[as.character(res@object_list[[ik]]$membership_each)])
				dim(res@object_list[[ik]]$membership_each) = odim
	        }
	    }
	    lt[[i]] = res
	}

	res@list = lt
	names(res@list) = paste(comb[, 1], comb[, 2], sep = ":")

	return(res)
}


#' Print the consensus_partition_all_methods object
#'
#' @param x a `consensus_partition_all_methods` object
#' @param ... other arguments
#' 
#' @export
#' @import GetoptLong
setMethod("show",
	signature = "ConsensusPartitionList",
	definition = function(object, ...) {
	qqcat("Top rows are extracted by '@{paste(object@top_method, collapse = ', ')}' methods.\n")
	qqcat("Subgroups are detected by '@{paste(object@partition_method, collapse = ', ')}' method.\n")
	qqcat("Number of partitions are tried for k = @{paste(object@list[[1]]@k, collapse = ', ')}\n")
}

#' Get object for a single combination of top method and partition method
#'
#' @param res_list object returned from [run_all()].
#' @param top_method a single string which is used in [run_all()]
#' @param partition_method a single string which is used in [run_all()]
#'
#' @return a `consensus_partition` class object.
#' @export
setMethod("get_single_run",
	signature = "ConsensusPartitionList",
	definition = function(object, top_method = object@top_method[1], 
		partition_method = object@partition_method[1]) {
	nm = paste0(top_method, ":", partition_method)
	object@list[[nm]]
}
