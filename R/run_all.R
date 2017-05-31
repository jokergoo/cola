
# == title
# Run subgroup classification for all methods
#
# == param
# -data a numeric matrix where subgroups are found by columns.
# -top_method method which are used to extract top n rows. Allowed methods
#        are in `all_top_value_methods` and can be self-added by `register_top_value_fun`.
# -partition_method method which are used to do partition on data columns. 
#        Allowed methods are in `all_partition_methods` and can be self-added 
#        by `register_partition_fun`.
# -mc.cores number of cores to use`.
# -... other arguments passed to `consensus_partition`.
#
# == details
# The function run consensus partitions for all combination of top methods and parittion methods.
#
# == return 
# a `ConsensusPartitionList-class` object. 
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
run_all_consensus_partition_methods = function(data, top_method = all_top_value_methods(), 
	partition_method = all_partition_methods(), 
	mc.cores = 1, ...) {
		
	.env = new.env()
	
	if(is.data.frame(data)) data = as.matrix(data)
	if(is.null(rownames(data))) rownames(data) = seq_len(nrow(data))

	l = rowSds(data) == 0
	data = data[!l, , drop = FALSE]
	if(sum(l)) qqcat("removed @{sum(l)}/@{length(l)} rows with sd = 0\n")

	all_value_list = lapply(top_method, function(tm) {
		qqcat("calculate @{tm} score for all rows\n")
		all_value = get_top_value_fun(tm)(data)
		all_value[is.na(all_value)] = -Inf
		return(all_value)
	})
	names(all_value_list) = top_method
	.env$all_value_list = all_value_list

	.env$data = data
	res_list = ConsensusPartitionList(
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
		for(k in res@k) {
			try(get_signatures(res, k = k, plot = FALSE))
		}
		return(res)
	}, mc.cores = mc.cores)
	names(lt) = paste(comb[, 1], comb[, 2], sep = ":")
browser()
	i_error = which(sapply(lt, inherits, "try-error"))
	if(length(i_error)) {
		sapply(lt[[i_error]], cat)
		stop("There are errors when doing mclapply.")
	}
	
	res_list@list = lt
	# cat("adjust class labels according to the consensus classification\n")
	reference_class = lapply(lt[[1]]@k, function(k) {
		get_class(res_list, k)$class
	})
	for(i in seq_along(lt)) {
		res = lt[[i]]
		for(k in res@k) {
			# relabel the class according to the class in the first object
	        ik = which(res@k == k)
	        
        	# following elements need to be relabeled
        	# - res$object_list[[ik]]$classification$class
        	# - column order of res$object_list[[ik]]$membership
        	# - res$object_list[[ik]]$membership_each
        	class_df = get_class(res, k)
        	class = class_df[, "class"]
        	map = relabel_class(class, reference_class[[ik]])
        	map2 = structure(names(map), names = map)
        	res@object_list[[ik]]$class_df$class = as.numeric(map[as.character(class)])
        	
        	res@object_list[[ik]]$membership = res@object_list[[ik]]$membership[, as.numeric(map2[as.character(1:k)]) ]
			colnames(res@object_list[[ik]]$membership) = paste0("p", 1:k)
			
			odim = dim(res@object_list[[ik]]$membership_each)
			res@object_list[[ik]]$membership_each = as.numeric(map[as.character(res@object_list[[ik]]$membership_each)])
			dim(res@object_list[[ik]]$membership_each) = odim
	        
	    }
	    lt[[i]] = res
	}

	return(res_list)
}

# == title
# Print the ConsensusPartitionList object
#
# == param
# -object a `ConsensusPartitionList-class` object
# -... other arguments
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "show",
	signature = "ConsensusPartitionList",
	definition = function(object) {
	qqcat("Top rows are extracted by '@{paste(object@top_method, collapse = ', ')}' methods.\n")
	qqcat("Subgroups are detected by '@{paste(object@partition_method, collapse = ', ')}' method.\n")
	qqcat("Number of partitions are tried for k = @{paste(object@list[[1]]@k, collapse = ', ')}\n")
})

# == title
# Get result for a single top method and partition method
#
# == param
# -object a `ConsensusPartitionList-class` object
# -top_method a single string which is used in `run_all_consensus_partition_methods`
# -partition_method a single string which is used in `run_all_consensus_partition_methods`
#
# == return 
# A `ConsensusPartition-class` object.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "get_single_run",
	signature = "ConsensusPartitionList",
	definition = function(object, top_method = object@top_method[1], 
		partition_method = object@partition_method[1]) {
	if(!top_method %in% object@top_method) {
		stop(qq("@{top_method} was not used in `run_all_consensus_partition_methods`."))
	}
	if(!partition_method %in% object@partition_method) {
		stop(qq("@{partition_method} was not used in `run_all_consensus_partition_methods`."))
	}
	nm = paste0(top_method, ":", partition_method)
	object@list[[nm]]
})
