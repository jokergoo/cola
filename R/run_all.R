
try_and_trace = function(expr) {
	env = new.env()
	try(withCallingHandlers(expr, error = function(e) {
		stack = sys.calls()
		env$stack = stack
	}), silent = TRUE)
	return(env$stack)
}



# == title
# Run subgroup classification from all methods
#
# == param
# -data a numeric matrix where subgroups are found by samples
# -top_method method which are used to extract top n rows. Allowed methods
#        are in `all_top_value_methods` and can be self-added by `register_top_value_fun`.
# -partition_method method which are used to do partition on samples. 
#        Allowed methods are in `all_partition_methods` and can be self-added 
#        by `register_partition_fun`.
# -k a list number of partitions.
# -mc.cores number of cores to use.
# -known_anno a data frame with known annotation of samples
# -known_col a list of colors for the annotations in ``known_anno``.
# -... other arguments passed to `consensus_partition`.
#
# == details
# The function run consensus partitions for all combinations of top methods and parittion methods.
#
# == return 
# A `ConsensusPartitionList-class` object. 
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
run_all_consensus_partition_methods = function(data, top_method = all_top_value_methods(), 
	partition_method = all_partition_methods(), k = 2:6,
	mc.cores = 1, known_anno = NULL, known_col = NULL, ...) {
		
	.env = new.env()
	
	if(is.data.frame(data)) data = as.matrix(data)
	if(is.null(rownames(data))) rownames(data) = seq_len(nrow(data))

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
		.env = .env
	)

	if(!is.null(known_anno)) {
		if(is.atomic(known_anno)) {
			known_nm = deparse(substitute(known_anno))
			known_anno = data.frame(known_anno)
			colnames(known_anno) = known_nm
			if(!is.null(known_col)) {
				known_col = list(known_col)
				names(known_col) = known_nm
			}
		}
	}

	if(is.null(known_col)) {
		known_col = lapply(known_anno, ComplexHeatmap:::default_col)
	} else {
		if(is.null(names(known_col))) {
			if(length(known_col) == ncol(known_anno)) {
				names(known_col) = colnames(known_anno)
			} else {
				known_col = lapply(known_anno, ComplexHeatmap:::default_col)
			}
		}
		for(nm in names(known_anno)) {
			if(is.null(known_col[[nm]])) {
				known_col[[nm]] = ComplexHeatmap:::default_col(known_anno[[nm]])
			}
		}
	}
	if(is.null(known_anno)) {
		known_col = NULL
	}

	comb = expand.grid(top_method, partition_method, stringsAsFactors = FALSE)
	# comb = comb[sample(nrow(comb), nrow(comb)), ]
	od = order(rep(sapply(partition_method, function(x) attr(get_partition_fun(x), "execution_scale")), each = length(top_method)), decreasing = TRUE)
	comb = comb[od, ]
	lt = mclapply(seq_len(nrow(comb)), function(i) {
		tm = comb[i, 1]
		pm = comb[i, 2]
		qqcat("running classification for @{tm}:@{pm}. @{i}/@{nrow(comb)}\n")
		x = try_and_trace(res <- consensus_partition(top_method = tm, partition_method = pm, k = k,
			known_anno = known_anno, known_col = known_col, .env = .env, verbose = FALSE, ...))
		
		if(inherits(x, "pairlist")) {
			print(x)
			stop(qq("You have an error when doing partition for @{tm}:@{pm}."))
		}
		for(k in res@k) {
			try(get_signatures(res, k = k, plot = FALSE, verbose = FALSE))
		}
		return(res)
	}, mc.cores = mc.cores)
	names(lt) = paste(comb[, 1], comb[, 2], sep = ":")

	i_error = which(sapply(lt, inherits, "try-error"))
	if(length(i_error)) {
		for(i in i_error) {
			cat(names(lt)[i], ": ", lt[[i]], "\n", sep = "")
		}
		stop("There are errors when doing mclapply.")
	}
	
	res_list@list = lt
	cat("adjust class labels according to the consensus classification from all methods.\n")
	reference_class = lapply(lt[[1]]@k, function(k) {
		class_ids = get_class(res_list, k)$class
		mean_dist = tapply(seq_len(ncol(data)), class_ids, function(ind) {
			n = length(ind)
			if(n == 1) {
				return(Inf)
			}
			sum(dist(t(data[, ind, drop = FALSE]))^2)/(n*(n-1)/2)
		})
		map = structure(names = names(mean_dist)[order(mean_dist)], names(mean_dist))
		class_ids = as.numeric(map[as.character(class_ids)])
		return(class_ids)
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
	qqcat("\n")
	qqcat("Following methods can be applied to this 'ConsensusPartitionList' object:\n")
	txt = showMethods(classes = "ConsensusPartitionList", where = topenv(), printTo = FALSE)
	txt = grep("Function", txt, value = TRUE)
	fname = gsub("Function: (.*?) \\(package.*$", "\\1", txt)
	print(fname)
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



# == title
# Get the best number of partitions
#
# == param
# -object a `ConsensusPartitionList-class` object
#
# == value
# A data frame with best k for each combination of methods
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "get_best_k",
	signature = "ConsensusPartitionList",
	definition = function(object) {

	best_k = NULL
	for(tm in object@top_method) {
		for(pm in object@partition_method) {
			nm = paste0(tm, ":", pm)
			obj = object@list[[nm]]
			best_k[nm] = get_best_k(obj)
		}
	}
	return(data.frame(best_k = best_k))
})
