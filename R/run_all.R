
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
# -data a numeric matrix where subgroups are found by columns.
# -top_value_method method which are used to extract top n rows. Allowed methods
#        are in `all_top_value_methods` and can be self-added by `register_top_value_method`.
# -partition_method method which are used to do partition on samples. 
#        Allowed methods are in `all_partition_methods` and can be self-added 
#        by `register_partition_method`.
# -max_k maximum number of partitions to try. It will use ``2:max_k`` partitions.
# -mc.cores number of cores to use.
# -anno a data frame with known annotation of samples
# -anno_col a list of colors for the annotations in ``known_anno``.
# -... other arguments passed to `consensus_partition`.
#
# == details
# The function run consensus partitions for all combinations of top methods and parittion methods.
#
# == return 
# A `ConsensusPartitionList-class` object. Simply type the name of the object in the R interactive session
# to see which functions can be applied on it.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# \dontrun{
# set.seed(123)
# m = cbind(rbind(matrix(rnorm(20*20, mean = 1), nr = 20),
# 			    matrix(rnorm(20*20, mean = -1), nr = 20)),
# 	rbind(matrix(rnorm(20*20, mean = -1), nr = 20),
# 			    matrix(rnorm(20*20, mean = 1), nr = 20))
# ) + matrix(rnorm(40*40), nr = 40)
# rl = run_all_consensus_partition_methods(data = m, top_n = c(20, 30, 40))
# }
# data(cola_rl)
# cola_rl
run_all_consensus_partition_methods = function(data, 
	top_value_method = all_top_value_methods(), 
	partition_method = all_partition_methods(), max_k = 6,
	mc.cores = 1, anno = NULL, anno_col = NULL, ...) {
	
	cl = match.call()

	.env = new.env()
	
	if(is.data.frame(data)) data = as.matrix(data)
	if(is.null(rownames(data))) rownames(data) = seq_len(nrow(data))

	all_top_value_list = lapply(top_value_method, function(tm) {
		qqcat("calculate @{tm} score for all rows\n")
		all_top_value = get_top_value_method(tm)(data)
		all_top_value[is.na(all_top_value)] = -Inf
		return(all_top_value)
	})
	names(all_top_value_list) = top_value_method
	.env$all_top_value_list = all_top_value_list

	.env$data = data
	res_list = ConsensusPartitionList(
		list = list(), 
		top_value_method = top_value_method, 
		partition_method = partition_method, 
		consensus_class = NULL,
		.env = .env
	)

	if(!is.null(anno)) {
		if(is.atomic(anno)) {
			known_nm = deparse(substitute(anno))
			anno = data.frame(anno)
			colnames(anno) = known_nm
			if(!is.null(anno_col)) {
				anno_col = list(anno_col)
				names(anno_col) = known_nm
			}
		}
	}

	if(is.null(anno_col)) {
		anno_col = lapply(anno, ComplexHeatmap:::default_col)
	} else {
		if(is.null(names(anno_col))) {
			if(length(anno_col) == ncol(anno)) {
				names(anno_col) = colnames(anno)
			} else {
				anno_col = lapply(anno, ComplexHeatmap:::default_col)
			}
		}
		for(nm in names(anno)) {
			if(is.null(anno_col[[nm]])) {
				anno_col[[nm]] = ComplexHeatmap:::default_col(anno[[nm]])
			}
		}
	}
	if(is.null(anno)) {
		anno_col = NULL
	}

	comb = expand.grid(top_value_method, partition_method, stringsAsFactors = FALSE)
	# comb = comb[sample(nrow(comb), nrow(comb)), ]
	od = order(rep(sapply(partition_method, function(x) attr(get_partition_method(x), "execution_time")), each = length(top_value_method)), decreasing = TRUE)
	comb = comb[od, , drop = FALSE]
	lt = mclapply(seq_len(nrow(comb)), function(i) {
		tm = comb[i, 1]
		pm = comb[i, 2]
		qqcat("running classification for @{tm}:@{pm}. @{i}/@{nrow(comb)}\n")
		x = try_and_trace(res <- consensus_partition(top_value_method = tm, partition_method = pm, max_k = max_k,
			anno = anno, anno_col = anno_col, .env = .env, verbose = interactive() & mc.cores == 1, ...))
		
		if(inherits(x, "pairlist")) {
			print(x)
			qqcat("You have an error when doing partition for @{tm}:@{pm}.\n")
			stop("You have an error.")
		}
		# for(k in res@k) {
		# 	try(get_signatures(res, k = k, plot = FALSE, verbose = FALSE))
		# }
		return(res)
	}, mc.cores = mc.cores)
	names(lt) = paste(comb[, 1], comb[, 2], sep = ":")

	for(i in seq_along(lt)) {
		if(!identical(.env, lt[[i]]@.env)) {
			lt[[i]]@.env = .env
		}
	}

	i_error = which(sapply(lt, inherits, "try-error"))
	if(length(i_error)) {
		for(i in i_error) {
			cat(names(lt)[i], ": ", lt[[i]], "\n", sep = "")
		}
		stop("There are errors when doing mclapply.")
	}
	
	res_list@list = lt
	data2 = t(scale(t(data)))
	cat("adjust class labels according to the consensus classification from all methods.\n")
	reference_class = lapply(lt[[1]]@k, function(k) {
		cl_df = get_consensus_from_multiple_methods(res_list, k)
		class_ids = cl_df$class
		mean_dist = tapply(seq_len(ncol(data2)), class_ids, function(ind) {
			n = length(ind)
			if(n == 1) {
				return(Inf)
			}
			sum(dist(t(data2[, ind, drop = FALSE]))^2)/(n*(n-1)/2)
		})
		map = structure(names = names(mean_dist)[order(mean_dist)], names(mean_dist))
		class_ids = as.numeric(map[as.character(class_ids)])
		cl_df$class = class_ids
		return(cl_df)
	})
	names(reference_class) = as.character(lt[[1]]@k)
	res_list@consensus_class = reference_class

	for(i in seq_along(lt)) {
		res = lt[[i]]
		for(k in res@k) {
			# relabel the class according to the class in the first object
	        ik = which(res@k == k)
	        
        	# following elements need to be relabeled
        	# - res$object_list[[ik]]$classification$class
        	# - column order of res$object_list[[ik]]$membership
        	# - res$object_list[[ik]]$membership_each
        	class_df = get_classes(res, k)
        	class = class_df[, "class"]
        	map = relabel_class(class, reference_class[[ik]]$class)
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
	res_list@list = lt
	res_list@comb = comb
	res_list@calling = cl

	return(res_list)
}

get_consensus_from_multiple_methods = function(object, k) {

	res = object
	partition_list = NULL
	mean_cophcor = NULL
	mean_silhouette = NULL
	reference_class = NULL
	for(tm in object@top_value_method) {
		for(pm in object@partition_method) {
			nm = paste0(tm, ":", pm)
			obj = object@list[[nm]]
			ik = which(obj@k == k)

			membership = get_membership(obj, k)
			if(is.null(reference_class)) {
	        	reference_class = get_classes(obj, k)[, "class"]
	        } else {
	        	map = relabel_class(get_classes(obj, k)[, "class"], reference_class)
	        	map2 = structure(names(map), names = map)
	        	membership = membership[, as.numeric(map2[as.character(1:k)]) ]
				colnames(membership) = paste0("p", 1:k)
			}

			partition_list = c(partition_list, list(as.cl_partition(membership)))
			mean_silhouette = c(mean_silhouette, mean(get_classes(obj, k)[, "silhouette"]))
		}
	}

	mean_silhouette[mean_silhouette < 0] = 0
	consensus = cl_consensus(cl_ensemble(list = partition_list), weights = mean_silhouette)
	m = cl_membership(consensus)
	class(m) = "matrix"
	colnames(m) = paste0("p", 1:k)
	class = as.vector(cl_class_ids(consensus))
	df = cbind(class = class, as.data.frame(m))
	df$entropy = apply(m, 1, entropy)

	membership_each = do.call("cbind", lapply(partition_list, function(x) {
		as.vector(cl_class_ids(x))
	}))

	consensus_mat = matrix(1, nrow = nrow(m), ncol = nrow(m))
	for(i in 1:(nrow(membership_each)-1)) {
		for(j in (i+1):nrow(membership_each)) {
			consensus_mat[i, j] = sum(membership_each[i, ] == membership_each[j, ])/ncol(membership_each)
			consensus_mat[j, i] = consensus_mat[i, j]
		}
 	}
 
	df$silhouette = silhouette(class, dist(t(consensus_mat)))[, "sil_width"]

	return(df)
}


# == title
# Print the ConsensusPartitionList object
#
# == param
# -object a `ConsensusPartitionList-class` object.
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
	
	# obj_name = deparse(substitute(object, env = env))
	obj_name = "object"

	qqcat("A 'ConsensusPartitionList' object with @{length(object@top_value_method)*length(object@partition_method)} methods.\n")
	qqcat("  On a matrix with @{nrow(object@.env$data)} rows and @{ncol(object@.env$data)} columns.\n")
	qqcat("  Top rows are extracted by '@{paste(object@top_value_method, collapse = ', ')}' methods.\n")
	qqcat("  Subgroups are detected by '@{paste(object@partition_method, collapse = ', ')}' method.\n")
	qqcat("  Number of partitions are tried for k = @{paste(object@list[[1]]@k, collapse = ', ')}\n")
	qqcat("\n")
	qqcat("Following methods can be applied to this 'ConsensusPartitionList' object:\n")
	txt = showMethods(classes = "ConsensusPartitionList", where = topenv(), printTo = FALSE)
	txt = grep("Function", txt, value = TRUE)
	fname = gsub("Function: (.*?) \\(package.*$", "\\1", txt)
	print(fname)
	cat("\n")

	qqcat("You can get result for a single method by e.g. @{obj_name}[\"@{object@top_value_method[1]}\", \"@{object@partition_method[1]}\"] or @{obj_name}[\"@{object@top_value_method[1]}:@{object@partition_method[1]}\"]\n")
	if(length(object@top_value_method) == 1) {
		ri = qq("\"@{object@top_value_method[1]}\"")
	} else {
		ri = qq("c(\"@{object@top_value_method[1]}\", \"@{object@top_value_method[2]}\")")
	}
	if(length(object@partition_method) == 1) {
		ci = qq("\"@{object@partition_method[1]}\"")
	} else {
		ci = qq("c(\"@{object@partition_method[1]}\", \"@{object@partition_method[2]}\")")
	}
	if(length(object@top_value_method) > 1 | length(object@partition_method) > 1) {
		qqcat("or a subset of methods by @{obj_name}[@{ri}], @{ci}]\n")
	}
})


# == title
# Subset a ConsensusPartitionList object
#
# == param
# -x a `ConsensusPartition-class` object.
# -i index for top value methods, character or nummeric.
# -j index for partition methods, character or nummeric.
# -drop whether drop class
#
# == value
# A `ConsensusPartitionList-class` object or a `ConsensusPartition-class` object.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# cola_rl[c("sd", "MAD"), c("hclust", "kmeans")]
# cola_rl["sd", "kmeans"]
# cola_rl["sd", "kmeans", drop = FALSE]
# cola_rl["sd", ]
# cola_rl[, "hclust"]
# cola_rl[1:2, 1:2]
"[.ConsensusPartitionList" = function (x, i, j, drop = TRUE) {

	cl = as.list(match.call())
	called_args = names(cl)[-1]

	all_top_value_methods = x@top_value_method
	all_partition_methods = x@partition_method
	n_top_value_methods = length(all_top_value_methods)
	n_partition_methods = length(all_partition_methods)

	if(nargs() == 1) {
		return(x)
	}
    if("i" %in% called_args & "j" %in% called_args) {
        if(is.numeric(i)) {
        	i = all_top_value_methods[i]
        }
        if(is.numeric(j)) {
        	j = all_partition_methods[j]
        }
        i = intersect(i, all_top_value_methods)
        j = intersect(j, all_partition_methods)
        l = x@comb[, 1] %in% i & x@comb[, 2] %in% j
        l[is.na(l)] = FALSE
        x@comb = x@comb[l, , drop = FALSE]
        x@list = x@list[l]
        x@top_value_method = i
        x@partition_method = j
        if(length(x@list) == 0) {
        	return(NULL)
        }
        if(length(x@list) == 1 && drop) {
        	x = x@list[[1]]
        }
        return(x)
    }
    if(!"i" %in% called_args & "j" %in% called_args) {
        if(is.numeric(j)) {
        	j = all_partition_methods[j]
        }
        j = intersect(j, all_partition_methods)
        l = x@comb[, 2] %in% j
        l[is.na(l)] = FALSE
        x@comb = x@comb[l, , drop = FALSE]
        x@list = x@list[l]
        x@top_value_method = all_top_value_methods
        x@partition_method = j
        if(length(x@list) == 0) {
        	return(NULL)
        }
        if(length(x@list) == 1 && drop) {
        	x = x@list[[1]]
        } 
        return(x)
    }
    if("i" %in% called_args & !"j" %in% called_args) {
    	if(nargs() == 3 & "drop" %in% called_args) {
    		if(length(i) > 1) {
    			stop("index can only be length 1.")
    		}
    		if(is.numeric(i)) {
	    		i = names(x)[i]
	    	}
	    	a = strsplit(i, ":")[[1]]
	    	return(x[a[1], a[2], cl$drop])
    	}
    	if(nargs() == 2) {
    		if(length(i) > 1) {
    			stop("index can only be length 1.")
    		}
	    	if(is.numeric(i)) {
	    		i = names(x)[i]
	    	}
	    	a = strsplit(i, ":")[[1]]
	    	return(x[a[1], a[2]])
	    }

        if(is.numeric(i)) {
        	i = all_top_value_methods[j]
        }
        i = intersect(i, all_top_value_methods)
        l = x@comb[, 1] %in% i
        l[is.na(l)] = FALSE
        x@comb = x@comb[l, , drop = FALSE]
        x@list = x@list[l]
        x@top_value_method = i
        x@partition_method = all_partition_methods
        if(length(x@list) == 0) {
        	return(NULL)
        }
        if(length(x@list) == 1 && drop) {
        	x = x@list[[1]]
        }
        return(x)
    }
    
    return(x)
}


