
############################
# merge k partitions into 2 partitions
#
merge_into_two_groups = function(object, verbose = TRUE) {

	object2 = new("ConsensusPartition")

	object2@object_list = list()
	object2@k = 2
	object2@n_partition = object@n_partition
	object2@partition_method = object@partition_method
	object2@top_value_method = object@top_value_method
	object2@top_n = object@top_n
	object2@anno = object@anno
	object2@anno_col = object@anno_col
	object2@scale_rows = object@scale_rows
	object2@sample_by = object@sample_by
	object2@column_index = object@column_index
	object2@cache = object@cache
	object2@others = object@others
	object2@hash = ""
	object2@.env = new.env()
	object2@.env$all_top_value_list = object@.env$all_top_value_list
	object2@.env$column_index = object@.env$column_index
	object2@.env$data = object@.env$data

	sample_by = object@sample_by
	object_list = lapply(object@object_list, merge_single, verbose = verbose)

	stat_df = data.frame(
		k = sapply(object_list, function(x) x$param$k[1]),
		'1-PAC' = sapply(object_list, function(x) x$stat$`1-PAC`),
		"mean_silhouette" = sapply(object_list, function(x) x$stat$mean_silhouette),
		"concordance" = sapply(object_list, function(x) x$stat$concordance),
		"Jaccard" = sapply(object_list, function(x) x$stat$Jaccard),
		check.names = FALSE
	)
	stat_df2 = stat_df[stat_df$Jaccard > 0.95, , drop = FALSE]

	if(nrow(stat_df2)) {
		select = order(stat_df2$`1-PAC`, stat_df2$mean_silhouette, stat_df2$concordance, -stat_df2$k, decreasing = TRUE)[1]
		select = stat_df2$k[select] - 1
	} else {
		select = order(stat_df$`1-PAC`, stat_df$mean_silhouette, stat_df$concordance, -stat_df$k, decreasing = TRUE)[1]
		select = stat_df$k[select] - 1
	}

	if(verbose) qqcat("select to merge k = @{select+1}\n")

	reference_class = object@object_list[[1]]$class_df[, "class"]

	object2@object_list = list("2" = object_list[[select]])

	# adjust labels according to the original 2-group classification
    class = object2@object_list[[1]]$class_df[, "class"]
	map = relabel_class(class, reference_class, full_set = 1:2)
	map2 = structure(names(map), names = map)

	object2@object_list[[1]]$class_df$class = as.numeric(map[as.character(class)])
	object2@object_list[[1]]$membership = object2@object_list[[1]]$membership[, as.numeric(map2[as.character(1:2)]) ]
	colnames(object2@object_list[[1]]$membership) = paste0("p", 1:2)
	
	odim = dim(object2@object_list[[1]]$membership_each)
	object2@object_list[[1]]$membership_each = as.numeric(map[as.character(object2@object_list[[1]]$membership_each)])
	dim(object2@object_list[[1]]$membership_each) = odim
	        
	object2@hash =  digest(object2)
	object2@running_time = object@running_time
	
	return(object2)
}

# merge k groups into 2 groups
merge_single = function(obj, verbose = TRUE, sample_by = "row") {

	k = obj$param$k[1]

	if(verbose) qqcat("- merge @{k} groups to 2 groups\n")
	if(k == 2) {
		if(verbose) cat("  - skip k = 2\n")
		return(obj)
	}

	class_df = obj$class_df
	membership = obj$membership
	consensus_mat = obj$consensus
	param = obj$param
	membership_each = obj$membership_each
	stat = obj$stat

	cl_mean = tapply(1:nrow(class_df), class_df$class, function(ind) {
		if(length(ind) == 1) {
			return(0)
		} else {
			mm = consensus_mat[ind, ind, drop = FALSE]
			mean(mm[lower.tri(mm)])
		}
	})
	group1 = as.numeric(names(which.max(cl_mean)))
	if(verbose) qqcat("  - group @{group1} is the most stable group\n")

	membership_each[membership_each != group1] = 0
	membership_each[membership_each == group1] = 1
	membership_each[membership_each == 0] = 2

	partition_list = data.frame(membership_each)
	partition_list = lapply(partition_list, as.cl_hard_partition)

	partition_list = cl_ensemble(list = partition_list)
	if(sample_by == "row") {
		partition_consensus = cl_consensus(partition_list)
	} else {
		partition_consensus = cl_consensus2(partition_list, k)
	}

	# note: number of class_ids may be less than k
	class_ids = as.vector(cl_class_ids(partition_consensus))
	
	membership_mat = cl_membership(partition_consensus)
	class(membership_mat) = "matrix"
	
	colnames(membership_mat) = paste0("p", 1:ncol(membership_mat))
	rownames(membership_mat) = rownames(membership_each)
	attr(membership_mat, "n_of_classes") = NULL
	attr(membership_mat, "is_cl_hard_partition") = NULL

	class_ids_by_top_n = tapply(seq_along(partition_list), param$top_n, function(ind) {
		if(sample_by == "row") {
			partition_consensus = cl_consensus(cl_ensemble(list = partition_list[ind]))
		} else {
			partition_consensus = cl_consensus2(cl_ensemble(list = partition_list[ind]), 2)
		}
		ci = as.vector(cl_class_ids(partition_consensus))
		map = relabel_class(ci, class_ids, full_set = 1:2)
		as.numeric(map[as.character(ci)])
	})

	# adjust class labels in each membership matrix to fit to the consensus class labels
	membership_each = do.call("cbind", lapply(seq_along(partition_list), function(i) {
		x = partition_list[[i]]
		class_ids = class_ids_by_top_n[[as.character(param$top_n[i])]]
		class = as.vector(cl_class_ids(x))
		map = relabel_class(class, class_ids, full_set = 1:2)
		class = as.numeric(map[as.character(class)])
		as.integer(class)
	}))
	rownames(membership_each) = rownames(membership_mat)

	membership_each2 = membership_each
	membership_each2[is.na(membership_each2)] = as.integer(0)
	consensus_mat = get_consensus_matrix(membership_each2)
 	rownames(consensus_mat) = rownames(membership_mat)
 	colnames(consensus_mat) = rownames(membership_mat)

 	class_df = data.frame(
 		class = class_ids,
 		entropy = apply(membership_mat, 1, entropy),
 		stringsAsFactors = FALSE
 	)
 	rownames(class_df) = rownames(membership_each)

 	if(length(unique(class_ids)) == 1) {
 		class_df$silhouette = rep(0, length(class_ids))
 	} else {
		class_df$silhouette = silhouette(class_ids, dist(t(consensus_mat)))[, "sil_width"]
	}

	if(length(class_df$silhouette)*0.05 > 1) {
		l = class_df$silhouette >= quantile(class_df$silhouette, 0.05)
	} else {
		l = rep(TRUE, length(class_df$silhouette))
	}
	stat = list(
		ecdf = ecdf(consensus_mat[lower.tri(consensus_mat)]),
		"1-PAC" = 1 - PAC(consensus_mat[l, l, drop = FALSE]),
		mean_silhouette = mean(class_df$silhouette),
		concordance = concordance(membership_each, class_ids),
		cophcor = cophcor(consensus_mat),
		aPAC = aPAC(consensus_mat),
		FCC = FCC(consensus_mat[l, l, drop = FALSE])
	)

	if(verbose) qqcat("  - 1-PAC for the merged 2-group is @{stat$'1-PAC'}\n")
	# an additional metric for determine "best k"
	f = stat$ecdf
	x = seq(0, 1, length = 1000)
	n = length(x)
	stat$area_increased = sum((x[2:n] - x[1:(n-1)])*f(x[2:n]))

	n_sample = nrow(class_df)
	for(method in c("Rand", "Jaccard")) {
		cl1 = rep(1, n_sample)
		cl2 = class_df$class
		stat[[method]] = cl_agreement(as.cl_hard_partition(cl1), as.cl_hard_partition(cl2), method = method)[[1]]
	}

	param$k = 2

	return(list(
		class_df = class_df, 
		membership = membership_mat, 
		consensus = consensus_mat, 
		param = param, 
		membership_each = membership_each,
		stat = stat
	))
}
