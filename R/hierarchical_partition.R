
HierarchicalPartition = setClass("HierarchicalPartition",
    slots = list(
        list = "list",
        hierarchy = "matrix",
        subgroup = "character"
    )
)

hierarchical_partition = function(data, column_index = seq_len(ncol(data)), PAC_cutoff = 0.2, min_samples = 6, ...) {
	
	.h_obj = new.env()
	.h_obj$hierarchy = matrix(nrow = 0, ncol = 2)
	.h_obj$list = list()
	.h_obj$subgroup = rep("", length(column_index))

	.hierarchical_partition = function(data, column_index = seq_len(ncol(data)), PAC_cutoff = 0.2, 
		node_id = '0', min_samples = 6, ...) {

		qqcat("\nsubmatrix with @{length(column_index)} columns, node_id: @{node_id}.\n")
	    part = consensus_partition(data, column_index = column_index, ...)

	    .h_obj$list[[node_id]] = part
	    .h_obj$subgroup[column_index] = node_id

	    k = get_best_k(part)
	    cl = get_class(part, k = k)

		oe = try(sig <- get_signatures(part, k = k, plot = FALSE))
		if(inherits(oe, "try-error")) {
			return(NULL)
		}
		if(length(sig$group) < 100) {
			qqcat("Number of signatures are too small (< 100), stop.\n")
			return(NULL)
		}

	    mat = part@.env$data[, part@.env$column_index, drop = FALSE]
	    mat = t(scale(t(mat)))
	    mat = mat[order(part@.env$all_value_list[[1]], decreasing = TRUE)[1:max(part@top_n)], , drop = FALSE]
	    s = most_tight_subgroup(mat, cl$class)
	    set1 = which(cl$class == s)
	    set2 = which(cl$class != s)

	    qqcat("best k = @{k}, most tight subgroup: @{s}\n")
		PAC_score = get_stat(part, k = k)[, "PAC"]
	    if(PAC_score > PAC_cutoff) {
	    	qqcat("PAC score too big @{PAC_score}, stop.\n")
	    	return(NULL)
	    }

	    # dist_decrease = mean_dist_decrease(mat, set1, set2)
	    # if(dist_decrease < reduce) {
	    # 	qqcat("mean distance does not decrease too much @{dist_decrease}, stop.\n")
	    # 	return(NULL)
	    # } else {

	    	if(length(set1) <= 1 || length(set2) <= 1) {
	    		return(NULL)
	    	}

	    	qqcat("partitioned into two subgroups with @{length(set1)} and @{length(set2)} columns.\n")
	    	# insert the two subgroups into the hierarchy
	    	sub_node_1 = paste0(node_id, s)
	    	sub_node_2 = paste0(node_id, "0")

	    	.h_obj$hierarchy = rbind(.h_obj$hierarchy, c(node_id, sub_node_1))
	    	.h_obj$hierarchy = rbind(.h_obj$hierarchy, c(node_id, sub_node_2))

	    	if(length(set1) > min_samples) {
	    		.hierarchical_partition(data, column_index = column_index[set1], node_id = sub_node_1, ...)
	    	}

	    	if(length(set2) > min_samples) {
	    		.hierarchical_partition(data, column_index = column_index[set2], node_id = sub_node_2, ...)
	    	}
	    # }

	    return(NULL)
	}

	.hierarchical_partition(data, column_index = column_index, PAC_cutoff = PAC_cutoff, 
		min_samples = min_samples, node_id = "0", ...)

	hp = new("HierarchicalPartition")
	hp@hierarchy = .h_obj$hierarchy
	hp@list = .h_obj$list
	hp@subgroup = .h_obj$subgroup
	names(hp@subgroup) = colnames(data)[column_index]

	return(hp)
}

most_tight_subgroup = function(mat, subgroup) {
	
	mean_dist = tapply(seq_len(ncol(mat)), subgroup, function(ind) {
		n = length(ind)
		if(n == 1) {
			return(Inf)
		}
		sum(dist(t(mat[, ind, drop = FALSE]))^2)/(n*(n-1)/2)
	})
	as.numeric(names(mean_dist[which.min(mean_dist)]))
}


# mat = matrix(rnorm(100), 10)
# mat2 = cbind(matrix(rnorm(50, mean = -1), nrow = 10), matrix(rnorm(50, mean = 1), nrow = 10))
# subset = sample(10, 10)
# mean_dist_decrease(mat, subset[1:5], subset[6:10])
mean_dist_decrease = function(mat, subset1, subset2) {
	n = ncol(mat)
	global_mean_dist = sum(dist(t(mat))^2) / (n*(n-1)/2)

	n1 = length(subset1)
	n2 = length(subset2)
	mean_dist = (sum(dist(t(mat[, subset1, drop = FALSE]))^2) +  sum(dist(t(mat[, subset2, drop = FALSE]))^2)) / ( n1*(n1-1)/2 + n2*(n2-1)/2)

	(global_mean_dist - mean_dist)/global_mean_dist
}


