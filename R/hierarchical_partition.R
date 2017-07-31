
HierarchicalPartition = setClass("HierarchicalPartition",
    slots = list(
        list = "list",
        hierarchy = "matrix",
        subgroup = "character"
    )
)

hierarchical_partition = function(data, column_index = seq_len(ncol(data)), 
	top_method = "MAD", partition_method = "kmeans",
	PAC_cutoff = 0.2, min_samples = 6, k = 2:4, ...) {
	
	.h_obj = new.env()
	.h_obj$hierarchy = matrix(nrow = 0, ncol = 2)
	.h_obj$list = list()
	.h_obj$subgroup = rep("", length(column_index))

	.hierarchical_partition = function(.env, column_index, PAC_cutoff = 0.2, node_id = '0', min_samples = 6, k = 2:4, ...) {

		if(node_id != "0") cat("\n")
		qqcat("submatrix with @{length(column_index)} columns and @{nrow(data)} rows, node_id: @{node_id}.\n")
	   	.env$all_value_list = NULL
	   	.env$column_index = column_index
		part = consensus_partition(verbose = FALSE, .env = .env, k = k, 
			top_method = top_method, partition_method = partition_method, ...)

	    .h_obj$list[[node_id]] = part
	    .h_obj$subgroup[column_index] = node_id

	    best_k = get_best_k(part)
	    cl = get_class(part, k = best_k)

		oe = try(sig <- get_signatures(part, k = best_k, plot = FALSE, verbose = FALSE))
		if(inherits(oe, "try-error")) {
			cat("caught an error, stop.\n")
			return(NULL)
		}
		if(length(sig$group) < 100) {
			qqcat("Number of signatures are too small (< 100), stop.\n")
			return(NULL)
		}

	    mat = .env$data[, column_index, drop = FALSE]
	    mat = t(scale(t(mat)))
	    mat = mat[order(.env$all_value_list[[1]], decreasing = TRUE)[1:max(part@top_n)], , drop = FALSE]
	    s = most_tight_subgroup(mat, cl$class)
	    set1 = which(cl$class == s)
	    set2 = which(cl$class != s)

	    qqcat("best k = @{best_k}, most tight subgroup: @{s}\n")
		PAC_score = get_stat(part, k = best_k)[, "PAC"]
	    if(PAC_score > PAC_cutoff) {
	    	qqcat("PAC score too big @{PAC_score}, stop.\n")
	    	return(NULL)
	    }

	    # dist_decrease = mean_dist_decrease(mat, set1, set2)
	    # if(dist_decrease < reduce) {
	    # 	qqcat("mean distance does not decrease too much @{dist_decrease}, stop.\n")
	    # 	return(NULL)
	    # } else {

	    	if(length(set1) <= min_samples || length(set2) <= min_samples) {
	    		cat("subgroups have too few columns, stop.\n")
	    		return(NULL)
	    	}

	    	qqcat("partitioned into two subgroups with @{length(set1)} and @{length(set2)} columns.\n")
	    	# insert the two subgroups into the hierarchy
	    	sub_node_1 = paste0(node_id, s)
	    	sub_node_2 = paste0(node_id, "0")

	    	.h_obj$hierarchy = rbind(.h_obj$hierarchy, c(node_id, sub_node_1))
	    	.h_obj$hierarchy = rbind(.h_obj$hierarchy, c(node_id, sub_node_2))

	    	if(length(set1) > min_samples) {
	    		.hierarchical_partition(.env, column_index = column_index[set1], node_id = sub_node_1, 
	    			PAC_cutoff = PAC_cutoff, min_samples = min_samples, k = k, ...)
	    	}

	    	if(length(set2) > min_samples) {
	    		.hierarchical_partition(.env, column_index = column_index[set2], node_id = sub_node_2,
	    			PAC_cutoff = PAC_cutoff, min_samples = min_samples, k = k, ...)
	    	}
	    # }

	    return(NULL)
	}

	.env = new.env()
	.env$data = data
	.hierarchical_partition(.env = .env, column_index = seq_len(ncol(data)), PAC_cutoff = PAC_cutoff, min_samples = min_samples, node_id = "0", k = k, ...)

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

setMethod("get_class",
	signature = "HierarchicalPartition",
	definition = function(object) {

	object@subgroup
})

setMethod("show",
	signature = "HierarchicalPartition",
	definition = function(object) {

	qqcat("A 'HierarchicalPartition' object with '@{object@list[[1]]@top_method}:@{object@list[[1]]@partition_method}' method.\n")
})

setMethod("get_signatures",
	signature = "HierarchicalPartition",
	definition = function(object) {

	all_subgroup_labels = unique(object@subgroup)
	sig_lt = list()
	for(nm in all_subgroup_labels) {
		nc = nchar(nm)
		which_group = as.numeric(substr(nm, nc, nc))
		parent = substr(nm, 1, nc - 1)
		obj = object@list[[parent]]
		best_k = get_best_k(obj)
		qqcat("get signatures at node @{parent} with @{best_k} subgroups.\n")
		sig = get_signatures(obj, k = best_k, verbose = FALSE, plot = FALSE)$group
		if(which_group != 0) {
			sig_lt[[nm]] = names(sig[sig == which_group])
		} else {
			l = grepl(qq("^@{parent}[1-9]$"), all_subgroup_labels)
			if(any(l)) {
				px = as.numeric(substr(all_subgroup_labels[l], nc, nc))
				sig_lt[[nm]] = names(sig[! sig %in% px])
			}
		}
	}
	sig_lt
})

setMethod(f = "collect_classes",
	signature = "HierarchicalPartition",
	definition = function(object, anno = object@list[[1]]@known_anno, 
	anno_col = if(missing(anno)) object@list[[1]]ject@known_col else NULL,
	...) {

	data = object@list[[1]]@.env$data
	cl = get_class(object)

	rm = do.call("cbind", tapply(1:ncol(data), cl, function(ind) rowMeans(data[, ind])))
	fake_mat = matrix(nrow = nrow(data), ncol = ncol(data))
	colnames(fake_mat) = colnames(data)
	rownames(fake_mat) = rownames(data)
	for(nm in colnames(rm)) {
		l = cl == nm
		fake_mat[, l] = rm[, nm]
	}
	# signature_pool = unique(unlist(signature_list))
	# fake_mat = fake_mat[signature_pool, ]
	hc = as.dendrogram(hclust(dist(t(fake_mat))))
	hc = reorder(hc, colMeans(fake_mat))

	ht_list = Heatmap(cl, name = "subgroups", width = unit(1, "cm"), 
		row_title_rot = 0, cluster_rows = hc, row_dend_width = unit(2, "cm"))
	if(!is.null(anno)) {
		if(is.atomic(anno)) {
			anno_nm = deparse(substitute(anno))
			anno = data.frame(anno)
			colnames(anno) = anno_nm
			if(!is.null(anno_col)) {
				anno_col = list(anno_col)
				names(anno_col) = anno_nm
			}
		}
		if(is.null(anno_col))
			ht_list = ht_list + rowAnnotation(df = anno, show_annotation_name = TRUE,
				annotation_name_side = "bottom", width = unit(ncol(anno)*5, "mm"))
		else {
			ht_list = ht_list + rowAnnotation(df = anno, col = anno_col, show_annotation_name = TRUE,
				annotation_name_side = "bottom", width = unit(ncol(anno)*5, "mm"))
		}
	}
	draw(ht_list)
})

