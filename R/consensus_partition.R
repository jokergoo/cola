
consensus_partition = function(data,
	top_method = "sd",
	top_n = seq(min(2000, round(nrow(data)*0.2)), min(c(6000, round(nrow(data)*0.6))), length.out = 5),
	partition_method = "kmeans",
	k = 2:6, p_sampling = 0.8,
	partition_repeat = 50,
	partition_param = list(),
	known = NULL, .env) {

	if(missing(.env)) {
		if(is.data.frame(data)) data = as.matrix(data)
		if(is.null(rownames(data))) rownames(data) = seq_len(nrow(data))

		l = rowSds(data) == 0
		data = data[!l, , drop = FALSE]
		qqcat("removed @{sum(l)} rows with sd = 0\n")

		.env = new.env()
		.env$data = data
	} else if(is.null(.env$data)) {
		if(is.data.frame(data)) data = as.matrix(data)
		if(is.null(rownames(data))) rownames(data) = seq_len(nrow(data))

		l = rowSds(data) == 0
		data = data[!l, , drop = FALSE]
		qqcat("removed @{sum(l)} rows with sd = 0\n")

		.env$data = data
	} else {
		data = .env$data
	}
	
	top_n = round(top_n)

	partition_fun = get_partition_fun(partition_method, partition_param)

	get_value_fun = get_top_value_fun(top_method)

	param = data.frame(top_n = numeric(0), k = numeric(0), n_row = numeric(0))
	partition_list = list()
	
	# .env is already defined
	if(is.null(.env$all_value_list)) {
		all_value = get_value_fun(data)
		all_value[is.na(all_value)] = -Inf
		.env$all_value_list = list()
		.env$all_value_list[[top_method]] = all_value
	} else if(is.null(.env$all_value_list[[top_method]])) {
		all_value = get_value_fun(data)
		all_value[is.na(all_value)] = -Inf
		.env$all_value_list[[top_method]] = all_value
	} else {
		all_value = .env$all_value_list[[top_method]]
	}

	scale_method = attr(partition_fun, "scale_method")
	if("normal" %in% scale_method) {
		cat("rows are scaled before sent to partitioning\n")
		data = t(scale(t(data)))
	} else if("positive" %in% scale_method) {
		cat("rows are scaled before sent to partitioning\n")
		row_min = rowMeans(data)
		row_max = rowMaxs(data)
		row_range = row_max - row_min
		data = apply(data, 2, function(x) (x - row_min)/row_range)
	}

	for(i in seq_along(top_n)) {
		qqcat("get top @{top_n[i]} by @{top_method} method\n")
		ind = order(all_value, decreasing = TRUE)[1:top_n[i]]

		for(j in 1:partition_repeat) {
			ind_sub = sample(ind, round(p_sampling*length(ind)))
			mat = data[ind_sub, , drop = FALSE]

			for(y in k) {
				if(interactive()) cat(strrep("\b", 100))
				if(interactive()) qqcat("  [k = @{y}] @{partition_method} repeated for @{j}th sampling from top @{top_n[i]} rows.")
				partition_list = c(partition_list, list(list(partition_fun(mat, y))))
				param = rbind(param, data.frame(top_n = top_n[i], k = y, n_row = nrow(mat)))
			}
		}
		if(interactive()) cat("\n")
	}

	construct_consensus_object = function(param, partition_list, k, prefix = "  ") {

		partition_list = do.call("c", partition_list)
		partition_list = cl_ensemble(list = partition_list)

		qqcat("@{prefix}merging @{length(partition_list)} partitions into an ensemble object.\n")
		partition_consensus = cl_consensus(partition_list)

		membership_mat = cl_membership(partition_consensus)

		class(membership_mat) = "matrix"
		colnames(membership_mat) = paste0("p", 1:ncol(membership_mat))

		qqcat("@{prefix}calculate consensus matrix for samples clustered in a same group.\n")
		group_mat = do.call("cbind", lapply(partition_list, function(x) {
			membership = cl_membership(x)
			v = numeric(nrow(membership))
			for(i in seq_len(ncol(membership))) {
				v[as.logical(membership[, i])] = i
			}
			return(v)
		}))
		consensus_mat = matrix(1, nrow = nrow(membership_mat), ncol = nrow(membership_mat))
		for(i in 1:(nrow(group_mat)-1)) {
			for(j in (i+1):nrow(group_mat)) {
				consensus_mat[i, j] = sum(group_mat[i, ] == group_mat[j, ])/ncol(group_mat)
				consensus_mat[j, i] = consensus_mat[i, j]
			}
	 	}
	 	rownames(consensus_mat) = rownames(membership_mat)
	 	colnames(consensus_mat) = rownames(membership_mat)

	 	class_ids = as.vector(cl_class_ids(partition_consensus))
	 	if(length(unique(class_ids)) == 1) {
	 		classification = data.frame(class = class_ids, entropy = apply(membership_mat, 1, entropy), 
				silhouette = rep(1, length(class_ids)))
	 	} else {
			classification = data.frame(class = class_ids, entropy = apply(membership_mat, 1, entropy), 
				silhouette = silhouette(class_ids, dist(t(consensus_mat)))[, "sil_width"])
		}
		rownames(classification) = colnames(data)
		unique_class = sort(unique(class_ids))
		n_class = length(unique_class)

		suppressWarnings(class_color <- structure(brewer.pal(k, "Set2")[1:k], names = 1:k))

		concordance_to_known = NULL
		if(!is.null(known)) {
			concordance_to_known = cl_dissimilarity(cl_ensemble(as.cl_partition(classification$class), 
				                                                as.cl_partition(known)))
		}

		return(list(classification = classification, class_color = class_color, membership = membership_mat, 
			consensus = consensus_mat, param = param, PAC = PAC(membership_mat), membership_each = group_mat,
			concordance_to_known = concordance_to_known,
			ecdf = ecdf(membership_mat[lower.tri(membership_mat)])))
	}

	object_list = lapply(k, function(y) {
		l = param$k == y
		qqcat("wrapping results for k = @{y}\n")
		construct_consensus_object(param[l, ], partition_list[l], y)
	})
	names(object_list) = as.character(k)

	rm(partition_list)
	gc(verbose = FALSE)

	if(!is.null(known)) {
		known_color = structure(seq_along(unique(known)), names = unique(known))
	} else {
		known_color = NULL
	}
	res = list(object_list = object_list, k = k, partition_method = partition_method, top_method = top_method, .env = .env,
		known = known, known_color = known_color)
	class(res) = c("consensus_partition", "list")
	return(invisible(res))
}

print.consensus_partition = function(x, ...) {
	qqcat("top rows are extracted by '@{x$top_method}' method.\n")
	qqcat("Subgroups are detected by '@{x$partition_method}' method.\n")
	qqcat("Number of partitionings are tried for k = @{paste(x$k, collapse = ', ')}\n")
}

plot_ecdf = function(res, ...) {
	plot(NULL, xlim = c(0, 1), ylim = c(0, 1), main = "ECDF of consensus matrix", xlab = "p")
	for(i in 1:length(res$object_list)) {
		f = res$object_list[[i]]$ecdf
		x = seq(0, 1, length = 100)
		lines(x, f(x), col = i)
	}
	legend("bottomright", pch = 15, legend = paste0("k = ", res$k), col = res$k-1)
}

select_k = function(res, plot = 0) {
	op = par(no.readonly = TRUE)

	if(plot == 0) {
		par(mfrow = c(2, 2))
	}

	if(plot %in% c(0, 1)) {
		plot_ecdf(res)
	}

	if(plot %in% c(0, 2)) {
		# plot(res$k, sapply(res$object_list, function(x) x$PAC),type = "b", xlab = "k, number of clusters", ylab = "PAC")
		plot(res$k, sapply(res$object_list, function(x) cophcor(x$consensus)), type = "b", xlab = "k, number of clusters", ylab = "Cophenetic Correlation Coefficient")
	}

	if(plot %in% c(0, 3)) {
		mean_silhouette = sapply(res$object_list, function(x) {
			y = x$classification$silhouette
			mean(y[y > 0])
		})
		plot(res$k, mean_silhouette, type = "b", xlab = "k, number of clusters", ylab = "mean silhouette")
	}

	if(plot %in% c(0, 4)) {
		ak = sapply(res$object_list, function(obj) {
			consensus_mat = obj$consensus
			x = sort(consensus_mat[lower.tri(consensus_mat)])
			Fn = ecdf(x)
			x = c(0, x, 1)
			n = length(x)
			sum((x[2:n] - x[1:(n-1)])*Fn(x[2:n]))
		})
		delta_k = ak
		for(i in 2:length(res$k)) {
			delta_k[i] = (ak[i] - ak[i-1])/ak[i-1]
		}
		plot(res$k, delta_k, type = "b", xlab = "k, number of clusters", ylab = "Proportion increase")
	}
	par(op)
}

consensus_heatmap = function(res, k, show_legend = TRUE,
	annotation = data.frame(known = res$known),
	annotation_color = list(known = res$known_color)) {

	i = which(res$k == k)

	obj = res$object_list[[i]]

	class_ids = obj$classification$class
	class_color = obj$class_color

	consensus_mat = obj$consensus

	mat_col_od = column_order_by_group(class_ids, consensus_mat)

	membership_mat = obj$membership

	ht_list = Heatmap(membership_mat, name = "membership", cluster_columns = FALSE, show_row_names = FALSE, 
		width = unit(5, "mm")*k, col = colorRamp2(c(0, 1), c("white", "red"))) + 
	Heatmap(obj$classification$entropy, name = "entropy", width = unit(5, "mm"), 
		show_row_names = FALSE, col = colorRamp2(c(0, 1), c("white", "orange")),
		heatmap_legend_param = list(at = c(0, 0.5, 1), labels = c("hetergeneous", "", "homogeneous"),
			color_bar = "continuous", legend_height = unit(2, "cm"))) +
	Heatmap(obj$classification$silhouette, name = "silhouette", width = unit(5, "mm"), 
		show_row_names = FALSE, col = colorRamp2(c(0, 1), c("white", "purple"))) +
	Heatmap(class_ids, name = "class", col = class_color, show_row_names = FALSE, width = unit(5, "mm"))
	
	ht_list = ht_list +	Heatmap(consensus_mat, name = "consensus", show_row_names = FALSE, show_row_dend = FALSE,
		col = colorRamp2(c(0, 1), c("white", "blue")), row_order = mat_col_od, column_order = mat_col_od,
		cluster_rows = FALSE, cluster_columns = FALSE, show_column_name = FALSE)
	
	if(nrow(annotation)) {
		ht_list = ht_list + rowAnnotation(df = annotation, col = annotation_color, show_annotation_name = TRUE,
			annotation_name_side = "bottom", width = unit(ncol(annotation)*5, "mm"))
	}
	draw(ht_list, main_heatmap = "consensus", column_title = qq("consensus @{res$partition_method} with @{k} groups from @{nrow(obj$param)} partitions"),
		show_heatmap_legend = show_legend, show_annotation_legend = show_legend)

	return(invisible(consensus_mat))
}

get_class = function(x, ...) {
	UseMethod("get_class", x)
}

get_class.consensus_partition = function(x, k, ...) {
	res = x
	i = which(res$k == k)
	cbind(as.data.frame(res$object_list[[i]]$membership), 
		silhouette = res$object_list[[i]]$classification$silhouette,
		class = res$object_list[[i]]$classification$class)
}

get_class.run_all = function(x, k, ...) {
	res = x
	partition_list = NULL
	mean_cophcor = NULL
	mean_silhouette = NULL
	reference_class = NULL
	for(tm in res$top_method) {
		for(pm in res$partition_method) {
			nm = paste0(tm, ":", pm)
			obj = res$list[[nm]]
			ik = which(obj$k == k)

			membership = obj$object_list[[ik]]$membership
			if(is.null(reference_class)) {
	        	reference_class = obj$object_list[[ik]]$classification$class
	        } else {
	        	map = relabel_class(obj$object_list[[ik]]$classification$class, reference_class, 1:k)
	        	map2 = structure(names(map), names = map)
	        	membership = membership[, as.numeric(map2[as.character(1:k)]) ]
				colnames(membership) = paste0("p", 1:k)
			}

			partition_list = c(partition_list, list(as.cl_partition(membership)))
			mean_cophcor = c(mean_cophcor, cophcor(obj$object_list[[ik]]$consensus))
			mean_silhouette = c(mean_silhouette, mean(obj$object_list[[ik]]$classification$silhouette))
		}
	}
	consensus = cl_consensus(cl_ensemble(list = partition_list), weights = 1)
	m = cl_membership(consensus)
	class(m) = "matrix"
	colnames(m) = paste0("p", 1:k)
	class = as.vector(cl_class_ids(consensus))
	cbind(as.data.frame(m), class = class)
}

membership_heatmap = function(res, k, show_legend = TRUE, 
	annotation = data.frame(known = res$known),
	annotation_color = list(known = res$known_color)) {

	ik = which(res$k == k)
	obj = res$object_list[[ik]]
	
	class_ids = obj$classification$class
	class_color = obj$class_color

	col_list = rep(list(colorRamp2(c(0, 1), c("white", "red"))), k)
	names(col_list) = colnames(obj$membership)

	m = obj$membership_each
	m2 = matrix(nrow = nrow(m), ncol = ncol(m))
	for(i in seq_len(ncol(m))) {
		map = relabel_class(m[, i], class_ids)
		m2[, i] = map[as.character(m[, i])]
	}

	m2 = as.numeric(m2)
	dim(m2) = dim(m)

	m3 = t(m2)
	suppressWarnings(col <- structure(brewer.pal(k, "Set1"), names = 1:k))

	mat_col_od = do.call("c", lapply(sort(unique(class_ids)), function(le) {
		m = m3[, class_ids == le, drop = FALSE]
		if(ncol(m) == 1) {
			which(class_ids == le)
		} else {
			hc1 = hclust(dist(t(m)))
			oe = try({ 
				hc1 = as.hclust(reorder(as.dendrogram(hc1), colSums(m)))
			}, silent = TRUE)
			col_od1 = hc1$order
			which(class_ids == le)[col_od1]
		}
	}))

	if(nrow(annotation) == 0) {
		bottom_anno = NULL
	} else {
		bottom_anno = HeatmapAnnotation(df = annotation, col = annotation_color,
			show_annotation_name = TRUE, annotation_name_side = "right")
	}

	n_row_level = unique(obj$param$n_row)
	n_row_col = structure(brewer.pal(length(n_row_level), "Set2"), names = n_row_level)
	ht = Heatmap(m3, name = "cluster", show_row_dend = FALSE, show_column_dend = FALSE, col = col,
		column_title = qq("membership heatmap, k = @{k}"), column_order = mat_col_od, cluster_columns = FALSE,
		row_title = qq("@{nrow(m2)} samplings"),
		split = obj$param$n_row,
		top_annotation = HeatmapAnnotation(df = as.data.frame(obj$membership),
			class = class_ids, col = c(list(class = class_color), col_list),
			show_annotation_name = TRUE, annotation_name_side = "right",
			show_legend = c(TRUE, rep(FALSE, k - 1), TRUE)),
		bottom_annotation = bottom_anno,
		combined_name_fun = function(x) paste0(x, " rows")
		) + 
	Heatmap(as.character(obj$param$n_row), name = "n_row", col = n_row_col,
		width = unit(5, "mm"), show_row_names = FALSE)

	draw(ht, row_title = qq("@{round(ncol(m2)/length(n_row_level))} x @{length(n_row_level)} random samplings"),
		show_heatmap_legend = show_legend, show_annotation_legend = show_legend)

	return(invisible(m3))
}

top_rows_overlap = function(res_list, top_n = 2000, type = c("venn", "correspondance"), show_heatmap = FALSE) {

	all_value_list = res_list$list[[1]]$.env$all_value_list

	type = match.arg(type)

    lt = lapply(all_value_list, function(x) order(x, decreasing = TRUE)[1:top_n])
    
    if(type == "venn") {
   		venn_euler(lt, main = qq("top @{top_n} rows"))
	} else if(type == "correspondance") {
		correspond_between_rankings(all_value_list, top_n = top_n)
	}

    if(show_heatmap) {
    	for(i in seq_along(lt)) {
			if(dev.interactive()) {
				readline(prompt = "press enter to load next plot: ")
			}
    		mat = res_list$.env$data[lt[[i]], ]
    		draw(Heatmap(t(scale(t(mat))), show_row_names = FALSE, column_title = qq("top @{length(lt[[i]])} rows of @{res_list$top_method[i]}")))
    	}
    }
}

venn_euler = function(lt, ...) {

	if(!requireNamespace("venneuler")) {
		stop("You need to install venneuler package.")
	}
    foo = venn(lt, show.plot = FALSE)
    foo = foo[-1, ]
    set = foo[, "num"]
    category = foo[, -1]
    names(set) = apply(category, 1, function(x) {
        paste(colnames(category)[as.logical(x)], collapse = "&")
    })
    plot(getFromNamespace("venneuler", "venneuler")(set), ...)
}
