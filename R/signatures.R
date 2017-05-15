
#' Get signature rows
#'
#' @param res A `consensus_partition` class object. The object can be returned
#'        from `get_single_run()`.
#' @param k number of partitions
#' @param silhouette_cutoff cutoff for silhouette values. Columns with values 
#'        less than it are not used for finding signature rows.
#' @param fdr_cutoff cutoff for fdr of the difference between subgroups.
#' @param scale_rows whether apply row scaling when making the heatmap.
#' @param annotation a data frame which contains annotation of columns.
#' @param annotation_color colors for the annotations.
#' @param show_legend whether draw the legends on the heatmap.
#' @param ... other arguments
#' 
#' @details 
#' Basically the function apply test for the difference of subgroups in every
#' row. Also, to call it a signature for a given subgroup, the values in the
#' corresponding subgroup should have the highest mean value compared to all
#' other subgroups. The minimal p-value compared to all other subgroups is taken
#' as the p-value of the row and used for FDR calculation.
#'
#' @return 
#' A list of three elements:
#' -`mat` the matrix for the signatures
#' -`fdr` FDR for rows
#' -`gropu` subgroups that the rows are significant for
#' 
#' @export
#' @import GetoptLong
#' @import stats
#' @import ComplexHeatmap
#' @import circlize
get_signatures = function(res, k,
	silhouette_cutoff = 0.5, 
	fdr_cutoff = 0.05, 
	scale_rows = TRUE,
	annotation = data.frame(known = res$known),
	annotation_color = list(known = res$known_color),
	show_legend = TRUE,
	...) {

	i = which(res$k == k)

	obj = res$object_list[[i]]

	data = res$.env$data
	class_ids = res$classification$class
	class_color = obj$class_color

	l = obj$classification$silhouette >= silhouette_cutoff
	data2 = data[, l]
	class = obj$classification$class[l]
	column_used_index = which(l)
	tb = table(class)
	l = as.character(class) %in% names(which(tb <= 1))
	data2 = data2[, !l]
	class = class[!l]
	column_used_index = column_used_index[!l]
	column_used_logical = rep(FALSE, ncol(data))
	column_used_logical[column_used_index] = TRUE
	has_ambiguous = sum(!column_used_logical)

	qqcat("@{length(class)}/@{length(obj$classification$class)} samples (in @{length(unique(class))} classes) remain after filtering by silhouette (>= @{silhouette_cutoff}).\n")

	nm = paste0(res$partition_method, "_", res$top_method, "_signatures_", k)
	if(is.null(res$.env[[nm]])) {
		res$.env[[nm]] = list()	
	}

	if(sum(tb > 1) <= 1) {
		stop("not enough samples.")
	}
	if(length(unique(class)) <= 1) {
		stop("not enough classes.")
	}

	find_signatures = TRUE
	if(!is.null(res$.env[[nm]]$silhouette_cutoff)) {
		if(abs(res$.env[[nm]]$silhouette_cutoff - silhouette_cutoff) < 1e-6) {
			p = res$.env[[nm]]$p
			group = res$.env[[nm]]$group
			find_signatures = FALSE
		}
	}

	if(find_signatures) {
		p = numeric(nrow(data2))
		group = NULL
		for(i in seq_len(nrow(data2))) {
			if(interactive()) cat(strrep("\b", 100))
			if(interactive()) qqcat("pairwise t-test: @{i}/@{nrow(data2)}")
			x = data2[i, ]
			# df = data.frame(value = x, class = class)
			# oneway.test(value ~ class, data = df)$p.value

			group_mean = tapply(x, class, mean)
			max_group = names(which.max(group_mean))
			fit = pairwise.t.test(x, class)
			pmat = fit$p.value
			all_p = c(pmat[, which(colnames(pmat) == max_group)], pmat[which(rownames(pmat) == max_group), ])
			all_p = all_p[!is.na(all_p)]
			if(length(all_p)) {
				# if(all(all_p < 0.05)) {
					p[i] = min(all_p)
				# } else {
				# 	p[i] = NA
				# }
			} else {
				p[i] = NA
			}
			group[i] = max_group
		}
		if(interactive()) cat("\n")
	}

	fdr = p.adjust(p, "BH")
	fdr[is.na(fdr)] = Inf

	res$.env[[nm]]$silhouette_cutoff = silhouette_cutoff
	res$.env[[nm]]$p = p
	res$.env[[nm]]$group = group

	col_list = rep(list(colorRamp2(c(0, 1), c("white", "red"))), k)
	names(col_list) = colnames(obj$membership)

	if(is.null(fdr_cutoff)) {
		fdr_cutoff_v = c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)
		i = which.min(abs(500 - sapply(fdr_cutoff_v, function(x) sum(fdr < x))))
		fdr_cutoff = fdr_cutoff_v[i]
	}
	mat = data[fdr < fdr_cutoff, ]
	fdr2 = fdr[fdr < fdr_cutoff]
	group2 = group[fdr < fdr_cutoff]
	names(group2) = rownames(mat)
	mat_return = list(mat = mat, fdr = fdr2, group = group2)
	qqcat("@{nrow(mat)} signatures under fdr < @{fdr_cutoff}\n")

	p = sort(p)
	f = ecdf(p)
	p = c(0, p, 1)
	n2 = length(p)
	res$.env[[nm]]$AUC = sum((p[2:n2] - p[1:(n2-1)])*f(p[-n2]))

	more_than_5k = FALSE
	if(nrow(mat) > 5000) {
		more_than_5k = TRUE
		mat1 = mat[order(fdr2)[1:5000], column_used_logical, drop = FALSE]
		mat2 = mat[order(fdr2)[1:5000], !column_used_logical, drop = FALSE]
		group2 = group2[order(fdr2)[1:5000]]
		cat("Only take top 5000 signatures with highest fdr\n")
	} else {
		mat1 = mat[, column_used_logical, drop = FALSE]
		mat2 = mat[, !column_used_logical, drop = FALSE]
		
	}
	base_mean = rowMeans(mat1)

	if(nrow(annotation) == 0) {
		bottom_anno1 = NULL
	} else {
		bottom_anno1 = HeatmapAnnotation(df = annotation[column_used_logical, ,drop = FALSE], col = annotation_color,
			show_annotation_name = !has_ambiguous, annotation_name_side = "right")
	}

	if(scale_rows) {
		scaled_mean = base_mean
		scaled_sd = rowSds(mat1)
		scaled_mat1 = t(scale(t(mat1)))
		scaled_mat2 = mat2
		if(has_ambiguous) {
			for(i in seq_len(nrow(mat2))) {
				scaled_mat2[i, ] = (scaled_mat2[i, ] - scaled_mean[i])/scaled_sd[i]
			}
		}

		use_mat1 = scaled_mat1
		use_mat2 = scaled_mat2
		mat_range = quantile(abs(scaled_mat1), 0.95)
		col_fun = colorRamp2(c(-mat_range, 0, mat_range), c("green", "white", "red"))
		heatmap_name = "scaled_expr"
	} else {
		use_mat1 = mat1
		use_mat2 = mat2
		mat_range = quantile(mat1, c(0.05, 0.95))
		col_fun = colorRamp2(c(mat_range[1], mean(mat_range), mat_range[2]), c("blue", "white", "red"))
		heatmap_name = "expr"
	}


	mat_col_od1 = column_order_by_group(class, use_mat1)

	if(has_ambiguous) {
		class2 = obj$classification$class[!column_used_logical]
		mat_col_od2 = column_order_by_group(class2, use_mat2)

		if(nrow(annotation) == 0) {
			bottom_anno2 = NULL
		} else { 
			bottom_anno2 = HeatmapAnnotation(df = annotation[!column_used_logical, ,drop = FALSE], col = annotation_color, 
				show_annotation_name = TRUE, annotation_name_side = "right", show_legend = FALSE)
		}
	}
	silhouette_range = range(obj$classification$silhouette)
	silhouette_range[2] = 1

	group2 = factor(group2, levels = sort(unique(group2)))
	ht_list = Heatmap(group2, name = "group", show_row_names = FALSE, width = unit(5, "mm"), col = class_color)

	ht_list = ht_list + Heatmap(use_mat1, name = heatmap_name, col = col_fun,
		top_annotation = HeatmapAnnotation(as.data.frame(obj$membership)[column_used_logical, ],
			class = obj$classification$class[column_used_logical],
			silhouette = anno_barplot(obj$classification$silhouette[column_used_logical], ylim = silhouette_range,
				gp = gpar(fill = ifelse(obj$classification$silhouette[column_used_logical] >= silhouette_cutoff, "grey", "#EEEEEE"),
					      col = ifelse(obj$classification$silhouette[column_used_logical] >= silhouette_cutoff, "black", NA)),
				baseline = 0, axis = !has_ambiguous, axis_side = "right"),
			col = c(list(class = class_color), col_list),
			annotation_height = unit(c(rep(4, k+1), 15), "mm"),
			show_annotation_name = if(has_ambiguous) FALSE else c(rep(TRUE, k+1), FALSE),
			annotation_name_side = "right",
			show_legend = c(TRUE, rep(FALSE, k - 1), TRUE)),
		cluster_columns = FALSE, column_order = mat_col_od1,
		show_row_names = FALSE, show_row_dend = FALSE, column_title = "confident samples",
		bottom_annotation = bottom_anno1, show_column_names = FALSE, split = group2, combined_name_fun = NULL)
 	if(scale_rows) {
		ht_list = ht_list + Heatmap(base_mean, show_row_names = FALSE, name = "base_mean", width = unit(5, "mm"))
	}
	if(has_ambiguous) {
		ht_list = ht_list + Heatmap(use_mat2, name = paste0(heatmap_name, 2), col = col_fun,
			top_annotation = HeatmapAnnotation(as.data.frame(obj$membership)[!column_used_logical, ],
				class = obj$classification$class[!column_used_logical],
				silhouette2 = anno_barplot(obj$classification$silhouette[!column_used_logical], ylim = silhouette_range,
					gp = gpar(fill = ifelse(obj$classification$silhouette[!column_used_logical] >= silhouette_cutoff, "grey", "grey"),
					      col = ifelse(obj$classification$silhouette[!column_used_logical] >= silhouette_cutoff, "black", NA)),
					baseline = 0, axis = TRUE, axis_side = "right"), 
				col = c(list(class = class_color), col_list),
				annotation_height = unit(c(rep(4, k+1), 15), "mm"),
				show_annotation_name = c(rep(TRUE, k+1), FALSE),
				annotation_name_side = "right",
				show_legend = FALSE),
			cluster_columns = FALSE, column_order = mat_col_od2,
			show_row_names = FALSE, show_row_dend = FALSE, show_heatmap_legend = FALSE,
			column_title = "ambiguous samples", bottom_annotation = bottom_anno2, show_column_names = FALSE)
	}

	draw(ht_list, main_heatmap = heatmap_name, column_title = qq("@{k} subgroups, @{nrow(mat)} signatures with fdr < @{fdr_cutoff}@{ifelse(more_than_5k, ', top 5k signatures', '')}"),
		show_heatmap_legend = show_legend, show_annotation_legend = show_legend)
	# https://www.stat.berkeley.edu/~s133/Cluster2a.html
	decorate_annotation("silhouette", {
		grid.rect(gp = gpar(fill = "transparent"))
		grid.lines(c(0, 1), unit(c(silhouette_cutoff, silhouette_cutoff), "native"), gp = gpar(lty = 2, col = "#CCCCCC"))
		if(!has_ambiguous) grid.text("silhouette", x = unit(1, "npc") + unit(10, "mm"), just = "left")
	})
	if(has_ambiguous) {
		decorate_annotation("silhouette2", {
			grid.rect(gp = gpar(fill = "transparent"))
			grid.lines(c(0, 1), unit(c(silhouette_cutoff, silhouette_cutoff), "native"), gp = gpar(lty = 2, col = "#CCCCCC"))
			if(has_ambiguous) grid.text("silhouette", x = unit(1, "npc") + unit(10, "mm"), just = "left")
		})
	}

	return(invisible(mat_return))
}
