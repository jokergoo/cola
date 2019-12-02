
# == title
# Get signature rows
#
# == param
# -object A `ConsensusPartition-class` object.
# -k Number of partitions.
# -silhouette_cutoff Cutoff for silhouette scores. Samples with values 
#        less than it are not used for finding signature rows. For selecting a 
#        proper silhouette cutoff, please refer to https://www.stat.berkeley.edu/~s133/Cluster2a.html#tth_tAb1.
# -fdr_cutoff Cutoff for FDR of the difference test between subgroups.
# -group_diff Cutoff for the maximal difference between group means.
# -scale_rows Whether apply row scaling when making the heatmap.
# -row_km Number of groups for performing k-means clustering on rows. By default it is automatically selected.
# -diff_method Methods to get rows which are significantly different between subgroups, see 'Details' section.
# -anno A data frame of annotations for the original matrix columns. 
#       By default it uses the annotations specified in `consensus_partition` or `run_all_consensus_partition_methods`.
# -anno_col A list of colors (color is defined as a named vector) for the annotations. If ``anno`` is a data frame,
#       ``anno_col`` should be a named list where names correspond to the column names in ``anno``.
# -internal Used internally.
# -show_row_dend Whether show row dendrogram.
# -show_column_names Whether show column names in the heatmap.
# -use_raster Internally used.
# -plot Whether to make the plot.
# -verbose Whether to print messages.
# -seed Random seed.
# -left_annotation Annotation put on the left of the heatmap. It should be a `ComplexHeatmap::HeatmapAnnotation-class` object. 
#              The number of items should be the same as the number of the original matrix rows. The subsetting to the significant 
#              rows are automatically performed on the annotation object.
# -right_annotation Annotation put on the right of the heatmap. Same format as ``left_annotation``.
# -col Colors.
# -simplify Only use internally.
# -... Other arguments.
# 
# == details 
# Basically the function applies statistical test for the difference in subgroups for every
# row. There are following methods which test significance of the difference:
#
# -ttest First it looks for the subgroup with highest mean value, compare to each of the 
#        other subgroups with t-test and take the maximum p-value. Second it looks
#        for the subgroup with lowest mean value, compare to each of the other subgroups
#        again with t-test and take the maximum p-values. Later for these two list of p-values
#        take the minimal p-value as the final p-value. 
# -samr/pamr use SAM (from samr package)/PAM (from pamr package) method to find significantly different rows between subgroups.
# -Ftest use F-test to find significantly different rows between subgroups.
# -one_vs_others For each subgroup i in each row, it uses t-test to compare samples in current 
#        subgroup to all other samples, denoted as p_i. The p-value for current row is selected as min(p_i).
#
# ``diff_method`` can also be a self-defined function. The function needs two arguments which are the matrix for the analysis
# and the predicted classes. The function should returns a vector of FDR from the difference test.
#
# == return 
# A data frame with more than two columns:
#
# -``which_row``: row index corresponding to the original matrix.
# -``fdr``: the FDR.
# -``km``: the k-means groups if ``row_km`` is set.
# -other_columns: the mean value (depending rows are scaled or not) in each subgroup.
# 
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "get_signatures",
	signature = "ConsensusPartition",
	definition = function(object, k,
	silhouette_cutoff = 0.5, 
	fdr_cutoff = cola_opt$fdr_cutoff, 
	group_diff = cola_opt$group_diff,
	scale_rows = object@scale_rows,
	row_km = NULL,
	diff_method = c("Ftest", "ttest", "samr", "pamr", "one_vs_others"),
	anno = get_anno(object), 
	anno_col = get_anno_col(object),
	internal = FALSE,
	show_row_dend = FALSE,
	show_column_names = FALSE, use_raster = TRUE,
	plot = TRUE, verbose = TRUE, seed = 888,
	left_annotation = NULL, right_annotation = NULL,
	col = if(scale_rows) c("green", "white", "red") else c("blue", "white", "red"),
	simplify = FALSE,
	...) {

	if(missing(k)) stop_wrap("k needs to be provided.")
	
	raster_resize = cola_opt$raster_resize

	class_df = get_classes(object, k)
	class_ids = class_df$class

	data = object@.env$data[, object@column_index, drop = FALSE]

	l = class_df$silhouette >= silhouette_cutoff
	data2 = data[, l, drop = FALSE]
	class = class_df$class[l]
	column_used_index = which(l)
	tb = table(class)
	l = as.character(class) %in% names(which(tb <= 1))
	data2 = data2[, !l, drop = FALSE]
	class = class[!l]
	column_used_index = column_used_index[!l]
	column_used_logical = rep(FALSE, ncol(data))
	column_used_logical[column_used_index] = TRUE
	has_ambiguous = sum(!column_used_logical)
	n_sample_used = length(class)

	if(verbose) qqcat("* @{n_sample_used}/@{nrow(class_df)} samples (in @{length(unique(class))} classes) remain after filtering by silhouette (>= @{silhouette_cutoff}).\n")

	tb = table(class)
	if(sum(tb > 1) <= 1) {
		if(plot) {
			grid.newpage()
			fontsize = convertUnit(unit(0.1, "npc"), "char", valueOnly = TRUE)*get.gpar("fontsize")$fontsize
			grid.text("not enough samples", gp = gpar(fontsize = fontsize))
		}
		return(invisible(data.frame(which_row = integer(0))))
	}
	if(length(unique(class)) <= 1) {
		if(plot) {
			grid.newpage()
			fontsize = convertUnit(unit(0.1, "npc"), "char", valueOnly = TRUE)*get.gpar("fontsize")$fontsize
			grid.text("not enough classes", gp = gpar(fontsize = fontsize))
		}
		return(invisible(data.frame(which_row = integer(0))))
	}

	do_row_clustering = TRUE
	if(inherits(diff_method, "function")) {
		if(verbose) qqcat("* calculate row difference between subgroups by user-defined function.\n")
		diff_method_fun = diff_method
		diff_method = digest(diff_method)
	} else {
		diff_method = match.arg(diff_method)
	}

	hash = digest(list(used_samples = which(l), 
		               class = class,
		               n_group = k, 
		               diff_method = diff_method,
		               column_index = object@column_index,
		               fdr_cutoff = fdr_cutoff,
		               group_diff = group_diff,
		               seed = seed),
				algo = "md5")
	nm = paste0("signature_fdr_", hash)
	if(verbose) qqcat("* cache hash: @{hash} (seed @{seed}).\n")

	find_signature = TRUE
	if(!is.null(object@.env[[nm]])) {
		if(diff_method == "samr") {
			if(object@.env[[nm]]$diff_method == "samr" && 
			   object@.env[[nm]]$n_sample_used == n_sample_used && 
			   abs(object@.env[[nm]]$fdr_cutoff - fdr_cutoff) < 1e-6) {
				fdr = object@.env[[nm]]$fdr
				find_signature = FALSE
			}
		} else if(diff_method == "pamr") {
			if(object@.env[[nm]]$diff_method == "pamr" && 
			   object@.env[[nm]]$n_sample_used == n_sample_used && 
			   abs(object@.env[[nm]]$fdr_cutoff - fdr_cutoff) < 1e-6) {
				fdr = object@.env[[nm]]$fdr
				find_signature = FALSE
			}
		} else {
			if(object@.env[[nm]]$diff_method == diff_method &&
			   object@.env[[nm]]$n_sample_used == n_sample_used) {
				fdr = object@.env[[nm]]$fdr
				find_signature = FALSE
			}
		}
	}

	if(verbose) qqcat("* calculating row difference between subgroups by @{diff_method}.\n")
	if(find_signature) {
		if(diff_method == "ttest") {
			fdr = ttest(data2, class)
		} else if(diff_method == "samr") {
			fdr = samr(data2, class, fdr.output = fdr_cutoff)
		} else if(diff_method == "Ftest") {
			fdr = Ftest(data2, class)
		} else if(diff_method == "pamr") {
			fdr = pamr(data2, class, fdr.ouput = fdr_cutoff)
		} else {
			fdr = diff_method_fun(data2, class)
		}
	} else {
		if(verbose) qqcat("  - row difference is extracted from cache.\n")
	}

	object@.env[[nm]]$diff_method = diff_method
	object@.env[[nm]]$fdr_cutoff = fdr_cutoff
	object@.env[[nm]]$fdr = fdr
	object@.env[[nm]]$n_sample_used = n_sample_used
	object@.env[[nm]]$group_diff = group_diff

	if(scale_rows && !is.null(object@.env[[nm]]$row_order_scaled)) {
		row_order = object@.env[[nm]]$row_order_scaled
		if(verbose) qqcat("  - row order for the scaled matrix is extracted from cache.\n")
		do_row_clustering = FALSE
	} else if(!scale_rows && !is.null(object@.env[[nm]]$row_order_unscaled)) {
		row_order = object@.env[[nm]]$row_order_unscaled
		if(verbose) qqcat("  - row order for the unscaled matrix is extracted from cache.\n")
		do_row_clustering = FALSE
	}

	# filter by fdr
	fdr[is.na(fdr)] = 1

	l_fdr = fdr < fdr_cutoff
	mat = data[l_fdr, , drop = FALSE]
	fdr2 = fdr[l_fdr]
	
	if(!is.null(left_annotation)) left_annotation = left_annotation[l_fdr, ]
	if(!is.null(right_annotation)) right_annotation = right_annotation[l_fdr, ]

	returned_df = data.frame(which_row = which(l_fdr), fdr = fdr2)

	# filter by group_diff
	mat1 = mat[, column_used_logical, drop = FALSE]
	if(nrow(mat) == 1) {
		group_mean = rbind(tapply(mat1, class, mean))
	} else {
		group_mean = do.call("cbind", tapply(seq_len(ncol(mat1)), class, function(ind) {
			rowMeans(mat1[, ind, drop = FALSE])
		}))
	}
	colnames(group_mean) = paste0("mean_", colnames(group_mean))
	returned_df = cbind(returned_df, group_mean)
	returned_df$group_diff = apply(group_mean, 1, function(x) max(x) - min(x))

	if(group_diff > 0) {
		l_diff = returned_df$group_diff >= group_diff
		mat = mat[l_diff, , drop = FALSE]
		mat1 = mat1[l_diff, , drop = FALSE]
		returned_df = returned_df[l_diff, , drop = FALSE]
	}

	if(scale_rows) {
		mat1_scaled = t(scale(t(mat[, column_used_logical, drop = FALSE])))
		if(nrow(mat) == 1) {
			group_mean_scaled = rbind(tapply(mat1_scaled, class, mean))
		} else {
			group_mean_scaled = do.call("cbind", tapply(seq_len(ncol(mat1_scaled)), class, function(ind) {
				rowMeans(mat1_scaled[, ind, drop = FALSE])
			}))
		}
		colnames(group_mean_scaled) = paste0("scaled_mean_", colnames(group_mean_scaled))
		returned_df = cbind(returned_df, group_mean_scaled)
	}

	returned_obj = returned_df
	rownames(returned_obj) = NULL

	attr(returned_obj, "sample_used") = column_used_logical

	## add k-means
	row_km_fit = NULL
	if(!internal) {
		if(nrow(mat1) > 10) {
			do_kmeans = TRUE
			if(scale_rows) {
				mat_for_km = t(scale(t(mat1)))
				row_km_fit = object@.env[[nm]]$row_km_fit_scaled
			} else {
				mat_for_km = mat1
				row_km_fit = object@.env[[nm]]$row_km_fit_unscaled
			}

			if(nrow(mat_for_km) > 5000) {
				set.seed(seed)
				mat_for_km2 = mat_for_km[sample(nrow(mat_for_km), 5000), , drop = FALSE]
			} else {
				mat_for_km2 = mat_for_km
			}

			if(!is.null(row_km_fit)) {
				if(is.null(row_km) || identical(as.integer(row_km), length(row_km_fit$size))) {
					returned_obj$km = apply(pdist(row_km_fit$centers, mat_for_km), 2, which.min)
					do_kmeans = FALSE
					if(verbose) qqcat("* use k-means partition that are already calculated in previous runs.\n")
				}
			}
			if(do_kmeans) {
				set.seed(seed)
				if(is.null(row_km)) {
					wss = (nrow(mat_for_km2)-1)*sum(apply(mat_for_km2,2,var))
					max_km = min(c(nrow(mat_for_km) - 1, 15))
					# if(verbose) qqcat("* apply k-means on rows with 2~@{max_km} clusters.\n")
					for (i in 2:max_km) {
						# if(verbose) qqcat("  - applying k-means with @{i} clusters.\n")
						wss[i] = sum(kmeans(mat_for_km2, centers = i, iter.max = 50)$withinss)
					}
					row_km = min(elbow_finder(1:max_km, wss)[1], knee_finder(1:max_km, wss)[1])
					if(length(unique(class)) == 1) row_km = 1
					if(length(unique(class)) == 2) row_km = min(row_km, 2)
				}
				if(row_km > 1) {
					row_km_fit = kmeans(mat_for_km2, centers = row_km)
					returned_obj$km = apply(pdist(row_km_fit$centers, mat_for_km), 2, which.min)
					if(scale_rows) {
						object@.env[[nm]]$row_km_fit_scaled = row_km_fit
					} else {
						object@.env[[nm]]$row_km_fit_unscaled = row_km_fit
					}
				}
				if(verbose) qqcat("* split rows into @{row_km} groups by k-means clustering.\n")
			}
		}
	}

	if(verbose) qqcat("* @{nrow(mat)} signatures (@{sprintf('%.1f',nrow(mat)/nrow(object)*100)}%) under fdr < @{fdr_cutoff}, group_diff > @{group_diff}.\n")

	if(nrow(mat) == 0) {
		if(plot) {
			grid.newpage()
			fontsize = convertUnit(unit(0.1, "npc"), "char", valueOnly = TRUE)*get.gpar("fontsize")$fontsize
			grid.text("no sigatures", gp = gpar(fontsize = fontsize))
		}
		return(invisible(NULL))
	}

	if(!plot) {
		return(invisible(returned_obj))
	}

	set.seed(seed)
	more_than_5k = FALSE
	if(!is.null(object@.env[[nm]]$row_index)) {
		if(verbose) qqcat("  - use the 2000 signatures what are already generated in previous runs.\n")
		row_index = object@.env[[nm]]$row_index
		mat1 = mat[row_index, column_used_logical, drop = FALSE]
		mat2 = mat[row_index, !column_used_logical, drop = FALSE]
		more_than_5k = TRUE
		if(!is.null(left_annotation)) left_annotation = left_annotation[row_index, ]
		if(!is.null(right_annotation)) right_annotation = right_annotation[row_index, ]
	} else if(nrow(mat) > 2000) {
		more_than_5k = TRUE
		row_index = sample(1:nrow(mat), 2000)
		object@.env[[nm]]$row_index = row_index
		# mat1 = mat[order(fdr2)[1:top_k_genes], column_used_logical, drop = FALSE]
		# mat2 = mat[order(fdr2)[1:top_k_genes], !column_used_logical, drop = FALSE]
		mat1 = mat[row_index, column_used_logical, drop = FALSE]
		mat2 = mat[row_index, !column_used_logical, drop = FALSE]
		# group2 = group2[order(fdr2)[1:top_k_genes]]
		if(verbose) cat(paste0("  - randomly sample 2000 signatures.\n"))
		if(!is.null(left_annotation)) left_annotation = left_annotation[row_index, ]
		if(!is.null(right_annotation)) right_annotation = right_annotation[row_index, ]
	} else {
		row_index = seq_len(nrow(mat))
		mat1 = mat[, column_used_logical, drop = FALSE]
		mat2 = mat[, !column_used_logical, drop = FALSE]
		
	}
	base_mean = rowMeans(mat1)
	if(nrow(mat) == 1) {
		group_mean = matrix(tapply(mat1, class, mean), nrow = 1)
	} else {
		group_mean = do.call("cbind", tapply(seq_len(ncol(mat1)), class, function(ind) {
			rowMeans(mat1[, ind, drop = FALSE])
		}))
	}
	rel_diff = (rowMaxs(group_mean) - rowMins(group_mean))/base_mean/2

	if(is.null(anno)) {
		bottom_anno1 = NULL
	} else {
		if(is.atomic(anno)) {
			anno_nm = deparse(substitute(anno))
			anno = data.frame(anno)
			colnames(anno) = anno_nm
			if(!is.null(anno_col)) {
				if(is.atomic(anno_col)) {
					anno_col = list(anno_col)
					names(anno_col) = anno_nm
				}
			}
		} else if(ncol(anno) == 1) {
			if(!is.null(anno_col)) {
				if(is.atomic(anno_col)) {
					anno_col = list(anno_col)
					names(anno_col) = colnames(anno)
				}
			}
		}

		if(is.null(anno_col)) {
			bottom_anno1 = HeatmapAnnotation(df = anno[column_used_logical, , drop = FALSE],
				show_annotation_name = !has_ambiguous & !internal, annotation_name_side = "right")
		} else {
			bottom_anno1 = HeatmapAnnotation(df = anno[column_used_logical, , drop = FALSE], col = anno_col,
				show_annotation_name = !has_ambiguous & !internal, annotation_name_side = "right")
		}
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
		use_mat1[is.infinite(use_mat1)] = 0
		use_mat1[is.na(use_mat1)] = 0
		use_mat2[is.infinite(use_mat2)] = 0
		use_mat2[is.na(use_mat2)] = 0
		mat_range = quantile(abs(scaled_mat1), 0.95, na.rm = TRUE)
		col_fun = colorRamp2(c(-mat_range, 0, mat_range), col)
		heatmap_name = "z-score"
	} else {
		use_mat1 = mat1
		use_mat2 = mat2
		mat_range = quantile(mat1, c(0.05, 0.95))
		col_fun = colorRamp2(c(mat_range[1], mean(mat_range), mat_range[2]), col)
		heatmap_name = "Value"
	}

	if(has_ambiguous) {
		class2 = class_df$class[!column_used_logical]

		if(is.null(anno)) {
			bottom_anno2 = NULL
		} else {
			anno_col = lapply(bottom_anno1@anno_list, function(anno) {
				if(is.null(anno@color_mapping)) {
					return(NULL)
				} else {
					if(anno@color_mapping@type == "discrete") {
						anno@color_mapping@colors
					} else {
						anno@color_mapping@col_fun
					}
				}
			})
			names(anno_col) = names(bottom_anno1@anno_list)
			anno_col = anno_col[!sapply(anno_col, is.null)]

			if(!is.null(object@anno_col)) {
				nmd = setdiff(names(object@anno_col), names(anno_col))
				if(length(nmd)) {
					anno_col[nmd] = object@anno_col[nmd]
				}
			}

			bottom_anno2 = HeatmapAnnotation(df = anno[!column_used_logical, , drop = FALSE], col = anno_col,
				show_annotation_name = !internal, annotation_name_side = "right")	
		}
	}
	silhouette_range = range(class_df$silhouette)
	silhouette_range[2] = 1

	if(verbose) qqcat("* making heatmaps for signatures.\n")

	row_split = NULL
	if(!internal) {
		if(scale_rows) {
			row_km_fit = object@.env[[nm]]$row_km_fit_scaled
		} else {
			row_km_fit = object@.env[[nm]]$row_km_fit_unscaled
		}
		if(!is.null(row_km_fit)) {
			row_split = factor(returned_obj$km[row_index], levels = sort(unique(returned_obj$km[row_index])))
		}
	}

	# group2 = factor(group2, levels = sort(unique(group2)))
	# ht_list = Heatmap(group2, name = "Group", show_row_names = FALSE, width = unit(5, "mm"), col = cola_opt$color_set_2)
	ht_list = NULL

	membership_mat = get_membership(object, k)
	prop_col_fun = colorRamp2(c(0, 1), c("white", "red"))

	if(internal) {
		ha1 = HeatmapAnnotation(Prob = membership_mat[column_used_logical, ],
				Class = class_df$class[column_used_logical],
				col = list(Class = cola_opt$color_set_2, Prob = prop_col_fun),
				show_annotation_name = !has_ambiguous & !internal,
				annotation_name_side = "right",
				show_legend = TRUE)
	} else {
		ha1 = HeatmapAnnotation(Prob = membership_mat[column_used_logical, ],
			Class = class_df$class[column_used_logical],
			silhouette = anno_barplot(class_df$silhouette[column_used_logical], ylim = silhouette_range,
				gp = gpar(fill = ifelse(class_df$silhouette[column_used_logical] >= silhouette_cutoff, "black", "#EEEEEE"),
					      col = NA),
				bar_width = 1, baseline = 0, axis = !has_ambiguous, axis_param = list(side= "right"),
				height = unit(15, "mm")),
			col = list(Class = cola_opt$color_set_2, Prob = prop_col_fun),
			show_annotation_name = !has_ambiguous & !internal,
			annotation_name_side = "right",
			show_legend = TRUE)
	}
	ht_list = ht_list + Heatmap(use_mat1, name = heatmap_name, col = col_fun,
		top_annotation = ha1, row_split = row_split,
		cluster_columns = TRUE, cluster_column_slices = FALSE, cluster_row_slices = FALSE,
		column_split = factor(class_df$class[column_used_logical], levels = sort(unique(class_df$class[column_used_logical]))), 
		show_column_dend = FALSE,
		show_row_names = FALSE, show_row_dend = show_row_dend, column_title = {if(internal) NULL else qq("@{ncol(use_mat1)} confident samples")},
		use_raster = use_raster, raster_resize = raster_resize,
		bottom_annotation = bottom_anno1, show_column_names = show_column_names, 
		left_annotation = left_annotation, right_annotation = {if(has_ambiguous) NULL else right_annotation})
 	
	all_value_positive = !any(data < 0)
 	if(scale_rows && all_value_positive && !simplify) {
		ht_list = ht_list + Heatmap(base_mean, show_row_names = FALSE, name = "base_mean", width = unit(5, "mm"), show_column_names = !internal) +
			Heatmap(rel_diff, col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), 
				show_row_names = FALSE, show_column_names = !internal, name = "rel_diff", width = unit(5, "mm"))
	}

	if(has_ambiguous) {
		if(internal) {
			ha2 = HeatmapAnnotation(Prob = membership_mat[!column_used_logical, ,drop = FALSE],
				Class = class_df$class[!column_used_logical],
				col = list(Class = cola_opt$color_set_2, Prob = prop_col_fun),
				show_annotation_name = !internal,
				annotation_name_side = "right",
				show_legend = FALSE)
		} else {
			ha2 = HeatmapAnnotation(Prob = membership_mat[!column_used_logical, ,drop = FALSE],
				Class = class_df$class[!column_used_logical],
				silhouette2 = anno_barplot(class_df$silhouette[!column_used_logical], ylim = silhouette_range,
					gp = gpar(fill = ifelse(class_df$silhouette[!column_used_logical] >= silhouette_cutoff, "grey", "grey"),
					      col = ifelse(class_df$silhouette[!column_used_logical] >= silhouette_cutoff, "black", NA)),
					bar_width = 1, baseline = 0, axis = TRUE, axis_param = list(side = "right"),
					height = unit(15, "mm")), 
				col = list(Class = cola_opt$color_set_2, Prob = prop_col_fun),
				show_annotation_name = c(TRUE, TRUE, FALSE) & !internal,
				annotation_name_side = "right",
				show_legend = FALSE)
		}
		ht_list = ht_list + Heatmap(use_mat2, name = paste0(heatmap_name, 2), col = col_fun,
			top_annotation = ha2,
			cluster_columns = TRUE, show_column_dend = FALSE,
			show_row_names = FALSE, show_row_dend = FALSE, show_heatmap_legend = FALSE,
			use_raster = use_raster, raster_resize = raster_resize,
			bottom_annotation = bottom_anno2, show_column_names = show_column_names,
			right_annotation = right_annotation)
	}

	if(has_ambiguous) {
		lgd = Legend(title = "Status (barplots)", labels = c("confident", "ambiguous"), legend_gp = gpar(fill = c("black", "grey")))
		heatmap_legend_list = list(lgd)
	} else {
		heatmap_legend_list = NULL
	}

	if(do_row_clustering) {
		ht_list = draw(ht_list, main_heatmap = heatmap_name, column_title = ifelse(internal, "", qq("@{k} subgroups, @{nrow(mat)} signatures (@{sprintf('%.1f',nrow(mat)/nrow(object)*100)}%) with fdr < @{fdr_cutoff}@{ifelse(group_diff > 0, paste0(', group_diff > ', group_diff), '')}")),
			show_heatmap_legend = !internal, show_annotation_legend = !internal,
			heatmap_legend_list = heatmap_legend_list,
			row_title = {if(length(unique(row_split)) <= 1) NULL else qq("k-means with @{length(unique(row_split))} groups")}
		)
		
		row_order = row_order(ht_list)
		if(!is.list(row_order)) row_order = list(row_order)
		if(scale_rows) {
			object@.env[[nm]]$row_order_scaled = do.call("c", row_order)
		} else {
			object@.env[[nm]]$row_order_unscaled = do.call("c", row_order)
		}
		
	} else {
		if(verbose) cat("  - use row order from cache.\n")
		draw(ht_list, main_heatmap = heatmap_name, column_title = ifelse(internal, "", qq("@{k} subgroups, @{nrow(mat)} signatures (@{sprintf('%.1f',nrow(mat)/nrow(object)*100)}%) with fdr < @{fdr_cutoff}@{ifelse(group_diff > 0, paste0(', group_diff > ', group_diff), '')}")),
			show_heatmap_legend = !internal, show_annotation_legend = !internal,
			cluster_rows = FALSE, row_order = row_order, heatmap_legend_list = heatmap_legend_list,
			row_title = {if(length(unique(row_split)) <= 1) NULL else qq("k-means with @{length(unique(row_split))} groups")}
		)
	}
	# the cutoff
	# https://www.stat.berkeley.edu/~s133/Cluster2a.html
	if(!internal) {
		for(i in seq_along(unique(class_df$class[column_used_logical]))) {
			decorate_annotation("silhouette",  slice = i, {
				grid.rect(gp = gpar(fill = "transparent"))
				grid.lines(c(0, 1), unit(c(silhouette_cutoff, silhouette_cutoff), "native"), gp = gpar(lty = 2, col = "#CCCCCC"))
				if(!has_ambiguous) grid.text("Silhouette\nscore", x = unit(1, "npc") + unit(10, "mm"), just = "top", rot = 90, gp = gpar(fontsize = 8))
			})
		}
		if(has_ambiguous) {
			decorate_annotation("silhouette2", {
				grid.rect(gp = gpar(fill = "transparent"))
				grid.lines(c(0, 1), unit(c(silhouette_cutoff, silhouette_cutoff), "native"), gp = gpar(lty = 2, col = "#CCCCCC"))
				if(has_ambiguous) grid.text("Silhouette\nscore", x = unit(1, "npc") + unit(10, "mm"), just = "top", rot = 90, gp = gpar(fontsize = 8))
			})
		}
	}

	return(invisible(returned_obj))
})

# https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
elbow_finder <- function(x_values, y_values) {
  # Max values to create line
  max_x_x <- max(x_values)
  max_x_y <- y_values[which.max(x_values)]
  max_y_y <- max(y_values)
  max_y_x <- x_values[which.max(y_values)]
  max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))

  # Creating straight line between the max values
  fit <- lm(max_df$y ~ max_df$x)

  # Distance from point to line
  distances <- c()
  for(i in 1:length(x_values)) {
    distances <- c(distances, abs(coef(fit)[2]*x_values[i] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2))
  }

  # Max distance point
  x_max_dist <- x_values[which.max(distances)]
  y_max_dist <- y_values[which.max(distances)]

  return(c(x_max_dist, y_max_dist))
}

# https://raghavan.usc.edu//papers/kneedle-simplex11.pdf
knee_finder = function(x, y) {
	n = length(x)
	a = (y[n] - y[1])/(x[n] - x[1])
	b = y[1] - a*x[1]
	d = a*x - y
	x[which.max(d)]
}

compare_to_subgroup = function(mat, class, which = "highest") {

	od = order(class)
	class = class[od]
	mat = mat[, od, drop = FALSE]
	class = as.numeric(factor(class))

	group = apply(mat, 1, function(x) {
		group_mean = tapply(x, class, mean)
		if(which == "highest") {
			which.max(group_mean)
		} else {
			which.min(group_mean)
		}
	})

	# class and subgroup_index are all numeric
	compare_to_one_subgroup = function(mat, class, subgroup_index) {
		
		oc = setdiff(class, subgroup_index)
		pmat = matrix(NA, nrow = nrow(mat), ncol = length(oc))
		for(i in seq_along(oc)) {
			l = class == subgroup_index | class == oc[i]
			m2 = mat[, l, drop = FALSE]
			fa = factor(class[l])
			pmat[, i] = genefilter::rowttests(m2, fa)[, "p.value"]
		}
		rowMaxs(pmat, na.rm = TRUE)
	}

	p = tapply(seq_len(nrow(mat)), group, function(ind) {
		compare_to_one_subgroup(mat[ind, , drop = FALSE], class, group[ind][1])
	})
	p2 = numeric(length(group))
	for(i in seq_along(p)) {
		p2[group == i] = p[[i]]
	}

	return(p2)
}

ttest = function(mat, class) {
	p1 = compare_to_subgroup(mat, class, "highest")
	p2 = compare_to_subgroup(mat, class, "lowest")
	p = pmin(p1, p2, na.rm = TRUE)
	fdr = p.adjust(p, method = "BH")
	fdr[is.na(fdr)] = Inf
	fdr[is.infinite(fdr)] = Inf
	fdr
}

one_vs_others = function(mat, class) {
	le = unique(class)
	dfl = list()
	for(x in le) {
		fa = as.vector(class)
		fa[class == le] = "a"
		fa[class != le] = "b"
		fa = factor(fa, levels = c("a", "b"))
		dfl[[x]] = genefilter::rowttests(mat, fa)[, "p.value"]
	}
	df = do.call("cbind", dfl)
	p = rowMins(df)
	p.adjust(p, "BH")
}

samr = function(mat, class, ...) {
	on.exit(if(sink.number()) sink(NULL))
	class = as.numeric(factor(class))
	n_class = length(unique(class))
	
	tempf = tempfile()
	sink(tempf)
	if(n_class == 2) {
		samfit = samr::SAM(mat, class, resp.type = "Two class unpaired", nperms = 1000, ...)
	} else {
		samfit = samr::SAM(mat, class, resp.type = "Multiclass", nperms = 1000, ...)
	}
	sink(NULL)
	file.remove(tempf)

	sig_index = NULL
	if(!is.null(samfit$siggenes.table$genes.up)) {
		id = samfit$siggenes.table$genes.up
		if(is.null(dim(id))) id = matrix(id, nrow = 1)
		sig_index = c(sig_index, id[, 2])
	}
	if(!is.null(samfit$siggenes.table$genes.lo)) {
		id = samfit$siggenes.table$genes.lo
		if(is.null(dim(id))) id = matrix(id, nrow = 1)
		sig_index = c(sig_index, id[, 2])
	}
	sig_index = as.numeric(sig_index)
	fdr = rep(1, nrow(mat))
	fdr[sig_index] = 0

	return(fdr)
}

pamr = function(mat, class, fdr.cutoff = 0.1, ...) {
	on.exit(if(sink.number()) sink(NULL))

	class = as.numeric(factor(class))
	
	tempf = tempfile()
	sink(tempf)
	mydata <- list(x=mat, y=class, geneid = rownames(mat))
	mydata.fit <- pamr::pamr.train(mydata)
	mydata.cv <- pamr::pamr.cv(mydata.fit, mydata)
	mydata.fdr <- pamr::pamr.fdr(mydata.fit, mydata)
	threshold = min(mydata.fdr$results[mydata.fdr$results[,"Median FDR"] < fdr.cutoff, "Threshold"])
	mydata.genelist <- pamr::pamr.listgenes(mydata.fit, mydata, threshold = threshold, fitcv=mydata.cv)
	sink(NULL)
	file.remove(tempf)

	fdr = rep(1, nrow(mat))
	fdr[rownames(mat) %in% mydata.genelist[,"id"]] = 0
	
	return(fdr)
}


Ftest = function(mat, class) {
	if(requireNamespace("genefilter")) {
		p = getFromNamespace("rowFtests", "genefilter")(mat, factor(class))[, "p.value"]
		fdr = p.adjust(p, "BH")
		fdr[is.na(fdr)] = Inf
		return(fdr)
	} else {
		stop_wrap("Cannot find 'genefilter' package.")
	}
}

# test_row_diff_fun = function(fun, fdr_cutoff = 0.1) {
# 	set.seed(100)
# 	x = matrix(rnorm(1000 * 20), ncol = 20)
# 	rownames(x) = rep(paste0("gene1", 1:1000))
# 	dd = sample(1:1000, size = 100)
# 	u = matrix(2 * rnorm(100), ncol = 10, nrow = 100)
# 	x[dd, 11:20] = x[dd, 11:20] + u
# 	row_diff = rep("no", 1000)
# 	row_diff[dd] = "yes"
# 	y = c(rep(1, 10), rep(2, 10))
# 	fdr = fun(x, y)

# 	ht = Heatmap(x, top_annotation = HeatmapAnnotation(foo = as.character(y), col = list(foo = c("1" = "blue", "2" = "red"))), show_row_names = FALSE) +
# 	Heatmap(row_diff, name = "diff", col = c("yes" = "red", "no" = "white"), width = unit(5, "mm")) +
# 	Heatmap(fdr, name = "fdr", width = unit(5, "mm"), show_row_names = FALSE)
# 	draw(ht, split = fdr < fdr_cutoff)
# }




# title
# Density for the signatures
#
# == param
# -object A `ConsensusPartition-class` object. 
# -k number of partitions
# -... pass to `get_signatures,ConsensusPartition-method`
#
# == details
# The function makes density distributio nf of signatures in all columns.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# setMethod(f = "signature_density",
# 	signature = "ConsensusPartition",
# 	definition = function(object, k, ...) {

# 	cl = get_class(object, k = k)$class
# 	data = object@.env$data[, object@column_index, drop = FALSE]

# 	all_den_list = lapply(seq_len(ncol(data)), function(i) {
# 		x = data[, i]
# 		density(x)
# 	})
# 	x_range = range(unlist(lapply(all_den_list, function(x) x$x)))
# 	y_range = range(unlist(lapply(all_den_list, function(x) x$y)))

# 	x = get_signatures(object, k = k, plot = FALSE, verbose = FALSE, ...)
# 	gp_tb = table(x$df$group)
# 	n_gp = sum(gp_tb > 5)
# 	gp_tb = gp_tb[gp_tb > 5]

# 	op = par(no.readonly = TRUE)
# 	par(mfrow = c(n_gp + 1, 1), mar = c(2, 4, 1, 3))
# 	plot(NULL, type = "n", xlim = x_range, ylim = y_range, ylab = "density", xlab = NULL)
# 	for(i in 1:ncol(data)) {
# 		lines(all_den_list[[i]], col = cola_opt$color_set_2[cl[i]], lwd = 1)
# 	}
# 	mtext("all rows", side = 4, line = 1)

# 	gp = x$df$group
# 	for(j in as.numeric(names(gp_tb))) {
# 		gp2 = gp[gp == as.character(j)]
# 		all_den_list = lapply(seq_len(ncol(data)), function(i) density(data[names(gp2), i]))
# 		# x_range = range(unlist(lapply(all_den_list, function(x) x$x)))
# 		y_range = range(unlist(lapply(all_den_list, function(x) x$y)))

# 		plot(NULL, type = "n", xlim = x_range, ylim = y_range, ylab = "density", xlab = NULL)
# 		for(i in 1:ncol(data)) {
# 			lines(all_den_list[[i]], col = cola_opt$color_set_2[cl[i]], lwd = ifelse(cl[i] == j, 2, 0.5))
# 		}
# 		mtext(qq("subgroup @{j}/@{k}"), side = 4, line = 1)
# 	}
# 	par(op)
# })


# == title
# Compare Signatures from Different k
#
# == param
# -object A `ConsensusPartition-class` object. 
# -k Number of partitions. Value should be a vector.
# -... Other arguments passed to `get_signatures,ConsensusPartition-method`.
#
# == details
# It plots an Euler diagram showing the overlap of signatures from different k.
#
setMethod(f = "compare_signatures",
	signature = "ConsensusPartition",
	definition = function(object, k = object@k, ...) {

	sig_list = sapply(k, function(x) {
		tb = get_signatures(object, k = x, ..., plot = FALSE)
		if(is.null(tb)) {
			return(integer(0))
		} else {
			return(tb$which_row)
		}
	})

	names(sig_list) = paste(k, "-group", sep = "")

	plot(eulerr::euler(sig_list), legend = TRUE, quantities = TRUE, main = "Signatures from different k")

})


# == title
# Find a best k for the k-means clustering
#
# == param
# -mat A matrix where k-means clustering is executed by rows.
# -max_km Maximal k to try.
#
# == details
# The best k is determined by looking for the knee/elbow of the WSS curve (within-cluster sum of square).
#
# Note this function is only for a rough and quick determination of the best k.
#
find_best_km = function(mat, max_km = 15) {
	wss = (nrow(mat)-1)*sum(apply(mat,2,var))
	max_km = min(c(nrow(mat) - 1, max_km))
	for (i in 2:max_km) wss[i] = sum(kmeans(mat, centers = i, iter.max = 50)$withinss)
	row_km = min(elbow_finder(1:max_km, wss)[1], knee_finder(1:max_km, wss)[1])
	return(row_km)
}
