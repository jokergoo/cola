
# == title
# Get signature rows
#
# == param
# -object A `ConsensusPartition-class` object. The object can be returned
#        from `get_single_run`.
# -k number of partitions
# -silhouette_cutoff cutoff for silhouette values. Columns with values 
#        less than it are not used for finding signature rows. For selecting a 
#        proper silhouette value, please refer to https://www.stat.berkeley.edu/~s133/Cluster2a.html#tth_tAb1.
# -fdr_cutoff cutoff for FDR of the difference between subgroups.
# -scale_rows whether apply row scaling when making the heatmap.
# -row_diff_by methods to get rows which are significantly different between subgroups.
# -anno a data frame with known annotation of columns.
# -anno_col a list of colors for the annotations in ``anno``.
# -show_legend whether draw the legends on the heatmap.
# -show_column_names whether show column names on the heatmap.
# -use_raster internally used
# -... other arguments
# 
# == details 
# Basically the function apply test for the difference of subgroups in every
# row. Also, to call it a signature for a given subgroup, the values in the
# corresponding subgroup should have the highest mean value compared to all
# other subgroups. The minimal p-value compared to all other subgroups is taken
# as the p-value of the row and used for FDR calculation.
#
# == return 
# A list of three elements:
# 
# -``mat`` the matrix for the signatures
# -``fdr`` FDR for rows
# -``group`` subgroups that the rows are significantly high
# 
setMethod(f = "get_signatures",
	signature = "ConsensusPartition",
	definition = function(object, k,
	silhouette_cutoff = 0.5, 
	fdr_cutoff = ifelse(row_diff_by == "samr", 0.1, 0.05), 
	scale_rows = TRUE,
	row_diff_by = c("compare_to_highest_subgroup", "Ftest", "samr"),
	anno = object@known_anno, 
	anno_col = if(missing(anno)) object@known_col else NULL,
	show_legend = TRUE,
	show_column_names = TRUE, use_raster = TRUE,
	...) {

	class_df = get_class(object, k)
	class_ids = class_df$class
	class_col = class_df$class_col
	class_col = structure(unique(class_col), names = unique(class_ids))
	class_col = class_col[order(names(class_col))]

	l = class_df$silhouette >= silhouette_cutoff
	data2 = data[, l]
	class = class_df$class[l]
	column_used_index = which(l)
	tb = table(class)
	l = as.character(class) %in% names(which(tb <= 1))
	data2 = data2[, !l]
	class = class[!l]
	column_used_index = column_used_index[!l]
	column_used_logical = rep(FALSE, ncol(data))
	column_used_logical[column_used_index] = TRUE
	has_ambiguous = sum(!column_used_logical)
	n_sample_used = length(class)

	qqcat("@{n_sample_used}/@{nrow(class_df)} samples (in @{length(unique(class))} classes) remain after filtering by silhouette (>= @{silhouette_cutoff}).\n")

	if(sum(tb > 1) <= 1) {
		stop("not enough samples.")
	}
	if(length(unique(class)) <= 1) {
		stop("not enough classes.")
	}

	if(inherits(row_diff_by, "function")) {
		qqcat("calculate row difference between subgroups by user-defined function\n")
		fdr = row_diff_by(data2, class)
	} else {
		row_diff_by = match.arg(row_diff_by)

		nm = paste0("signature_fdr_", object@top_method, "_", object@partition_method, "_", k, "groups")

		find_signature = TRUE
		if(!is.null(object@.env[[nm]])) {
			if(row_diff_by == "samr") {
				if(object@.env[[nm]]$row_diff_by == "samr" && 
				   object@.env[[nm]]$n_sample_used == n_sample_used && 
				   abs(object@.env[[nm]]$fdr_cutoff - fdr_cutoff) < 1e-6) {
					fdr = object@.env[[nm]]$fdr
					find_signature = FALSE
				}
			} else {
				if(object@.env[[nm]]$row_diff_by == row_diff_by &&
				   object@.env[[nm]]$n_sample_used == n_sample_used) {
					fdr = object@.env[[nm]]$fdr
					find_signature = FALSE
				}
			}
		}

		if(find_signature) {
			qqcat("calculate row difference between subgroups by @{row_diff_by}\n")
			if(row_diff_by == "compare_to_highest_subgroup") {
				fdr = compare_to_highest_subgroup(data2, class)
			} else if(row_diff_by == "samr") {
				fdr = samr(data2, class, fdr.output = fdr_cutoff)
			} else if(row_diff_by == "Ftest") {
				fdr = Ftest(data2, class)
			}
		}

		object@.env[[nm]]$row_diff_by = row_diff_by
		object@.env[[nm]]$fdr_cutoff = fdr_cutoff
		object@.env[[nm]]$fdr = fdr
		object@.env[[nm]]$n_sample_used = n_sample_used
	}

	group = character(nrow(data2))
	for(i in seq_len(nrow(data2))) {
		x = data2[i, ]
		group_mean = tapply(x, class, mean)
		max_group = names(which.max(group_mean))
		group[i] = max_group
	}

	mat = data[fdr < fdr_cutoff, ]
	fdr2 = fdr[fdr < fdr_cutoff]
	group2 = group[fdr < fdr_cutoff]
	names(group2) = rownames(mat)
	mat_return = list(mat = mat, fdr = fdr2, group = group2)
	qqcat("@{nrow(mat)} signatures under fdr < @{fdr_cutoff}\n")

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

	if(is.null(anno)) {
		bottom_anno1 = NULL
	} else {
		if(is.atomic(anno)) {
			anno_nm = deparse(substitute(anno))
			anno = data.frame(anno)
			colnames(anno) = anno_nm
			if(!is.null(anno_col)) {
				anno_col = list(anno_col)
				names(anno_col) = anno_nm
			}
		}

		if(is.null(anno_col)) {
			bottom_anno1 = HeatmapAnnotation(df = anno[column_used_logical, , drop = FALSE],
				show_annotation_name = !has_ambiguous, annotation_name_side = "right")
		} else {
			bottom_anno1 = HeatmapAnnotation(df = anno[column_used_logical, , drop = FALSE], col = anno_col,
				show_annotation_name = !has_ambiguous, annotation_name_side = "right")
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
		class2 = class_df$class[!column_used_logical]
		mat_col_od2 = column_order_by_group(class2, use_mat2)

		if(is.null(anno)) {
			bottom_anno2 = NULL
		} else {

			if(is.null(anno_col)) {
				bottom_anno2 = HeatmapAnnotation(df = anno[!column_used_logical, , drop = FALSE],
					show_annotation_name = TRUE, annotation_name_side = "right")
			} else {
				bottom_anno2 = HeatmapAnnotation(df = anno[!column_used_logical, , drop = FALSE], col = anno_col,
					show_annotation_name = TRUE, annotation_name_side = "right")
				for(i in seq_along(bottom_anno2$anno_list)) {
					bottom_anno2@anno_list[[i]]@color_mapping = bottom_anno1@anno_list[[i]]@color_mapping
				}
			}
		}
	}
	silhouette_range = range(class_df$silhouette)
	silhouette_range[2] = 1

	group2 = factor(group2, levels = sort(unique(group2)))
	ht_list = Heatmap(group2, name = "group", show_row_names = FALSE, width = unit(5, "mm"), col = class_col)

	membership_mat = as.data.frame(get_membership(object, k))
	col_list = rep(list(colorRamp2(c(0, 1), c("white", "red"))), k)
	names(col_list) = colnames(membership_mat)

	ht_list = ht_list + Heatmap(use_mat1, name = heatmap_name, col = col_fun,
		top_annotation = HeatmapAnnotation(df = membership_mat[column_used_logical, ],
			class = class_df$class[column_used_logical],
			silhouette = anno_barplot(class_df$silhouette[column_used_logical], ylim = silhouette_range,
				gp = gpar(fill = ifelse(class_df$silhouette[column_used_logical] >= silhouette_cutoff, "grey", "#EEEEEE"),
					      col = ifelse(class_df$silhouette[column_used_logical] >= silhouette_cutoff, "black", NA)),
				baseline = 0, axis = !has_ambiguous, axis_side = "right"),
			col = c(list(class = class_col), col_list),
			annotation_height = unit(c(rep(4, k+1), 15), "mm"),
			show_annotation_name = if(has_ambiguous) FALSE else c(rep(TRUE, k+1), FALSE),
			annotation_name_side = "right",
			show_legend = c(TRUE, rep(FALSE, k - 1), TRUE)),
		cluster_columns = FALSE, column_order = mat_col_od1,
		show_row_names = FALSE, show_row_dend = FALSE, column_title = "confident samples",
		use_raster = use_raster, raster_quality = 2,
		bottom_annotation = bottom_anno1, show_column_names = show_column_names, split = group2, combined_name_fun = NULL)
 	if(scale_rows) {
		ht_list = ht_list + Heatmap(base_mean, show_row_names = FALSE, name = "base_mean", width = unit(5, "mm"))
	}
	if(has_ambiguous) {
		ht_list = ht_list + Heatmap(use_mat2, name = paste0(heatmap_name, 2), col = col_fun,
			top_annotation = HeatmapAnnotation(df = membership_mat[!column_used_logical, ],
				class = class_df$class[!column_used_logical],
				silhouette2 = anno_barplot(class_df$silhouette[!column_used_logical], ylim = silhouette_range,
					gp = gpar(fill = ifelse(class_df$silhouette[!column_used_logical] >= silhouette_cutoff, "grey", "grey"),
					      col = ifelse(class_df$silhouette[!column_used_logical] >= silhouette_cutoff, "black", NA)),
					baseline = 0, axis = TRUE, axis_side = "right"), 
				col = c(list(class = class_col), col_list),
				annotation_height = unit(c(rep(4, k+1), 15), "mm"),
				show_annotation_name = c(rep(TRUE, k+1), FALSE),
				annotation_name_side = "right",
				show_legend = FALSE),
			cluster_columns = FALSE, column_order = mat_col_od2,
			show_row_names = FALSE, show_row_dend = FALSE, show_heatmap_legend = FALSE,
			use_raster = use_raster, raster_quality = 2,
			column_title = "ambiguous samples", bottom_annotation = bottom_anno2, show_column_names = show_column_names)
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
})

compare_to_highest_subgroup = function(mat, class) {

	od = order(class)
	class = class[od]
	mat = mat[, od, drop = FALSE]
	class = as.numeric(factor(class))

	group = apply(mat, 1, function(x) {
		group_mean = tapply(x, class, mean)
		which.max(group_mean)
	})

	# class and subgroup_index are all numeric
	compare_to_one_subgroup = function(mat, class, subgroup_index) {
		
		oc = setdiff(class, subgroup_index)
		pmat = matrix(NA, nrow = nrow(mat), ncol = length(oc))
		for(i in seq_along(oc)) {
			l = class == subgroup_index | class == oc[i]
			m2 = mat[, l]
			fa = factor(class[l])
			pmat[, i] = rowttests(m2, fa)[, "p.value"]
		}
		rowMins(pmat, na.rm = TRUE)
	}

	p = tapply(seq_len(nrow(mat)), group, function(ind) {
		compare_to_one_subgroup(mat[ind, , drop = FALSE], class, group[ind][1])
	})
	p2 = numeric(length(group))
	for(i in seq_along(p)) {
		p2[group == i] = p[[i]]
	}

	fdr = p.adjust(p2, "BH")
	fdr[is.na(fdr)] = Inf
	return(fdr)
}

samr = function(mat, class, ...) {
	class = as.numeric(factor(class))
	n_class = length(unique(class))
	
	tempf = tempfile()
	sink(tempf)
	if(n_class == 2) {
		samfit = SAM(mat, class, resp.type = "Two class unpaired", nperms = 1000, ...)
	} else {
		samfit = SAM(mat, class, resp.type = "Multiclass", nperms = 1000, ...)
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

Ftest = function(mat, class) {
	p = rowFtests(mat, factor(class))[, "p.value"]
	fdr = p.adjust(p, "BH")
	fdr[is.na(fdr)] = Inf
	return(fdr)
}

test_row_diff_fun = function(fun, fdr_cutoff = 0.1) {
	set.seed(100)
	x = matrix(rnorm(1000 * 20), ncol = 20)
	dd = sample(1:1000, size = 100)
	u = matrix(2 * rnorm(100), ncol = 10, nrow = 100)
	x[dd, 11:20] = x[dd, 11:20] + u
	row_diff = rep("no", 1000)
	row_diff[dd] = "yes"
	y = c(rep(1, 10), rep(2, 10))
	fdr = fun(x, y)

	ht = Heatmap(x, top_annotation = HeatmapAnnotation(foo = as.character(y), col = list(foo = c("1" = "blue", "2" = "red")))) +
	Heatmap(row_diff, name = "diff", col = c("yes" = "red", "no" = "white"), width = unit(5, "mm")) +
	Heatmap(fdr, name = "fdr", width = unit(5, "mm"), show_row_names = FALSE)
	draw(ht, split = fdr < fdr_cutoff)
}



# lines = strsplit(readLines("h.all.v6.0.symbols.gmt"), "\t")
# genesets = lapply(lines, function(x) x[-(1:2)])
# names(genesets) = sapply(lines, function(x) x[1])
BHI = function(x, map = NULL, genesets, min_count = 50, max_count = 5000) {
	match_mat = matrix(0, nrow = nrow(x$mat), ncol = length(genesets))

	gl = sapply(genesets, length)
	l = gl >= min_count & gl <= max_count
	genesets = genesets[l]

	g = rownames(x$mat)
	if(!is.null(map)) {
		g2 = map[g]
		l = is.na(g2)
		g2[l] = g[l]
		g = g2
	}
	rownames(match_mat) = g
	colnames(match_mat) = names(genesets)

	for(i in seq_along(genesets)) {
		match_mat[genesets[[i]], i] = 1
	}

	mean(tapply(seq_len(nrow(x$mat)), x$group, function(ind) {
		g2 = g[ind]
		submat = matrix(0, nrow = length(g2), ncol = length(g2))
		rownames(submat) = g2
		colnames(submat) = g2

		for(i in 1:(length(g2) - 1)) {
			for(j in i:length(g2)) {
				submat[i, j] = colSums(match_mat[g2[c(i, j)], ]) == 2
				submat[j, i] = submat[i, j]
			}
		}

		mean(submat[lower.tri(submat)])
	}))

}
