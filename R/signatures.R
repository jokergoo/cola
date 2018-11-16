
# == title
# Get signature rows
#
# == param
# -object A `ConsensusPartition-class` object.
# -k number of partitions.
# -silhouette_cutoff cutoff for silhouette values. Columns with values 
#        less than it are not used for finding signature rows. For selecting a 
#        proper silhouette value, please refer to https://www.stat.berkeley.edu/~s133/Cluster2a.html#tth_tAb1.
# -fdr_cutoff cutoff for FDR of the difference between subgroups.
# -scale_rows whether apply row scaling when making the heatmap.
# -diff_method methods to get rows which are significantly different between subgroups, see 'Details' section.
# -anno a data frame with known annotation of samples.
# -anno_col a list of colors for the annotations in ``anno``.
# -internal used internally.
# -show_column_names whether show column names on the heatmap.
# -use_raster internally used.
# -plot whether to make the plot.
# -verbose whether to print messages.
# -top_k_genes top k genes to show on the heatmap if the number of signatures exceed it.
# -... other arguments.
# 
# == details 
# Basically the function applies test for the difference of subgroups for every
# row. There are following methods which test significance of the difference:
#
# -ttest First it looks for the subgroup with highest mean value, compare to each of the 
#        other subgroups with t-test and take the maximum p-value. Second it looks
#        for the subgroup with lowest mean value, compare to each of the other subgroups
#        again with t-test and take the maximum p-values. Later for these two list of p-values
#        take the minimal p-value as the final p-value. 
# -samr/pamr use SAM/PAM method to find significantly different rows between subgroups.
# -Ftest use F-test to find significantly different rows between subgroups.
#
# Also, to call it a signature for a given subgroup, the values in the
# corresponding subgroup should have the highest mean value compared to all
# other subgroups. The minimal p-value compared to all other subgroups is taken
# as the p-value of the row and used for FDR calculation.
#
# == return 
# A list of three elements:
# 
# -``df`` a data frame.
# -``sample_used`` sample index used.
# 
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "get_signatures",
	signature = "ConsensusPartition",
	definition = function(object, k,
	silhouette_cutoff = 0.5, 
	fdr_cutoff = ifelse(diff_method == "samr", 0.05, 0.1), 
	scale_rows = object@scale_rows,
	diff_method = c("ttest", "Ftest", "samr", "pamr"),
	anno = get_anno(object), 
	anno_col = get_anno_col(object),
	internal = FALSE,
	show_column_names = FALSE, use_raster = TRUE,
	plot = TRUE, verbose = TRUE,
	top_k_genes = 5000,
	...) {

	class_df = get_classes(object, k)
	class_ids = class_df$class

	data = object@.env$data[, object@column_index, drop = FALSE]

	l = class_df$silhouette >= silhouette_cutoff
	data2 = data[, l]
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
		return(invisible(NULL))
	}
	if(length(unique(class)) <= 1) {
		if(plot) {
			grid.newpage()
			fontsize = convertUnit(unit(0.1, "npc"), "char", valueOnly = TRUE)*get.gpar("fontsize")$fontsize
			grid.text("not enough classes", gp = gpar(fontsize = fontsize))
		}
		return(invisible(NULL))
	}

	do_row_clustering = TRUE
	if(inherits(diff_method, "function")) {
		if(verbose) qqcat("* calculate row difference between subgroups by user-defined function\n")
		diff_method_fun = diff_method
		diff_method = digest(diff_method)
	} else {
		diff_method = match.arg(diff_method)
	}

	hash = digest(list(used_samples = which(l), 
		               class = class,
		               n_group = k, 
		               diff_method = diff_method,
		               column_index = object@column_index),
				algo = "md5")
	nm = paste0("signature_fdr_", hash)

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

	if(verbose) qqcat("* calculating row difference between subgroups by @{diff_method}\n")
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
		if(verbose) qqcat("  - row difference is extracted from cache\n")
	}

	object@.env[[nm]]$diff_method = diff_method
	object@.env[[nm]]$fdr_cutoff = fdr_cutoff
	object@.env[[nm]]$fdr = fdr
	object@.env[[nm]]$n_sample_used = n_sample_used

	if(scale_rows && !is.null(object@.env[[nm]]$row_order_scaled)) {
		row_order = object@.env[[nm]]$row_order_scaled
		if(verbose) qqcat("  - row order for the scaled matrix is extracted from cache\n")
		do_row_clustering = FALSE
	} else if(!scale_rows && !is.null(object@.env[[nm]]$row_order_unscaled)) {
		row_order = object@.env[[nm]]$row_order_unscaled
		if(verbose) qqcat("  - row order for the unscaled matrix is extracted from cache\n")
		do_row_clustering = FALSE
	}

	fdr[is.na(fdr)] = 1

	group = character(nrow(data2))
	for(i in seq_len(nrow(data2))) {
		x = data2[i, ]
		group_mean = tapply(x, class, mean)
		od = order(group_mean)

		if(group_mean[od[2]] - group_mean[od[1]] > group_mean[od[length(od)]] - group_mean[od[length(od)-1]]) {
			group[i] = names(group_mean)[od[1]]
		} else {
			group[i] = names(group_mean)[od[length(od)]]
		}
	}

	l_fdr = fdr < fdr_cutoff
	mat = data[l_fdr, , drop = FALSE]
	fdr2 = fdr[l_fdr]
	group2 = group[l_fdr]
	names(group2) = rownames(mat)

	returned_df = data.frame(which_row = which(l_fdr), fdr = fdr2, group = group2)
	returned_df = returned_df[order(returned_df[, "fdr"]), ]
	returned_obj = returned_df
	attr(returned_obj, "sample_used") = column_used_logical

	if(verbose) qqcat("* @{nrow(mat)} signatures under fdr < @{fdr_cutoff}\n")

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

	more_than_5k = FALSE
	if(nrow(mat) > top_k_genes) {
		more_than_5k = TRUE
		mat1 = mat[order(fdr2)[1:top_k_genes], column_used_logical, drop = FALSE]
		mat2 = mat[order(fdr2)[1:top_k_genes], !column_used_logical, drop = FALSE]
		group2 = group2[order(fdr2)[1:top_k_genes]]
		if(verbose) cat(paste0("* Only take top ", top_k_genes, " signatures with highest fdr for the plot\n"))
	} else {
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
		mat_range = quantile(abs(scaled_mat1), 0.95, na.rm = TRUE)
		col_fun = colorRamp2(c(-mat_range, 0, mat_range), c("green", "white", "red"))
		heatmap_name = "z-score"
	} else {
		use_mat1 = mat1
		use_mat2 = mat2
		mat_range = quantile(mat1, c(0.05, 0.95))
		col_fun = colorRamp2(c(mat_range[1], mean(mat_range), mat_range[2]), c("blue", "white", "red"))
		heatmap_name = "expr"
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
				anno_col[names(object@anno_col)] = object@anno_col
			}

			bottom_anno2 = HeatmapAnnotation(df = anno[!column_used_logical, , drop = FALSE], col = anno_col,
				show_annotation_name = TRUE, annotation_name_side = "right")	
		}
	}
	silhouette_range = range(class_df$silhouette)
	silhouette_range[2] = 1

	if(verbose) qqcat("* making heatmaps for signatures\n")

	group2 = factor(group2, levels = sort(unique(group2)))
	ht_list = Heatmap(group2, name = "Group", show_row_names = FALSE, width = unit(5, "mm"), col = brewer_pal_set2_col)

	membership_mat = get_membership(object, k)
	prop_col_fun = colorRamp2(c(0, 1), c("white", "red"))

	if(internal) {
		ha1 = HeatmapAnnotation(Prob = membership_mat[column_used_logical, ],
				Class = class_df$class[column_used_logical],
				col = list(Class = brewer_pal_set2_col, Prob = prop_col_fun),
				show_annotation_name = !has_ambiguous,
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
			col = list(Class = brewer_pal_set2_col, Prob = prop_col_fun),
			show_annotation_name = !has_ambiguous,
			annotation_name_side = "right",
			show_legend = TRUE)
	}
	ht_list = ht_list + Heatmap(use_mat1, name = heatmap_name, col = col_fun,
		top_annotation = ha1,
		cluster_columns = TRUE, column_split = class_df$class[column_used_logical], show_column_dend = FALSE,
		show_row_names = FALSE, show_row_dend = FALSE, column_title = "confident samples",
		use_raster = use_raster,
		bottom_annotation = bottom_anno1, show_column_names = show_column_names, row_split = group2, row_title = NULL)
 	
	all_value_positive = !any(data < 0)
 	if(scale_rows && all_value_positive) {
		ht_list = ht_list + Heatmap(base_mean, show_row_names = FALSE, name = "base_mean", width = unit(5, "mm")) +
			Heatmap(rel_diff, col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), 
				show_row_names = FALSE, name = "rel_diff", width = unit(5, "mm"))
	}

	if(has_ambiguous) {
		if(internal) {
			ha2 = HeatmapAnnotation(Prob = membership_mat[!column_used_logical, ,drop = FALSE],
				Class = class_df$class[!column_used_logical],
				col = list(Class = brewer_pal_set2_col, Prob = prop_col_fun),
				show_annotation_name = TRUE,
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
				col = list(Class = brewer_pal_set2_col, Prob = prop_col_fun),
				show_annotation_name = c(TRUE, TRUE, FALSE),
				annotation_name_side = "right",
				show_legend = FALSE)
		}
		ht_list = ht_list + Heatmap(use_mat2, name = paste0(heatmap_name, 2), col = col_fun,
			top_annotation = ha2,
			cluster_columns = TRUE, show_column_dend = FALSE,
			show_row_names = FALSE, show_row_dend = FALSE, show_heatmap_legend = FALSE,
			use_raster = use_raster,
			bottom_annotation = bottom_anno2, show_column_names = show_column_names)
	}

	if(do_row_clustering) {
		ht_list = draw(ht_list, main_heatmap = heatmap_name, column_title = qq("@{k} subgroups, @{nrow(mat)} signatures with fdr < @{fdr_cutoff}@{ifelse(more_than_5k, paste0(', top ', top_k_genes,' signatures'), '')}"),
			show_heatmap_legend = !internal, show_annotation_legend = !internal)
		
		row_order = row_order(ht_list)
		if(!is.list(row_order)) row_order = list(row_order)
		if(scale_rows) {
			object@.env[[nm]]$row_order_scaled = do.call("c", row_order)
		} else {
			object@.env[[nm]]$row_order_unscaled = do.call("c", row_order)
		}
		
	} else {
		if(verbose) cat("  - using row order from cache\n")
		draw(ht_list, main_heatmap = heatmap_name, column_title = qq("@{k} subgroups, @{nrow(mat)} signatures with fdr < @{fdr_cutoff}@{ifelse(more_than_5k, paste0(', top ', top_k_genes,' signatures'), '')}"),
			show_heatmap_legend = !internal, show_annotation_legend = !internal,
			cluster_rows = FALSE, row_order = row_order)
	}
	# the cutoff
	# https://www.stat.berkeley.edu/~s133/Cluster2a.html
	if(!internal) {
		decorate_annotation("silhouette", {
			grid.rect(gp = gpar(fill = "transparent"))
			grid.lines(c(0, 1), unit(c(silhouette_cutoff, silhouette_cutoff), "native"), gp = gpar(lty = 2, col = "#CCCCCC"))
			if(!has_ambiguous) grid.text("Silhouette\nscore", x = unit(1, "npc") + unit(10, "mm"), just = "top", rot = 90, gp = gpar(fontsize = 8))
		})
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

get_signature_index = function(sig, group) {
	which(sig$df$group == group)
}

compare_to_highest_subgroup = function(mat, class) {

	od = order(class)
	class = class[od]
	mat = mat[, od, drop = FALSE]
	class = as.numeric(factor(class))

	group = apply(mat, 1, function(x) {
		group_mean = tapply(x, class, mean)
		which.max(group_mean)
	})

	if(requireNamespace("genefilter")) {
		rowttests = getFromNamespace("rowttests", "genefilter")
	} else {
		stop("Cannot find 'genefilter' package.")
	}

	# class and subgroup_index are all numeric
	compare_to_one_subgroup = function(mat, class, subgroup_index) {
		
		oc = setdiff(class, subgroup_index)
		pmat = matrix(NA, nrow = nrow(mat), ncol = length(oc))
		for(i in seq_along(oc)) {
			l = class == subgroup_index | class == oc[i]
			m2 = mat[, l, drop = FALSE]
			fa = factor(class[l])
			pmat[, i] = rowttests(m2, fa)[, "p.value"]
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

compare_to_lowest_subgroup = function(mat, class) {

	od = order(class)
	class = class[od]
	mat = mat[, od, drop = FALSE]
	class = as.numeric(factor(class))

	group = apply(mat, 1, function(x) {
		group_mean = tapply(x, class, mean)
		which.min(group_mean)
	})

	if(requireNamespace("genefilter")) {
		rowttests = getFromNamespace("rowttests", "genefilter")
	} else {
		stop("Cannot find 'genefilter' package.")
	}

	# class and subgroup_index are all numeric
	compare_to_one_subgroup = function(mat, class, subgroup_index) {
		
		oc = setdiff(class, subgroup_index)
		pmat = matrix(NA, nrow = nrow(mat), ncol = length(oc))
		for(i in seq_along(oc)) {
			l = class == subgroup_index | class == oc[i]
			m2 = mat[, l, drop = FALSE]
			fa = factor(class[l])
			pmat[, i] = rowttests(m2, fa)[, "p.value"]
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
	p1 = compare_to_highest_subgroup(mat, class)
	p2 = compare_to_lowest_subgroup(mat, class)
	p = pmin(p1, p2, na.rm = TRUE)
	fdr = p.adjust(p, method = "BH")
	fdr[is.na(fdr)] = Inf
	fdr
}

samr = function(mat, class, ...) {
	on.exit(if(sink.number()) sink(NULL))
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

pamr = function(mat, class, fdr.cutoff = 0.1, ...) {
	on.exit(if(sink.number()) sink(NULL))

	class = as.numeric(factor(class))
	
	tempf = tempfile()
	sink(tempf)
	mydata <- list(x=mat, y=class, geneid = rownames(mat))
	mydata.fit <- pamr.train(mydata)
	mydata.cv <- pamr.cv(mydata.fit, mydata)
	mydata.fdr <- pamr.fdr(mydata.fit, mydata)
	threshold = min(mydata.fdr$results[mydata.fdr$results[,"Median FDR"] < fdr.cutoff, "Threshold"])
	mydata.genelist <- pamr.listgenes(mydata.fit, mydata, threshold = threshold, fitcv=mydata.cv)
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
		stop("Cannot find 'genefilter' package.")
	}
}

test_row_diff_fun = function(fun, fdr_cutoff = 0.1) {
	set.seed(100)
	x = matrix(rnorm(1000 * 20), ncol = 20)
	rownames(x) = rep(paste0("gene1", 1:1000))
	dd = sample(1:1000, size = 100)
	u = matrix(2 * rnorm(100), ncol = 10, nrow = 100)
	x[dd, 11:20] = x[dd, 11:20] + u
	row_diff = rep("no", 1000)
	row_diff[dd] = "yes"
	y = c(rep(1, 10), rep(2, 10))
	fdr = fun(x, y)

	ht = Heatmap(x, top_annotation = HeatmapAnnotation(foo = as.character(y), col = list(foo = c("1" = "blue", "2" = "red"))), show_row_names = FALSE) +
	Heatmap(row_diff, name = "diff", col = c("yes" = "red", "no" = "white"), width = unit(5, "mm")) +
	Heatmap(fdr, name = "fdr", width = unit(5, "mm"), show_row_names = FALSE)
	draw(ht, split = fdr < fdr_cutoff)
}




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
# 		lines(all_den_list[[i]], col = brewer_pal_set2_col[cl[i]], lwd = 1)
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
# 			lines(all_den_list[[i]], col = brewer_pal_set2_col[cl[i]], lwd = ifelse(cl[i] == j, 2, 0.5))
# 		}
# 		mtext(qq("subgroup @{j}/@{k}"), side = 4, line = 1)
# 	}
# 	par(op)
# })
