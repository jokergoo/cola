
# == title
# The DownSamplingConsensusPartition class
#
# == alias
# DownSamplingConsensusPartition
#
# == details
# The ``DownSamplingConsensusPartition`` performs consensus partitioning only with a small subset
# of columns and the class of other columns are predicted by `predict_classes,ConsensusPartition-method`.
#
# The ``DownSamplingConsensusPartition-class`` is a child class of `ConsensusPartition-class`. It inherits
# all methods of `ConsensusPartition-class`.
#
# == seealso
# The constructor function `consensus_partition_by_down_sampling`.
#
DownSamplingConsensusPartition = setClass("DownSamplingConsensusPartition",
	slots = list(predict = "ANY",
		         full_column_index = "ANY",
		         full_anno = "ANY"),
	contains = "ConsensusPartition"
)

# == title
# Consensus partitioning only with a subset of columns
#
# == param
# -data A numeric matrix where subgroups are found by columns.
# -top_value_method A single top-value method. Available methods are in `all_top_value_methods`.
#                   Use `register_top_value_methods` to add a new top-value method.
# -top_n Number of rows with top values. The value can be a vector with length > 1. When n > 5000, 
#        the function only randomly sample 5000 rows from top n rows. If ``top_n`` is a vector, paritition
#        will be applied to every values in ``top_n`` and consensus partition is summarized from all partitions.
# -partition_method A single partitioning method. Available methods are in `all_partition_methods`.
#                   Use `register_partition_methods` to add a new partition method.
# -max_k Maximal number of subgroups to try. The function will try for ``2:max_k`` subgroups
# -subset Number of columns to randomly sample, or a vector of selected indices.
# -verbose Whether to print messages.
# -prefix Internally used.
# -anno Annotation data frame.
# -anno_col Annotation colors.
# -dist_method Method for predict the class for other columns.
# -.env An environment, internally used.
# -.predict Internally used.
# -mc.cores Number of cores.
# -... All pass to `consensus_partition`.
#
# == details
# The function performs consensus partitioning only with a small subset
# of columns and the class of other columns are predicted by `predict_classes,ConsensusPartition-method`.
#
# == example
# \dontrun{
# data(golub_cola)
# m = get_matrix(golub_cola)
#
# set.seed(123)
# golub_cola_ds = consensus_partition_by_down_sampling(m, subset = 50,
# 	anno = get_anno(golub_cola), anno_col = get_anno_col(golub_cola),
# 	top_value_method = "SD", partition_method = "kmeans")
# }
consensus_partition_by_down_sampling = function(data, 
	top_value_method = "ATC",
	top_n = seq(min(1000, round(nrow(data)*0.1)), 
		        min(3000, round(nrow(data)*0.3)), 
		        length.out = 3),
	partition_method = "skmeans",
	max_k = 6, 
	subset = min(round(ncol(data)*0.2), 250),
	verbose = TRUE, prefix = "", anno = NULL, anno_col = NULL,
	dist_method = c("euclidean", "correlation", "cosine"),
	.env = NULL, .predict = TRUE, mc.cores = 1, ...) {

	if(is.null(.env)) {
		if(is.data.frame(data)) data = as.matrix(data)

		.env = new.env(parent = emptyenv())
		.env$data = data
		.env$column_index = seq_len(ncol(data))
	} else if(is.null(.env$data)) {
		data = as.matrix(data)

		.env$data = data
		.env$column_index = seq_len(ncol(data))
	} else if(is.null(.env$column_index)) {
		data = .env$data
		.env$column_index = seq_len(ncol(data))
	} else {
		data = .env$data
	}

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

		if(nrow(anno) != length(.env$column_index)) {
			stop_wrap("nrow of `anno` should be the same as ncol of the matrix.")
		}
	}

	if(is.null(anno_col)) {
		anno_col = lapply(anno, ComplexHeatmap:::default_col)
	} else {
		if(ncol(anno) == 1 && is.atomic(anno_col)) {
			anno_col = list(anno_col)
			names(anno_col) = colnames(anno)
		} else if(is.null(names(anno_col))) {
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

	qqcat("@{prefix}* @{ifelse(length(subset) == 1, subset, length(subset))} columns are randomly sampled from @{length(.env$column_index)} columns.\n")

	column_index = .env$column_index
	if(length(subset) == 1) {
		if(is.null(.env$node_0_top_value_list)) {
			subset = sample(length(column_index), subset)
		} else {
			if(is.null(.env$node_0_top_value_list[[top_value_method]])) {
				subset = sample(length(column_index), subset)
			} else {
				qqcat("@{prefix}* assign sampling probability to columns by a pre-partitioning.\n")
				km = kmeans(t(.env$data[order(.env$node_0_top_value_list[[top_value_method]], decreasing = TRUE)[1:top_n[1]], .env$column_index, drop = FALSE]), centers = max_k)$cluster
				tb = table(km)
				p = km/tb[as.character(km)]
				subset = unique(sample(seq_along(column_index), subset, prob = p))
			}
		}
		subset_index = column_index[subset]
	} else {
		if(!is.numeric(subset)) {
			stop_wrap("If `subset` is specified as an indices vector, it should be in numeric.")
		}
		subset_index = column_index[subset]
	}
	.env$column_index = subset_index

	if(is.null(anno)) {
		anno2 = NULL
	} else {
		anno2 = anno[subset, , drop = FALSE]
	}
	
	# top_value_list cannot be repetitively used here
	.env$all_top_value_list = NULL
	cp = consensus_partition(.env = .env, top_value_method = top_value_method, partition_method = partition_method, 
		top_n = top_n, max_k = max_k, prefix = prefix, anno = anno2, anno_col = anno_col, mc.cores = mc.cores, ...)

	attr(cp, "full_anno") = anno
	
	if(.predict) {
		obj = convert_to_DownSamplingConsensusPartition(cp, column_index, dist_method, verbose, prefix, mc.cores)
		return(obj)
	} else {
		return(cp)
	}
}

convert_to_DownSamplingConsensusPartition = function(cp, column_index, dist_method, verbose, prefix, mc.cores) {

	data = cp@.env$data[, column_index, drop = FALSE]

	obj = new("DownSamplingConsensusPartition")
	for(nm in slotNames(cp)) {
		slot(obj, nm) = slot(cp, nm)
	}

	cl = list()
	if(cp@scale_rows) {
		if("z-score" %in% attr(cp@scale_rows, "scale_method")) {
			data2 = t(scale(t(data)))
		} else if("min-max" %in% attr(cp@scale_rows, "scale_method")) {
			row_min = rowMins(data)
			row_max = rowMaxs(data)
			row_range = row_max - row_min
			data2 = (data - row_min)/row_range
		}
	}

	for(k in cp@k) {
		qqcat("@{prefix}* predict class for @{ncol(data)} samples with k = @{k}\n")
		if(cp@scale_rows) {
			cl[[as.character(k)]] = predict_classes(cp, k = k, data2, p_cutoff = Inf, dist_method = dist_method, 
				plot = FALSE, verbose = verbose, force = TRUE, help = FALSE, prefix = qq("@{prefix}  "), mc.cores = mc.cores)
		} else {
			cl[[as.character(k)]] = predict_classes(cp, k = k, data, p_cutoff = Inf, dist_method = dist_method, 
				plot = FALSE, verbose = verbose, force = TRUE, help = FALSE, prefix = qq("@{prefix}  "), mc.cores = mc.cores)
		}
	}

	obj@full_anno = attr(cp, "full_anno")
	obj@full_column_index = column_index
	obj@predict$class = cl

	return(obj)

}

# == title
# Print the DownSamplingConsensusPartition object
#
# == param
# -object A `DownSamplingConsensusPartition-class` object.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola_ds)
# golub_cola_ds
setMethod(f = "show",
	signature = "DownSamplingConsensusPartition",
	definition = function(object) {

	# fix older version where there was no sample_by slot
	error = try(object@sample_by, silent = TRUE)
	if(inherits(error, "try-error")) {
		object@sample_by = "row"
	}
	qqcat("A 'DownSamplingConsensusPartition' object with k = @{paste(object@k, collapse = ', ')}.\n")
	qqcat("  On a matrix with @{nrow(object@.env$data)} rows and @{length(object@column_index)} columns, randomly sampled from @{length(object@full_column_index)} columns.\n")
	top_n_str = object@top_n
	qqcat("  Top rows (@{paste(top_n_str, collapse = ', ')}) are extracted by '@{object@top_value_method}' method.\n")
	qqcat("  Subgroups are detected by '@{object@partition_method}' method.\n")
	qqcat("  Performed in total @{object@n_partition} partitions by @{object@sample_by} resampling.\n")
	best_k = suggest_best_k(object)
	if(is.na(best_k)) {
		qqcat("  There is no best k.\n")
	} else {
		qqcat("  Best k for subgroups seems to be @{best_k}.\n")
	}
	qqcat("\n")
	qqcat("Following methods can be applied to this 'DownSamplingConsensusPartition' object:\n")
	txt = showMethods(classes = c("DownSamplingConsensusPartition", "ConsensusPartition"), where = topenv(), printTo = FALSE)
	txt = grep("Function", txt, value = TRUE)
	fname = gsub("Function: (.*?) \\(package.*$", "\\1", txt)
	print(unique(fname))

})

# == title
# Get annotations
#
# == param
# -object A `DownSamplingConsensusPartition-class` object.
# -reduce Used internally.
#
# == value
# A data frame if ``anno`` was specified in `consensus_partition_by_down_sampling`, or else ``NULL``.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola_ds)
# get_anno(golub_cola_ds)
setMethod(f = "get_anno",
	signature = "DownSamplingConsensusPartition",
	definition = function(object, reduce = FALSE) {
	if(reduce) {
		object@anno
	} else {
		object@full_anno
	}
})

# == title
# Get subgroup labels
#
# == param
# -object A `DownSamplingConsensusPartition-class` object.
# -k Number of subgroups.
# -p_cutoff Cutoff of p-values of class label prediction. It is only used when ``k`` is a vector.
# -reduce Used internally.
#
# == return
# If ``k`` is a scalar, it returns a data frame with two columns:
#
# - the class labels
# - the p-value for the prediction of class labels.
#
# If ``k`` is a vector, it returns a data frame of class labels for each k. The class
# label with prediction p-value > ``p_cutoff`` is set to ``NA``.
#
# == example
# data(golub_cola_ds)
# get_classes(golub_cola_ds, k = 3)
# get_classes(golub_cola_ds)
setMethod(f = "get_classes",
	signature = "DownSamplingConsensusPartition",
	definition = function(object, k = object@k, p_cutoff = 0.05, reduce = FALSE) {

	if(reduce) {
		return(selectMethod("get_classes", signature = "ConsensusPartition")(object, k = k))
	}
	if(length(k) == 1) {
		i = which(object@k == k)
		object@predict$class[[i]]
	} else {
		df = do.call("cbind", lapply(k, function(i) {
			df = object@predict$class[[as.character(i)]]
			cl = df[, "class"]
			cl[df[, "p"] > p_cutoff] = NA
			cl
		}))
		colnames(df) = paste0("k=", k)
		rownames(df) = rownames(object@predict$class[[1]])
		df
	}
})

# == title
# Test correspondance between predicted subgroups and known factors
#
# == param
# -object A `DownSamplingConsensusPartition-class` object.
# -k Number of subgroups. It uses all ``k`` if it is not specified.
# -known A vector or a data frame with known factors. By default it is the annotation table set in `consensus_partition_by_down_sampling`.
# -p_cutoff Cutoff for p-values for the class prediction. Samples with p-value higher than it are omit.
# -verbose Whether to print messages.
#
# == details
# The test is performed by `test_between_factors` between the predicted classes and user's annotation table.
# 
# == value
# A data frame with the following columns:
#
# - number of samples used to test after filtered by ``p_cutoff``,
# - p-values from the tests,
# - number of subgroups.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola_ds)
# test_to_known_factors(golub_cola_ds, k = 3)
# test_to_known_factors(golub_cola_ds)
setMethod(f = "test_to_known_factors",
	signature = "DownSamplingConsensusPartition",
	definition = function(object, k, known = get_anno(object), 
	p_cutoff = 0.05, verbose = FALSE) {

	if(missing(k)) {
		all_k = object@k
		m = do.call("rbind", lapply(all_k, function(k) test_to_known_factors(object, k, known = known, 
			p_cutoff = p_cutoff, verbose = verbose)))
		return(m)
	}

	class_df = get_classes(object, k)
	l = class_df$p <= p_cutoff
	class = as.character(class_df$class)[l]

	if(is.null(known)) {
		stop_wrap("`known` should not be NULL.")
	}

	if(is.data.frame(known)) {
		known = known[l, , drop = FALSE]
	} else if(is.matrix(known)) {
		known = known[l, ,drop = FALSE]
	} else {
		known = known[l]
	}

	m = test_between_factors(class, known, verbose = verbose)
	rownames(m) = paste(object@top_value_method, object@partition_method, sep = ":")
	colnames(m) = paste0(colnames(m), "(p-value)")
	m = cbind(n_sample = sum(l), m, k = k)
	return(m)
})

# == title
# Visualize column after dimension reduction
#
# == description
# Visualize samples (the matrix columns) after dimension reduction
#
# == param
# -object A `DownSamplingConsensusPartition-class` object.
# -k Number of subgroups.
# -top_n Top n rows to use. By default it uses all rows in the original matrix.
# -method Which method to reduce the dimension of the data. ``MDS`` uses `stats::cmdscale`,
#         ``PCA`` uses `stats::prcomp`. ``t-SNE`` uses `Rtsne::Rtsne`. ``UMAP`` uses
#         `umap::umap`.
# -color_by If annotation table is set, an annotation name can be set here.
# -control A list of parameters for `Rtsne::Rtsne` or `umap::umap`.
# -internal Internally used.
# -nr If number of matrix rows is larger than this value, random ``nr`` rows are used.
# -p_cutoff Cutoff of p-value of class label prediction. Data points with values higher
#        than it will be mapped with cross symbols.
# -remove Whether to remove columns which have high p-values than
#        the cutoff.
# -scale_rows Whether to perform scaling on matrix rows.
# -verbose Whether print messages.
# -... Other arguments.
#
# == details
# This function is basically very similar as `dimension_reduction,ConsensusPartition-method`.
#
# == value
# No value is returned.
#
# == example
# data(golub_cola_ds)
# dimension_reduction(golub_cola_ds, k = 2)
# dimension_reduction(golub_cola_ds, k = 3)
setMethod(f = "dimension_reduction",
	signature = "DownSamplingConsensusPartition",
	definition = function(object, k, top_n = NULL,
	method = c("PCA", "MDS", "t-SNE", "UMAP"), 
	control = list(), color_by = NULL,
	internal = FALSE, nr = 5000,
	p_cutoff = 0.05, remove = FALSE,
	scale_rows = TRUE, verbose = TRUE, ...) {

	# the following code is basically the same as dimension_reduction,ConsensusPartition-method

	if(missing(k)) stop_wrap("k needs to be provided.")

	cl = as.list(match.call())
	# default value
	if(! "method" %in% names(cl)) {
		method = NULL
		if(is.null(method)) {
			oe = try(loadNamespace("umap"), silent = TRUE)
			if(inherits(oe, "try-error")) {
				if(verbose) cat("umap package is not installed.\n")
			} else {
				if(verbose) cat("use UMAP\n")
				method = "UMAP"
			}
		}
		if(is.null(method)) {
			oe = try(loadNamespace("Rtsne"), silent = TRUE)
			if(inherits(oe, "try-error")) {
				if(verbose) cat("Rtsne package is not installed.\n")
			} else {
				if(verbose) cat("use t-SNE\n")
				method = "t-SNE"
			}
		}
		if(is.null(method)) {
			if(verbose) cat("use PCA\n")
			method = "PCA"
		}
	}

	method = match.arg(method)
	data = object@.env$data[, object@full_column_index, drop = FALSE]

	if(!is.null(top_n)) {
		top_n = min(c(top_n, nrow(data)))
		all_value = object@top_value_list
		ind = order(all_value)[1:top_n]
		if(length(ind) > nr) ind = sample(ind, nr)
		data = data[ind, , drop = FALSE]
	} else {
		top_n = nrow(data)
		if(nrow(data) > nr) data = data[sample(1:nrow(data), nr), , drop = FALSE]
	}

	class_df = get_classes(object, k)

	l = class_df$p <= p_cutoff
	
	op = par(c("mar", "xpd"))
	par(mar = c(4.1, 4.1, 4.1, 6), xpd = NA)
	on.exit(par(op))

	class_level = sort(as.character(unique(class_df$class)))
	n_class = length(class_level)
	
	if(is.null(color_by)) {
		col = cola_opt$color_set_2[as.character(class_df$class)]
	} else {
		if(!color_by %in% colnames(object@anno)) {
			stop_wrap("`color_by` should only contain the annotation names.")
		}
		col = object@anno_col[[color_by]][ object@anno[, color_by] ]
	}
	if(remove) {
		dimension_reduction(data[, l], pch = 16, col = col[l],
			cex = ifelse(which(l) %in% object@column_index, 1, 0.5), 
			main = qq("@{method} on @{top_n} rows with highest @{object@top_value_method} scores@{ifelse(scale_rows, ', rows are scaled', '')}\n@{sum(l)}/@{length(l)} confident samples (p < @{p_cutoff})"),
			method = method, control = control, scale_rows = scale_rows, nr = nr, internal = internal, verbose = verbose, ...)
		if(!internal) {
			legend(x = par("usr")[2], y = mean(par("usr")[3:4]), 
				legend = c(paste0("group", class_level), "predicted"), 
				pch = c(rep(16, n_class), 16),
				pt.cex = c(rep(1, n_class), 0.5),
				col = c(cola_opt$color_set_2[class_level], "black"), 
				xjust = 0, yjust = 0.5,
				title = "Class", title.adj = 0.1, bty = "n",
				text.col = c(rep("black", n_class), "black"))
		}
	} else {
		dimension_reduction(data, pch = ifelse(l, 16, 4), col = col,
			cex = ifelse(1:ncol(data) %in% object@column_index, 1, 0.5),
			main = qq("@{method} on @{top_n} rows with highest @{object@top_value_method} scores@{ifelse(scale_rows, ', rows are scaled', '')}\n@{sum(l)}/@{length(l)} confident samples (p < @{p_cutoff})"),
			method = method, control = control, scale_rows = scale_rows, nr = nr, internal = internal, verbose = verbose, ...)
		if(!internal) {
			if(any(!l)) {
				legend(x = par("usr")[2], y = mean(par("usr")[3:4]), 
					legend = c(paste0("group", class_level), "ambiguous", "predicted", "predicted"), 
						pch = c(rep(16, n_class), 4, 16, 4),
						pt.cex = c(rep(1, n_class), 1, 0.5, 0.5),
						col = c(cola_opt$color_set_2[class_level], "black", "black", "black"), 
						xjust = 0, yjust = 0.5,
						title = "Class", title.adj = 0.1, bty = "n")
			} else {
				legend(x = par("usr")[2], y = mean(par("usr")[3:4]), 
					legend = c(paste0("group", class_level), "predicted"), 
					pch = c(rep(16, n_class), 16),
					pt.cex = c(rep(1, n_class), 0.5),
					col = c(cola_opt$color_set_2[class_level], "black"), 
					xjust = 0, yjust = 0.5,
					title = "Class", title.adj = 0.1, bty = "n",
					text.col = c(rep("black", n_class), "black"))
			}
		}
	}
})

# == title
# Get signature rows
#
# == param
# -object A `DownSamplingConsensusPartition-class` object.
# -k Number of subgroups.
# -p_cutoff Cutoff for p-values of class label prediction. Samples with values 
#        higher than it are not used for finding signature rows.
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
# This function is very similar as `get_signatures,ConsensusPartition-method`.
#
# == example
# \donttest{
# data(golub_cola_ds)
# get_signatures(golub_cola_ds, k = 2)
# get_signatures(golub_cola_ds, k = 3)
# }
setMethod(f = "get_signatures",
	signature = "DownSamplingConsensusPartition",
	definition = function(object, k,
	p_cutoff = 0.05, 
	fdr_cutoff = cola_opt$fdr_cutoff, 
	group_diff = cola_opt$group_diff,
	scale_rows = object@scale_rows,
	row_km = NULL,
	diff_method = c("Ftest", "ttest", "samr", "pamr", "one_vs_others", "uniquely_high_in_one_group"),
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
	
	class_df = get_classes(object, k)
	class_ids = class_df$class

	data = object@.env$data[, object@full_column_index, drop = FALSE]

	l = class_df$p <= p_cutoff
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

	if(verbose) qqcat("* @{n_sample_used}/@{nrow(class_df)} samples (in @{length(unique(class))} classes) remain after filtering by p-value (<= @{p_cutoff}).\n")

	tb = table(class)
	if(sum(tb > 1) <= 1) {
		if(plot) {
			grid.newpage()
			fontsize = convertUnit(unit(0.1, "npc"), "char", valueOnly = TRUE)*get.gpar("fontsize")$fontsize
			grid.text("not enough samples", gp = gpar(fontsize = fontsize))
		}
		if(verbose) cat("* Not enough samples.\n")
		return(invisible(data.frame(which_row = integer(0))))
	}
	if(length(unique(class)) <= 1) {
		if(plot) {
			grid.newpage()
			fontsize = convertUnit(unit(0.1, "npc"), "char", valueOnly = TRUE)*get.gpar("fontsize")$fontsize
			grid.text("not enough classes", gp = gpar(fontsize = fontsize))
		}
		if(verbose) cat("* Not enough classes.\n")
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

	if(diff_method %in% c("samr", "pamr")) {
		hash = digest(list(used_samples = which(l), 
			               class = class,
			               n_group = k, 
			               diff_method = diff_method,
			               column_index = object@column_index,
			               fdr_cutoff = fdr_cutoff,
			               group_diff = group_diff,
			               seed = seed),
					algo = "md5")
	} else {
		hash = digest(list(used_samples = which(l), 
			               class = class,
			               n_group = k, 
			               diff_method = diff_method,
			               column_index = object@column_index,
			               group_diff = group_diff,
			               seed = seed),
					algo = "md5")
	}
	nm = paste0("signature_fdr_", hash)
	if(verbose) qqcat("* cache hash: @{hash} (seed @{seed}).\n")

	find_signature = TRUE
	if(!is.null(object@.env[[nm]])) {
		if(diff_method == "samr") {
			if(object@.env[[nm]]$diff_method == "samr" && 
			   object@.env[[nm]]$n_sample_used == n_sample_used && 
			   abs(object@.env[[nm]]$fdr_cutoff - fdr_cutoff) < 1e-10) {
				fdr = object@.env[[nm]]$fdr
				find_signature = FALSE
			}
		} else if(diff_method == "pamr") {
			if(object@.env[[nm]]$diff_method == "pamr" && 
			   object@.env[[nm]]$n_sample_used == n_sample_used && 
			   abs(object@.env[[nm]]$fdr_cutoff - fdr_cutoff) < 1e-10) {
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
					returned_obj$km = apply(pdist(row_km_fit$centers, mat_for_km, as.integer(1)), 2, which.min)
					do_kmeans = FALSE
					if(verbose) qqcat("* use k-means partition that are already calculated in previous runs.\n")
				}
			}
			if(do_kmeans) {
				set.seed(seed)
				if(is.null(row_km)) {
					wss = (nrow(mat_for_km2)-1)*sum(apply(mat_for_km2,1,var))
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
					returned_obj$km = apply(pdist(row_km_fit$centers, mat_for_km, as.integer(1)), 2, which.min)
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

	if(internal) {
		ha1 = HeatmapAnnotation(
				Class = class_df$class[column_used_logical],
				col = list(Class = cola_opt$color_set_2),
				show_annotation_name = !has_ambiguous & !internal,
				annotation_name_side = "right",
				show_legend = TRUE)
	} else {
		if(simplify) {
			ha1 = HeatmapAnnotation(
				Class = class_df$class[column_used_logical],
				col = list(Class = cola_opt$color_set_2),
				show_annotation_name = !has_ambiguous & !internal,
				annotation_name_side = "right",
				show_legend = TRUE)
		} else {
			ha1 = HeatmapAnnotation(
				Class = class_df$class[column_used_logical],
				col = list(Class = cola_opt$color_set_2),
				show_annotation_name = !has_ambiguous & !internal,
				annotation_name_side = "right",
				show_legend = TRUE)
		}
	}
	ht_list = ht_list + Heatmap(use_mat1, name = heatmap_name, col = col_fun,
		top_annotation = ha1, row_split = row_split,
		cluster_columns = TRUE, cluster_column_slices = FALSE, cluster_row_slices = FALSE,
		column_split = factor(class_df$class[column_used_logical], levels = sort(unique(class_df$class[column_used_logical]))), 
		show_column_dend = FALSE,
		show_row_names = FALSE, show_row_dend = show_row_dend, column_title = {if(internal) NULL else qq("@{ncol(use_mat1)} confident samples")},
		use_raster = use_raster, raster_by_magick = requireNamespace("magick", quietly = TRUE),
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
			ha2 = HeatmapAnnotation(
				Class = class_df$class[!column_used_logical],
				col = list(Class = cola_opt$color_set_2),
				show_annotation_name = !internal,
				annotation_name_side = "right",
				show_legend = FALSE)
		} else {
			if(simplify) {
				ha2 = HeatmapAnnotation(
					Class = class_df$class[!column_used_logical],
					col = list(Class = cola_opt$color_set_2),
					show_annotation_name = c(TRUE) & !internal,
					annotation_name_side = "right",
					show_legend = FALSE)
			} else {
				ha2 = HeatmapAnnotation(
					Class = class_df$class[!column_used_logical],
					col = list(Class = cola_opt$color_set_2),
					show_annotation_name = c(TRUE, TRUE) & !internal,
					annotation_name_side = "right",
					show_legend = FALSE)
			}
			
		}
		ht_list = ht_list + Heatmap(use_mat2, name = paste0(heatmap_name, 2), col = col_fun,
			top_annotation = ha2,
			cluster_columns = TRUE, show_column_dend = FALSE,
			show_row_names = FALSE, show_row_dend = FALSE, show_heatmap_legend = FALSE,
			use_raster = use_raster, raster_by_magick = requireNamespace("magick", quietly = TRUE),
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

	return(invisible(returned_obj))
})



# == title
# Number of columns in the matrix
#
# == param
# -x A `DownSamplingConsensusPartition-class` object.
#
setMethod(f = "ncol",
	signature = "DownSamplingConsensusPartition",
	definition = function(x) {
	length(x@full_column_index)
})


# == title
# Column names of the matrix
#
# == param
# -x A `DownSamplingConsensusPartition-class` object.
#
setMethod(f = "colnames",
	signature = "DownSamplingConsensusPartition",
	definition = function(x) {
	colnames(x@.env$data)[x@full_column_index]
})


# == title
# Dimension of the matrix
#
# == param
# -x A `DownSamplingConsensusPartition-class` object.
#
dim.DownSamplingConsensusPartition = function(x) {
	c(nrow(x), ncol(x))
}


