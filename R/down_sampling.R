
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
# -k Alternatively, you can specify a vector k.
# -subset Number of columns to randomly sample, or a vector of selected indices.
# -pre_select Whether to pre-select by k-means.
# -verbose Whether to print messages.
# -prefix Internally used.
# -anno Annotation data frame.
# -anno_col Annotation colors.
# -predict_method Method for predicting class labels. Possible values are "centroid", "svm" and "randomForest".
# -dist_method Method for predict the class for other columns.
# -.env An environment, internally used.
# -.predict Internally used.
# -mc.cores Number of cores. This argument will be removed in future versions.
# -cores Number of cores, or a ``cluster`` object returned by `parallel::makeCluster`.
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
	top_n = NULL,
	partition_method = "skmeans",
	max_k = 6, k = NULL,
	subset = min(round(ncol(data)*0.2), 250), pre_select = TRUE,
	verbose = TRUE, prefix = "", anno = NULL, anno_col = NULL,
	predict_method = "centroid",
	dist_method = c("euclidean", "correlation", "cosine"),
	.env = NULL, .predict = TRUE, mc.cores = 1, cores = mc.cores, ...) {

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

	
	column_index = .env$column_index

	if(is.null(.env$row_index)) {
		.env$row_index = seq_len(nrow(data))
	}

	if(length(subset) == 1) {
		# top n smallest knn mean distance
		subset = sample(length(column_index), subset)
		# qqcat("@{prefix}* @{ifelse(length(subset) == 1, subset, length(subset))} columns are picked from @{length(.env$column_index)} columns based on distance to nearest columns.\n")
		# qqcat("@{prefix}* calculate @{top_value_method} scores based on the complete matrix.\n")
		# top_value = get_top_value_method(top_value_method)(data[.env$row_index, .env$column_index])
		# if(is.null(top_n)) {
		# 	m = data[.env$row_index, .env$column_index][order(top_value, decreasing = TRUE)[1:round(length(.env$row_index)*0.1)], ]
		# } else {
		# 	m = data[.env$row_index, .env$column_index][order(top_value, decreasing = TRUE)[1:min(top_n)], ]
		# }
		# kn = min(10, round(ncol(dm)*0.1))
		# qqcat("@{prefix}* calculate mean distance to the nearest @{kn} columns.\n")
		# nearest_dm = ATC(t(m), k_neighbours = kn, cores = cores)
		# subset = order(nearest_dm, decreasing = TRUE)[1:subset]

		subset_index = column_index[subset]
	} else {
		if(!is.numeric(subset)) {
			stop_wrap("If `subset` is specified as an indices vector, it should be numeric.")
		}
		subset_index = column_index[subset]
	}
	.env$column_index = subset_index

	if(is.null(anno)) {
		anno2 = NULL
	} else {
		anno2 = anno[subset, , drop = FALSE]
	}

	if(verbose) qqcat("@{prefix}* apply consensus_partition_by_down_sampling() with @{length(subset_index)} columns.\n")

	# top_value_list cannot be repetitively used here
	.env$all_top_value_list = NULL
	cp = consensus_partition(.env = .env, top_value_method = top_value_method, partition_method = partition_method, 
		top_n = top_n, max_k = max_k, k = k, prefix = prefix, anno = anno2, anno_col = anno_col, cores = cores, ...)

	attr(cp, "full_anno") = anno
	
	if(.predict) {
		obj = convert_to_DownSamplingConsensusPartition(cp, column_index, predict_method, dist_method, verbose, prefix, cores)
		return(obj)
	} else {
		return(cp)
	}
}

convert_to_DownSamplingConsensusPartition = function(cp, column_index, predict_method, dist_method, verbose, prefix, cores) {

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
		if(verbose) qqcat("@{prefix}* predict class for @{ncol(data)} samples with k = @{k}\n")
		if(cp@scale_rows) {
			cl[[as.character(k)]] = predict_classes(cp, k = k, data2, p_cutoff = Inf, method = predict_method, dist_method = dist_method, 
				plot = FALSE, verbose = verbose, force = TRUE, help = FALSE, prefix = qq("@{prefix}  "), cores = cores)
		} else {
			cl[[as.character(k)]] = predict_classes(cp, k = k, data, p_cutoff = Inf, method = predict_method, dist_method = dist_method, 
				plot = FALSE, verbose = verbose, force = TRUE, help = FALSE, prefix = qq("@{prefix}  "), cores = cores)
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
	best_k = suggest_best_k(object, help = FALSE)
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
		ind = order(all_value, decreasing = TRUE)[1:top_n]
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
		col = object@anno_col[[color_by]][ object@full_anno[, color_by] ]
	}
	if(remove) {
		dimension_reduction(data[, l], pch = 16, col = col[l],
			cex = ifelse(which(l) %in% object@column_index, 1, 0.5), 
			main = qq("@{method} on @{top_n} rows with highest @{object@top_value_method} scores@{ifelse(scale_rows, ', rows are scaled', '')}\n@{sum(l)}/@{length(l)} confident samples (p < @{p_cutoff})"),
			method = method, control = control, scale_rows = scale_rows, nr = nr, internal = internal, verbose = verbose, ...)
		if(!internal) {
			if(is.null(color_by)) {
				legend(x = par("usr")[2], y = mean(par("usr")[3:4]), 
					legend = c(paste0("group", class_level), "predicted"), 
					pch = c(rep(16, n_class), 16),
					pt.cex = c(rep(1, n_class), 0.5),
					col = c(cola_opt$color_set_2[class_level], "black"), 
					xjust = 0, yjust = 0.5,
					title = "Class", title.adj = 0.1, bty = "n",
					text.col = c(rep("black", n_class), "black"))
			} else {
				legend(x = par("usr")[2], y = mean(par("usr")[3:4]), legend = names(object@anno_col[[color_by]]), 
					pch = 16,
					col = object@anno_col[[color_by]], xjust = 0, yjust = 0.5,
					title = color_by, title.adj = 0.1, bty = "n")
			}
		}
	} else {
		dimension_reduction(data, pch = ifelse(l, 16, 4), col = col,
			cex = ifelse(1:ncol(data) %in% object@column_index, 1, 0.5),
			main = qq("@{method} on @{top_n} rows with highest @{object@top_value_method} scores@{ifelse(scale_rows, ', rows are scaled', '')}\n@{sum(l)}/@{length(l)} confident samples (p < @{p_cutoff})"),
			method = method, control = control, scale_rows = scale_rows, nr = nr, internal = internal, verbose = verbose, ...)
		if(!internal) {
			if(any(!l)) {
				if(is.null(color_by)) {
					legend(x = par("usr")[2], y = mean(par("usr")[3:4]), 
						legend = c(paste0("group", class_level), "ambiguous", "predicted", "predicted"), 
							pch = c(rep(16, n_class), 4, 16, 4),
							pt.cex = c(rep(1, n_class), 1, 0.5, 0.5),
							col = c(cola_opt$color_set_2[class_level], "black", "black", "black"), 
							xjust = 0, yjust = 0.5,
							title = "Class", title.adj = 0.1, bty = "n")
				} else {
					legend(x = par("usr")[2], y = mean(par("usr")[3:4]), legend = names(object@anno_col[[color_by]]), 
						pch = 16,
						col = object@anno_col[[color_by]], xjust = 0, yjust = 0.5,
						title = color_by, title.adj = 0.1, bty = "n")
				}
			} else {
				if(is.null(color_by)) {
					legend(x = par("usr")[2], y = mean(par("usr")[3:4]), 
						legend = c(paste0("group", class_level), "predicted"), 
						pch = c(rep(16, n_class), 16),
						pt.cex = c(rep(1, n_class), 0.5),
						col = c(cola_opt$color_set_2[class_level], "black"), 
						xjust = 0, yjust = 0.5,
						title = "Class", title.adj = 0.1, bty = "n",
						text.col = c(rep("black", n_class), "black"))
				} else {
					legend(x = par("usr")[2], y = mean(par("usr")[3:4]), legend = names(object@anno_col[[color_by]]), 
						pch = 16,
						col = object@anno_col[[color_by]], xjust = 0, yjust = 0.5,
						title = color_by, title.adj = 0.1, bty = "n")
				}
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
# -... Other arguments passed to `get_signatures,ConsensusPartition-method`.
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
	p_cutoff = 1, ...) {

	callNextMethod(object, k, p_cutoff = p_cutoff, ...)
	
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
# Get the original matrix
#
# == param
# -object A `DownSamplingConsensusPartition-class` object.
# -reduce Whether to return the reduced matrix where columns are randomly sampled.
#
# == value
# A numeric matrix
setMethod(f = "get_matrix",
	signature = "DownSamplingConsensusPartition",
	definition = function(object, reduce = FALSE) {

	if(reduce) {
		object@.env$data[, object@column_index, drop = FALSE]
	} else {
		object@.env$data[, object@full_column_index, drop = FALSE]
	}
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


