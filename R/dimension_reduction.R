
# == title
# Visualize column after dimension reduction
#
# == description
# Visualize samples (the matrix columns) after dimension reduction
#
# == param
# -object A `ConsensusPartition-class` object.
# -k Number of subgroups.
# -top_n Top n rows to use. By default it uses all rows in the original matrix.
# -method Which method to reduce the dimension of the data. ``MDS`` uses `stats::cmdscale`,
#         ``PCA`` uses `stats::prcomp`. ``t-SNE`` uses `Rtsne::Rtsne`. ``UMAP`` uses
#         `umap::umap`.
# -color_by If annotation table is set, an annotation name can be set here.
# -control A list of parameters for `Rtsne::Rtsne` or `umap::umap`.
# -internal Internally used.
# -nr If number of matrix rows is larger than this value, random ``nr`` rows are used.
# -silhouette_cutoff Cutoff of silhouette score. Data points with values less
#        than it will be mapped with cross symbols.
# -remove Whether to remove columns which have less silhouette scores than
#        the cutoff.
# -scale_rows Whether to perform scaling on matrix rows.
# -verbose Whether print messages.
# -... Other arguments.
#
# == value
# Locations of the points.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# dimension_reduction(golub_cola["ATC", "skmeans"], k = 3)
setMethod(f = "dimension_reduction",
	signature = "ConsensusPartition",
	definition = function(object, k, top_n = NULL,
	method = c("PCA", "MDS", "t-SNE", "UMAP"), 
	control = list(), color_by = NULL,
	internal = FALSE, nr = 5000,
	silhouette_cutoff = 0.5, remove = FALSE,
	scale_rows = object@scale_rows, verbose = TRUE, ...) {

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
	data = object@.env$data[object@row_index, object@column_index, drop = FALSE]

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

	l = class_df$silhouette >= silhouette_cutoff
	
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
		loc = dimension_reduction(data[, l], pch = 16, col = col[l],
			cex = 1, main = qq("@{method} on @{top_n} rows with highest @{object@top_value_method} scores@{ifelse(scale_rows, ', rows are scaled', '')}\n@{sum(l)}/@{length(l)} confident samples (silhouette > @{silhouette_cutoff})"),
			method = method, control = control, scale_rows = scale_rows, nr = nr, internal = internal, verbose = verbose, ...)
		if(!internal) {
			if(is.null(color_by)) {
				legend(x = par("usr")[2], y = mean(par("usr")[3:4]), legend = c(paste0("group", class_level), "ambiguous"), 
					pch = c(rep(16, n_class), 0),
					col = c(cola_opt$color_set_2[class_level], "white"), xjust = 0, yjust = 0.5,
					title = "Class", title.adj = 0.1, bty = "n",
					text.col = c(rep("black", n_class), "white"))
			} else {
				legend(x = par("usr")[2], y = mean(par("usr")[3:4]), legend = names(object@anno_col[[color_by]]), 
					pch = 16,
					col = object@anno_col[[color_by]], xjust = 0, yjust = 0.5,
					title = color_by, title.adj = 0.1, bty = "n")
			}
		}
	} else {
		loc = dimension_reduction(data, pch = ifelse(l, 16, 4), col = col,
			cex = 1, main = qq("@{method} on @{top_n} rows with highest @{object@top_value_method} scores@{ifelse(scale_rows, ', rows are scaled', '')}\n@{sum(l)}/@{length(l)} confident samples (silhouette > @{silhouette_cutoff})"),
			method = method, control = control, scale_rows = scale_rows, nr = nr, internal = internal, verbose = verbose, ...)
		if(!internal) {
			if(any(!l)) {
				if(is.null(color_by)) {
					legend(x = par("usr")[2], y = mean(par("usr")[3:4]), legend = c(paste0("group", class_level), "ambiguous"), 
							pch = c(rep(16, n_class), 4),
							col = c(cola_opt$color_set_2[class_level], "black"), xjust = 0, yjust = 0.5,
							title = "Class", title.adj = 0.1, bty = "n")
				} else {
					legend(x = par("usr")[2], y = mean(par("usr")[3:4]), legend = names(object@anno_col[[color_by]]), 
						pch = 16,
						col = object@anno_col[[color_by]], xjust = 0, yjust = 0.5,
						title = color_by, title.adj = 0.1, bty = "n")
				}
			} else {
				if(is.null(color_by)) {
					legend(x = par("usr")[2], y = mean(par("usr")[3:4]), legend = c(paste0("group", class_level), "ambiguous"), 
						pch = c(rep(16, n_class), 0),
						col = c(cola_opt$color_set_2[class_level], "white"), xjust = 0, yjust = 0.5,
						title = "Class", title.adj = 0.1, bty = "n",
						text.col = c(rep("black", n_class), "white"))
				} else {
					legend(x = par("usr")[2], y = mean(par("usr")[3:4]), legend = names(object@anno_col[[color_by]]), 
						pch = 16,
						col = object@anno_col[[color_by]], xjust = 0, yjust = 0.5,
						title = color_by, title.adj = 0.1, bty = "n")
				}
			}
		}
	}
	return(invisible(loc))
})


# == title
# Visualize columns after dimension reduction
#
# == param
# -object A numeric matrix.
# -method Which method to reduce the dimension of the data. ``MDS`` uses `stats::cmdscale`,
#         ``PCA`` uses `stats::prcomp`. ``t-SNE`` uses `Rtsne::Rtsne`. ``UMAP`` uses
#         `umap::umap`.
# -pc Which two principle components to visualize
# -control A list of parameters for `Rtsne::Rtsne` or `umap::umap`.
# -pch Ahape of points.
# -col Color of points.
# -cex Aize of points.
# -main Title of the plot.
# -scale_rows Whether perform scaling on matrix rows.
# -nr If number of matrix rows is larger than this value, random ``nr`` rows are used.
# -internal Internally used.
# -verbose Whether print messages.
#
# == value
# Locations of the points.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "dimension_reduction",
	signature = "matrix",
	definition = function(object, 
	pch = 16, col = "black", cex = 1, main = NULL,
	method = c("PCA", "MDS", "t-SNE", "UMAP"),
	pc = NULL, control = list(), 
	scale_rows = FALSE, nr = 5000,
	internal = FALSE, verbose = TRUE) {

	data = object

	if(nrow(data) > nr) {
		data = data[sample(1:nrow(data), nr), , drop = FALSE]
	}
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

	method = match.arg(method)[1]

	if(is.null(main)) {
		main = qq("@{method} on a matrix with @{nrow(object)} rows @{ifelse(scale_rows, ', rows are scaled', '')}")
	}

	if(is.null(pc)) {
		if(method == "PCA") {
			pc = 1:2
		} else if(method %in% c("UMAP", "t-SNE")) {
			pc = seq_len(min(c(10, ncol(data))))
		}
	} else if(length(pc) == 1) {
		pc = 1:pc
	}
	if(method %in% c("UMAP", "t-SNE")) {
		pc = 1:max(pc)
	}

	if(method %in% c("UMAP", "t-SNE")) {
		main = paste0(main, qq(", with @{length(pc)} PCs"))
	}

	if(scale_rows) data = t(scale(t(data)))
	l = apply(data, 1, function(x) any(is.na(x)))
	data = data[!l, ]

	if(internal) {
		omar = par("mar")
		par(mar = c(4.1, 4.1, 1, 1))
		on.exit(par(mar = omar))

		main = NULL
	}

	if(length(col) == 1) col = rep(col, ncol(data))

	## format pch
	n_col = length(unique(col))
	if(length(pch) == 1) {
		pch = rep(pch, ncol(data))
	} else if(length(pch) == n_col) {
		col2index = split(1:ncol(data), col)
		pch2 = numeric(ncol(data))
		for(i in seq_along(col2index)) {
			pch2[col2index[[i]]] = pch[i]
		}
		pch = pch2
	}

	if(method == "MDS") {
		loc = cmdscale(dist(t(data)))
		plot(loc, pch = pch, col = col, cex = cex, main = main, xlab = "Coordinate 1", ylab = "Coordinate 2")
	} else if(method == "PCA") {
		fit = prcomp_irlba(t(data), n = max(pc))
		sm = summary(fit)
		prop = sm$importance[2, pc]
		loc = fit$x[, pc]

		plot(loc, pch = pch, col = col, cex = cex, main = main, xlab = qq("PC@{pc[1]} (@{round(prop[1]*100)}%)"), ylab = qq("PC@{pc[2]} (@{round(prop[2]*100)}%)"))
	} else if(method == "t-SNE") {

		check_pkg("Rtsne", bioc = FALSE)
		fit = prcomp_irlba(t(data), n = max(pc))
		sm = summary(fit)
		loc = fit$x[, pc]

		param = list(X = loc)
		if(! "perplexity" %in% names(control)) {
			control$perplexity = min(30, floor(ncol(data)-1)/3)
		}
		control$pca = FALSE
		param = c(param, control)
		param$pca = FALSE
		fit = do.call(Rtsne::Rtsne, param)
		loc = fit$Y
		plot(loc, pch = pch, col = col, cex = cex, main = main, xlab = "t-SNE 1", ylab = "t-SNE 2")
	} else if(method == "UMAP") {

		check_pkg("umap", bioc = FALSE)

		fit = prcomp_irlba(t(data), n = max(pc))
		sm = summary(fit)
		loc = fit$x[, pc]

		param = list(d = loc)
		if(!"config" %in% names(control)) {
			# reset n_neighbors for small dataset
			control$config = umap::umap.defaults
			control$config$n_neighbors = min(umap::umap.defaults$n_neighbors, 10)
		} else {
			if(!"n_neighbors" %in% names(control$config)) {
				# reset n_neighbors for small dataset
				control$config$n_neighbors = min(umap::umap.defaults$n_neighbors, 10)
			}
		}
		param = c(param, control)
		fit = do.call(umap::umap, param)
		loc = fit$layout
		plot(loc, pch = pch, col = col, cex = cex, main = main, xlab = "UMAP 1", ylab = "UMAP 2")
	}

	return(invisible(loc))
})
