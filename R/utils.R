
if(!exists("strrep")) {
	strrep = function(x, times) {
		x = rep(x, times = times)
		return(paste0(x, collapse = ""))
	}
}

# == title
# Relabel class labels according to the reference labels
#
# == param
# -class A vector of class labels.
# -ref A vector of reference labels.
# -full_set The full set of labels. 	
# -return_map Whether return the mapping or the adjusted labels.
#
# == details
# In partitions, the exact value of the class label is not of importance. E.g. for two partitions
# ``a, a, a, b, b, b, b`` and ``b, b, b, a, a, a, a``, they are the same partitions although the labels
# of ``a`` and ``b`` are switched in the two partitions. Even the partition ``c, c, c, d, d, d, d``
# is the same as the previous two although it uses a different set of labels. Here `relabel_class` function relabels
# ``class`` vector according to the labels in ``ref`` vector by looking for a mapping ``m()`` to maximize ``sum(m(class) == ref)``.
#
# Mathematically, this is called linear sum assignment problem and it is solved by `clue::solve_LSAP`.
#
# == value
# A named vector where names correspond to the labels in ``class`` and values correspond to ``ref``,
# which means ``map = relabel_class(class, ref); map[class]`` returns the relabelled labels.
#
# The returned object attaches a data frame with three columns:
#
# - original labels. in ``class``
# - adjusted labels. according to ``ref``
# - reference labels. in ``ref``
#
# If ``return_map`` in the `relabel_class` is set to `FALSE`, the function simply returns
# a vector of adjusted class labels.
#
# If the function returns the mapping vector (when ``return_map = TRUE``), the mapping variable
# is always character, which means, if your ``class`` and ``ref`` are numeric, you need to convert
# them back to numeric explicitely. If ``return_map = FALSE``, the returned relabelled vector has
# the same mode as ``class``.
#
# == example
# class = c(rep("a", 10), rep("b", 3))
# ref = c(rep("b", 4), rep("a", 9))
# relabel_class(class, ref)
# relabel_class(class, ref, return_map = FALSE)
# # if class and ref are from completely different sets
# class = c(rep("A", 10), rep("B", 3))
# relabel_class(class, ref)
#
# # class labels are numeric
# class = c(rep(1, 10), rep(2, 3))
# ref = c(rep(2, 4), rep(1, 9))
# relabel_class(class, ref)
# relabel_class(class, ref, return_map = FALSE)
relabel_class = function(class, ref, full_set = union(class, ref), return_map = TRUE) {
	md = mode(class)

	class = as.character(class)
	ref = as.character(ref)
	full_set = as.character(full_set)

	# if(length(intersect(class, ref)) == 0) {
	# 	stop("class and ref must be from same set of labels.")
	# }
	all_class = union(class, ref)
	n = length(all_class)

	m = matrix(0, nrow = n, ncol = n, dimnames = list(all_class, all_class))
	tb = table(class, ref)
	m[rownames(tb), colnames(tb)] = tb

	imap = clue::solve_LSAP(m, maximum = TRUE)
	map = structure(rownames(m)[imap], names = rownames(m))
	map = map[unique(class)]
	unmapped = setdiff(setdiff(full_set, class), names(map))
	if(length(unmapped)) {
		map = c(map, structure(unmapped, names = unmapped))
	}
	    	
	df = data.frame(
		class = class,
		adjusted = map[class], 
		ref = ref,
		stringsAsFactors = FALSE)
	attr(map, "df") = df

	if(return_map) {
		return(map)
	} else {
		return(as(df$adjusted, md))
	}
}

# columns only in one same level are clustered
column_order_by_group = function(factor, mat) {
	if(!is.factor(factor)) {
		factor = factor(factor, levels = unique(factor))
	}
	do.call("c", lapply(levels(factor), function(le) {
		m = mat[, factor == le, drop = FALSE]
		if(ncol(m) == 1) {
			which(factor == le)
		} else {
			hc1 = hclust(dist(t(m)))
			oe = try({ 
				hc1 = as.hclust(reorder(as.dendrogram(hc1), colSums(m)))
			}, silent = TRUE)
			col_od1 = hc1$order
			which(factor == le)[col_od1]
		}
	}))
}


# == title
# Adjust outliers
#
# == param
# -x A numeric vector.
# -q Percentile to adjust.
# 
# == details
# Vaules larger than percentile ``1 - q`` are adjusted to the ``1 - q`` percentile and 
# values smaller than percentile ``q`` are adjusted to the ``q`` percentile
#
# == value
# A numeric vector with same length as the original one.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# set.seed(123)
# x = rnorm(40)
# x[1] = 100
# adjust_outlier(x)
adjust_outlier = function(x, q = 0.05) {
	qu = quantile(x, c(q, 1 - q), na.rm = TRUE)
	x[x < qu[1]] = qu[1]
	x[x > qu[2]] = qu[2]
	x
}

# == title
# Remove rows with low variance and impute missing values
#
# == param
# -m A numeric matrix.
# -sd_quantile Cutoff of the quantile of standard deviation. Rows with standard deviation less than it are removed.
# -max_na Maximum NA fraction in each row. Rows with NA fraction larger than it are removed.
#
# == details
# The function uses `impute::impute.knn` to impute missing values, then
# uses `adjust_outlier` to adjust outliers and 
# removes rows with low standard deviations.
#
# == value
# A numeric matrix.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# set.seed(123)
# m = matrix(rnorm(100), nrow = 10)
# m[sample(length(m), 5)] = NA
# m[1, ] = 0
# m
# m2 = adjust_matrix(m)
# m2
adjust_matrix = function(m, sd_quantile = 0.05, max_na = 0.25) {
	if(is.data.frame(m)) {
		m = as.matrix(m)
	}

	l = apply(m, 1, function(x) sum(is.na(x))/length(x)) < max_na
	if(sum(!l)) {
		message_wrap(qq("removed @{sum(!l)} rows where more than @{round(max_na*100)}% of samples have NA values."))
	}
	m = m[l, , drop = FALSE]

	

	n_na = sum(is.na(m))
	if(n_na > 0) {
		message_wrap(qq("There are NA values in the data, now impute missing data."))
		
		on.exit(if(sink.number()) sink(NULL))	
		tempf = tempfile()
		sink(tempf)
		
		m = impute.knn(m)$data

		sink(NULL)
		file.remove(tempf)
	}
	m = t(apply(m, 1, adjust_outlier))
	row_sd = rowSds(m, na.rm = TRUE)
	l = abs(row_sd) == 0
	m2 = m[!l, , drop = FALSE]
	row_sd = row_sd[!l]
	if(sum(l)) {
		message_wrap(qq("@{sum(l)} rows have been removed with zero variance."))
	}

	qa = quantile(unique(row_sd), sd_quantile, na.rm = TRUE)
	l = row_sd > qa
	if(sum(!l)) {
		message_wrap(qq("@{sum(!l)} rows have been removed with too low variance (sd <= @{sd_quantile} quantile)"))
	}
	m2[l, , drop = FALSE]
}

dev.off2 = ComplexHeatmap:::dev.off2
dev.null = ComplexHeatmap:::dev.null

# https://stackoverflow.com/questions/15282471/get-stack-trace-on-trycatched-error-in-r
try_and_trace = function(expr, message) {
	withCallingHandlers(expr, 
		error = function(e) {
			cat(message, "\n")
			cat("Following is the call stack:\n")
			print(sys.calls())
		}
	)
}


stop_wrap = function (...) {
    x = paste0(...)
    x = paste(strwrap(x), collapse = "\n")
    stop(x, call. = FALSE)
}

warning_wrap = function (...) {
    x = paste0(...)
    x = paste(strwrap(x), collapse = "\n")
    warning(x, call. = FALSE)
}

message_wrap = function (...) {
    x = paste0(...)
    x = paste(strwrap(x), collapse = "\n")
    message(x)
}

# http://conjugateprior.org/2015/06/identifying-the-os-from-r/
os_type <- function() {
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "OSX"
  } else {
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "OSX"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

multicore_supported = function() {
	!(os_type() %in% c("osx", "windows"))
}


# == title
# Number of rows in the matrix
#
# == param
# -x A `ConsensusPartition-class` object.
#
setMethod(f = "nrow",
	signature = "ConsensusPartition",
	definition = function(x) {
	nrow(x@.env$data)
})


# == title
# Number of rows in the matrix
#
# == param
# -x A `ConsensusPartitionList-class` object.
#
setMethod(f = "nrow",
	signature = "ConsensusPartitionList",
	definition = function(x) {
	nrow(x@.env$data)
})

# == title
# Number of columns in the matrix
#
# == param
# -x A `ConsensusPartition-class` object.
#
setMethod(f = "ncol",
	signature = "ConsensusPartition",
	definition = function(x) {
	length(x@column_index)
})

# == title
# Number of columns in the matrix
#
# == param
# -x A `ConsensusPartitionList-class` object.
#
setMethod(f = "ncol",
	signature = "ConsensusPartitionList",
	definition = function(x) {
	ncol(x@.env$data)
})


# == title
# Row names of the matrix
#
# == param
# -x A `ConsensusPartition-class` object.
#
setMethod(f = "rownames",
	signature = "ConsensusPartition",
	definition = function(x) {
	rownames(x@.env$data)
})

# == title
# Row names of the matrix
#
# == param
# -x A `ConsensusPartitionList-class` object.
#
setMethod(f = "rownames",
	signature = "ConsensusPartitionList",
	definition = function(x) {
	rownames(x@.env$data)
})


# == title
# Column names of the matrix
#
# == param
# -x A `ConsensusPartition-class` object.
#
setMethod(f = "colnames",
	signature = "ConsensusPartition",
	definition = function(x) {
	colnames(x@.env$data[, x@column_index, drop = FALSE])
})


# == title
# Column names of the matrix
#
# == param
# -x A `ConsensusPartitionList-class` object.
#
setMethod(f = "colnames",
	signature = "ConsensusPartitionList",
	definition = function(x) {
	colnames(x@.env$data)
})

# == title
# Dimension of the matrix
#
# == param
# -x A `ConsensusPartition-class` object.
#
dim.ConsensusPartition = function(x) {
	c(nrow(x), ncol(x))
}

# == title
# Dimension of the matrix
#
# == param
# -x A `ConsensusPartitionList-class` object.
#
dim.ConsensusPartitionList = function(x) {
	dim(x@.env$data)
}

check_pkg = function(pkg, bioc = FALSE) {
	if(requireNamespace(pkg, quietly = TRUE)) {
		return(NULL)
	} else {

		if(!interactive()) {
			if(bioc) {
				stop_wrap(qq("You need to manually install package '@{pkg}' from Bioconductor."))
			} else {
				stop_wrap(qq("You need to manually install package '@{pkg}' from CRAN."))
			}
		}

		if(bioc) {
			answer = readline(qq("Package '@{pkg}' is required but not installed. Do you want to install it from Bioconductor? [y|n] "))
		} else {
			answer = readline(qq("Package '@{pkg}' is required but not installed. Do you want to install it from CRAN? [y|n] "))
		}

		if(bioc) {
			if(tolower(answer) %in% c("y", "yes")) {
				if(!requireNamespace("BiocManager", quietly = TRUE)) {
					install.packages("BiocManager")
				}
				BiocManager::install(pkg)
			} else {
				stop_wrap(qq("You need to manually install package '@{pkg}' from Bioconductor."))
			}
		} else {
			if(tolower(answer) %in% c("y", "yes")) {
				install.packages(pkg)
			} else {
				stop_wrap(qq("You need to manually install package '@{pkg}' from CRAN."))
			}
		}
	}
}

guess_best_km = function(mat, max_km = 15) {
	wss = NA
	max_km = min(c(nrow(mat) - 1, max_km))
	for (i in 2:max_km) {
		oe = try(fit <- kmeans(mat, centers = i, iter.max = 50), silent = TRUE)
		if(inherits(oe, "error")) {
			break
		}
		wss[i] = fit$tot.withinss
		if(is.na(wss[1])) wss[1] = fit$totss
	}
	k = 1:max_km
	k = k[seq_along(wss)]
	if(length(k) == 1) {
		return(1)
	} else if(length(k) == 2) {
		return(2)
	} else {
		min(elbow_finder(k, wss)[1], knee_finder(k, wss)[1])
	}
}






