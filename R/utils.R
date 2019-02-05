
if(!exists("strrep")) {
	strrep = function(x, times) {
		x = rep(x, times = times)
		return(paste0(x, collapse = ""))
	}
}

# change the label of `class` to let `class` fits `ref` more, maximize sum(class ==ref)
# class = c(rep("a", 10), rep("b", 3))
# ref = c(rep("b", 4), rep("a", 9))
# relabel_class(class, ref)
relabel_class = function(class, ref, full_set = union(class, ref)) {
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
	unmapped = setdiff(full_set, map)
	if(length(unmapped)) {
		map = c(map, structure(unmapped, names = unmapped))
	}
	    	
	df = data.frame(class = class, adjusted = map[class], ref = ref)
	attr(map, "df") = df
	return(map)
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
# -x a numeric vector.
# -q quantile to adjust.
# 
# == details
# Vaules larger than quantile ``1 - q`` are adjusted to the ``1 - q`` quantile and 
# values smaller than quantile ``q`` are adjusted to the ``q`` quantile.
#
# == value
# A numeric vector with same length as the original one.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# x = rnorm(10)
# x[1] = 100
# adjust_outlier(x)
adjust_outlier = function(x, q = 0.05) {
	qu = quantile(x, c(q, 1 - q), na.rm = TRUE)
	x[x < qu[1]] = qu[1]
	x[x > qu[2]] = qu[2]
	x
}

# == title
# Remove rows with low variance and impute missing data
#
# == param
# -m a numeric matrix.
# -sd_quantile cutoff the quantile of standard variation. Rows with variance less than it are removed.
# -max_na maximum NA rate in each row. Rows with NA rate larger than it are removed.
#
# == details
# The function uses `impute::impute.knn` to impute missing data, then
# uses `adjust_outlier` to adjust outliers and 
# removes rows with low standard variation.
#
# == value
# A numeric matrix.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# m = matrix(rnorm(200), 10)
# rownames(m) = letters[1:10]
# m[1, 1] = 1000
# range(m)
# m2 = adjust_matrix(m)
# range(m2)
adjust_matrix = function(m, sd_quantile = 0.05, max_na = 0.25) {
	if(is.data.frame(m)) {
		m = as.matrix(m)
	}

	l = apply(m, 1, function(x) sum(is.na(x))/length(x)) < max_na
	if(sum(!l)) {
		qqcat("removed @{sum(!l)} rows where more than @{round(max_na*100)}% of samples have NA values.\n")
	}
	m = m[l, , drop = FALSE]

	n_na = sum(is.na(m))
	if(n_na > 0) {
		qqcat("There are NA values in the data, now impute missing data.\n")
		m = impute.knn(m)$data
	}
	m = t(apply(m, 1, adjust_outlier))
	row_sd = rowSds(m, na.rm = TRUE)
	l = abs(row_sd) == 0
	m2 = m[!l, , drop = FALSE]
	row_sd = row_sd[!l]
	if(sum(l)) {
		qqcat("@{sum(l)} rows have been removed with zero variance.\n")
	}

	qa = quantile(row_sd, sd_quantile, na.rm = TRUE)
	l = row_sd >= qa
	if(sum(!l)) {
		qqcat("@{sum(!l)} rows have been removed with too low variance (sd < @{sd_quantile} quantile)\n")
	}
	m2[l, , drop = FALSE]
}

dev.off2 = ComplexHeatmap:::dev.off2



############### tracer
try_and_trace = function(expr) {
	oe = try(expr, silent = TRUE)
	if(inherits(oe, "try-error")) {
		cat("Caught an error. Following is the trace:\n")
		format_trace(getFromNamespace("data", ns = "tracer")$last_dump)
	}
}

format_trace <- function(data, style = NULL) {

  if (is.null(style)) {
    style <- getOption("tracer.style", getFromNamespace("tracer_default_style", ns = "tracer")())
  }

  if (length(data$calls) == 0) {

  } else {

    envs <- data$envs
    envs[envs == "R_GlobalEnv"] <- ""

    arg_cols <- (
      getOption("width", 80)
      - max(nchar(data$nums)) - 2
      - nchar(envs) - (nchar(envs) > 0) * 2
      - nchar(data$fnams) - 3
    )

    args <- ifelse(
      nchar(data$fargs) <= arg_cols,
      data$fargs,
      paste(
        substring(data$fargs, 1, arg_cols - nchar(clisymbols::symbol$ellipsis, "width")),
        crayon::cyan(clisymbols::symbol$ellipsis)
      )
    )

    locs <- paste0(
      "\n   at ", data$dirs, .Platform$file.sep,
      data$files, ":", data$lines, ":", data$cols
    )

    col_nums  <- style$num(data$nums)
    col_envs  <- vapply(envs, FUN.VALUE = "", getFromNamespace("style_env", ns = "tracer"), style = style)
    col_args  <- style$arg(args)
    col_fnams <- style$fnam(data$fnams)
    col_locs  <- style$location(
      ifelse(is.na(data$files) | data$files == "", "", locs)
    )

    str_calls <- paste(
      paste0(col_nums, " ", col_envs, col_fnams, " ", col_args),
      col_locs
    )

    cat("", str_calls, sep = "\n")

  }

  if (!is.null(data$error)) {
    cat("\n", style$error(data$error), sep = "")
  }

  invisible(NULL)
}

