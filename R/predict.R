# == title
# Predict classes for new samples based on cola classification
#
# == param
# -object A `ConsensusPartition-class` object.
# -k Number of subgroups to get the classifications.
# -mat The new matrix where the sample classes are going to be predicted. The number of rows should be
#      the same as the original matrix for cola analysis (also make sure the row orders are the same).
#      Be careful that the scaling of ``mat`` should be the same as that applied in cola analysis.
# -silhouette_cutoff Send to `get_signatures,ConsensusPartition-method` for determining signatures.
# -fdr_cutoff Send to `get_signatures,ConsensusPartition-method` for determining signatures.
# -group_diff Send to `get_signatures,ConsensusPartition-method` for determining signatures.
# -scale_rows Send to `get_signatures,ConsensusPartition-method` for determining signatures.
# -diff_method Send to `get_signatures,ConsensusPartition-method` for determining signatures.
# -dist_method Distance method. Value should be "euclidean", "correlation" or "cosine". Send to `predict_classes,matrix-method`.
# -nperm Number of permutatinos. It is used when ``dist_method`` is set to "euclidean" or "cosine". Send to `predict_classes,matrix-method`.
# -p_cutoff Cutoff for the p-values for determining class assignment. Send to `predict_classes,matrix-method`.
# -plot Whether to draw the plot that visualizes the process of prediction. Send to `predict_classes,matrix-method`.
# -force If the value is ``TRUE`` and when `get_signatures,ConsensusPartition-method` internally failed, top 1000 rows
#        with the highest between-group mean difference are used for constructing the signature centroid matrix.
#        It is basically used internally.
# -verbose Whether to print messages. Send to `predict_classes,matrix-method`.
# -help Whether to print help messages.
# -prefix Used internally.
#
# == details
# The prediction is based on the signature centroid matrix from cola classification. The processes are as follows:
#
# 1. For the provided `ConsensusPartition-class` object and a selected k, the signatures that discriminate classes
#    are extracted by `get_signatures,ConsensusPartition-method`. If number of signatures is more than 2000, only 2000 signatures are randomly sampled.
# 2. The signature centroid matrix is a k-column matrix where each column is the centroid of samples in the corresponding
#    class, i.e. the mean across samples. If rows were scaled in cola analysis, the signature centroid matrix is the mean of scaled
#    values and vise versa. Please note the samples with silhouette score less than ``silhouette_cutoff`` are removed
#    for calculating the centroids.
# 3. With the signature centroid matrix and the new matrix, it calls `predict_classes,matrix-method` to perform the prediction.
#    Please see more details of the prediction on that help page. Please note, the scales of the new matrix should be the same as the matrix
#    used for cola analysis.
#
# == value
# A data frame with two columns: the class labels (in numeric) and the corresponding p-values.
#
# == seealso
# `predict_classes,matrix-method` that predicts the classes for new samples.
#
# == example
# \donttest{
# data(golub_cola)
# res = golub_cola["ATC:skmeans"]
# mat = get_matrix(res)
# # note scaling should be applied here because the matrix was scaled in the cola analysis
# mat2 = t(scale(t(mat)))
# cl = predict_classes(res, k = 3, mat2)
# # compare the real classification and the predicted classification
# data.frame(cola_class = get_classes(res, k = 3)[, "class"],
#            predicted = cl[, "class"])
# # change to correlation method
# cl = predict_classes(res, k = 3, mat2, dist_method = "correlation")
# # compare the real classification and the predicted classification
# data.frame(cola_class = get_classes(res, k = 3)[, "class"],
#            predicted = cl[, "class"]) 
# }
setMethod(f = "predict_classes",
	signature = "ConsensusPartition",
	definition = function(object, k, mat, 
	silhouette_cutoff = 0.5, 
	fdr_cutoff = cola_opt$fdr_cutoff, 
	group_diff = cola_opt$group_diff, 
	scale_rows = object@scale_rows, 
	diff_method = "Ftest",
	dist_method = c("euclidean", "correlation", "cosine"), nperm = 1000,
	p_cutoff = 0.05, plot = TRUE, force = FALSE, 
	verbose = TRUE, help = TRUE, prefix = "") {

	if(help) {
		if(object@scale_rows) {
			if(object@partition_method %in% all_partition_methods()) {
				pm = get_partition_method(object@partition_method)
				scale_method = attr(pm, "scale_method")
			} else {
				scale_method = NULL
			}
			msg = strwrap(qq("The matrix has been scaled in cola analysis, thus the new matrix should also be scaled with the same method ('@{scale_method}'). Please double check."))
			
		} else {
			msg = strwrap(qq("The matrix is not scaled in cola analysis, thus the new matrix should not be scaled either. Please double check."))
		}

		msg = c(msg, "Set `help = FALSE` to suppress this message.")

		cat(paste(msg, collapse = "\n"), "\n\n")
	}

	if(!is.matrix(mat) && is.atomic(mat)) mat = matrix(mat, ncol = 1)

	if(nrow(object) != nrow(mat)) {
		stop_wrap("Number of rows in the new matrix should be the same as the matrix for cola analysis (also the row order).")
	}

	tb = get_signatures(object, k = k, plot = FALSE, silhouette_cutoff = silhouette_cutoff, 
		fdr_cutoff = fdr_cutoff, group_diff = group_diff, scale_rows = scale_rows, 
		diff_method = diff_method, verbose = verbose, prefix = prefix)

	if(nrow(tb) < 20) {
		if(force) {
			hash = attr(tb, "hash")
			data = object@.env$data
			if(object@scale_rows) {
				data = t(scale(t(data)))
			}

			class_df = get_classes(object, k = k)
			sig_mat = do.call(cbind, tapply(1:nrow(class_df), class_df$class, function(ind) {
				rowMeans(data[, ind, drop = FALSE])
			}))

			# if no hash attribute, randomly select one sample from each group
			if(is.null(hash)) {
				ind = order(apply(sig_mat, 1, function(x) max(x) - min(x)), decreasing = TRUE)[1:min(1000, nrow(data))]
				sig_mat = sig_mat[ind, , drop = FALSE]
				if(verbose) qqcat("@{prefix}* simply take top @{min(1000, nrow(data))} rows with the highest row range.\n")
			} else {
				# if there is hash attached, adjust fdr_cutoff
				hash_nm = paste0("signature_fdr_", hash)
				fdr = object@.env[[hash_nm]]$fdr
				ind = order(-fdr, apply(sig_mat, 1, function(x) max(x) - min(x)), decreasing = TRUE)[1:min(1000, nrow(data))]
				sig_mat = sig_mat[ind, , drop = FALSE]
				if(verbose) qqcat("@{prefix}* simply take top @{min(1000, nrow(data))} rows with the most significant FDRs.\n")
			}

		} else {
			stop_wrap("Number of signatures is too small.")
		}
	}
	
	if(object@scale_rows) {
		sig_mat = tb[, grepl("^scaled_mean_\\d+$", colnames(tb))]
	} else {
		sig_mat = tb[, grepl("^mean_\\d+$", colnames(tb))]
	}

	sig_mat = as.matrix(sig_mat)
	colnames(sig_mat) = NULL

	mat = mat[tb$which_row, , drop = FALSE]

	if(nrow(mat) > 1000) {
		ind = sample(nrow(mat), 1000)
		mat = mat[ind, , drop = FALSE]
		sig_mat = sig_mat[ind, , drop = FALSE]
	}
	
	predict_classes(sig_mat, mat, nperm = nperm, dist_method = dist_method, p_cutoff = p_cutoff, 
		plot = plot, verbose = verbose, prefix = prefix)
})


# == title
# Predict classes for new samples based on signature centroid matrix
#
# == param
# -object The signature centroid matrix. See the Details section.
# -mat The new matrix where the classes are going to be predicted. The number of rows should be
#      the same as the signature centroid matrix (also make sure the row orders are the same).
#      Be careful that ``mat`` should be in the same scale as the centroid matrix.
# -dist_method Distance method. Value should be "euclidean", "correlation" or "cosine".
# -nperm Number of permutatinos. It is used when ``dist_method`` is set to "euclidean" or "cosine".
# -p_cutoff Cutoff for the p-values for determining class assignment.
# -plot Whether to draw the plot that visualizes the process of prediction.
# -verbose Whether to print messages.
# -prefix Used internally.
#
# == details
# The signature centroid matrix is a k-column matrix where each column is the centroid of samples 
# in the corresponding class (k-group classification).
# 
# For each sample in the new matrix, the task is basically to test which signature centroid the 
# current sample is the closest to. There are two methods: the Euclidean distance and the 
# correlation (Spearman) distance.
#
# For the Euclidean/cosine distance method, for the vector denoted as x which corresponds to sample i 
# in the new matrix, to test which class should be assigned to sample i, the distance between 
# sample i and all k signature centroids are calculated and denoted as d_1, d_2, ..., d_k. The class with the smallest distance is assigned to sample i.
# The distances for k centroids are sorted increasingly, and we design a statistic named "difference ratio", denoted as r
# and calculated as: (|d_(1) - d_(2)|)/mean(d), which is the difference between the smallest distance and the second
# smallest distance, normalized by the mean distance. 
# To test the statistical significance of r, we randomly permute rows of the signature centroid matrix and calculate r_rand. 
# The random permutation is performed ``n_perm`` times and the p-value is calculated as the proportion of r_rand being
# larger than r.
#
# For the correlation method, the distance is calculated as the Spearman correlation between sample i and signature
# centroid k. The label for the class with the maximal correlation value is assigned to sample i. The 
# p-value is simply calculated by `stats::cor.test` between sample i and centroid k.
#
# If a sample is tested with a p-value higher than ``p_cutoff``, the corresponding class label is set to ``NA``.
#
# == value
# A data frame with two columns: the class labels (the column names of the signature centroid matrix are treated as class labels) and the corresponding p-values.
#
# == example
# \donttest{
# data(golub_cola)
# res = golub_cola["ATC:skmeans"]
# mat = get_matrix(res)
# # note scaling should be applied here because the matrix was scaled in the cola analysis
# mat2 = t(scale(t(mat)))
#
# tb = get_signatures(res, k = 3, plot = FALSE)
# sig_mat = tb[, grepl("scaled_mean", colnames(tb))]
# sig_mat = as.matrix(sig_mat)
# colnames(sig_mat) = paste0("class", seq_len(ncol(sig_mat)))
# # this is how the signature centroid matrix looks like:
# head(sig_mat)
#
# mat2 = mat2[tb$which_row, , drop = FALSE]
# 
# # now we predict the class for `mat2` based on `sig_mat`
# predict_classes(sig_mat, mat2)
# }
setMethod(f = "predict_classes",
	signature = "matrix",
	definition = function(object, mat, dist_method = c("euclidean", "correlation", "cosine"), 
	nperm = 1000, p_cutoff = 0.05, plot = TRUE, verbose = TRUE, prefix = "") {

	sig_mat = object

	if(nrow(mat) != nrow(sig_mat)) {
		stop_wrap("nrow of the matrix and the signature matrix should be the same.")
	}

	dist_method = match.arg(dist_method)[1]
	n_sig = ncol(sig_mat)

	if(verbose) qqcat("@{prefix}* Predict classes based on @{ncol(sig_mat)}-group classification (@{dist_method} method) on a @{ncol(mat)}-column matrix.\n")

	if(dist_method %in% c("euclidean", "cosine")) {

		if(dist_method == "euclidean") {
			dm = as.integer(1)
		} else if(dist_method == "cosine") {
			dm = as.integer(2)
		}

		dist_to_signatures = as.matrix(pdist(t(mat), t(sig_mat), dm))
		diff_ratio = apply(dist_to_signatures, 1, function(x) { 
			x = sort(x)
			abs(x[1] - x[2])/mean(x)
		})

# 		diff_ratio_r = NULL
# 		counter = set_counter(nperm, fmt = qq("@{prefix}Permute rows of the signature centroid matrix, run %s..."))
# 		for(i in 1:nperm) {
# 			dist_to_signatures_r = as.matrix(pdist(t(mat), t(sig_mat[sample(nrow(sig_mat)), , drop = FALSE]), dm))
# 			diff_ratio_r = cbind(diff_ratio_r, apply(dist_to_signatures_r, 1, function(x) { 
# 				x = sort(x)
# 				abs(x[1] - x[2])/mean(x)
# 			}))
# 			if(verbose) counter()
# 		}

		diff_ratio_r = cal_diff_ratio_r(t(mat), t(sig_mat), nperm, dm)

		p = rowSums(diff_ratio_r - diff_ratio > 0)/nperm

		predicted_class = apply(dist_to_signatures, 1, which.min)
	} else {
		dist_to_signatures = as.matrix(cor(mat, sig_mat, method = "spearman"))
		predicted_class = apply(dist_to_signatures, 1, which.max)
		p = sapply(seq_along(predicted_class), function(i) {
			cor.test(mat[, i], sig_mat[, predicted_class[i]])$p.value
		})
	}

	if(!is.null(colnames(sig_mat))) {
		predicted_class = colnames(sig_mat)[predicted_class]
		predicted_class2 = predicted_class
		level = colnames(sig_mat)
	} else {
		predicted_class2 = paste0("class", predicted_class)
		level = paste0("class", seq_len(ncol(sig_mat)))
	}

	predicted_class2[p > p_cutoff] = "unclear"

	if(plot) {
		predicted_class2 = factor(predicted_class2, levels = c(level, "unclear"))
		predicted_col = structure(1:ncol(sig_mat)+1, names = level)
		
		group_col = structure(seq_len(n_sig) + 1, names = level)
		group_col = c(group_col, c("unclear" = "grey"))
		row_split = level[apply(sig_mat, 1, which.max)]
		row_split = factor(row_split, levels = level)

		if(dist_method != "correlation") {
			dist_col_fun = colorRamp2(range(dist_to_signatures), c("white", "purple"))
		} else {
			mabs = max(abs(dist_to_signatures))
			dist_col_fun = colorRamp2(c(-mabs, 0, mabs), c("green", "white", "red"))
		}
		ha = HeatmapAnnotation(
				"Dist to centroid" = dist_to_signatures,
				"Predicted classes" = predicted_class2, 
				col = list(
					"Dist to centroid" = dist_col_fun,
					"Predicted classes" = group_col),
				show_annotation_name = TRUE,
				simple_anno_size = unit(4, "mm"))
		if(dist_method == "correlation") {
			names(ha@anno_list)[1] = "Corr to centroid"
			ha@anno_list[[1]]@color_mapping@name = "Corr to centroid"
		}

		wss = (nrow(mat)-1)*sum(apply(mat,1,var))
		max_km = min(c(nrow(mat) - 1, 15))
		# if(verbose) qqcat("* apply k-means on rows with 2~@{max_km} clusters.\n")
		for (i in 2:max_km) {
			# if(verbose) qqcat("  - applying k-means with @{i} clusters.\n")
			wss[i] = sum(kmeans(mat, centers = i, iter.max = 50)$withinss)
		}
		row_km = min(elbow_finder(1:max_km, wss)[1], knee_finder(1:max_km, wss)[1])
			
		ht_list = Heatmap(mat, name = "New matrix",
			top_annotation = ha, row_km = row_km,
			# row_split = row_split, cluster_row_slices = FALSE,
			show_column_names = FALSE, row_title = NULL,
			cluster_columns = TRUE, cluster_column_slices = FALSE, show_column_dend = FALSE,
			column_split = predicted_class2,
			show_row_dend = FALSE,
			column_title = qq("Based on @{nrow(sig_mat)} signatures")
		) + Heatmap(sig_mat, cluster_columns = FALSE, width = unit(4*ncol(sig_mat), "mm"),
			heatmap_legend_param = list(title = "Signature centroid"))

		draw(ht_list, merge_legend = TRUE)
	}

	df = data.frame(class = predicted_class, p = p)
	rownames(df) = colnames(mat)
	return(df)
})


set_counter = function(n, fmt = "%s") {

	n = as.integer(n)
	i = 1

	f = function() {
		if(interactive()) {
			pct = round(i/n*100, 1)
			str = paste0(i, "/", n, " (", pct, "%)")
			str = sprintf(fmt, str)

			cat(strrep("\r", nchar(str)))
			cat(str)
			if(i == n) cat("\n")

			i = i + 1
			assign("i", i, envir = parent.env(environment()))
			return(invisible(i))
		}
	}
}
