
# == title
# Co-occurrence of genes in gene sets
#
# == param
# -x the object returned from `get_signatures`
# -map mapping between rows of ``x$mat`` and genes in ``genesets``
# -genesets a object constructed from `msigdb_catalogue`
# -min_count minimal number of genes in genesets
# -max_count maximal number of genes in genesets
#
# == details
# For genes in each row group, the co-occurence of every gene pair to be in a same gene set
# is calculated. The mean co-occurence of all genes is used as the final statistic which can
# be understanded as the mean number of gene sets that two genes co-exist.
#
# == value
# The mean co-occurrence in each subgroup.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
gene_co_occurrence = function(x, genesets, map = NULL, min_count = 50, max_count = 5000) {

	genesets = genesets$list

	gl = sapply(genesets, length)
	l = gl >= min_count & gl <= max_count
	if(sum(l) == 0) {
		stop("no genesets left with cutoff of `min_count` and `max_count`.")
	}
	qqcat("@{sum(l)}/@{length(l)} gene sets used\n")
	genesets = genesets[l]

	match_mat = matrix(FALSE, nrow = nrow(x$mat), ncol = length(genesets))

	g = rownames(x$mat)
	if(!is.null(map)) {
		g2 = map[g]
		l = is.na(g2)
		g2[l] = g[l]
		g = g2
	}
	rownames(match_mat) = g

	for(i in seq_along(genesets)) {
		match_mat[intersect(genesets[[i]], g), i] = TRUE
	}

	if(sum(match_mat) == 0) {
		stop("gene names in `x` have no overlap to `genesets`.")
	}

	unique_group = sort(unique(x$group))

	cooccurrence = NULL
	for(i in seq_along(unique_group)) {
		ind = which(x$group == unique_group[i])
		qqcat("gene co-occurrence in group @{unique_group[i]}, @{length(ind)} rows\n")
		submat = gene_cooccurrence_in_geneset(match_mat[ind, ])
		cooccurrence[i] = mean(submat[lower.tri(submat)])
	}
	names(cooccurrence) = unique_group

	return(cooccurrence)
}


# == title
# Enrich signature genes to genesets
#
# == param
# -x the object returned from `get_signatures`
# -map mapping between rows of ``x$mat`` and genes in ``genesets``
# -bg background gene list
# -genesets a object constructed from `msigdb_catalogue`
# -min_count minimal number of genes in genesets
# -max_count maximal number of genes in genesets
# -fdr_cutoff1 cutoff of FDR for the geneset to be significantly enriched
# -fdr_cutoff2 cutoff of RDR for the geneset to be not enriched
#
# == details
# The function tries to find significantly enriched genesets which at the same time
# are also subgroup specific.
#
# == value
# A list with significant genesets in each subgroup
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
enrich_signatures_to_genesets = function(x, genesets, map = NULL, bg, min_count = 50, max_count = 5000,
	fdr_cutoff1 = 0.05, fdr_cutoff2 = 0.5) {
	
	genesets0 = genesets

	genesets = lapply(genesets0$list, function(x) intersect(x, bg))
	gl = sapply(genesets, length)
	if(all(gl == 0)) {
		stop("gene names in `bg` have no overlap to `genesets`.")
	}
	l = gl >= min_count & gl <= max_count
	if(sum(l) == 0) {
		stop("no genesets left with cutoff of `min_count` and `max_count`.")
	}
	qqcat("@{sum(l)}/@{length(l)} gene sets used\n")
	genesets = genesets[l]
	n_genesets = length(genesets)

	match_mat = matrix(0, nrow = nrow(x$mat), ncol = length(genesets))

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
		match_mat[intersect(genesets[[i]], g), i] = 1
	}

	unique_group = sort(unique(x$group))
	n_groups = length(unique_group)
	p_mat = matrix(nrow = length(genesets), ncol = n_groups)
	stat_list = lapply(1:n_groups, function(i) {
		data.frame(geneset = names(genesets),
			       gene_in_set = numeric(n_genesets),
			       geneset_size = numeric(n_genesets),
			       row_group_size = numeric(n_genesets),
			       p_value = numeric(n_genesets),
			       fdr = numeric(n_genesets),
			       stringsAsFactors = FALSE)
	})
	names(stat_list) = sort(unique(x$group))

	colnames(p_mat) = unique_group
	for(gi in unique(x$group)) {
		gp = unique(g[x$group == gi])
		for(i in seq_along(genesets)) {
			v1 = length(intersect(gp, genesets[[i]]))
			v2 = length(gp) - v1
			v3 = length(genesets[[i]]) - v1
			v4 = length(bg) - v1 - v2 - v3
			p_mat[i, gi] = fisher.test(matrix(c(v1, v2, v3, v4), nrow = 2), alternative = "greater")$p.value
			stat_list[[gi]][i, "gene_in_set"] = v1
			stat_list[[gi]][i, "geneset_size"] = v3 + v1
			stat_list[[gi]][i, "row_group_size"] = v2 + v1
			stat_list[[gi]][i, "p_value"] = p_mat[i, gi]
		}
	}

	fdr_mat = p.adjust(p_mat, "BH")
	dim(fdr_mat) = dim(p_mat)
	colnames(fdr_mat) = unique_group

	for(i in seq_along(stat_list)) {
		stat_list[[i]]$fdr = fdr_mat[, i]
	}
	l = apply(fdr_mat, 1, function(x) {
		ind = x < fdr_cutoff1
		if(sum(ind) != 1) {
			return(FALSE)
		} else {
			all(x[!ind] > fdr_cutoff2)
		}
	})

	if(sum(l)) {
		match_mat2 = match_mat[, l, drop = FALSE]
		fdr_mat2 = fdr_mat[l, , drop = FALSE]
		min_fdr= rowMins(fdr_mat2)
		stat_list = lapply(stat_list, function(df) {
			df[l, , drop = FALSE]
			df[order(df$fdr), , drop = FALSE]
		})

		l_hit = rowSums(match_mat2) > 0
		ht_list = Heatmap(x$group[l_hit], name = "group", show_row_names = FALSE, width = unit(5, "mm"), col = x$class_col)
		
		column_anno = as.character(apply(fdr_mat2, 1, function(x) which(x < fdr_cutoff1)))
		
		fdr_col_fun = replicate(ncol(fdr_mat2), colorRamp2(c(0, -log10(fdr_cutoff1), -log10(fdr_cutoff1)*2), c("green", "white", "red")))
		names(fdr_col_fun) = colnames(fdr_mat2)
		ht_list = ht_list + Heatmap(match_mat2[l_hit, , drop = FALSE], name = "in geneset", col = c("1" = "purple", "0" = "#FFFFFFFF"), 
			top_annotation = HeatmapAnnotation(group = column_anno, 
				df = -log10(fdr_mat2),
				col = c(list(group = x$class_col), fdr_col_fun), show_legend = TRUE, show_annotation_name = TRUE),
			show_row_names = FALSE, show_row_dend = FALSE, split = x$group[l_hit], combined_name_fun = NULL,
			cluster_columns = FALSE, column_order = order(column_anno, min_fdr), show_column_dend = TRUE, clustering_method_rows = "ward.D",
			show_column_names = FALSE, column_names_gp = gpar(fontsize = 8))
		draw(ht_list, main_heatmap = "in geneset", column_title = qq("@{sum(l)}/@{length(l)} genesets which are subgroup specific"))

		for(i in 1:n_groups) {
			decorate_heatmap_body("in geneset", slice = i, {
				grid.rect(gp = gpar(fill = "transparent", col = "black"))
			})
		}

		sig_geneset_list = list()
		for(ug in sort(unique(column_anno))) {
			sig_geneset_list[[ug]] = stat_list[[ug]][column_anno == ug, , drop = FALSE]
		}
		attr(sig_geneset_list, "genesets") = genesets0
		return(invisible(sig_geneset_list))
	} else {
		cat("no significant geneset has enrichment for row groups.\n")
		return(invisible(NULL))
	}
}

# == title
# Word cloud for the enriched functions
#
# == param
# -x the object returned from `enrich_signatures_to_genesets`
# -stopwords words that are not take into account for the word cloud
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
enriched_functions_word_cloud = function(x, stopwords = GS_STOPWORDS) {

	genesets = attr(x, "genesets")

	for(i in seq_along(x)) {
		df = x[[i]]
		gs = df$geneset
		gs_desc = genesets$meta[which(genesets$meta$id %in% gs), "desc"]

		docs = Corpus(VectorSource(gs_desc))

		docs = tm_map(docs, content_transformer(tolower))
		docs = tm_map(docs, removeNumbers)
		docs = tm_map(docs, removeWords, c(stopwords("SMART"), stopwords("english")))
		docs = tm_map(docs, removeWords, GS_STOPWORDS)
		docs = tm_map(docs, removePunctuation)
		docs = tm_map(docs, stripWhitespace)

		tdm.bigram = TermDocumentMatrix(docs)

		freq = sort(rowSums(as.matrix(tdm.bigram)), decreasing = TRUE)
		freq.df = data.frame(word = names(freq), freq = freq)

		if(dev.interactive() && interactive()) {
			readline(prompt = "press enter to load next plot: ")
		}

		wordcloud(words = freq.df$word, freq = freq.df$freq, min.freq = 2,
		          max.words = 200, random.order = FALSE, rot.per = 0.35, 
		          colors = brewer.pal(8, "Dark2"))
		
		text(0.5, 1, qq("group @{i}, @{nrow(df)} genesets"), adj = c(0.5, 0))
	}
}

GS_STOPWORDS = c("gene", "genes", "geneid", "regulated", "cells", "cell", "involved", "compared")

# == title
# Load from MSigDB
#
# == param
# -f path of the xml file.
#
# == details
# The xml file can be downloaded from http://software.broadinstitute.org/gsea/downloads.jsp .
#
# == value
# A ``msigdb`` class object with two elements:
#
# -``meta`` a data frame with geneset ids, organisms, description of the gene sets and categories.
# -``list`` a list of gene sets of gene symbols.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
load_msigdb = function(f) {
	x = read_xml(f)
	genesets = xml_children(x)
	geneset_id = xml_attr(genesets, "STANDARD_NAME")
	geneset_organism = xml_attr(genesets, "ORGANISM")
	geneset_desc = xml_attr(genesets, "DESCRIPTION_BRIEF")
	geneset_category = xml_attr(genesets, "CATEGORY_CODE")
	geneset_subcategory = xml_attr(genesets, "SUB_CATEGORY_CODE")
	geneset_member = xml_attr(genesets, "MEMBERS_SYMBOLIZED")

	geneset_meta = data.frame(id = geneset_id,
		organism = geneset_organism,
		desc = geneset_desc,
		category = geneset_category,
		sub_category = geneset_subcategory
	)
	geneset_list = strsplit(geneset_member, ",")
	names(geneset_list) = geneset_id

	l = geneset_meta$category == "ARCHIVED" | geneset_meta$organism %in% c("Danio rerio", "Macaca mulatta", "Human", "Mouse")
	geneset_meta = geneset_meta[!l, ]
	geneset_list = geneset_list[!l]

	env = new.env()
	env$meta = geneset_meta
	env$list = geneset_list
	class(env) = c("msigdb", "environment")
	return(env)
}

# == title
# Print the msigdb class object
#
# == param
# -x the ``msigdb`` object
# -... other arguments
#
# == value
# No value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
print.msigdb = function(x, ...) {
	all_organisms = sort(unique(x$meta$organism))
	for(organism in all_organisms) {
		l1 = x$meta$organism == organism
		qqcat("@{organism}: @{sum(l1)} gene sets\n")
		all_categories = sort(unique(x$meta$category[l1]))
		for(category in all_categories) {
			l2 = l1 & x$meta$category == category
			all_sub_categories = sort(unique(x$meta$sub_category[l2]))
			qqcat("  @{category}: @{sum(l2)} gene sets\n")
			if(!identical(all_sub_categories, "")) {
				for(sub_category in all_sub_categories) {
					l3 = l2 & x$meta$sub_category == sub_category
					qqcat("    @{sub_category}: @{sum(l3)} gene sets\n")
				}
			}
		}
	}
}

# == title
# Get a catelogue from the whole msigdb database
#
# == param
# -x the msigdb object from `load_msigdb`
# -category category
# -sub_category sub category
# -organism organism
#
# == value
# a ``msigdb`` object
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
msigdb_catalogue = function(x, category = "H", sub_category, organism = "Homo sapiens") {
	if(missing(sub_category)) {
		l = x$meta$organism == organism & x$meta$category == category
	} else {
		l = x$meta$organism == organism & x$meta$category == category & x$meta$sub_category == sub_category
	}

	if(sum(l) == 0) {
		stop("No gene set found.")
	}

	x2 = new.env()
	x2$meta = x$meta[l, , drop = FALSE]
	x2$list = x$list[l]

	class(x2) = c("msigdb", "environment")
	return(x2)
}

# == title
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
setMethod(f = "signature_density",
	signature = "ConsensusPartition",
	definition = function(object, k, ...) {

	cl = get_class(object, k = k)$class
	data = object@.env$data[, object@.env$column_index, drop = FALSE]

	all_den_list = lapply(seq_len(ncol(data)), function(i) density(data[, i]))
	x_range = range(unlist(lapply(all_den_list, function(x) x$x)))
	y_range = range(unlist(lapply(all_den_list, function(x) x$y)))

	op = par(no.readonly = TRUE)
	par(mfrow = c(k + 1, 1), mar = c(3, 4, 4, 1))
	plot(NULL, type = "n", xlim = x_range, ylim = y_range, main = "all rows", ylab = "density", xlab = NULL)
	for(i in 1:ncol(data)) {
		lines(all_den_list[[i]], col = brewer_pal_set2_col[cl[i]], lwd = 1)
	}

	x = get_signatures(object, k = k, plot = FALSE, ...)
	gp = x$group
	for(j in 1:k) {
		gp2 = gp[gp == as.character(j)]
		all_den_list = lapply(seq_len(ncol(data)), function(i) density(data[names(gp2), i]))
		# x_range = range(unlist(lapply(all_den_list, function(x) x$x)))
		y_range = range(unlist(lapply(all_den_list, function(x) x$y)))

		plot(NULL, type = "n", xlim = x_range, ylim = y_range, main = qq("signatures in subgroup @{j}/@{k}"), ylab = "density", xlab = NULL)
		for(i in 1:ncol(data)) {
			lines(all_den_list[[i]], col = brewer_pal_set2_col[cl[i]], lwd = ifelse(cl[i] == j, 2, 0.5))
		}
	}
	par(op)
})
