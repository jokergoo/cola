
# == title
# Compare two partitionings
#
# == param
# -object A `ConsensusPartition` object.
# -object2  A `ConsensusPartition` object.
# -output_file The path of the output HTML file. If it is not specified, the report will be opened in the web browser.
# -k1 Number of subgroups in ``object``.
# -k2 Number of subgroups in ``object2``.
# -dimension_reduction_method Which dimension reduction method to use.
# -id_mapping Pass to `functional_enrichment,ConsensusPartition-method`.
# -row_km1 Number of k-means groups, see Details.
# -row_km2 Number of k-means groups, see Details.
# -row_km3 Number of k-means groups, see Details.
#
# == details
# The function produces a HTML report which includes comparisons between two partitioning results.
#
# In the report, there are three heatmaps which visualize A) the signature genes specific in the first partition, B) the signature genes
# both in the two partitionings and C) the signatures genes specific in the second partition. Argument ``row_km1``, ``row_km2`` and 
# ``row_km3`` control how many k-means groups should be applied on the three heatmaps.
#
# == example
# \dontrun{
# data(golub_cola)
# require(hu6800.db)
# x = hu6800ENTREZID
# mapped_probes = mappedkeys(x)
# id_mapping = unlist(as.list(x[mapped_probes]))
# compare_partitions(golub_cola["ATC:skmeans"], golub_cola["SD:kmeans"], 
#     id_mapping = id_mapping)
# }
setMethod(f = "compare_partitions",
	signature = "ConsensusPartition",
	definition = function(object, object2, output_file, k1 = 2, k2 = 2, 
	dimension_reduction_method = "UMAP",
	id_mapping = guess_id_mapping(rownames(object), "org.Hs.eg.db", FALSE),
	row_km1 = ifelse(k1 == 2, 2, 1),
	row_km2 = ifelse(k1 ==2 && k2 == 2, 2, 1),
	row_km3 = ifelse(k2 == 2, 2, 1)) {

	check_pkg("cowplot")
	check_pkg("rmarkdown")

	if(identical(topenv(), .GlobalEnv)) {
		template_file = "~/project/cola/inst/extdata/compare_partitions.Rmd"
	} else {
		template_file = system.file("extdata", "compare_partitions.Rmd", package = "cola")
	}

	res1 = object
	res2 = object2

	k1 = k1
	k2 = k2
	dimension_reduction_method = dimension_reduction_method
	id_mapping = id_mapping
	row_km1 = row_km1
	row_km2 = row_km2
	row_km3 = row_km3

	if(!inherits(res2, "ConsensusPartition")) {
		stop_wrap("The first two arguments should all be `ConsensusPartition` objects.")
	}

	method1 = paste0(res1@top_value_method, ":", res1@partition_method)
	method2 = paste0(res2@top_value_method, ":", res2@partition_method)

	if(method1 == method2) {
		do_top_row = FALSE
	} else {
		do_top_row = TRUE
	}

	feature = "row"
	if(is.null(rownames(object))) {
		do_go_enrichment = FALSE
	} else {
		if(is.function(id_mapping)) {
			mapped = id_mapping(rownames(res1))
			mapped = mapped[!is.na(mapped)]
			if(length(mapped) == 0) {
				do_go_enrichment = FALSE
			} else {
				do_go_enrichment = TRUE
			}
		} else if(is.atomic(id_mapping)) {
			if(is.null(names(id_mapping))) {
				do_go_enrichment = FALSE
			} else {
				mapped = id_mapping[rownames(res1)]
				mapped = mapped[!is.na(mapped)]
				if(length(mapped) == 0) {
					do_go_enrichment = FALSE
				} else {
					do_go_enrichment = TRUE
				}
			}
		}
	}

	if(do_go_enrichment) feature = "gene"

	tmpdir = tempdir()
	rmd_file = tempfile(tmpdir = tmpdir, fileext = ".Rmd")
	brew(template_file, output = rmd_file)

	if(missing(output_file)) {
		output_file = tempfile(tmpdir = tmpdir, fileext = ".html")
		rmarkdown::render(rmd_file, output_file = output_file, output_dir = tmpdir, quiet = TRUE)
		browseURL(output_file)
	} else {
		rmarkdown::render(rmd_file, output_file = output_file, output_dir = dirname(output_file), quiet = TRUE)
	}

	invisible(NULL)
})
	