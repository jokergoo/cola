 


DAVID_ALL_ID_TYPES = "AFFYMETRIX_3PRIME_IVT_ID,AFFYMETRIX_EXON_GENE_ID,AFFYMETRIX_SNP_ID,AGILENT_CHIP_ID,AGILENT_ID,AGILENT_OLIGO_ID,ENSEMBL_GENE_ID,ENSEMBL_TRANSCRIPT_ID,ENTREZ_GENE_ID,FLYBASE_GENE_ID,FLYBASE_TRANSCRIPT_ID,GENBANK_ACCESSION,GENOMIC_GI_ACCESSION,GENPEPT_ACCESSION,ILLUMINA_ID,IPI_ID,MGI_ID,PFAM_ID,PIR_ID,PROTEIN_GI_ACCESSION,REFSEQ_GENOMIC,REFSEQ_MRNA,REFSEQ_PROTEIN,REFSEQ_RNA,RGD_ID,SGD_ID,TAIR_ID,UCSC_GENE_ID,UNIGENE,UNIPROT_ACCESSION,UNIPROT_ID,UNIREF100_ID,WORMBASE_GENE_ID,WORMPEP_ID,ZFIN_ID"
DAVID_ALL_ID_TYPES = strsplit(DAVID_ALL_ID_TYPES, ",")[[1]]
DAVID_ALL_CATALOGS = "BBID,BIND,BIOCARTA,BLOCKS,CGAP_EST_QUARTILE,CGAP_SAGE_QUARTILE,CHROMOSOME,COG_NAME,COG_ONTOLOGY,CYTOBAND,DIP,EC_NUMBER,ENSEMBL_GENE_ID,ENTREZ_GENE_ID,ENTREZ_GENE_SUMMARY,GENETIC_ASSOCIATION_DB_DISEASE,GENERIF_SUMMARY,GNF_U133A_QUARTILE,GENETIC_ASSOCIATION_DB_DISEASE_CLASS,GOTERM_BP_2,GOTERM_BP_1,GOTERM_BP_4,GOTERM_BP_3,GOTERM_BP_FAT,GOTERM_BP_5,GOTERM_CC_1,GOTERM_BP_ALL,GOTERM_CC_3,GOTERM_CC_2,GOTERM_CC_5,GOTERM_CC_4,GOTERM_MF_1,GOTERM_MF_2,GOTERM_CC_FAT,GOTERM_CC_ALL,GOTERM_MF_5,GOTERM_MF_FAT,GOTERM_MF_3,GOTERM_MF_4,HIV_INTERACTION_CATEGORY,HOMOLOGOUS_GENE,GOTERM_MF_ALL,HIV_INTERACTION,MINT,NCICB_CAPATHWAY_INTERACTION,INTERPRO,KEGG_PATHWAY,PANTHER_FAMILY,PANTHER_BP_ALL,OMIM_DISEASE,OFFICIAL_GENE_SYMBOL,PANTHER_SUBFAMILY,PANTHER_PATHWAY,PANTHER_MF_ALL,PIR_SUMMARY,PIR_SEQ_FEATURE,PFAM,PRODOM,PRINTS,PIR_TISSUE_SPECIFICITY,PIR_SUPERFAMILY,SMART,SP_COMMENT,SP_COMMENT_TYPE,SP_PIR_KEYWORDS,PROSITE,PUBMED_ID,REACTOME_INTERACTION,REACTOME_PATHWAY,UNIGENE_EST_QUARTILE,UP_SEQ_FEATURE,UP_TISSUE,ZFIN_ANATOMY,SSF,TIGRFAMS,UCSC_TFBS"
DAVID_ALL_CATALOGS = strsplit(DAVID_ALL_CATALOGS, ",")[[1]]

# == title
# Perform DAVID analysis
#
# == param
# -genes A vector of gene identifiers.
# -email The email that user registered on DAVID web service (https://david.ncifcrf.gov/content.jsp?file=WS.html ).
# -catalog A vector of function catalogs. Valid values should be in ``cola:::DAVID_ALL_CATALOGS``.
# -idtype ID types for the input gene list. Valid values should be in ``cola:::DAVID_ALL_ID_TYPES``.
# -species Full species name if the ID type is not uniquely mapped to one single species.
#
# == details
# This function directly sends the HTTP request to DAVID web service (https://david.ncifcrf.gov/content.jsp?file=WS.html )
# and parses the returned XML. The reason of writing this function is I have problems with other
# R packages doing DAVID analysis (e.g. RDAVIDWebService, https://bioconductor.org/packages/devel/bioc/html/RDAVIDWebService.html )
# because the rJava package RDAVIDWebService depends on can not be installed on my machine.
#
# Users are encouraged to use more advanced
# gene set enrichment tools such as clusterProfiler (http://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html ), 
# or fgsea (http://www.bioconductor.org/packages/release/bioc/html/fgsea.html ).
#
# If you want to run this function multiple times, please set time intervals between runs.
# 
# == value
# A data frame with functional enrichment results.
#
# == seealso
# https://david.ncifcrf.gov
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
submit_to_david = function(genes, email, 
	catalog = c("GOTERM_CC_FAT", "GOTERM_BP_FAT", "GOTERM_MF_FAT", "KEGG_PATHWAY"),
	idtype = "ENSEMBL_GENE_ID", species = "Homo sapiens") {

	if(missing(email)) {
		stop_wrap("You need to register to DAVID web service.")
	}

	if(!idtype %in% DAVID_ALL_ID_TYPES) {
		cat("idtype is wrong, it should be in:\n")
		print(DAVID_ALL_ID_TYPES)
		stop_wrap("You have an error.")
	}
	if(!all(catalog %in% DAVID_ALL_CATALOGS)) {
		cat("catalog is wrong, it should be in:\n")
		print(DAVID_ALL_CATALOGS)
		stop_wrap("You have an error.")
	}

	if(grepl("ENSEMBL", idtype)) {
		map = structure(genes, names = gsub("\\.\\d+$", "", genes))
		genes = names(map)
	}

	message_wrap(qq("Idtype: @{idtype}"))
	message_wrap(qq("Catalog: @{paste(catalog, collapse = ' ')}"))
	
	DAVID_DWS = "https://david-d.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint"

	# login
	message_wrap(qq("log to DAVID web service with @{email}"))
	response = GET(qq("@{DAVID_DWS}/authenticate"),
		query = list("args0" = email)
	)
	if(xml_text(httr::content(response)) != "true") {
		stop_wrap(qq("@{email} has not been registed."))
	}

	# add gene list
	message_wrap(qq("add gene list (@{length(genes)} genes)"))
	if(length(genes) > 2500) {
		genes = sample(genes, 2500)
		message_wrap("There are more than 2500 genes, only randomly sample 2500 from them.")
	}
	response = POST(qq("@{DAVID_DWS}/addList"),
		body = list("args0" = paste(genes, collapse = ","),  # inputIds
			         "args1" = idtype,             # idType
			         "args2" = as.character(Sys.time()),                    # listName
			         "args3" = 0))                           # listType

	response = GET(qq("@{DAVID_DWS}/getSpecies"))
	all_species = sapply(xml_children(httr::content(response)), xml_text)
	if(length(all_species) > 1) {
		i = grep(species, all_species)
		if(length(i) != 1) {
			cat("check your species, mapped species are:\n")
			print(all_species)
			stop_wrap("you have an error.")
		}
		GET(qq("@{DAVID_DWS}/getSpecies"),
			query = list("arg0" = i))
	} else {
		message_wrap(qq("There is one unique species (@{all_species}) mapped, no need to check species."))
	}

	message_wrap("set catalogs")
	response = GET(qq("@{DAVID_DWS}/setCategories"),
		query = list("args0" = paste(catalog, collapse = ","))
	)

	message_wrap(qq("doing enrichment"))
	response = GET(qq("@{DAVID_DWS}/getTermClusterReport"),
		query = list("args0" = 3,         # overlap, int
			         "args1" = 3,         # initialSeed, int
			         "args2" = 3,         # finalSeed, int
			         "args3" = 0.5,       # linkage, double
			         "args4" = 1))        # kappa, int

	message_wrap(qq("formatting results"))
	xml = httr::content(response)
	clusters = xml_children(xml)
	lt = lapply(clusters, function(x) {
		terms = xml_children(x)[-(1:2)]
		lt = lapply(terms, function(t) {
			fileds = xml_children(t)
			field_name = sapply(fileds, xml_name)
			field_value = sapply(fileds, xml_text)
			l = !field_name %in% c("scores", "listName")
			field_name = field_name[l]
			field_value = field_value[l]
			names(field_value) = field_name
			return(field_value)
		})
		do.call("rbind", lt)
	})
	for(i in seq_along(lt)) {
		lt[[i]] = cbind(lt[[i]], cluster = i)
	}
	tb = do.call("rbind", lt)
	tb = as.data.frame(tb, stringsAsFactors = FALSE)
	for(i in c(1, 2, 3, 4, 6, 7, 8, 11, 12, 13, 14, 15, 16, 18)) {
		tb[[i]] = as.numeric(tb[[i]])
	}

	gene_ids = lapply(strsplit(tb$geneIds, ", "), function(x) map[x])
	tb$geneIds = gene_ids
	return(tb)
}



# == title
# Perform functional enrichment on signature genes
#
# == param
# -object A `ConsensusPartitionList-class` object from `run_all_consensus_partition_methods`.
# -gene_fdr_cutoff Cutoff of FDR to define significant signature genes.
# -id_mapping If the gene IDs which are row names of the original matrix are not Entrez IDs, a
#       named vector should be provided where the names are the gene IDs in the matrix and values
#       are correspoinding Entrez IDs. The value can also be a function that converts gene IDs.
# -org_db Annotation database.
# -ontology See corresponding argumnet in `functional_enrichment,ANY-method`.
# -min_set_size The minimal size of the gene sets.
# -max_set_size The maximal size of the gene sets.
# -... Pass to `functional_enrichment,ANY-method`.
#
# == details
# For each method, the signature genes are extracted based on the best k.
#
# It calls `functional_enrichment,ConsensusPartition-method` on the consensus partitioning results for each method.
#
# For how to control the parameters of functional enrichment, see help page of `functional_enrichment,ANY-method`.
#
# == values
# A list where each element in the list corresponds to enrichment results from a single method.
#
setMethod(f = "functional_enrichment",
    signature = "ConsensusPartitionList",
    definition = function(object, gene_fdr_cutoff = cola_opt$fdr_cutoff,
    id_mapping = guess_id_mapping(rownames(object), org_db, FALSE), 
    org_db = "org.Hs.eg.db", ontology = "BP",
    min_set_size = 10, max_set_size = 1000, ...) {

    if(!grepl("\\.db$", org_db)) org_db = paste0(org_db, ".db")

    id_mapping = id_mapping
    mc.cores = 1

    if(mc.cores == 1) {
        lt = list()
        for(i in seq_along(object@list)) {
            nm = names(object@list)[i]

            best_k = suggest_best_k(object@list[[i]])
            cat("-----------------------------------------------------------\n")
            qqcat("* enrich signature genes (k = @{best_k}) to @{ontology} terms for @{nm} on @{org_db}, @{i}/@{length(object@list)}\n")
            lt[[nm]] = functional_enrichment(object@list[[i]], gene_fdr_cutoff = gene_fdr_cutoff, id_mapping = id_mapping, org_db = org_db,
                min_set_size = min_set_size, max_set_size = max_set_size, prefix = "  ", ontology = ontology, ...)
        }
    } else {
        lt = mclapply(seq_along(object@list), function(i) {
            nm = names(object@list)[i]

            best_k = suggest_best_k(object@list[[i]])
            qqcat("* enrich signature genes (k = @{best_k}) to @{ontology} terms for @{nm} on @{org_db}, @{i}/@{length(object@list)}\n")
            functional_enrichment(object@list[[i]], gene_fdr_cutoff = gene_fdr_cutoff, id_mapping = id_mapping, org_db = org_db,
                min_set_size = min_set_size, max_set_size = max_set_size, prefix = "  ", verbose = FALSE, ontology = ontology, ...)
        }, mc.cores = mc.cores)
        names(lt) = names(object@list)
    }

    return(lt)
})

# == title
# Perform functional enrichment on signature genes
#
# == param
# -object a `ConsensusPartition-class` object from `run_all_consensus_partition_methods`.
# -gene_fdr_cutoff Cutoff of FDR to define significant signature genes.
# -k Number of subgroups.
# -row_km Number of row clusterings by k-means to separate the matrix that only contains signatures.
# -id_mapping If the gene IDs which are row names of the original matrix are not Entrez IDs, a
#       named vector should be provided where the names are the gene IDs in the matrix and values
#       are correspoinding Entrez IDs. The value can also be a function that converts gene IDs.
# -org_db Annotation database.
# -ontology See corresponding argumnet in `functional_enrichment,ANY-method`.
# -min_set_size The minimal size of the gene sets.
# -max_set_size The maximal size of the gene sets.
# -verbose Whether to print messages.
# -... Pass to `functional_enrichment,ANY-method`.
#
# == details
# For how to control the parameters of functional enrichment, see help page of `functional_enrichment,ANY-method`.
#
# == value
# A list of data frames which correspond to results for the functional ontologies:
#
setMethod(f = "functional_enrichment",
    signature = "ConsensusPartition",
    definition = function(object, gene_fdr_cutoff = cola_opt$fdr_cutoff, k = suggest_best_k(object),
    row_km = NULL, id_mapping = guess_id_mapping(rownames(object), org_db, verbose), 
    org_db = "org.Hs.eg.db", ontology = "BP",
    min_set_size = 10, max_set_size = 1000, 
    verbose = TRUE, ...) {

    if(!grepl("\\.db$", org_db)) org_db = paste0(org_db, ".db")
	arg_lt = list(...)
	if("prefix" %in% names(arg_lt)) {
		prefix = arg_lt$prefix
	} else {
		prefix = ""
	}

    if(is.na(k)) {
        if(verbose) qqcat("@{prefix}- no best k was defined.\n")
        return(list(BP = NULL, MF = NULL, CC = NULL))
    }

    id_mapping = id_mapping
    
    if(is.na(k)) {
    	if(verbose) qqcat("@{prefix}- no proper number of groups found.\n")
    	return(lt)
    }
    sig_df = get_signatures(object, k = k, fdr_cutoff = gene_fdr_cutoff, plot = FALSE, verbose = FALSE, row_km = row_km)
    m = get_matrix(object)
    if(is.null(sig_df)) {
        sig_gene = NULL
    } else {
        sig_gene = rownames(m)[ sig_df[, "which_row"]]
    }
    if(verbose) qqcat("@{prefix}- @{length(sig_gene)}/@{nrow(m)} significant genes are taken from @{k}-group comparisons\n")
    
    if(length(sig_gene) == 0) {
        if(verbose) qqcat("@{prefix}- no significant genes (fdr < @{gene_fdr_cutoff}) found.\n")
        return(list(BP = NULL, MF = NULL, CC = NULL))
    }

    if("km" %in% colnames(sig_df)) {
        lt = list()
        for(km in sort(unique(sig_df$km))) {
            l = sig_df$km == km
            if(verbose) qqcat("@{prefix}- on k-means group @{km}/@{max(sig_df$km)}, @{sum(l)} genes\n")
            lt2 = functional_enrichment(sig_gene[l], id_mapping = id_mapping, org_db = org_db, ontology = ontology,
                min_set_size = min_set_size, max_set_size = max_set_size, verbose = verbose, prefix = paste0(prefix, "  "), ...)
            names(lt2) = paste0(names(lt2), "_km", km)
            lt = c(lt, lt2)
        }
    } else {
        lt = functional_enrichment(sig_gene, id_mapping = id_mapping, org_db = org_db, ontology = ontology,
            min_set_size = min_set_size, max_set_size = max_set_size, verbose = verbose, prefix = prefix, ...)
    }
    return(lt)
})

# == title
# Perform functional enrichment on signature genes
#
# == param
# -object A vector of gene IDs.
# -id_mapping If the gene IDs which are row names of the original matrix are not Entrez IDs, a
#       named vector should be provided where the names are the gene IDs in the matrix and values
#       are correspoinding Entrez IDs. The value can also be a function that converts gene IDs.
# -org_db Annotation database.
# -ontology Following ontologies are allowed: ``BP``, ``CC``, ``MF``, ``KEGG``, ``Reactome``.
#           ``MSigDb`` with the gmt file set by ``gmt_file`` argument, or ``gmt`` for general gmt gene sets.
# -min_set_size The minimal size of the gene sets.
# -max_set_size The maximal size of the gene sets.
# -verbose Whether to print messages.
# -prefix Used internally.
# -... Pass to `clusterProfiler::enrichGO`, `clusterProfiler::enrichKEGG`, `clusterProfiler::enricher`, `DOSE::enrichDO` or `ReactomePA::enrichPathway`.
#
# == details
# The function enrichment is applied by clusterProfiler, DOSE or ReactomePA packages.
#
# == value
# A data frame.
#
setMethod(f = "functional_enrichment",
    signature = "ANY",
    definition = function(object, 
    id_mapping = guess_id_mapping(object, org_db, verbose), 
    org_db = "org.Hs.eg.db", ontology = "BP",
    min_set_size = 10, max_set_size = 1000, 
    verbose = TRUE, prefix = "", ...) {

    if(!grepl("\\.db$", org_db)) org_db = paste0(org_db, ".db")

    arg_lt = list(...)

    # if(length(setdiff(ontology, c("BP", "MF", "CC", "KEGG")))) {
    #     stop_wrap("ontology can only be in 'BP', 'MF' and 'CC'")
    # }

    # lt = list(BP = NULL, MF = NULL, CC = NULL, KEGG = NULL)[ontology]
    lt = list()
    sig_gene = object
    if(length(sig_gene)) {
        if(!is.null(id_mapping)) {
            if(is.function(id_mapping)) {
                sig_gene = id_mapping(sig_gene)
            } else {
                sig_gene = id_mapping[sig_gene]
            }
        }
        sig_gene = sig_gene[!is.na(sig_gene)]
        sig_gene = unique(sig_gene)

        if(!is.null(id_mapping)) {
            if(length(sig_gene) == 0) {
                warning_wrap("Cannot match to any gene by the id mapping that user provided.")
            }
        }

        if(!is.null(id_mapping)) {
            if(verbose) qqcat("@{prefix}- @{length(sig_gene)} genes left after id mapping\n")
        }

        if(length(sig_gene)) {
            if("BP" %in% ontology) {
                if(verbose) qqcat("@{prefix}- gene set enrichment, GO:BP\n")
                res = clusterProfiler::enrichGO(
                    gene = sig_gene,
                    OrgDb         = org_db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    minGSSize = min_set_size,
                    maxGSSize = max_set_size,
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    readable      = TRUE,
                    ...)
                res = as.data.frame(res)
                # res$geneID = NULL
                lt$BP = res
            }

            if("MF" %in% ontology) {
                if(verbose) qqcat("@{prefix}- gene set enrichment, GO:MF\n")
                res = clusterProfiler::enrichGO(
                    gene = sig_gene,
                    OrgDb         = org_db,
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    minGSSize = min_set_size,
                    maxGSSize = max_set_size,
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    readable      = TRUE,
                    ...)
                res = as.data.frame(res)
                # res$geneID = NULL
                lt$MF = res
            }

            if("CC" %in% ontology) {
                if(verbose) qqcat("@{prefix}- gene set enrichment, GO:CC\n")
                res = clusterProfiler::enrichGO(
                    gene = sig_gene,
                    OrgDb         = org_db,
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    minGSSize = min_set_size,
                    maxGSSize = max_set_size,
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    readable      = TRUE,
                    ...)
                res = as.data.frame(res)
                # res$geneID = NULL
                lt$CC = res
            }

            if("KEGG" %in% ontology) {
                if(verbose) qqcat("@{prefix}- gene set enrichment, KEGG\n")
                if(is.null(arg_lt$organism)) {
                    stop_wrap("Please specify 'organism' argument for KEGG enrichment.")
                }

                res = clusterProfiler::enrichKEGG(
                    gene = sig_gene,
                    pAdjustMethod = "BH",
                    minGSSize = min_set_size,
                    maxGSSize = max_set_size,
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    ...)
                res = as.data.frame(res)
                # res$geneID = NULL
                lt$KEGG = res
            }

            if("DO" %in% ontology) {
                if(verbose) qqcat("@{prefix}- gene set enrichment, DO\n")
                res = DOSE::enrichDO(
                    gene = sig_gene,
                    ont = "DO",
                    pAdjustMethod = "BH",
                    minGSSize = min_set_size,
                    maxGSSize = max_set_size,
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    readable = TRUE,
                    ...)
                res = as.data.frame(res)
                # res$geneID = NULL
                lt$DO = res
            }

            if("MSigDb" %in% ontology) {
                if(verbose) qqcat("@{prefix}- gene set enrichment, MSigDb\n")

                if(is.null(arg_lt$gmt_file)) {
                    stop_wrap("Please specify 'gmt_file' argument for MSigDb enrichment.")
                }
                if(!grepl(".entrez.gmt$", basename(arg_lt$gmt_file))) {
                    stop_wrap("Please use the original file (*.entrez.gmt) from MsgDb.")
                }
                msigdb_ontology = gsub(".entrez.gmt$", "", basename(arg_lt$gmt_file))
                gmt = clusterProfiler::read.gmt(arg_lt$gmt_file)

                res = clusterProfiler::enricher(
                    gene = sig_gene,
                    pAdjustMethod = "BH",
                    minGSSize = min_set_size,
                    maxGSSize = max_set_size,
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1, 
                    TERM2GENE = gmt,
                    ...)
                res = as.data.frame(res)
                # res$geneID = NULL
                lt[[msigdb_ontology]] = res
            }

            if("Reactome" %in% ontology) {
                if(verbose) qqcat("@{prefix}- gene set enrichment, Reactome\n")
                if(is.null(arg_lt$organism)) {
                    stop_wrap("Please specify 'organism' argument for Reactome enrichment.")
                }
                res = ReactomePA::enrichPathway(
                    gene = sig_gene,
                    pAdjustMethod = "BH",
                    minGSSize = min_set_size,
                    maxGSSize = max_set_size,
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    ...)
                res = as.data.frame(res)
                # res$geneID = NULL
                lt$Reactome = res
            }

            if("gmt" %in% ontology) {
                if(verbose) qqcat("@{prefix}- gene set enrichment, custom gmt file\n")

                if(is.null(arg_lt$gmt_file)) {
                    stop_wrap("Please specify 'gmt_file' argument for enrichment.")
                }
                gmt_name = gsub(".gmt$", "", basename(arg_lt$gmt_file))
                gmt = clusterProfiler::read.gmt(arg_lt$gmt_file)

                res = clusterProfiler::enricher(
                    gene = sig_gene,
                    pAdjustMethod = "BH",
                    minGSSize = min_set_size,
                    maxGSSize = max_set_size,
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1, 
                    TERM2GENE = gmt,
                    ...)
                res = as.data.frame(res)
                # res$geneID = NULL
                lt[[gmt_name]] = res
            }

            # lt = lt[ontology]
        }
    }
    return(lt)
})

guess_id_type = function(id, org_db = "org.Hs.eg.db", verbose = TRUE) {
    l = grepl("^ENS.*G", id)
    if(sum(l)/length(l) > 0.5) {
        return("ENSEMBL")
    }
    l = grepl("^ENS.*T", id)
    if(sum(l)/length(l) > 0.5) {
        return("ENSEMBLTRANS")
    }
    l = grepl("^(NC|NG|NM|NR|NP|XM|XR|XP|WP)_\\d+", id)
    if(sum(l)/length(l) > 0.5) {
        return("REFSEQ")
    }

    if(!grepl("\\.db$", org_db)) {
        org_db = paste0(org_db, ".db")
    }

    all_var = getNamespaceExports(org_db)
    col = eval(parse(text = qq("AnnotationDbi::columns(@{org_db}::@{org_db})")))

    l = all_var %in% paste0(gsub(".db$", "", org_db), col)
    all_var = all_var[l]
    col = gsub(gsub(".db$", "", org_db), "", all_var)

    od_ensembl = which(col == "ENSEMBL")
    od_refseq = which(col == "REFSEQ")
    od_symbol = which(col == "SYMBOL")

    od = c(od_ensembl, od_refseq, od_symbol, setdiff(seq_along(col), c(od_ensembl, od_refseq, od_symbol)))
    col = col[od]
    all_var = all_var[od]

    p_match = numeric(length(all_var))
    for(i in seq_along(all_var)) {
        oe = try(all_ids <- AnnotationDbi::contents(getFromNamespace(all_var[i], ns = org_db)), silent = TRUE)
        if(inherits(oe, "try-error")) {
            next
        }
        if(verbose) qqcat("testing @{all_var[i]}...\n")
        all_ids = unlist(all_ids)
        l = sample(id, min(100, length(id))) %in% all_ids
        p_match[i] = sum(l)/length(l)

        if(p_match[i] > 0.5) {
            if(verbose) qqcat("  gene id might be @{col[i]} (p = @{sprintf('%.3f', p_match[i])})\n")
            return(col[i])
        }
    }
    if(verbose) qqcat("  cannot decide which gene id to use.\n")
    return(NULL)
}


guess_id_mapping = function(id, org_db = "org.Hs.eg.db", verbose = TRUE) {
    col = guess_id_type(id, org_db, verbose = verbose)
    if(is.null(col)) {
        return(NULL)
    }

    id_mapping = map_to_entrez_id(col, org_db)

    l = grepl("^ENS.*(G|T)", id) | grepl("^(NC|NG|NM|NR|NP|XM|XR|XP|WP)_\\d+", id)
    if(sum(l)/length(l) > 0.5) {
        fun = local({
            id_mapping = id_mapping
            function(x) {
                x = gsub("\\.\\d+$", "", x)
                id_mapping[x]
            }
        })
        return(fun)
    } else {
        return(id_mapping)
    }
}

# == title
# Map to Entrez IDs
#
# == param
# -from The input gene ID type. Valid values should be in, e.g. ``columns(org.Hs.eg.db::org.Hs.eg.db)``.
# -org_db The annotation database.
#
# == details
# If there are multiple mappings from the input ID type to an unique Entrez ID, randomly picked one.
#
# == value
# A named vectors where names are IDs with input ID type and values are the Entrez IDs.
#
# The returned object normally is used in `functional_enrichment`.
#
# == example
# \dontrun{
#     map_to_entrez_id("ENSEMBL")
# }
map_to_entrez_id = function(from, org_db = "org.Hs.eg.db") {

    prefix = gsub("\\.db$", "", org_db)
    if(!grepl("\\.db$", org_db)) org_db = paste0(org_db, ".db")

    x = getFromNamespace(qq("@{prefix}@{from}"), ns = org_db)
    mapped_genes = AnnotationDbi::mappedkeys(x)
    xx = AnnotationDbi::as.list(x[mapped_genes])

    ENTREZID = rep(names(xx), times = sapply(xx, length))
    df = data.frame(ENTREZID, unlist(xx), stringsAsFactors = FALSE)
    names(df) = c("ENTREZID", from)
    df = df[sample(nrow(df), nrow(df)), ]
    df = df[!duplicated(df[, 2]), ]
    x = structure(df[, 1], names = df[, 2])
    x = x[!is.na(x)]
    return(x)
}
		  
