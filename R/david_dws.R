 


DAVID_ALL_ID_TYPES = "AFFYMETRIX_3PRIME_IVT_ID,AFFYMETRIX_EXON_GENE_ID,AFFYMETRIX_SNP_ID,AGILENT_CHIP_ID,AGILENT_ID,AGILENT_OLIGO_ID,ENSEMBL_GENE_ID,ENSEMBL_TRANSCRIPT_ID,ENTREZ_GENE_ID,FLYBASE_GENE_ID,FLYBASE_TRANSCRIPT_ID,GENBANK_ACCESSION,GENOMIC_GI_ACCESSION,GENPEPT_ACCESSION,ILLUMINA_ID,IPI_ID,MGI_ID,PFAM_ID,PIR_ID,PROTEIN_GI_ACCESSION,REFSEQ_GENOMIC,REFSEQ_MRNA,REFSEQ_PROTEIN,REFSEQ_RNA,RGD_ID,SGD_ID,TAIR_ID,UCSC_GENE_ID,UNIGENE,UNIPROT_ACCESSION,UNIPROT_ID,UNIREF100_ID,WORMBASE_GENE_ID,WORMPEP_ID,ZFIN_ID"
DAVID_ALL_ID_TYPES = strsplit(DAVID_ALL_ID_TYPES, ",")[[1]]
DAVID_ALL_CATALOGS = "BBID,BIND,BIOCARTA,BLOCKS,CGAP_EST_QUARTILE,CGAP_SAGE_QUARTILE,CHROMOSOME,COG_NAME,COG_ONTOLOGY,CYTOBAND,DIP,EC_NUMBER,ENSEMBL_GENE_ID,ENTREZ_GENE_ID,ENTREZ_GENE_SUMMARY,GENETIC_ASSOCIATION_DB_DISEASE,GENERIF_SUMMARY,GNF_U133A_QUARTILE,GENETIC_ASSOCIATION_DB_DISEASE_CLASS,GOTERM_BP_2,GOTERM_BP_1,GOTERM_BP_4,GOTERM_BP_3,GOTERM_BP_FAT,GOTERM_BP_5,GOTERM_CC_1,GOTERM_BP_ALL,GOTERM_CC_3,GOTERM_CC_2,GOTERM_CC_5,GOTERM_CC_4,GOTERM_MF_1,GOTERM_MF_2,GOTERM_CC_FAT,GOTERM_CC_ALL,GOTERM_MF_5,GOTERM_MF_FAT,GOTERM_MF_3,GOTERM_MF_4,HIV_INTERACTION_CATEGORY,HOMOLOGOUS_GENE,GOTERM_MF_ALL,HIV_INTERACTION,MINT,NCICB_CAPATHWAY_INTERACTION,INTERPRO,KEGG_PATHWAY,PANTHER_FAMILY,PANTHER_BP_ALL,OMIM_DISEASE,OFFICIAL_GENE_SYMBOL,PANTHER_SUBFAMILY,PANTHER_PATHWAY,PANTHER_MF_ALL,PIR_SUMMARY,PIR_SEQ_FEATURE,PFAM,PRODOM,PRINTS,PIR_TISSUE_SPECIFICITY,PIR_SUPERFAMILY,SMART,SP_COMMENT,SP_COMMENT_TYPE,SP_PIR_KEYWORDS,PROSITE,PUBMED_ID,REACTOME_INTERACTION,REACTOME_PATHWAY,UNIGENE_EST_QUARTILE,UP_SEQ_FEATURE,UP_TISSUE,ZFIN_ANATOMY,SSF,TIGRFAMS,UCSC_TFBS"
DAVID_ALL_CATALOGS = strsplit(DAVID_ALL_CATALOGS, ",")[[1]]

# == title
# Doing DAVID analysis
#
# == param
# -genes a vector of gene identifiers
# -email the email that user registered on DAVID web service (https://david.ncifcrf.gov/content.jsp?file=WS.html )
# -catalog a vector of function catalogs. Valid values are in ``cola:::DAVID_ALL_CATALOGS``.
# -idtype id types for the input gene list. Valid values are in ``cola:::DAVID_ALL_ID_TYPES``.
# -species full species name if the id type is not uniquely mapped to one single species
#
# == details
# If you want to run this function multiple times, please set time intervals between runs.
# 
# == value
# A data frame with functional enrichment results
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
		stop("You need to register to DAVID web service.")
	}

	if(!idtype %in% DAVID_ALL_ID_TYPES) {
		cat("idtype is wrong, it should be in:\n")
		print(DAVID_ALL_ID_TYPES)
		stop("You have an error.")
	}
	if(!all(catalog %in% DAVID_ALL_CATALOGS)) {
		cat("catalog is wrong, it should be in:\n")
		print(DAVID_ALL_CATALOGS)
		stop("You have an error.")
	}

	if(grepl("ENSEMBL", idtype)) {
		map = structure(genes, names = gsub("\\.\\d+$", "", genes))
		genes = names(map)
	}

	message(qq("Idtype: @{idtype}"))
	message(qq("Catalog: @{paste(catalog, collapse = ' ')}"))
	
	DAVID_DWS = "https://david-d.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint"

	# login
	message(qq("log to DAVID web service with @{email}"))
	response = GET(qq("@{DAVID_DWS}/authenticate"),
		query = list("args0" = email)
	)
	if(xml_text(httr::content(response)) != "true") {
		stop(qq("@{email} has not been registed."))
	}

	# add gene list
	message(qq("add gene list (@{length(genes)} genes)"))
	if(length(genes) > 2500) {
		genes = sample(genes, 2500)
		message("There are more than 2500 genes, only randomly sample 2500 from them.")
	}
	response = POST(qq("@{DAVID_DWS}/addList"),
		body = list("args0" = paste(genes, collapse = ","),  # inputIds
			         "args1" = idtype,             # idType
			         "args2" = Sys.time(),                    # listName
			         "args3" = 0))                           # listType

	response = GET(qq("@{DAVID_DWS}/getSpecies"))
	all_species = sapply(xml_children(httr::content(response)), xml_text)
	if(length(all_species) > 1) {
		i = grep(species, all_species)
		if(length(i) != 1) {
			cat("check your species, mapped species are:\n")
			print(all_species)
			stop("you have an error.")
		}
		GET(qq("@{DAVID_DWS}/getSpecies"),
			query = list("arg0" = i))
	} else {
		message(qq("There is one unique species (@{all_species}) mapped, no need to check species."))
	}

	message("set catalogs")
	response = GET(qq("@{DAVID_DWS}/setCategories"),
		query = list("args0" = paste(catalog, collapse = ","))
	)

	message(qq("doing enrichment"))
	response = GET(qq("@{DAVID_DWS}/getTermClusterReport"),
		query = list("args0" = 3,         # overlap, int
			         "args1" = 3,         # initialSeed, int
			         "args2" = 3,         # finalSeed, int
			         "args3" = 0.5,         # linkage, double
			         "args4" = 1))        # kappa, int

	message(qq("formatting results"))
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
# Reduce DAVID results
#
# == param
# -tb object from `submit_to_david`
# -fdr_cutoff cutoff for fdr (the ``benjamini`` column)
# -hit_cutoff cutoff for number of genes in a term
#
# == details
# For each cluster and each functional category, the function picks the most
# significant term.
#
# == value
# A subset of rows in ``tb``.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
reduce_david_results = function(tb, fdr_cutoff = 0.05, hit_cutoff = 5) {
	l = tb$benjamini <= fdr_cutoff & tb$listHits >= hit_cutoff & !duplicated(tb$termName)
	if(sum(l) < 1) {
		cat("No record left.\n")
		return(NULL)
	}
	tb = tb[l, , drop = FALSE]

	tb_reduced = do.call("rbind", lapply(split(tb, tb$cluster), function(x) {
		do.call("rbind", lapply(split(x, x$categoryName), function(y) {
			i = which.min(y$benjamini)
			y[i, , drop = FALSE]
		}))
	}))
	rownames(tb_reduced) = NULL

	attr(tb_reduced, "reduced") = TRUE
	return(unique(tb_reduced))
}
