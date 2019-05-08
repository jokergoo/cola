# root = "/home/guz"
root = "/desktop-home/guz"


library(cola)
library(GetoptLong)

set.seed(123)
m = cbind(rbind(matrix(rnorm(20*20, mean = 1,   sd = 0.5), nr = 20),
              matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20),
              matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20)),
        rbind(matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20),
              matrix(rnorm(20*20, mean = 1,   sd = 0.5), nr = 20),
              matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20)),
        rbind(matrix(rnorm(20*20, mean = 0.5, sd = 0.5), nr = 20),
              matrix(rnorm(20*20, mean = 0.5, sd = 0.5), nr = 20),
              matrix(rnorm(20*20, mean = 1,   sd = 0.5), nr = 20))
       ) + matrix(rnorm(60*60, sd = 0.5), nr = 60)
cola_rl = run_all_consensus_partition_methods(data = m, top_n = c(20, 30, 40))
save(cola_rl, file = qq("@{root}/project/development/cola/data/cola_rl.rda"), compress = "xz")

set.seed(123)
m = cbind(rbind(matrix(rnorm(20*20, mean = 2, sd = 0.3), nr = 20),
              matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20),
              matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20)),
        rbind(matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20),
              matrix(rnorm(20*20, mean = 1, sd = 0.3), nr = 20),
              matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20)),
        rbind(matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20),
              matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20),
              matrix(rnorm(20*20, mean = 1, sd = 0.3), nr = 20))
       ) + matrix(rnorm(60*60, sd = 0.5), nr = 60)
cola_rh = hierarchical_partition(m, top_n = c(20, 30, 40), PAC_cutoff = 0.3)
save(cola_rh, file = qq("@{root}/project/development/cola/data/cola_rh.rda"), compress = "xz")

library(GetoptLong)
run_script = function(script) {
    cmd = qq("module load R/3.3.1; Rscript @{script};")
    name = qq("cola_@{gsub('.R$', '', basename(script))}")
    cmd = qq("perl @{root}/project/development/ngspipeline2/bsub_single_line.pl --hour 50 --memory 20 --core 4 --name @{name} --dir @{root}/project/development/cola_examples/ --command '@{cmd}' --enforce")
    # cmd = qq("perl @{root}/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=4:00:00,mem=10G,nodes=1:ppn=4 -N @{name}' '@{cmd}'")
    system(cmd)
}

run_script(qq("@{root}/project/development/cola/inst/extdata/test_Golub.R"))
run_script(qq("@{root}/project/development/cola/inst/extdata/test_tcga_gbm.R"))
run_script(qq("@{root}/project/development/cola/inst/extdata/test_Ritz_ALL.R"))
run_script(qq("@{root}/project/development/cola/inst/extdata/test_HSMM_single_cell.R"))
run_script(qq("@{root}/project/development/cola/inst/extdata/test_MCF10CA_scRNAseq.R"))

run_script(qq("@{root}/project/development/cola/inst/extdata/test_cola.R"))



###### just geneate reports ###
root = "/desktop-home/guz"

## Golub
# rl = readRDS(qq("@{root}/project/development/cola_examples/Golub_leukemia/Golub_leukemia_subgroup.rds"))
# cola_report(rl, output_dir = qq("@{root}/project/development/cola_examples/Golub_leukemia/Golub_leukemia_subgroup_cola_report"), mc.cores = 4)

rh = readRDS(qq("@{root}/project/development/cola_examples/Golub_leukemia/Golub_leukemia_subgroup_hierarchical_partition.rds"))
cola_report(rh, output_dir = qq("@{root}/project/development/cola_examples/Golub_leukemia/Golub_leukemia_subgroup_hierarchical_partition_cola_report"), mc.cores = 4)

## HSMM single cell
# rl = readRDS(qq("@{root}/project/development/cola_examples/HSMM_single_cell/HSMM_single_cell_subgroup.rds"))
# cola_report(rl, output_dir = qq("@{root}/project/development/cola_examples/HSMM_single_cell/HSMM_single_cell_subgroup_cola_report"), mc.cores = 4)

rh = readRDS(qq("@{root}/project/development/cola_examples/HSMM_single_cell/HSMM_single_cell_subgroup_hierarchical_partition.rds"))
cola_report(rh, output_dir = qq("@{root}/project/development/cola_examples/HSMM_single_cell/HSMM_single_cell_subgroup_hierarchical_partition_cola_report"), mc.cores = 4)

## MCF10CA
# rl = readRDS(qq("@{root}/project/development/cola_examples/MCF10CA_scRNAseq/MCF10CA_scRNAseq_subgroup.rds"))
# cola_report(rl, output_dir = qq("@{root}/project/development/cola_examples/MCF10CA_scRNAseq/MCF10CA_scRNAseq_subgroup_cola_report"), mc.cores = 4)

rh = readRDS(qq("@{root}/project/development/cola_examples/MCF10CA_scRNAseq/MCF10CA_scRNAseq_subgroup_hierarchical_partition.rds"))
cola_report(rh, output_dir = qq("@{root}/project/development/cola_examples/MCF10CA_scRNAseq/MCF10CA_scRNAseq_subgroup_hierarchical_partition_cola_report"), mc.cores = 4)

## Ritz ALl
# rl = readRDS(qq("@{root}/project/development/cola_examples/Ritz_ALL/Ritz_ALL_subgroup.rds"))
# cola_report(rl, output_dir = qq("@{root}/project/development/cola_examples/Ritz_ALL/Ritz_ALL_subgroup_cola_report"), mc.cores = 4)

rh = readRDS(qq("@{root}/project/development/cola_examples/Ritz_ALL/Ritz_ALL_subgroup_hierarchical_partition.rds"))
cola_report(rh, output_dir = qq("@{root}/project/development/cola_examples/Ritz_ALL/Ritz_ALL_subgroup_hierarchical_partition_cola_report"), mc.cores = 4)

## TCGA GBM
# rl = readRDS(qq("@{root}/project/development/cola_examples/TCGA_GBM/TCGA_GBM_subgroup.rds"))
# cola_report(rl, output_dir = qq("@{root}/project/development/cola_examples/TCGA_GBM/TCGA_GBM_subgroup_cola_report"), mc.cores = 4)

rh = readRDS(qq("@{root}/project/development/cola_examples/TCGA_GBM/TCGA_GBM_subgroup_hierarchical_partition.rds"))
cola_report(rh, output_dir = qq("@{root}/project/development/cola_examples/TCGA_GBM/TCGA_GBM_subgroup_hierarchical_partition_cola_report"), mc.cores = 4)


### doing function enrichment
library() 
setMethod(f = "GO_enrichment",
    signature = "ConsensusPartitionList",
    definition = function(object, cutoff = 0.05,
    id_mapping = NULL, id_mapping_pre_fun = NULL, org_db = "org.Hs.eg.db",
    min_set_size = 10, max_set_size = 1000) {

    if(!grepl("\\.db$", org_db)) org_db = paste0(org_db, ".db")

    lt = list()
    for(i in seq_along(object@list)) {
        nm = names(object@list)[i]
        lt[[nm]] = list(BP = NULL, MF = NULL, CC = NULL)

        qqcat("* enrich signature genes to GO terms for @{nm} on @{org_db}\n")
        lt[[nm]] = GO_enrichment(object@list[[i]], cutoff = cutoff, id_mapping = id_mapping, 
            id_mapping_pre_fun = id_mapping_pre_fun, org_db = org_db,
            min_set_size = min_set_size, max_set_size = max_set_size)
    }

    return(lt)
}

setMethod(f = "GO_enrichment",
    signature = "ConsensusPartition",
    definition = function(object, cutoff = 0.05, k = guess_best_k(object),
    id_mapping = NULL, id_mapping_pre_fun = NULL, org_db = "org.Hs.eg.db",
    min_set_size = 10, max_set_size = 1000) {

    if(!grepl("\\.db$", org_db)) org_db = paste0(org_db, ".db")

    lt = list(BP = NULL, MF = NULL, CC = NULL)
    sig_df = get_signatures(object, k = k, plot = FALSE, verbose = FALSE)
    if(is.null(sig_df)) {
        sig_gene = NULL
    } else {
        sig_gene = rownames(sig_df[sig_df$fdr < cutoff, , drop = FALSE])
    }

    if(length(sig_gene)) {
        if(!is.null(id_mapping_pre_fun)) sig_gene = id_mapping_pre_fun(sig_gene)
        if(!is.null(id_mapping)) sig_gene = id_mapping[sig_gene]
        sig_gene = sig_gene[!is.na(sig_gene)]
        sig_gene = unique(sig_gene)

        if(length(sig_gene)) {
            ego = clusterProfiler::enrichGO(gene = sig_gene,
                OrgDb         = org_db,
                ont           = "BP",
                pAdjustMethod = "BH",
                minGSSize = min_set_size,
                maxGSSize = max_set_size,
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)
            ego = as.data.frame(ego)
            lt$BP = ego

            ego = clusterProfiler::enrichGO(gene = sig_gene,
                OrgDb         = org_db,
                ont           = "MF",
                pAdjustMethod = "BH",
                minGSSize = min_set_size,
                maxGSSize = max_set_size,
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)
            ego = as.data.frame(ego)
            lt$MF = ego

            ego = clusterProfiler::enrichGO(gene = sig_gene,
                OrgDb         = org_db,
                ont           = "CC",
                pAdjustMethod = "BH",
                minGSSize = min_set_size,
                maxGSSize = max_set_size,
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)
            ego = as.data.frame(ego)
            lt$CC = ego
        }
    }
    return(lt)
})


map_to_entrez_id = function(from, org_db = "org.Hs.eg.db") {

    prefix = gsub("\\.db$", "", org_db)
    if(!grepl("\\.db$", org_db)) org_db = paste0(org_db, ".db")

    x <- getFromNamespace(qq("@{prefix}@{from}"), ns = org_db)
    mapped_genes <- mappedkeys(x)
    xx <- as.list(x[mapped_genes])

    ENTREZID = rep(names(xx), times = sapply(xx, length))
    df = data.frame(ENTREZID, unlist(xx))
    names(df) = c("ENTREZID", from)
    df = df[sample(nrow(df), nrow(df)), ]
    df = df[!duplicated(df[, 2]), ]
    x = structure(df[, 1], names = df[, 2])
    return(x)
}

# Golub_leukemia: HU6800, hu6800.db
library(hu6800.db)
x <- hu6800ENTREZID
mapped_probes <- mappedkeys(x)
id_mapping <- unlist(as.list(x[mapped_probes]))
rl = readRDS(qq("@{root}/project/development/cola_examples/Golub_leukemia/Golub_leukemia_subgroup.rds"))
GO_enrichment(rl, id_mapping = id_mapping)

# HSMM single cell: Ensembl, no version numbers
id_mapping = map_to_entrez_id("ENSEMBL")
rl = readRDS(qq("@{root}/project/development/cola_examples/HSMM_single_cell/HSMM_single_cell_subgroup.rds"))
GO_enrichment(rl, id_mapping = id_mapping)

# MCF10CA, Ensembl, gencode 19
id_mapping = map_to_entrez_id("ENSEMBL")
rl = readRDS(qq("@{root}/project/development/cola_examples/MCF10CA_scRNAseq/MCF10CA_scRNAseq_subgroup.rds"))
GO_enrichment(rl, id_mapping = id_mapping, id_mapping_pre_fun = function(x) gsub("\\.\\d+$", "", x))

# Ritz ALL: HGU95aV2, hgu95av2.db
library(hgu95av2.db)
x <- hgu95av2ENTREZID
mapped_probes <- mappedkeys(x)
id_mapping <- unlist(as.list(x[mapped_probes]))
rl = readRDS(qq("@{root}/project/development/cola_examples/Ritz_ALL/Ritz_ALL_subgroup.rds"))
GO_enrichment(rl, id_mapping = id_mapping)

# TCGA GBM, gene symbol
library(org.Hs.eg.db)
id_mapping = map_to_entrez_id("SYMBOL")
rl = readRDS(qq("@{root}/project/development/cola_examples/TCGA_GBM/TCGA_GBM_subgroup.rds"))
GO_enrichment(rl, id_mapping = id_mapping)
