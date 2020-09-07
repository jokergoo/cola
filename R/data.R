# == title (data:golub_cola)
# Example ConsensusPartitionList object from Golub dataset
#
# == details
# Following code was used to generate ``golub_cola``:
#
#     library(cola)
#     
#     library(golubEsets)  # from bioc
#     data(Golub_Merge)
#     m = exprs(Golub_Merge)
#     colnames(m) = paste0("sample_", colnames(m))
#     anno = pData(Golub_Merge)
#     
#     m[m <= 1] = NA
#     m = log10(m)
#     
#     m = adjust_matrix(m)
#     
#     library(preprocessCore)  # from bioc
#     cn = colnames(m)
#     rn = rownames(m)
#     m = normalize.quantiles(m)
#     colnames(m) = cn
#     rownames(m) = rn
#     
#     register_NMF()
#     
#     set.seed(123)
#     golub_cola = run_all_consensus_partition_methods(
#         m,
#         mc.cores = 4, 
#         anno = anno[, c("ALL.AML"), drop = FALSE],
#         anno_col = c("ALL" = "red", "AML" = "blue")
#     )
#
# == seealso
# https://jokergoo.github.io/cola_examples/Golub_leukemia/
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# golub_cola

# == title (data:golub_cola_rh)
# Example HierarchicalPartition object from Golub dataset
#
# == details
# Following code was used to generate ``golub_cola_rh``:
#
#     library(cola)
#     
#     library(golubEsets)  # from bioc
#     data(Golub_Merge)
#     m = exprs(Golub_Merge)
#     colnames(m) = paste0("sample_", colnames(m))
#     anno = pData(Golub_Merge)
#     
#     m[m <= 1] = NA
#     m = log10(m)
#     
#     m = adjust_matrix(m)
#     
#     library(preprocessCore)  # from bioc
#     cn = colnames(m)
#     rn = rownames(m)
#     m = normalize.quantiles(m)
#     colnames(m) = cn
#     rownames(m) = rn
#     
#     register_NMF()
#     
#     set.seed(123)
#     golub_cola_rh = hierarchical_partition(
#         m,
#         PAC_cutoff = 0.3,
#         mc.cores = 4, 
#         anno = anno[, c("ALL.AML"), drop = FALSE],
#         anno_col = c("ALL" = "red", "AML" = "blue")
#     )
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola_rh)
# golub_cola_rh


# == title (data:golub_cola_ds)
# Example DownSamplingConsensusPartition object from Golub dataset
#
# == details
# Following code was used to generate ``golub_cola_ds``:
#
#     library(cola)
#     data(golub_cola)
#     m = get_matrix(golub_cola)
#     set.seed(123)
#     golub_cola_ds = consensus_partition_by_down_sampling(m, subset = 50,
#         anno = get_anno(golub_cola), anno_col = get_anno_col(golub_cola),
#         top_value_method = "SD", partition_method = "kmeans")
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola_ds)
# golub_cola_ds

# == title (data:cola_rl)
# Example ConsensusPartitionList object
#
# == details
# Following code was used to generate ``cola_rl``:
#
#   set.seed(123)
#   m = cbind(rbind(matrix(rnorm(20*20, mean = 1,   sd = 0.5), nr = 20),
#                   matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20),
#                   matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20)),
#             rbind(matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20),
#                   matrix(rnorm(20*20, mean = 1,   sd = 0.5), nr = 20),
#                   matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20)),
#             rbind(matrix(rnorm(20*20, mean = 0.5, sd = 0.5), nr = 20),
#                   matrix(rnorm(20*20, mean = 0.5, sd = 0.5), nr = 20),
#                   matrix(rnorm(20*20, mean = 1,   sd = 0.5), nr = 20))
#            ) + matrix(rnorm(60*60, sd = 0.5), nr = 60)
#   cola_rl = run_all_consensus_partition_methods(data = m, top_n = c(20, 30, 40))
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# cola_rl

