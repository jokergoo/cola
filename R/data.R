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
#     rl = run_all_consensus_partition_methods(
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

