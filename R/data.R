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


# == title (data:cola_rh)
# Example HierarchicalPartition object
#
# == details
# Following code was used to generate ``cola_rh``:
#
#   set.seed(123)
#   m = cbind(rbind(matrix(rnorm(20*20, mean = 2, sd = 0.3), nr = 20),
#                   matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20),
#                   matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20)),
#             rbind(matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20),
#                   matrix(rnorm(20*20, mean = 1, sd = 0.3), nr = 20),
#                   matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20)),
#             rbind(matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20),
#                   matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20),
#                   matrix(rnorm(20*20, mean = 1, sd = 0.3), nr = 20))
#            ) + matrix(rnorm(60*60, sd = 0.5), nr = 60)
#   cola_rh = hierarchical_partition(m, top_n = c(20, 30, 40))
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rh)
# cola_rh


