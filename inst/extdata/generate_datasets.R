
############# cola_rl ##############
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


############# golub_cola ##############
library(golubEsets)  # from bioc
data(Golub_Merge)
m = exprs(Golub_Merge)
colnames(m) = paste0("sample_", colnames(m))
anno = pData(Golub_Merge)

m[m <= 1] = NA
m = log10(m)

m = adjust_matrix(m)

library(preprocessCore)  # from bioc
cn = colnames(m)
rn = rownames(m)
m = normalize.quantiles(m)
colnames(m) = cn
rownames(m) = rn

set.seed(123)
golub_cola = run_all_consensus_partition_methods(
    m, mc.cores = 2, 
    anno = anno[, c("ALL.AML"), drop = FALSE],
    anno_col = c("ALL" = "red", "AML" = "blue")
)


############# golub_cola_rh ##############
m = get_matrix(golub_cola)
set.seed(123)
golub_cola_rh = hierarchical_partition(
    m, mc.cores = 2, 
    anno = get_anno(golub_cola), 
    anno_col = get_anno_col(golub_cola)
)


############# golub_cola_ds ##############
m = get_matrix(golub_cola)
set.seed(123)
golub_cola_ds = consensus_partition_by_down_sampling(
  m, subset = 50, mc.cores = 2,
  anno = get_anno(golub_cola), 
  anno_col = get_anno_col(golub_cola)
)


save(cola_rl, file = "~/project/cola/data/cola_rl.rda", compress = "xz")
save(golub_cola, file = "~/project/cola/data/golub_cola.rda", compress = "xz")
save(golub_cola_ds, file = "~/project/cola/data/golub_cola_ds.rda", compress = "xz")
save(golub_cola_rh, file = "~/project/cola/data/golub_cola_rh.rda", compress = "xz")


############################
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

