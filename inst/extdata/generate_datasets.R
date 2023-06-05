
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
cola_rl = run_all_consensus_partition_methods(data = m, top_n = 30, cores = 6)


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
    m, cores = 6, 
    anno = anno[, c("ALL.AML"), drop = FALSE],
    anno_col = c("ALL" = "red", "AML" = "blue")
)


############# golub_cola_rh ##############
m = get_matrix(golub_cola)
set.seed(123)
golub_cola_rh = hierarchical_partition(
    m, cores = 6, 
    anno = anno[, c("ALL.AML"), drop = FALSE],
    anno_col = c("ALL" = "red", "AML" = "blue")
)


############# golub_cola_ds ##############
m = get_matrix(golub_cola)
set.seed(123)
golub_cola_ds = consensus_partition_by_down_sampling(
    m, subset = 50, cores = 6,
    anno = anno[, c("ALL.AML"), drop = FALSE],
    anno_col = c("ALL" = "red", "AML" = "blue")
)


save(cola_rl, file = "~/project/development/cola/data/cola_rl.rda", compress = "xz")
save(golub_cola, file = "~/project/development/cola/data/golub_cola.rda", compress = "xz")
save(golub_cola_ds, file = "~/project/development/cola/data/golub_cola_ds.rda", compress = "xz")
save(golub_cola_rh, file = "~/project/development/cola/data/golub_cola_rh.rda", compress = "xz")



