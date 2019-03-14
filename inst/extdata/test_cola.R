
if(grepl("tbi", Sys.info()["nodename"])) {
	root = "/home/guz"
} else {
	root = "/desktop-home/guz"
}

library(cola)

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
cola_rl = run_all_consensus_partition_methods(data = m, top_n = c(20, 30, 40), 
	top_value_method = c("sd", "MAD"), partition_method = c("hclust", "kmeans"),
	mc.cores = 4)
save(cola_rl, file = qq("@{root}/project/development/cola/data/cola_rl.rda"), compress = "xz")
cola_report(cola_rl, qq("@{root}/test_report/cola_rl_report"), mc.cores = 4)

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
cola_report(cola_rh, qq("@{root}/test_report/cola_rh_report"))
