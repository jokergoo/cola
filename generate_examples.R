
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
    cmd = qq("perl @{root}/project/development/ngspipeline2/bsub_single_line.pl --hour 50 --memory 20 --core 4 --name @{name} --output @{root}/project/development/cola_examples/@{name}.out --command '@{cmd}'")
    system(cmd)
}

run_script(qq("@{root}/project/development/cola/inst/extdata/test_Golub.R"))
run_script(qq("@{root}/project/development/cola/inst/extdata/test_tcga_gbm.R"))
run_script(qq("@{root}/project/development/cola/inst/extdata/test_Ritz_ALL.R"))
run_script(qq("@{root}/project/development/cola/inst/extdata/test_HSMM_single_cell.R"))
run_script(qq("@{root}/project/development/cola/inst/extdata/test_MCF10CA_scRNAseq.R"))
run_script(qq("@{root}/project/development/cola/inst/extdata/test_cola.R"))
