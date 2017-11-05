library(methods)
library(GetoptLong)
ncore = 1
GetoptLong(
	"ncore=i", "mc.cores"
)

source(qq("@{dirname(get_scriptname())}/../load.R"))
# library(cola)
register_top_value_fun(AAC = function(mat) AAC(t(mat), mc.cores = ncore))

data = read.table("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/Glioblastoma_expressed_genes.txt", header = TRUE, row.names = 1)

data = adjust_matrix(data)
res = run_all_consensus_partition_methods(data, mc.cores = ncore)

saveRDS(res, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/gbm_scrnaseq_subgroup.rds"))

cola_report(res, output_dir = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test/gbm_scrnaseq_subgroup_report")

# cmd = qq("Rscript-3.3.1 /home/guz/project/development/cola/test_not_run/test_gbm_scrnaseq.R --ncore 4")
# cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=20G,nodes=1:ppn=4 -N gbm_scrnaseq_subgroup' '@{cmd}'")
# system(cmd)

