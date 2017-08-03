library(methods)
library(GetoptLong)
p = 0.8
ncore = 1
GetoptLong(
	"p=f", "0.8",
	"ncore=i", "mc.cores"
)

# library(cola)
source(qq("@{dirname(get_scriptname())}/../load.R"))
register_top_value_fun(AAC = function(mat) AAC(t(mat), mc.cores = ncore))


source("/home/guz/project/analysis/LiNA/asthma_cohort/script/asthma_head.R")

data = deseq2[gt[rownames(deseq2)] == "protein_coding", ]

data = adjust_matrix(data)
res = run_all_consensus_partition_methods(data, top_n = c(1000, 2000, 4000), k = 2:6, p_sampling = p, 
	known_anno = data.frame(phenotype = SAMPLE$phenotype), 
	known_col = list(phenotype = c(Asthma = "red", Control = "black")), mc.cores = ncore)

saveRDS(res, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/hipo35_asthma_rnaseq_subgroup_p@{p}.rds"))


# for(p in c(0.2, 0.4, 0.6, 0.8)) {
# 	cmd = qq("Rscript-3.3.1 /home/guz/project/development/cola/test_not_run/test_asthma.R --p @{p} --ncore 4")
# 	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=20G,nodes=1:ppn=4 -N lina_asthma_subgroup_p@{p}' '@{cmd}'")
# 	system(cmd)
# }

