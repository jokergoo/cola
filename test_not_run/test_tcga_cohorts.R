 
library(methods)
library(GetoptLong)
file = ""
ncore = 1
GetoptLong(
	"file=s", "file",
	"ncore=i", "mc.cores"
)

library(cola)
# source(qq("@{dirname(get_scriptname())}/../load.R"))
register_top_value_fun(AAC = function(mat) AAC(t(mat), mc.cores = ncore))

x = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/TCGA/@{file}.rds"))
data = x$dat

is_methylation = grepl("Methylation", file)
is_rnaseq = grepl("RNASeq2", file)
data = data[apply(data, 1, function(x) all(!is.na(x))), ]

if(is_rnaseq) {
	data = log2(data + 1)
}

if(is_methylation) {
	q(save = "no")
	res_list = run_all_consensus_partition_methods(data, k = 2:6, top_n = c(10000, 20000, 30000), mc.cores = ncore,
		scale_rows = FALSE)
} else {
	data = adjust_matrix(data)
	res_list = run_all_consensus_partition_methods(data, k = 2:6, top_n = c(2000, 4000, 6000), mc.cores = ncore)
}
saveRDS(res_list, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/TCGA_subgroup/@{file}_subgroups.rds"))
dir.create(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/TCGA_subgroup/@{file}_subgroups_report"))
cola_report(res_list, qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/TCGA_subgroup/@{file}_subgroups_report"))

res_hc = hierarchical_partition(data)
saveRDS(res_hc, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/TCGA_subgroup/@{file}_hc.rds"))
dir.create(qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/TCGA_subgroup/@{file}_hc_report"))
cola_report(res_hc, qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/TCGA_subgroup/@{file}_hc_report"))


# for(f in dir("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/TCGA")) {
# 	f = gsub("\\.rds$", "", f)
#   if(grepl("Methylation", f)) next
# 	cmd = qq("Rscript-3.3.1 /home/guz/project/development/cola/test_not_run/test_tcga_cohorts.R --file @{f} --ncore 4")
# 	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=20G,nodes=1:ppn=4 -N @{f}_subgroup' '@{cmd}'")
# 	system(cmd)
# }
