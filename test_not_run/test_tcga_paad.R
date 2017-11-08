

library(methods)
library(GetoptLong)
datatype = "mRNA"
ncore = 1
GetoptLong(
	"datatype=s", "mRNA|miRNA|methylation|protein",
	"ncore=i", "mc.cores"
)

source(qq("@{dirname(get_scriptname())}/../load.R"))
register_top_value_fun(AAC = function(mat) AAC(t(mat), mc.cores = ncore))


library(TCGA2STAT)

if(datatype == "mRNA") {
	data = getTCGA("PAAD", data.type = "RNASeq2")$dat
	data = log2(data + 1)
	data = adjust_matrix(data)
} else if(datatype == "miRNA") {
	data = getTCGA("PAAD", data.type = "miRNASeq", type = "rpmmm")$dat
	data = log2(data + 1)
	data = adjust_matrix(data)
} else if(datatype == "methylation") {
	data = getTCGA("PAAD", data.type = "methylation", type = "450K")$dat
	data = data[apply(data, 1, function(x) all(!is.na(x))), ]
} else if(datatype == "protein") {
	data = read.table("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/TCGA_subgroup/PAAD.protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data.data.txt", sep = "\t",
		header = TRUE, row.names = 1, check.names = FALSE)
	data = adjust_matrix(data)
}

if(datatype == "methylation") {
	res_list = run_all_consensus_partition_methods(data, top_n = c(10000, 20000, 30000), mc.cores = ncore,
		scale_rows = FALSE)
} else {
	data = adjust_matrix(data)
	res_list = run_all_consensus_partition_methods(data, mc.cores = ncore)
}

saveRDS(res_list, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/TCGA_subgroup/TCGA_PAAD_@{datatype}_subgroup.rds"))


# for(datatype in c("mRNA", "miRNA", "protein")) {
# 	cmd = qq("Rscript-3.3.1 /home/guz/project/development/cola/test_not_run/test_tcga_paad.R --datatype @{datatype} --ncore 4")
# 	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=@{ifelse(datatype==\"methylation\", 100, 40)}:00:00,mem=20G,nodes=1:ppn=4 -N TCGA_PAAD_@{datatype}_subgroup' '@{cmd}'")
# 	system(cmd)
# }




