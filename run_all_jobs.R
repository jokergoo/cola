library(methods)
library(GetoptLong)

for(p in c(0.2, 0.4, 0.6, 0.8)) {
	cmd = qq("Rscript-3.1.2 /home/guz/project/development/cola/test_not_run/test_scrnaseq.R --p @{p} --ncore 4")
	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=60:00:00,mem=10G,nodes=1:ppn=4 -N scrnaseq_subgroup_p@{p}' '@{cmd}'")
	system(cmd)
}

for(p in c(0.2, 0.4, 0.6, 0.8)) {
	cmd = qq("Rscript-3.1.2 /home/guz/project/development/cola/test_not_run/test_hipo16_rnaseq.R --p @{p} --ncore 4")
	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=60:00:00,mem=10G,nodes=1:ppn=4 -N hipo16_rnaseq_subgroup_p@{p}' '@{cmd}'")
	system(cmd)
}

for(p in c(0.2, 0.4, 0.6, 0.8)) {
	cmd = qq("Rscript-3.1.2 /home/guz/project/development/cola/test_not_run/test_tcga_gbm.R --p @{p} --ncore 4")
	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=60:00:00,mem=10G,nodes=1:ppn=4 -N TCGA_subgroup_p@{p}' '@{cmd}'")
	system(cmd)
}

for(p in c(0.2, 0.4, 0.6, 0.8)) {
	cmd = qq("Rscript-3.1.2 /home/guz/project/development/cola/test_not_run/test_asthma.R --p @{p} --ncore 4")
	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=60:00:00,mem=10G,nodes=1:ppn=4 -N LiNA_asthma_subgroup_p@{p}' '@{cmd}'")
	system(cmd)
}

for(datatype in c("cell01", "cell02", "cell03", "primary_tumor", "xenograft")) {
	cmd = qq("Rscript-3.1.2 /home/guz/project/development/cola/test_not_run/hipo15_subgroup.R --datatype @{datatype}")
 	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=10G,nodes=1:ppn=4 -N hipo15_subgroup_@{datatype}' '@{cmd}'")
 	system(cmd)
}
