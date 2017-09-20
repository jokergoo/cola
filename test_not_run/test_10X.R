library(GetoptLong)
method = "normal"

GetoptLong("method=s", "normal")

setwd("/home/parkj/projects/Single_Cell_RNA/results_per_pid/HD1495-10X/featureCounts")


df = read.table("Colon_HD1495-10X_featureCounts.count.tsv", row.names = 1)
count = as.matrix(df[, -1])

require(limma)
expr = voom(count)$E
		

l = apply(count, 1, function(x) sum(x > 0)/length(x) > 0.5)
count = count[l ,]
expr = expr[l, ]

source(qq("@{dirname(get_scriptname())}/../load.R"))
register_top_value_fun(AAC = function(mat) AAC(t(mat), mc.cores = 4))


mat = adjust_matrix(expr)

if(method == "normal") {
	res_list = run_all_consensus_partition_methods(mat, k = 2:8, mc.cores = 4)
	saveRDS(res_list, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/Colon_HD1495-10X_subgroups.rds"))
} else {
	res = hierarchical_partition(mat)
	saveRDS(res, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/Colon_HD1495-10X_hierarchical.rds"))
}

# cmd = qq("Rscript-3.3.1 /home/guz/project/development/cola/test_not_run/test_10X.R --method normal")
# cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=100:00:00,mem=20G,nodes=1:ppn=4 -N Colon_HD1495-10X_scrnaseq_subgroup' '@{cmd}'")
# system(cmd)
# cmd = qq("Rscript-3.3.1 /home/guz/project/development/cola/test_not_run/test_10X.R --method hierarchical")
# cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=100:00:00,mem=10G -N Colon_HD1495-10X_scrnaseq_hierarchical' '@{cmd}'")
# system(cmd)
