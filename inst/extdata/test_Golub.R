library(GetoptLong)
library(cola)

library(golubEsets)
data(Golub_Merge)
m = exprs(Golub_Merge)
colnames(m) = paste0("sample_", colnames(m))
anno = pData(Golub_Merge)

m[m <= 1] = NA
m = log10(m)

register_NMF()

m = adjust_matrix(m)

library(preprocessCore)
cn = colnames(m)
m = normalize.quantiles(m)
colnames(m) = cn
res = run_all_consensus_partition_methods(m, max_k = 6, mc.cores = 4, 
	anno = anno[, c("ALL.AML"), drop = FALSE],
	anno_col = c("ALL" = "red", "AML" = "blue")
)

saveRDS(res, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/Golub_subgroup.rds"))
cola_report(res, output_dir = "/icgc/dkfzlsdf/analysis/B080/guz/cola_test/Golub_cl_report")


# cmd = qq("module load R/3.3.1; Rscript /home/guz/project/development/cola/test_not_run/test_Golub.R")
# cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=20G,nodes=1:ppn=4 -N Golub_subgroup' '@{cmd}'")
# system(cmd)
# par(mfrow = c(5, 6))
ha = HeatmapAnnotation(df = anno[, c("ALL.AML"), drop = FALSE])
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 5, nc = 6)))
i = 0
for(power in c(1, 2, 3, 4, 5)) {
	i = i + 1
	j = 0
	for(min_cor in c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) {
		j = j + 1
		s = ATC(m, min_cor = min_cor, power = power)
		# plot(sort(s), type = "l", main = qq("min_cor = @{min_cor}, power = @{power}"))
		ind = order(s, decreasing = TRUE)[1:500]
		col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
		ht = Heatmap(t(scale(t(m[ind, ]))), show_row_names = FALSE, column_title = qq("min_cor = @{min_cor}, power = @{power}"),
			show_row_dend = FALSE, show_column_names = FALSE, show_column_dend = FALSE,
			show_heatmap_legend = FALSE, col = col_fun, top_annotation = ha)
		pushViewport(viewport(layout.pos.row = i, layout.pos.col = j))
		draw(ht, newpage = FALSE)
		popViewport()
	}
}
popViewport()
