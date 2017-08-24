
load("~/analysis/hipo_050_rnaseq_gencode_v19_expression.Rdata")

source("~/project/cola/load.R")
library(GetoptLong)

l = expression$anno$sample_type == "tumor"
count = as.matrix(expression$count)[, l]
expr = as.matrix(expression$deseq2)[, l]
anno = expression$anno[l, ]

expr_control = as.matrix(expression$deseq2)[, expression$anno$sample_type == "control"]
data = expr

data = adjust_matrix(data)

anno2 = anno[, c("age_at_dx", "sex", "SUBTYPE_1", "SUBTYPE_2", "ABCGCB_Hans", "purity",
	             "MYC.Status.Consensus", "BCL2.Consensus", "BCL6.Consensus", "IGH.Consensus")]
for(i in seq_len(ncol(anno2))) {
	if(is.character(anno2[[i]])) {
		anno2[[i]][anno2[[i]] == ""] = NA
	}
}

anno2 = anno2[, c("SUBTYPE_1", "SUBTYPE_2", "ABCGCB_Hans", "age_at_dx", "purity")]

unique = function(x) base::unique(x[!is.na(x)])
library(circlize)
set.seed(123)
SUBTYPE_1_color = c("DLBCL" = "black", "FL" = "red", "DH-BL" = "blue")
SUBTYPE_2_color = c("DLBCL" = "black", "FL_1/2" = "red", "FL_3A" = rand_color(1, "red"), 
	"FL_2/3A" = rand_color(1, "red"), "FL_1/2/3A" = rand_color(1, "red"), "FL_3B" = rand_color(1, "red"),
	"DH-BL" = "blue", "PCNSL" = "green", "SCNSL" = "darkgreen")
hans_color = structure(rand_color(3, luminosity = "bright"), names = c("ABC", "GCB", "TypeIII"))
anno_col = list(SUBTYPE_2 = SUBTYPE_2_color,
	            SUBTYPE_1 = SUBTYPE_1_color,
	            ABCGCB_Hans = hans_color,
	            age_at_dx = colorRamp2(range(anno2[, "age_at_dx"], na.rm = TRUE), c("white", "orange")),
	            purity = colorRamp2(range(anno2[, "purity"], na.rm = TRUE), c("white", "purple")))


set.seed(123)
res = consensus_partition(data, known_anno = anno2, known_col = anno_col,
	top_method = "MAD", partition_method = "kmeans")

set.seed(456)
res = hierarchical_partition(data, known_anno = anno2, known_col = anno_col,
	top_method = "MAD", partition_method = "kmeans")


