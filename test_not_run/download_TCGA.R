library(GetoptLong)
GetoptLong(
	"disease=s", "disease"
)

data.type = c("mRNA_Array", "RNASeq2", "Methylation", "CNA_SNP", "CNV_SNP")
type = list("mRNA_Array" = c("G450", "U133"), 
	"RNASeq2" = "", 
	"Methylation" = c("450K"), 
	"CNA_SNP" = "", 
	"CNV_SNP" = "")

library(TCGA2STAT)

setwd("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/TCGA/")

for(d in disease) {
	for(t in data.type) {
		for(p in type[[t]]) {
			obj = getTCGA(disease = d, data.type = t, clinical = TRUE, type = p, p = 1)
			if(!is.null(obj)) {
				saveRDS(obj, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/TCGA/TCGA_@{d}_@{t}_@{p}.rds"))
			}
		}
	}
}

# disease = c("ACC", "BLCA", "BRCA", "CESC", 'CHOL', "COAD", "COADREAD",
# 	"DLBC", "ESCA", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC",
# 	"KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD",
# 	"PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
# 	"UCEC", "UCS", "UVM")
# for(d in disease) {
# 	cmd = qq("Rscript-3.3.1 /home/guz/project/development/cola/test_not_run/download_TCGA.R --disease @{d}")
# 	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=50:00:00,mem=10G -N download_TCGA_@{d}' '@{cmd}'")
# 	system(cmd)
# }
