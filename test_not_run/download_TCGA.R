disease = c("ACC", "BLCA", "BRCA", "CESC", 'CHOL', "COAD", "COADREAD",
	"DLBC", "ESCA", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC",
	"KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD",
	"PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
	"UCEC", "UCS", "UVM")
data.type = c("mRNA_Array", "RNASeq2", "Methylation", "CNA_SNP", "CNV_SNP")
type = list("mRNA_Array" = c("G450", "U133"), 
	"RNASeq2" = "", 
	"Methylation" = c("450K"), 
	"CNA_SNP" = "", 
	"CNV_SNP" = "")

library(TCGA2STAT)
library(GetoptLong)

for(d in disease) {
	for(t in data.type) {
		for(p in type[[t]]) {
			obj = getTCGA(disease = d, data.type = t, clinical = TRUE, type = p, p = 2)
			if(!is.null(obj)) {
				saveRDS(obj, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/TCGA/TCGA_@{d}_@{t}_@{p}.rds"))
			}
		}
	}
}
