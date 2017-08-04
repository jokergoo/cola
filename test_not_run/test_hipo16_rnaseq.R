library(methods)
library(GetoptLong)
p = 0.8
ncore = 1
subtype = c("IDH", "MES", "RTK_I", "RTK_II")
GetoptLong(
	"p=f", "0.8",
	"ncore=i", "mc.cores",
	"subtype=s{2,}", "subtypes"
)

# library(cola)
source(qq("@{dirname(get_scriptname())}/../load.R"))

############################################################
anno_text =
"id    hipo_id     subtype batch
AK015  H016-WFRL   IDH   1
AK041   H016-7EN2   IDH  1
AK066   H016-DVZSMF IDH  1
AK068   H016-1VS79M IDH  1
AK076   H016-JGR2   IDH  1
AK085   H016-GA6F   IDH  1
AK102   H016-D3VY   IDH  1
AK103   H016-XEMY   IDH  1
AK124   H016-3N2CQY IDH  1
AK199   H016-SEMSBV IDH  1
AK213   H016-C9EF5G IDH  1
AK231   H016-K82Q   IDH  1
AK005   H016-6BP868 MES  2
AK006   H016-RKB4RB MES  2
AK030   H016-AZH7   MES  1
AK055   H016-U9HSNM MES  1
AK071   H016-TELN6S MES  1
AK072   H016-STUK   MES  1
AK079   H016-N6KNX8 MES  2
AK081   H016-KG8EA4 MES  2
AK088   H016-K48U   MES  1
AK091   H016-DD22   MES  1
AK134   H016-7CCGYW MES  2
AK139   H016-V41MG6 MES  1
AK153   H016-2BH85A MES  1
AK185   H016-1SDSG7 MES  1
AK188   H016-XXJFNP MES  1
AK195   H016-77FF   MES  1
AK218   H016-82ZL4S MES  2
AK227   H016-H1M4FV MES  1
AK235   H016-7JVLAT MES  2
AK236   H016-U676   MES  1
AK256   H016-V1A93Q MES  2
AK002   H016-8S2Z4Z RTK_I  2
AK003   H016-KBJ2J5 RTK_I  1
AK043   H016-F7WG7L RTK_I  2
AK049   H016-YKZ5   RTK_I  1
AK051   H016-AYDUQX RTK_I  1
AK142   H016-6L2VZW RTK_I  1
AK149   H016-9JGR   RTK_I  1
AK156   H016-3LK6   RTK_I  1
AK165   H016-59ND   RTK_I  1
AK173   H016-9GTQ8S RTK_I  1
AK183   H016-D3H7   RTK_I  1
AK217   H016-YN4UXM RTK_I  1
AK035   H016-1AG619 RTK_II  2
AK053   H016-AZJVFM RTK_II  1
AK074   H016-4F2ZQN RTK_II  1
AK089   H016-761V   RTK_II  1
AK098   H016-BRD3LH RTK_II  1
AK099   H016-8DFN1P RTK_II  2
AK100   H016-TTQXDW RTK_II  1
AK117   H016-TCS133 RTK_II  2
AK123   H016-PNVGQB RTK_II  2
AK132   H016-L5XL   RTK_II  1
AK133   H016-9WEKX1 RTK_II  2
AK158   H016-DUHE   RTK_II  1
AK167   H016-LEZR   RTK_II  1
AK178   H016-F6J1VJ RTK_II  1
AK205   H016-K1RYMM RTK_II  2
AK216   H016-VNDF   RTK_II  1
AK226   H016-3ZCL2Y RTK_II  2
"

SAMPLE = read.table(textConnection(anno_text), header = TRUE, stringsAsFactors = FALSE)
SAMPLE = SAMPLE[SAMPLE$subtype %in% subtype, ]
rownames(SAMPLE) = SAMPLE$id
SAMPLE_ID = SAMPLE$id
SUBTYPE_COLOR = RColorBrewer::brewer.pal(4, "Set1")
names(SUBTYPE_COLOR) = c("IDH", "MES", "RTK_I", "RTK_II")

PROJECT_DIR = "/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort"
library(GenomicFeatures)
TXDB = loadDb(qq("@{PROJECT_DIR}/txdb/gencode19_protein_coding_txdb.sqlite"))
GENE = genes(TXDB)
GTF = qq("@{PROJECT_DIR}/txdb/gencode.v19.annotation.gtf")

####################################################
## expression data
load(qq("@{PROJECT_DIR}/expression/hipo16_rnaseq_count_rpkm.RData"))
count = count[names(GENE), SAMPLE_ID]
rpkm = rpkm[names(GENE), SAMPLE_ID]
data = log2(rpkm + 1)
data = adjust_matrix(data)
res = run_all_consensus_partition_methods(data, top_n = c(1000, 2000, 4000), k = 2:6, p_sampling = p, 
	known_anno = data.frame(subtype = SAMPLE$subtype), 
	known_col = list(subtype = SUBTYPE_COLOR[unique(SAMPLE$subtype)]), 
	mc.cores = ncore)

saveRDS(res, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/hipo16_rnaseq_subgroup_p@{p}_@{paste(subtype, collapse = '')}.rds"))

res = hierarchical_partition(data, top_n = c(1000, 2000, 4000), 
	known_anno = data.frame(subtype = SAMPLE$subtype), 
	known_col = list(subtype = SUBTYPE_COLOR[unique(SAMPLE$subtype)]))

saveRDS(res, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/hipo16_rnaseq_subgroup_hierarchical_partition_@{paste(subtype, collapse = '')}.rds"))

# for(p in c(0.2, 0.4, 0.6, 0.8)) {
# 	cmd = qq("Rscript-3.1.2 /home/guz/project/development/subgroup/test_hipo16_rnaseq.R --p @{p} --ncore 4")
# 	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=20G,nodes=1:ppn=4 -N hipo16_rnaseq_subgroup_p@{p}' '@{cmd}'")
# 	system(cmd)
# }
