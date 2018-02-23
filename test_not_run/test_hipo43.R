library(methods)
library(GetoptLong)
ncore = 1
GetoptLong(
    "ncore=i", "mc.cores"
)

# anno = read.table(textConnection(
# "PID    purity_RNA  purity_WGS  Epic_subtype
# H043-5FM91P_tumor 0.746   0.74    IDH
# H043-5FM91P_tumor02 0.968   0.93    IDH
# H043-6FRXV9_tumor 0.968   0.86    IDH
# H043-6FRXV9_tumor02 0.870   0.3 IDH
# H043-CWJQEW_tumor 0.955   0.84    IDH
# H043-CWJQEW_tumor02 0.960   0.86    IDH
# H043-GG5L52_tumor 0.976   0.89    IDH
# H043-GG5L52_tumor02 0.891   0.79    IDH
# H043-KE3H42_tumor 0.883   0.82    IDH
# H043-KE3H42_tumor02 0.794   0.78    IDH
# H043-LQDGD4_tumor 0.945   0.74    IDH
# H043-LQDGD4_tumor02 0.984   0.82    IDH
# H043-X33N8G_tumor 0.914   0.8 IDH
# H043-X33N8G_tumor02 0.887   0.92    IDH
# H043-ZK2PKS_tumor 0.960   0.89    IDH
# H043-ZK2PKS_tumor02 0.924   0.87    IDH
# H043-4PGF_tumor   0.877   0.81    MES
# H043-5VWP_tumor   0.556   0.34    MES
# H043-5VWP_tumor02   0.822   0.75    MES
# H043-63R6_tumor   0.814   0.75    MES
# H043-6F91_tumor   0.971   0.84    MES
# H043-6F91_tumor02   0.821   0.97    MES
# H043-BU96_tumor   0.756   0.83    MES
# H043-BU96_tumor02   0.754   0.81    MES
# H043-GESMJV_tumor 0.858   0.54    MES
# H043-GESMJV_tumor02 0.800   0.51    MES
# H043-GKS176_tumor 0.717   0.92    MES
# H043-GKS176_tumor02 0.507   0.45    MES
# H043-LNWEGT_tumor02 0.683   0.73    MES
# H043-ULLV_tumor   0.506   0.35    MES
# H043-W99H5K_tumor 0.973   0.89    MES
# H043-W99H5K_tumor02 0.736   0.58    MES
# H043-XG4KY2_tumor 0.769   0.85    MES
# H043-XG4KY2_tumor02 0.865   0.92    MES
# H043-4PGF_tumor02   0.558   0.21    RTK1
# H043-63R6_tumor02   0.646   0.48    RTK1
# H043-LNWEGT_tumor 0.775   0.52    RTK1
# H043-MXE7Y8_tumor02 0.816   0.42    RTK1
# H043-U65X_tumor   0.549   0.35    RTK1
# H043-U65X_tumor02   0.828   0.22    RTK1
# H043-ULLV_tumor02   0.796   0.18    RTK1
# H043-ZMHY_tumor   0.878   0.65    RTK1
# H043-ZMHY_tumor02   0.947   0.62    RTK1
# H043-D9MRCY_tumor 0.864   0.85    RTK2
# H043-D9MRCY_tumor02 0.916   0.9 RTK2
# H043-DSX2_tumor   0.767   0.83    RTK2
# H043-DSX2_tumor02   0.858   0.79    RTK2
# H043-MXE7Y8_tumor 0.907   0.51    RTK2
# H043-N7LCPV_tumor 0.903   0.92    RTK2
# H043-N7LCPV_tumor02 0.893   0.81    RTK2
# H043-PWC258_tumor 0.800   0.84    RTK2
# H043-PWC258_tumor02 0.804   0.72    RTK2
# H043-XACH_tumor   0.938   0.97    RTK2
# H043-XACH_tumor02   0.940   0.93    RTK2 
# "), header = TRUE, row.names = 1)

anno = read.table(textConnection(
"pid               purity  subtype
H043-5VWP_tumor02  0.75    MES
H043-5VWP_tumor   0.34    MES
H043-63R6_tumor   0.75    RTKI
H043-4PGF_tumor02   0.48    MES
H043-6F91_tumor02   0.97    RTKII
H043-63R6_tumor02   0.48    MES
H043-BU96_tumor   0.83    RTKI
H043-4PGF_tumor   0.81    RTKI
H043-DSX2_tumor   0.84    RTKII
H043-BU96_tumor02   0.81    RTKI
H043-XACH_tumor   0.97    RTKII
H043-DSX2_tumor02   0.79    RTKII
H043-XACH_tumor02   0.93    RTKII
H043-ZMHY_tumor   0.65    MES
H043-ZMHY_tumor02   0.63    RTKII
H043-6F91_tumor   0.84    RTKII
H043-LNWEGT_tumor 0.53    RTKI
H043-LNWEGT_tumor02 0.73    RTKI
H043-GESMJV_tumor 0.54    MES
H043-GKS176_tumor 0.92    RTKII
H043-D9MRCY_tumor 0.85    RTKII
H043-D9MRCY_tumor02 0.9 RTKI
H043-GESMJV_tumor02 0.51    MES
H043-GKS176_tumor02 0.45    MES
H043-N7LCPV_tumor02 0.81    RTKII
H043-PWC258_tumor02 0.72    RTKII
H043-N7LCPV_tumor 0.92    RTKII
H043-PWC258_tumor 0.84    RTKII
H043-XG4KY2_tumor 0.85    RTKII
H043-W99H5K_tumor02 0.58    RTKI
H043-XG4KY2_tumor02 0.92    RTKI
H043-W99H5K_tumor 0.89    MES
"), header = TRUE, row.names = 1)


source(qq("@{dirname(get_scriptname())}/../load.R"))
# library(cola)
register_top_value_fun(AAC = function(mat) AAC(t(mat), mc.cores = ncore))

# load("/icgc/dkfzlsdf/analysis/cnag/cnag_MCF10CA_scRNAseq_gencode19_expression.RData")
load("/icgc/dkfzlsdf/analysis/hipo/hipo_043/RNAseq/hipo_043_rnaseq_gencode19_expression.RData")

library(GenomicFeatures)
TXDB = loadDb(qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/txdb/gencode19_protein_coding_txdb.sqlite"))
GENE = genes(TXDB)

rpkm = as.matrix(expression$rpkm)
count = as.matrix(expression$count)

rpkm = rpkm[names(GENE), rownames(anno)]
count = count[names(GENE), rownames(anno)]
l = apply(count, 1, function(x) sum(x > 0)/length(x) > 0.5)
data = log10(rpkm[l, ] + 1)

data = adjust_matrix(data)
res = run_all_consensus_partition_methods(data, mc.cores = ncore, known_anno = anno)

saveRDS(res, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/cola_test/hipo43_rnaseq_subgroup_2018_01_25.rds"))

    # cmd = qq("Rscript-3.3.1 /home/guz/project/development/cola/test_not_run/test_hipo43.R --ncore 4")
    # cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=20G,nodes=1:ppn=4 -N hipo43_rnaseq_subgroup' '@{cmd}'")
    # system(cmd)
