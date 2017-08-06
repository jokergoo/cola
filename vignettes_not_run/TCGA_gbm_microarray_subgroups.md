---
title: "cola"

output:
  html_document:
    css: "style.css"
---


```r
library(cola)
```

Data is from https://tcga-data.nci.nih.gov/docs/publications/gbm_exp/.


```r
data = read.table("~/analysis/unifiedScaled.txt", 
	header = TRUE, row.names = 1, check.names = FALSE)
data = as.matrix(data)

subtype = read.table("~/analysis/TCGA_unified_CORE_ClaNC840.txt", 
	sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
subtype = structure(unlist(subtype[1, -(1:2)]), names = colnames(subtype)[-(1:2)])

data = data[, names(subtype)]
dim(data)
```

```
## [1] 11861   173
```

```r
table(subtype)
```

```
## subtype
##   Classical Mesenchymal      Neural   Proneural 
##          38          56          26          53
```

Get all supported top methods and partition methods:


```r
all_top_value_methods()
```

```
## [1] "sd"  "vc"  "MAD" "AAC"
```

```r
all_partition_methods()
```

```
## [1] "hclust"  "kmeans"  "skmeans" "pam"     "mclust"  "som"
```


```r
register_top_value_fun(AAC = function(mat) AAC(t(mat), mc.cores = 4))
```

Run clustering for all combination of methods in batch:


```r
res = run_all_consensus_partition_methods(data, top_n = c(1000, 2000, 4000), k = 2:6, mc.cores = 4,
	known_anno = data.frame(subtype = subtype), 
	known_col = list(subtype = structure(seq_len(4), names = unique(subtype))))
```


```r
res_list = readRDS("~/analysis/TCGA_subgroup_p0.8.rds")
res_list
```

```
## Top rows are extracted by 'sd, vc, MAD, AAC' methods.
## Subgroups are detected by 'hclust, kmeans, skmeans, pam, mclust, som' method.
## Number of partitions are tried for k = 2, 3, 4, 5, 6
## 
## Following methods can be applied to this 'ConsensusPartitionList' object:
##  [1] "collect_classes"       "collect_plots"        
##  [3] "get_best_k"            "get_class"            
##  [5] "get_single_run"        "get_stat"             
##  [7] "show"                  "test_to_known_factors"
##  [9] "top_rows_heatmap"      "top_rows_overlap"
```


```r
get_best_k(res_list)
```

```
##             best_k
## sd:hclust        6
## sd:kmeans        4
## sd:skmeans       3
## sd:pam           2
## sd:mclust        2
## sd:som           5
## vc:hclust        5
## vc:kmeans        2
## vc:skmeans       2
## vc:pam           3
## vc:mclust        2
## vc:som           6
## MAD:hclust       6
## MAD:kmeans       4
## MAD:skmeans      3
## MAD:pam          2
## MAD:mclust       2
## MAD:som          5
## AAC:hclust       3
## AAC:kmeans       2
## AAC:skmeans      2
## AAC:pam          2
## AAC:mclust       4
## AAC:som          2
```

Collect all plots for a k:


```r
collect_plots(res_list, k = 4, fun = plot_ecdf)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

```r
collect_plots(res_list, k = 4, fun = consensus_heatmap)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-2.png)

```r
collect_plots(res_list, k = 4, fun = membership_heatmap)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-3.png)

```r
# collect_plots(res_list, k = 3, fun = get_signatures)
```


```r
get_stat(res_list, k = 4)
```

```
##               cophcor        PAC mean_silhouette tot_withinss
## sd:skmeans  0.9902906 0.03065622       0.9461292     515908.2
## vc:skmeans  0.8879438 0.21614825       0.5618254     574878.4
## MAD:skmeans 0.9917860 0.04217366       0.9408289     515334.7
## AAC:skmeans 0.9838451 0.04969188       0.9055829     480082.7
## sd:mclust   0.9504076 0.20431452       0.7747583     530099.1
## vc:mclust   0.8740280 0.33580814       0.4834314     578804.4
## MAD:mclust  0.9728579 0.13903732       0.8329616     524471.1
## AAC:mclust  0.9507134 0.20952895       0.7212143     487082.8
## sd:som      0.9134786 0.28883730       0.6241455     519812.2
## vc:som      0.8246908 0.40207343       0.4806463     580457.5
## MAD:som     0.9151461 0.28035670       0.6364929     515681.8
## AAC:som     0.8799678 0.29925462       0.5662937     477219.6
## sd:pam      0.9597767 0.09899116       0.8184894     531469.5
## vc:pam      0.9548823 0.16490125       0.7773808     583282.7
## MAD:pam     0.9617178 0.09382500       0.8261648     527140.2
## AAC:pam     0.9644622 0.08439025       0.8551892     485802.1
## sd:kmeans   0.9874768 0.04465861       0.9148519     514742.7
## vc:kmeans   0.8626394 0.28486387       0.5165131     575378.7
## MAD:kmeans  0.9939087 0.03172691       0.9466794     513463.6
## AAC:kmeans  0.9832365 0.07962207       0.8755590     472618.4
## sd:hclust   0.8902708 0.30583397       0.5901399     520961.7
## vc:hclust   0.8412445 0.37525236       0.5039156     585461.5
## MAD:hclust  0.8718519 0.32276054       0.5544088     522535.7
## AAC:hclust  0.8989310 0.30475330       0.4974731     494902.5
##             area_increased
## sd:skmeans       0.1227831
## vc:skmeans       0.1212277
## MAD:skmeans      0.1233362
## AAC:skmeans      0.1490848
## sd:mclust        0.1567642
## vc:mclust        0.3533911
## MAD:mclust       0.1869009
## AAC:mclust       0.1602564
## sd:som           0.1793262
## vc:som           0.2517942
## MAD:som          0.1780210
## AAC:som          0.1317615
## sd:pam           0.1514919
## vc:pam           0.1767540
## MAD:pam          0.1422933
## AAC:pam          0.1378385
## sd:kmeans        0.1329748
## vc:kmeans        0.1152486
## MAD:kmeans       0.1296491
## AAC:kmeans       0.1255296
## sd:hclust        0.1974301
## vc:hclust        0.2258603
## MAD:hclust       0.1638498
## AAC:hclust       0.1462913
```


```r
collect_classes(res_list, k = 4)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)

Overlap of top rows in different top methods:


```r
par(mfrow = c(1, 3))
top_rows_overlap(res_list, top_n = 1000)
```

```
## Loading required namespace: venneuler
```

```r
top_rows_overlap(res_list, top_n = 2000)
top_rows_overlap(res_list, top_n = 4000)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)

Also visualize the correspondance of rankings between different scoreing methods:


```r
top_rows_overlap(res_list, top_n = 1000, type = "correspondance")
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)

Heatmaps for the top rows:


```r
top_rows_heatmap(res_list, top_n = 1000)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)

Get clustering in a specified combination of top method and partition method:


```r
res = get_single_run(res_list, top_method = "MAD", partition_method = "kmeans")
res
```

```
## A 'ConsensusPartition' object with k = 2, 3, 4, 5, 6.
##   top rows (1000, 2000, 4000) are extracted by 'MAD' method.
##   subgroups are detected by 'kmeans' method.
##   best k for subgroups seems to be 4.
## 
## Following methods can be applied to this 'ConsensusPartition' object:
##  [1] "collect_classes"         "collect_plots"          
##  [3] "consensus_heatmap"       "dimension_reduction"    
##  [5] "get_best_k"              "get_class"              
##  [7] "get_consensus"           "get_membership"         
##  [9] "get_param"               "get_signatures"         
## [11] "get_stat"                "membership_heatmap"     
## [13] "plot_ecdf"               "select_partition_number"
## [15] "show"                    "signature_density"      
## [17] "test_to_known_factors"
```

Collect all plots


```r
collect_plots(res)
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15-1.png)

plots:


```r
select_partition_number(res)
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)


```r
get_best_k(res)
```

```
## [1] 4
```

```r
consensus_heatmap(res, k = 4)
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17-1.png)

```r
membership_heatmap(res, k = 4)
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17-2.png)

```r
# get_signatures(res, k = 4)
```

Get classifications


```r
get_class(res, k = 4)
```

```
##                     class    entropy silhouette
## TCGA-02-0003-01A-01     2 0.26468043  0.8602827
## TCGA-02-0010-01A-01     2 0.00000000  0.9881322
## TCGA-02-0011-01B-01     2 0.02888898  0.9846284
## TCGA-02-0014-01A-01     2 0.00000000  0.9881322
## TCGA-02-0024-01B-01     2 0.00000000  0.9881322
## TCGA-02-0026-01B-01     2 0.00000000  0.9881322
## TCGA-02-0028-01A-01     2 0.00000000  0.9881322
## TCGA-02-0046-01A-01     2 0.02888898  0.9832806
## TCGA-02-0047-01A-01     4 0.02888898  0.9652453
## TCGA-02-0048-01A-01     2 0.07072027  0.9707940
## TCGA-02-0060-01A-01     4 0.17667967  0.9178131
## TCGA-02-0069-01A-01     2 0.00000000  0.9881322
## TCGA-02-0074-01A-01     2 0.00000000  0.9881322
## TCGA-02-0080-01A-01     2 0.02888898  0.9846284
## TCGA-02-0084-01A-03     4 0.13603239  0.9278645
## TCGA-02-0087-01A-01     4 0.25494457  0.8708645
## TCGA-02-0104-01A-01     2 0.00000000  0.9881322
## TCGA-02-0114-01A-01     2 0.00000000  0.9881322
## TCGA-02-0281-01A-01     2 0.00000000  0.9881322
## TCGA-02-0321-01A-01     3 0.36096405  0.7797467
## TCGA-02-0325-01A-01     2 0.11745329  0.9608990
## TCGA-02-0338-01A-01     2 0.00000000  0.9881322
## TCGA-02-0339-01A-01     2 0.00000000  0.9881322
## TCGA-02-0432-01A-02     2 0.02888898  0.9846284
## TCGA-02-0439-01A-01     4 0.00000000  0.9618755
## TCGA-02-0440-01A-01     2 0.00000000  0.9881322
## TCGA-02-0446-01A-01     4 0.02888898  0.9652453
## TCGA-06-0128-01A-01     4 0.26468043  0.8631171
## TCGA-06-0129-01A-01     2 0.00000000  0.9881322
## TCGA-06-0146-01A-01     2 0.02888898  0.9846284
## TCGA-06-0156-01A-01     4 0.02888898  0.9652453
## TCGA-06-0166-01A-01     4 0.02888898  0.9652453
## TCGA-06-0174-01A-01     2 0.00000000  0.9881322
## TCGA-06-0177-01A-01     2 0.00000000  0.9881322
## TCGA-06-0238-01A-02     4 0.47902101  0.3934988
## TCGA-06-0241-01A-02     2 0.00000000  0.9881322
## TCGA-06-0410-01A-01     2 0.00000000  0.9881322
## TCGA-06-0413-01A-01     2 0.00000000  0.9881322
## TCGA-06-0414-01A-01     2 0.00000000  0.9881322
## TCGA-06-0646-01A-01     4 0.02888898  0.9652453
## TCGA-06-0648-01A-01     2 0.00000000  0.9881322
## TCGA-08-0245-01A-01     2 0.00000000  0.9881322
## TCGA-08-0344-01A-01     2 0.00000000  0.9881322
## TCGA-08-0347-01A-01     4 0.02888898  0.9652453
## TCGA-08-0348-01A-01     2 0.37390308  0.7220132
## TCGA-08-0350-01A-01     2 0.02888898  0.9846284
## TCGA-08-0353-01A-01     3 0.02888898  0.9524119
## TCGA-08-0359-01A-01     4 0.00000000  0.9618755
## TCGA-08-0385-01A-01     2 0.00000000  0.9881322
## TCGA-08-0517-01A-01     2 0.00000000  0.9881322
## TCGA-08-0524-01A-01     2 0.00000000  0.9881322
## TCGA-12-0616-01A-01     2 0.00000000  0.9881322
## TCGA-12-0618-01A-01     2 0.00000000  0.9881322
## TCGA-02-0089-01A-01     4 0.02888898  0.9652453
## TCGA-02-0113-01A-01     3 0.39188847  0.7276572
## TCGA-02-0115-01A-01     4 0.02888898  0.9652453
## TCGA-02-0451-01A-01     4 0.02888898  0.9652453
## TCGA-06-0132-01A-02     4 0.02888898  0.9652453
## TCGA-06-0133-01A-02     4 0.02888898  0.9652453
## TCGA-06-0138-01A-02     4 0.02888898  0.9652453
## TCGA-06-0160-01A-01     4 0.17667967  0.9178131
## TCGA-06-0162-01A-01     4 0.02888898  0.9652453
## TCGA-06-0167-01A-01     4 0.17667967  0.9178131
## TCGA-06-0171-01A-02     4 0.02888898  0.9652453
## TCGA-06-0173-01A-01     4 0.02888898  0.9652453
## TCGA-06-0179-01A-02     4 0.02888898  0.9652453
## TCGA-06-0182-01A-01     3 0.39188847  0.7276572
## TCGA-06-0185-01A-01     3 0.39188847  0.7276572
## TCGA-06-0195-01B-01     4 0.07990323  0.9583896
## TCGA-06-0208-01B-01     4 0.02888898  0.9652453
## TCGA-06-0214-01A-02     4 0.02888898  0.9652453
## TCGA-06-0219-01A-01     4 0.02888898  0.9652453
## TCGA-06-0221-01A-01     4 0.18912082  0.9119053
## TCGA-06-0237-01A-02     4 0.02888898  0.9652453
## TCGA-06-0240-01A-02     4 0.00000000  0.9618755
## TCGA-08-0349-01A-01     4 0.02888898  0.9652453
## TCGA-08-0380-01A-01     4 0.02888898  0.9652453
## TCGA-08-0386-01A-01     3 0.41832037  0.6771913
## TCGA-08-0520-01A-01     4 0.34003852  0.7616114
## TCGA-02-0007-01A-01     3 0.34722304  0.7871946
## TCGA-02-0009-01A-01     3 0.02888898  0.9552673
## TCGA-02-0016-01A-01     3 0.02888898  0.9552673
## TCGA-02-0021-01A-01     3 0.02888898  0.9552673
## TCGA-02-0023-01B-01     3 0.02888898  0.9552673
## TCGA-02-0027-01A-01     3 0.08869473  0.9424507
## TCGA-02-0038-01A-01     3 0.49479376  0.2948284
## TCGA-02-0043-01A-01     3 0.02888898  0.9537879
## TCGA-02-0070-01A-01     3 0.02888898  0.9537879
## TCGA-02-0102-01A-01     3 0.02888898  0.9537879
## TCGA-02-0260-01A-03     3 0.05774568  0.9492765
## TCGA-02-0269-01B-01     3 0.05107902  0.9527670
## TCGA-02-0285-01A-01     3 0.02888898  0.9537879
## TCGA-02-0289-01A-01     3 0.02888898  0.9552673
## TCGA-02-0290-01A-01     3 0.02888898  0.9537879
## TCGA-02-0317-01A-01     3 0.05107902  0.9527670
## TCGA-02-0333-01A-02     3 0.05107902  0.9527670
## TCGA-02-0422-01A-01     3 0.02888898  0.9552673
## TCGA-02-0430-01A-01     3 0.02888898  0.9552673
## TCGA-06-0125-01A-01     3 0.02888898  0.9552673
## TCGA-06-0126-01A-01     3 0.02888898  0.9552673
## TCGA-06-0137-01A-03     3 0.02888898  0.9552673
## TCGA-06-0145-01A-04     3 0.02888898  0.9552673
## TCGA-06-0148-01A-01     3 0.02888898  0.9537879
## TCGA-06-0187-01A-01     3 0.02888898  0.9537879
## TCGA-06-0211-01B-01     3 0.02888898  0.9552673
## TCGA-06-0402-01A-01     3 0.02888898  0.9537879
## TCGA-08-0246-01A-01     3 0.02888898  0.9537879
## TCGA-08-0354-01A-01     3 0.02888898  0.9537879
## TCGA-08-0355-01A-01     3 0.02888898  0.9552673
## TCGA-08-0357-01A-01     3 0.02888898  0.9552673
## TCGA-08-0358-01A-01     3 0.05774568  0.9492765
## TCGA-08-0375-01A-01     3 0.02888898  0.9552673
## TCGA-08-0511-01A-01     3 0.02888898  0.9552673
## TCGA-08-0514-01A-01     3 0.02888898  0.9537879
## TCGA-08-0518-01A-01     3 0.02888898  0.9537879
## TCGA-08-0529-01A-02     3 0.02888898  0.9552673
## TCGA-08-0531-01A-01     3 0.02888898  0.9552673
## TCGA-02-0057-01A-01     4 0.05774568  0.9604716
## TCGA-02-0004-01A-01     1 0.02888898  0.9809814
## TCGA-02-0006-01B-01     1 0.34003852  0.8003203
## TCGA-02-0025-01A-01     1 0.02888898  0.9809814
## TCGA-02-0033-01A-01     1 0.02888898  0.9809814
## TCGA-02-0034-01A-01     1 0.02888898  0.9809814
## TCGA-02-0039-01A-01     1 0.07072027  0.9699716
## TCGA-02-0051-01A-01     1 0.02888898  0.9809814
## TCGA-02-0054-01A-01     1 0.34520223  0.8162763
## TCGA-02-0055-01A-01     1 0.02888898  0.9809814
## TCGA-02-0059-01A-01     1 0.02888898  0.9809814
## TCGA-02-0064-01A-01     1 0.00000000  0.9815527
## TCGA-02-0075-01A-01     1 0.00000000  0.9815527
## TCGA-02-0079-01A-03     1 0.00000000  0.9815527
## TCGA-02-0085-01A-01     4 0.02888898  0.9606365
## TCGA-02-0086-01A-01     1 0.00000000  0.9815527
## TCGA-02-0099-01A-01     1 0.32884118  0.8309643
## TCGA-02-0106-01A-01     1 0.02888898  0.9809814
## TCGA-02-0107-01A-01     1 0.00000000  0.9815527
## TCGA-02-0111-01A-01     1 0.00000000  0.9815527
## TCGA-02-0326-01A-01     3 0.05107902  0.9527670
## TCGA-02-0337-01A-01     1 0.17667967  0.9325983
## TCGA-06-0122-01A-01     1 0.00000000  0.9815527
## TCGA-06-0124-01A-01     1 0.00000000  0.9815527
## TCGA-06-0130-01A-01     1 0.02888898  0.9809814
## TCGA-06-0139-01A-01     1 0.02888898  0.9809814
## TCGA-06-0143-01A-01     1 0.00000000  0.9815527
## TCGA-06-0147-01A-01     1 0.02888898  0.9809814
## TCGA-06-0149-01A-05     1 0.07072027  0.9699716
## TCGA-06-0152-01A-02     3 0.05774568  0.9506049
## TCGA-06-0154-01A-02     1 0.00000000  0.9815527
## TCGA-06-0164-01A-01     1 0.00000000  0.9815527
## TCGA-06-0175-01A-01     1 0.20523571  0.9269461
## TCGA-06-0176-01A-03     1 0.02888898  0.9809814
## TCGA-06-0184-01A-01     1 0.17667967  0.9325983
## TCGA-06-0189-01A-05     1 0.02888898  0.9809814
## TCGA-06-0190-01A-01     1 0.00000000  0.9815527
## TCGA-06-0194-01A-01     1 0.00000000  0.9815527
## TCGA-06-0197-01A-02     1 0.02888898  0.9809814
## TCGA-06-0210-01A-01     1 0.00000000  0.9815527
## TCGA-06-0397-01A-01     1 0.02888898  0.9809814
## TCGA-06-0409-01A-02     1 0.00000000  0.9815527
## TCGA-06-0412-01A-01     1 0.00000000  0.9815527
## TCGA-06-0644-01A-02     1 0.02888898  0.9809814
## TCGA-06-0645-01A-01     1 0.00000000  0.9815527
## TCGA-08-0346-01A-01     1 0.00000000  0.9815527
## TCGA-08-0352-01A-01     1 0.00000000  0.9815527
## TCGA-08-0360-01A-01     1 0.00000000  0.9815527
## TCGA-08-0390-01A-01     1 0.07072027  0.9699716
## TCGA-08-0392-01A-01     1 0.02888898  0.9809814
## TCGA-08-0509-01A-01     1 0.00000000  0.9815527
## TCGA-08-0510-01A-01     1 0.00000000  0.9815527
## TCGA-08-0512-01A-01     1 0.02888898  0.9809814
## TCGA-08-0522-01A-01     1 0.02888898  0.9809814
## TCGA-12-0619-01A-01     1 0.00000000  0.9815527
## TCGA-12-0620-01A-01     1 0.17667967  0.9233442
```

MDS or T-sne plots:


```r
dimension_reduction(res, k = 4)
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png)

```r
dimension_reduction(res, k = 4, method = "tsne")
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-2.png)

Consistency of classes.


```r
collect_classes(res_list, k = 4)
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20-1.png)

```r
collect_classes(res)
```


```r
res = hierarchical_partition(data, top_n = c(1000, 2000, 4000), 
	known_anno = data.frame(subtype = subtype), 
	known_col = list(subtype = structure(seq_len(4), names = unique(subtype))))
```


```r
res = readRDS("~/analysis/TCGA_subgroup_hierarchical_partition.rds")
res
```

```
## A 'HierarchicalPartition' object with 'MAD:kmeans' method.
## 
## +-- 01, 52 cols
## |   |-- 011, 37 cols
## |   +-- 010, 15 cols
## +-- 00, 121 cols
##     |-- 001, 46 cols
##     +-- 000, 75 cols
##         |-- 0001, 37 cols
##         +-- 0000, 38 cols
##             |-- 00001, 16 cols
##             +-- 00000, 22 cols
##                 |-- 000001, 9 cols
##                 +-- 000000, 13 cols
## 
## Following methods can be applied to this 'HierarchicalPartition' object:
## [1] "collect_classes"       "get_class"             "get_signatures"       
## [4] "get_single_run"        "show"                  "test_to_known_factors"
```


```r
collect_classes(res)
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23-1.png)


```r
get_class(res)
```

```
## TCGA-02-0003-01A-01 TCGA-02-0010-01A-01 TCGA-02-0011-01B-01 
##             "00001"            "000000"            "000001" 
## TCGA-02-0014-01A-01 TCGA-02-0024-01B-01 TCGA-02-0026-01B-01 
##            "000000"            "000000"            "000001" 
## TCGA-02-0028-01A-01 TCGA-02-0046-01A-01 TCGA-02-0047-01A-01 
##            "000000"             "00001"              "0001" 
## TCGA-02-0048-01A-01 TCGA-02-0060-01A-01 TCGA-02-0069-01A-01 
##             "00001"              "0001"            "000000" 
## TCGA-02-0074-01A-01 TCGA-02-0080-01A-01 TCGA-02-0084-01A-03 
##             "00001"            "000001"              "0001" 
## TCGA-02-0087-01A-01 TCGA-02-0104-01A-01 TCGA-02-0114-01A-01 
##              "0001"            "000000"            "000000" 
## TCGA-02-0281-01A-01 TCGA-02-0321-01A-01 TCGA-02-0325-01A-01 
##            "000000"               "001"            "000001" 
## TCGA-02-0338-01A-01 TCGA-02-0339-01A-01 TCGA-02-0432-01A-02 
##            "000000"            "000000"            "000001" 
## TCGA-02-0439-01A-01 TCGA-02-0440-01A-01 TCGA-02-0446-01A-01 
##              "0001"             "00001"              "0001" 
## TCGA-06-0128-01A-01 TCGA-06-0129-01A-01 TCGA-06-0146-01A-01 
##              "0001"            "000001"            "000001" 
## TCGA-06-0156-01A-01 TCGA-06-0166-01A-01 TCGA-06-0174-01A-01 
##              "0001"              "0001"             "00001" 
## TCGA-06-0177-01A-01 TCGA-06-0238-01A-02 TCGA-06-0241-01A-02 
##             "00001"              "0001"             "00001" 
## TCGA-06-0410-01A-01 TCGA-06-0413-01A-01 TCGA-06-0414-01A-01 
##             "00001"            "000000"             "00001" 
## TCGA-06-0646-01A-01 TCGA-06-0648-01A-01 TCGA-08-0245-01A-01 
##              "0001"             "00001"             "00001" 
## TCGA-08-0344-01A-01 TCGA-08-0347-01A-01 TCGA-08-0348-01A-01 
##            "000000"              "0001"             "00001" 
## TCGA-08-0350-01A-01 TCGA-08-0353-01A-01 TCGA-08-0359-01A-01 
##            "000001"               "001"              "0001" 
## TCGA-08-0385-01A-01 TCGA-08-0517-01A-01 TCGA-08-0524-01A-01 
##            "000001"             "00001"            "000000" 
## TCGA-12-0616-01A-01 TCGA-12-0618-01A-01 TCGA-02-0089-01A-01 
##             "00001"             "00001"              "0001" 
## TCGA-02-0113-01A-01 TCGA-02-0115-01A-01 TCGA-02-0451-01A-01 
##               "001"              "0001"              "0001" 
## TCGA-06-0132-01A-02 TCGA-06-0133-01A-02 TCGA-06-0138-01A-02 
##              "0001"              "0001"              "0001" 
## TCGA-06-0160-01A-01 TCGA-06-0162-01A-01 TCGA-06-0167-01A-01 
##              "0001"              "0001"              "0001" 
## TCGA-06-0171-01A-02 TCGA-06-0173-01A-01 TCGA-06-0179-01A-02 
##              "0001"              "0001"              "0001" 
## TCGA-06-0182-01A-01 TCGA-06-0185-01A-01 TCGA-06-0195-01B-01 
##               "001"               "001"              "0001" 
## TCGA-06-0208-01B-01 TCGA-06-0214-01A-02 TCGA-06-0219-01A-01 
##              "0001"              "0001"              "0001" 
## TCGA-06-0221-01A-01 TCGA-06-0237-01A-02 TCGA-06-0240-01A-02 
##              "0001"              "0001"              "0001" 
## TCGA-08-0349-01A-01 TCGA-08-0380-01A-01 TCGA-08-0386-01A-01 
##              "0001"              "0001"               "001" 
## TCGA-08-0520-01A-01 TCGA-02-0007-01A-01 TCGA-02-0009-01A-01 
##              "0001"               "001"               "001" 
## TCGA-02-0016-01A-01 TCGA-02-0021-01A-01 TCGA-02-0023-01B-01 
##               "001"               "001"               "001" 
## TCGA-02-0027-01A-01 TCGA-02-0038-01A-01 TCGA-02-0043-01A-01 
##               "001"               "001"               "001" 
## TCGA-02-0070-01A-01 TCGA-02-0102-01A-01 TCGA-02-0260-01A-03 
##               "001"               "001"               "001" 
## TCGA-02-0269-01B-01 TCGA-02-0285-01A-01 TCGA-02-0289-01A-01 
##               "001"               "001"               "001" 
## TCGA-02-0290-01A-01 TCGA-02-0317-01A-01 TCGA-02-0333-01A-02 
##               "001"               "001"               "001" 
## TCGA-02-0422-01A-01 TCGA-02-0430-01A-01 TCGA-06-0125-01A-01 
##               "001"               "001"               "001" 
## TCGA-06-0126-01A-01 TCGA-06-0137-01A-03 TCGA-06-0145-01A-04 
##               "001"               "001"               "001" 
## TCGA-06-0148-01A-01 TCGA-06-0187-01A-01 TCGA-06-0211-01B-01 
##               "001"               "001"               "001" 
## TCGA-06-0402-01A-01 TCGA-08-0246-01A-01 TCGA-08-0354-01A-01 
##               "001"               "001"               "001" 
## TCGA-08-0355-01A-01 TCGA-08-0357-01A-01 TCGA-08-0358-01A-01 
##               "001"               "001"               "001" 
## TCGA-08-0375-01A-01 TCGA-08-0511-01A-01 TCGA-08-0514-01A-01 
##               "001"               "001"               "001" 
## TCGA-08-0518-01A-01 TCGA-08-0529-01A-02 TCGA-08-0531-01A-01 
##               "001"               "001"               "001" 
## TCGA-02-0057-01A-01 TCGA-02-0004-01A-01 TCGA-02-0006-01B-01 
##              "0001"               "010"               "011" 
## TCGA-02-0025-01A-01 TCGA-02-0033-01A-01 TCGA-02-0034-01A-01 
##               "010"               "010"               "010" 
## TCGA-02-0039-01A-01 TCGA-02-0051-01A-01 TCGA-02-0054-01A-01 
##               "011"               "010"               "011" 
## TCGA-02-0055-01A-01 TCGA-02-0059-01A-01 TCGA-02-0064-01A-01 
##               "010"               "010"               "011" 
## TCGA-02-0075-01A-01 TCGA-02-0079-01A-03 TCGA-02-0085-01A-01 
##               "011"               "011"              "0001" 
## TCGA-02-0086-01A-01 TCGA-02-0099-01A-01 TCGA-02-0106-01A-01 
##               "011"               "011"               "010" 
## TCGA-02-0107-01A-01 TCGA-02-0111-01A-01 TCGA-02-0326-01A-01 
##               "011"               "011"               "001" 
## TCGA-02-0337-01A-01 TCGA-06-0122-01A-01 TCGA-06-0124-01A-01 
##               "011"               "011"               "011" 
## TCGA-06-0130-01A-01 TCGA-06-0139-01A-01 TCGA-06-0143-01A-01 
##               "010"               "010"               "011" 
## TCGA-06-0147-01A-01 TCGA-06-0149-01A-05 TCGA-06-0152-01A-02 
##               "011"               "011"               "001" 
## TCGA-06-0154-01A-02 TCGA-06-0164-01A-01 TCGA-06-0175-01A-01 
##               "011"               "011"               "011" 
## TCGA-06-0176-01A-03 TCGA-06-0184-01A-01 TCGA-06-0189-01A-05 
##               "010"               "011"               "010" 
## TCGA-06-0190-01A-01 TCGA-06-0194-01A-01 TCGA-06-0197-01A-02 
##               "011"               "011"               "011" 
## TCGA-06-0210-01A-01 TCGA-06-0397-01A-01 TCGA-06-0409-01A-02 
##               "011"               "011"               "011" 
## TCGA-06-0412-01A-01 TCGA-06-0644-01A-02 TCGA-06-0645-01A-01 
##               "011"               "010"               "011" 
## TCGA-08-0346-01A-01 TCGA-08-0352-01A-01 TCGA-08-0360-01A-01 
##               "011"               "011"               "011" 
## TCGA-08-0390-01A-01 TCGA-08-0392-01A-01 TCGA-08-0509-01A-01 
##               "011"               "010"               "011" 
## TCGA-08-0510-01A-01 TCGA-08-0512-01A-01 TCGA-08-0522-01A-01 
##               "011"               "011"               "010" 
## TCGA-12-0619-01A-01 TCGA-12-0620-01A-01 
##               "011"               "011"
```

