<!--
%\VignetteEngine{knitr}
%\VignetteIndexEntry{3. Gene Ontology Enrichment to Signature Genes}
-->

Gene Ontology Enrichment to Signature Genes
=============================================================

**Author**: Zuguang Gu ( z.gu@dkfz.de )

**Date**: `r Sys.Date()`

**Package version**: `r installed.packages()["cola", "Version"]`

-------------------------------------------------------------

```{r, echo = FALSE, message = FALSE}
library(markdown)
library(knitr)
knitr::opts_chunk$set(
    error = FALSE,
    tidy  = FALSE,
    message = FALSE,
    fig.align = "center")
options(width = 100)
```


```{r, results = "hide", eval = FALSE}
library(cola)
download.file("https://github.com/jokergoo/cola_examples/raw/master/TCGA_GBM/TCGA_GBM_subgroup.rds", 
    destfile = "TCGA_GBM_subgroup.rds", quiet = TRUE)
rl = readRDS("TCGA_GBM_subgroup.rds")
file.remove("TCGA_GBM_subgroup.rds")

res = rl["ATC:skmeans"]
```

```{r, eval = FALSE}
get_signatures(res, k = 4)
```

```{r, eval = FALSE}
GO_enrichment(res, k = 4, ontology = "BP")
```
