# Using fgsea package

`fgsea` is an R-package for fast preranked gene set enrichment analysis (GSEA). This package allows to quickly and accurately calculate arbitrarily low GSEA P-values for a collection of gene sets. P-value estimation is based on an adaptive multi-level split Monte-Carlo scheme. See the [preprint](https://www.biorxiv.org/content/10.1101/060012v2) for algorithmic details.

## Loading necessary libraries

```
library(fgsea)
library(data.table)
library(ggplot2)
```

## Quick run

Loading example pathways and gene-level statistics and setting random seed:

```
data(examplePathways)
data(exampleRanks)
set.seed(42)
```

Running fgsea:

```
fgseaRes <- fgsea(pathways = examplePathways,
                  stats    = exampleRanks,
                  minSize  = 15,
                  maxSize  = 500)
```

The resulting table contains enrichment scores and p-values:

```
head(fgseaRes[order(pval), ])
```

    ##                                            pathway         pval         padj
    ##                                             <char>        <num>        <num>
    ## 1:                     5990979_Cell_Cycle,_Mitotic 6.690481e-27 3.920622e-24
    ## 2:                              5990980_Cell_Cycle 3.312565e-26 9.705816e-24
    ## 3:                    5991851_Mitotic_Prometaphase 8.470173e-19 1.654507e-16
    ## 4: 5992217_Resolution_of_Sister_Chromatid_Cohesion 2.176649e-18 3.188791e-16
    ## 5:                                 5991454_M_Phase 1.873997e-14 2.196325e-12
    ## 6:         5991599_Separation_of_Sister_Chromatids 8.733223e-14 8.529448e-12
    ##      log2err        ES      NES  size
    ##        <num>     <num>    <num> <int>
    ## 1: 1.3422338 0.5594755 2.769070   317
    ## 2: 1.3267161 0.5388497 2.705894   369
    ## 3: 1.1239150 0.7253270 2.972690    82
    ## 4: 1.1053366 0.7347987 2.957518    74
    ## 5: 0.9759947 0.5576247 2.554076   173
    ## 6: 0.9545416 0.6164600 2.670030   116
    ##                                      leadingEdge
    ##                                           <list>
    ## 1: 66336,66977,12442,107995,66442,12571,...[141]
    ## 2: 66336,66977,12442,107995,66442,19361,...[155]
    ## 3:  66336,66977,12442,107995,66442,52276,...[44]
    ## 4:  66336,66977,12442,107995,66442,52276,...[40]
    ## 5:  66336,66977,12442,107995,66442,52276,...[67]
    ## 6:  66336,66977,107995,66442,52276,67629,...[50]

As you can see from the warning, `fgsea` has a default lower bound `eps=1e-10` for estimating P-values. If you need to estimate P-value more accurately, you can set the `eps` argument to zero in the `fgsea` function.

```
fgseaRes <- fgsea(pathways = examplePathways,
                  stats    = exampleRanks,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)

head(fgseaRes[order(pval), ])
```

    ##                                            pathway         pval         padj
    ##                                             <char>        <num>        <num>
    ## 1:                              5990980_Cell_Cycle 2.535645e-26 1.485888e-23
    ## 2:                     5990979_Cell_Cycle,_Mitotic 9.351994e-26 2.740134e-23
    ## 3:                    5991851_Mitotic_Prometaphase 3.633805e-19 7.098033e-17
    ## 4: 5992217_Resolution_of_Sister_Chromatid_Cohesion 2.077985e-17 3.044248e-15
    ## 5:                                 5991454_M_Phase 2.251818e-14 2.639131e-12
    ## 6:          5991502_Mitotic_Metaphase_and_Anaphase 3.196758e-14 3.122167e-12
    ##      log2err        ES      NES  size
    ##        <num>     <num>    <num> <int>
    ## 1: 1.3344975 0.5388497 2.664606   369
    ## 2: 1.3188888 0.5594755 2.740246   317
    ## 3: 1.1330899 0.7253270 2.926512    82
    ## 4: 1.0768682 0.7347987 2.920436    74
    ## 5: 0.9759947 0.5576247 2.547515   173
    ## 6: 0.9653278 0.6052907 2.639370   123
    ##                                      leadingEdge
    ##                                           <list>
    ## 1: 66336,66977,12442,107995,66442,19361,...[155]
    ## 2: 66336,66977,12442,107995,66442,12571,...[141]
    ## 3:  66336,66977,12442,107995,66442,52276,...[44]
    ## 4:  66336,66977,12442,107995,66442,52276,...[40]
    ## 5:  66336,66977,12442,107995,66442,52276,...[67]
    ## 6:  66336,66977,107995,66442,52276,67629,...[51]

One can make an enrichment plot for a pathway:

```
plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
               exampleRanks) + labs(title="Programmed Cell Death")
```

Or make a table plot for a bunch of selected pathways:

```
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes,
              gseaParam=0.5)
```

From the plot above one can see that there are very similar pathways in the table (for example `5991502_Mitotic_Metaphase_and_Anaphase` and `5991600_Mitotic_Anaphase`). To select only independent pathways one can use `collapsePathways` function:

```
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01],
                                      examplePathways, exampleRanks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
                         order(-NES), pathway]
plotGseaTable(examplePathways[mainPathways], exampleRanks, fgseaRes,
              gseaParam = 0.5)
```

To save the results in a text format `data:table::fwrite` function can be used:

```
fwrite(fgseaRes, file="fgseaRes.txt", sep="\t", sep2=c("", " ", ""))
```

To make leading edge more human-readable it can be converted using `mapIdsList` (similar to `AnnotationDbi::mapIds`) function and a corresponding database (here `org.Mm.eg.db` for mouse):

```
library(org.Mm.eg.db)
fgseaResMain <- fgseaRes[match(mainPathways, pathway)]
fgseaResMain[, leadingEdge := mapIdsList(
                                     x=org.Mm.eg.db,
                                     keys=leadingEdge,
                                     keytype="ENTREZID",
                                     column="SYMBOL")]
fwrite(fgseaResMain, file="fgseaResMain.txt", sep="\t", sep2=c("", " ", ""))
```

## Performance considerations

Also, `fgsea` is parallelized using `BiocParallel` package. By default the first registered backend returned by `bpparam()` is used. To tweak the parallelization one can either specify `BPPARAM` parameter used for `bplapply` of set `nproc` parameter, which is a shorthand for setting `BPPARAM=MulticoreParam(workers = nproc)`.

## Using Reactome pathways

For convenience there is `reactomePathways` function that obtains pathways from Reactome for given set of genes. Package `reactome.db` is required to be installed.

```
pathways <- reactomePathways(names(exampleRanks))
fgseaRes <- fgsea(pathways, exampleRanks, maxSize=500)
head(fgseaRes)
```

    ##                                       pathway       pval      padj    log2err
    ##                                        <char>      <num>     <num>      <num>
    ## 1: 5-Phosphoribose 1-diphosphate biosynthesis 0.87473461 0.9496503 0.05456806
    ## 2:      ABC transporters in lipid homeostasis 0.22093023 0.5704830 0.12750532
    ## 3:     ABC-family proteins mediated transport 0.47394137 0.7764357 0.07130530
    ## 4:               ABO blood group biosynthesis 0.97137014 0.9893377 0.04802204
    ## 5:  ADP signalling through P2Y purinoceptor 1 0.01983994 0.1533244 0.35248786
    ## 6: ADP signalling through P2Y purinoceptor 12 0.05108055 0.2777946 0.28201335
    ##            ES        NES  size                                 leadingEdge
    ##         <num>      <num> <int>                                      <list>
    ## 1:  0.4267378  0.6782286     2                                328099,19139
    ## 2: -0.4385385 -1.2280487    12 19299,27403,11307,11806,217265,27409,...[8]
    ## 3:  0.2512744  0.9890018    57 17463,26440,26444,19179,56325,26445,...[27]
    ## 4:  0.5120427  0.6927180     1                                       14344
    ## 5:  0.6097588  1.8256797    17  14696,14702,14700,14682,14676,66066,...[7]
    ## 6:  0.5746173  1.6217227    14  14696,14702,14700,66066,14697,14693,...[8]

## Starting from files

One can also start from `.rnk` and `.gmt` files as in original GSEA:

```
rnk.file <- system.file("extdata", "naive.vs.th1.rnk", package="fgsea")
gmt.file <- system.file("extdata", "mouse.reactome.gmt", package="fgsea")
```

Loading ranks:

```
ranks <- read.table(rnk.file,
                    header=TRUE, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$t, ranks$ID)
str(ranks)
```

    ##  Named num [1:12000] -63.3 -49.7 -43.6 -41.5 -33.3 ...
    ##  - attr(*, "names")= chr [1:12000] "170942" "109711" "18124" "12775" ...

Loading pathways:

```
pathways <- gmtPathways(gmt.file)
str(head(pathways))
```

    ## List of 6
    ##  $ 1221633_Meiotic_Synapsis                                                : chr [1:64] "12189" "13006" "15077" "15078" ...
    ##  $ 1368092_Rora_activates_gene_expression                                  : chr [1:9] "11865" "12753" "12894" "18143" ...
    ##  $ 1368110_Bmal1:Clock,Npas2_activates_circadian_gene_expression           : chr [1:16] "11865" "11998" "12753" "12952" ...
    ##  $ 1445146_Translocation_of_Glut4_to_the_Plasma_Membrane                   : chr [1:55] "11461" "11465" "11651" "11652" ...
    ##  $ 186574_Endocrine-committed_Ngn3+_progenitor_cells                       : chr [1:4] "18012" "18088" "18506" "53626"
    ##  $ 186589_Late_stage_branching_morphogenesis_pancreatic_bud_precursor_cells: chr [1:4] "11925" "15205" "21410" "246086"

And running fgsea:

```
fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500)
head(fgseaRes)
```

    ##                                                                                    pathway
    ##                                                                                     <char>
    ## 1:                                                                1221633_Meiotic_Synapsis
    ## 2:                                   1445146_Translocation_of_Glut4_to_the_Plasma_Membrane
    ## 3: 442533_Transcriptional_Regulation_of_Adipocyte_Differentiation_in_3T3-L1_Pre-adipocytes
    ## 4:                                                                  508751_Circadian_Clock
    ## 5:                                               5334727_Mus_musculus_biological_processes
    ## 6:                                        573389_NoRC_negatively_regulates_rRNA_expression
    ##          pval      padj    log2err         ES        NES  size
    ##         <num>     <num>      <num>      <num>      <num> <int>
    ## 1: 0.52321429 0.6952462 0.07096095  0.2885754  0.9510908    27
    ## 2: 0.68197279 0.8274038 0.05582647  0.2387284  0.8531383    39
    ## 3: 0.07964602 0.2044229 0.23779383 -0.3640706 -1.3498065    31
    ## 4: 0.80186916 0.9001826 0.05269731  0.2516324  0.7391697    17
    ## 5: 0.34032984 0.5390089 0.08431443  0.2469065  1.0678318   106
    ## 6: 0.39252336 0.5913077 0.08862611  0.3607407  1.0596751    17
    ##                                     leadingEdge
    ##                                          <list>
    ## 1:                            15270,12189,71846
    ## 2:  17918,19341,20336,22628,22627,20619,...[10]
    ## 3: 76199,19014,26896,229003,17977,17978,...[12]
    ## 4:                                  20893,59027
    ## 5:  60406,19361,15270,20893,12189,68240,...[12]
    ## 6:                                  60406,20018

## Over-representation test

`fgsea` package also contains a function called `fora` for over-representation analysis based enrichment tests, based on the hypergeometric test.

`fora` requires a foreground set of genes of interest, a background set consisting of all robustly detected genes (also called universe), and some pathways.

In the following example, the foreground (`fg`) consists of the 500 genes with the highest `stat` value.

```
fg <-  names(head(exampleRanks[order(exampleRanks, decreasing=TRUE)],500))
bg <- names(exampleRanks)
foraRes <- fora(genes=fg, universe=bg, pathways=examplePathways)
head(foraRes)
```

    ##                                            pathway         pval         padj
    ##                                             <char>        <num>        <num>
    ## 1:                              5990980_Cell_Cycle 2.498582e-22 3.547986e-19
    ## 2:                     5990979_Cell_Cycle,_Mitotic 6.527769e-20 4.634716e-17
    ## 3:                    5991851_Mitotic_Prometaphase 1.108012e-17 5.244589e-15
    ## 4: 5992217_Resolution_of_Sister_Chromatid_Cohesion 1.096637e-15 3.893062e-13
    ## 5:                                 5991454_M_Phase 8.051967e-13 2.286759e-10
    ## 6:            5991757_RHO_GTPases_Activate_Formins 4.101960e-12 9.707972e-10
    ##    foldEnrichment overlap  size                                overlapGenes
    ##             <num>   <int> <int>                                      <list>
    ## 1:       4.097561      63   369 11799,12189,12236,12442,12448,12534,...[63]
    ## 2:       4.164038      55   317 11799,12236,12442,12448,12534,12571,...[55]
    ## 3:       7.902439      27    82 11799,12236,12442,12534,12615,14211,...[27]
    ## 4:       7.783784      24    74 11799,12236,12442,12534,12615,18817,...[24]
    ## 5:       4.439306      32   173 11799,12236,12442,12534,12615,14211,...[32]
    ## 6:       6.461538      21    78 11799,12236,12615,18817,20877,52276,...[21]

## Session info

```
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ##
    ## Matrix products: default
    ## BLAS:   /home/biocbuild/bbs-3.22-bioc/R/lib/libRblas.so
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0
    ##
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
    ##  [3] LC_TIME=en_GB              LC_COLLATE=C
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
    ##
    ## time zone: America/New_York
    ## tzcode source: system (glibc)
    ##
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods
    ## [8] base
    ##
    ## other attached packages:
    ##  [1] org.Mm.eg.db_3.22.0  AnnotationDbi_1.72.0 IRanges_2.44.0
    ##  [4] S4Vectors_0.48.0     Biobase_2.70.0       BiocGenerics_0.56.0
    ##  [7] generics_0.1.4       BiocParallel_1.44.0  ggplot2_4.0.1
    ## [10] data.table_1.18.0    fgsea_1.36.2
    ##
    ## loaded via a namespace (and not attached):
    ##  [1] sass_0.4.10        RSQLite_2.4.5      lattice_0.22-7     digest_0.6.39
    ##  [5] magrittr_2.0.4     evaluate_1.0.5     grid_4.5.2         RColorBrewer_1.1-3
    ##  [9] blob_1.2.4         fastmap_1.2.0      jsonlite_2.0.0     Matrix_1.7-4
    ## [13] DBI_1.2.3          httr_1.4.7         scales_1.4.0       Biostrings_2.78.0
    ## [17] codetools_0.2-20   jquerylib_0.1.4    cli_3.6.5          crayon_1.5.3
    ## [21] rlang_1.1.6        XVector_0.50.0     bit64_4.6.0-1      cowplot_1.2.0
    ## [25] withr_3.0.2        cachem_1.1.0       yaml_2.3.12        otel_0.2.0
    ## [29] tools_4.5.2        parallel_4.5.2     memoise_2.0.1      dplyr_1.1.4
    ## [33] reactome.db_1.94.0 fastmatch_1.1-6    png_0.1-8          vctrs_0.6.5
    ## [37] R6_2.6.1           lifecycle_1.0.4    Seqinfo_1.0.0      KEGGREST_1.50.0
    ## [41] bit_4.6.0          pkgconfig_2.0.3    pillar_1.11.1      bslib_0.9.0
    ## [45] gtable_0.3.6       glue_1.8.0         Rcpp_1.1.0         xfun_0.55
    ## [49] tibble_3.3.0       tidyselect_1.2.1   knitr_1.51         dichromat_2.0-0.1
    ## [53] farver_2.1.2       htmltools_0.5.9    rmarkdown_2.30     labeling_0.4.3
    ## [57] compiler_4.5.2     S7_0.2.1
