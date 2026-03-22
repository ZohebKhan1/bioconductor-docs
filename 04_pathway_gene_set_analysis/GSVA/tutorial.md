
# GSVA on single-cell RNA-seq data

Robert Castelo[1\*], Axel Klenk[1\*\*] and Justin Guinney[2\*\*\*]

[1]Dept. of Medicine and Life Sciences, Universitat Pompeu Fabra, Barcelona, Spain\
[2]Tempus Labs, Inc.

[\*]<robert.castelo@upf.edu>\
[\*\*]<axel.klenk@gmail.com>\
[\*\*\*]<justin.guinney@tempus.com>

#### 18 March 2026

#### Abstract

Here we illustrate how to use GSVA with single-cell RNA sequencing (scRNA-seq) data.

#### Package

GSVA 2.4.8

**License**: Artistic-2.0

# [1] Introduction

GSVA provides now specific support for single-cell data in the algorithm that runs through the `gsvaParam()` parameter constructor, and originally described in the publication by [Hänzelmann, Castelo, and Guinney ([2013](#ref-haenzelmann2013gsva){role="doc-biblioref"})]. At the moment, this specific support consists of the following features:

-   The input expression data can be stored in different types of data containers prepared to store sparse single-cell data. These types of sparse data containers can be broadly categorized in those that only store the expression values, and those that may store additional row and column metadata. The currently available value-only containers for input are `dgCMatrix`, `SVT_SparseArray`, and `DelayedMatrix`. The currently available container for single-cell data that allows one to input additional row and column metadata is a `SingleCellExperiment` object.
-   While the input single-cell data is always sparse, the output of enrichment scores will be always dense, and therefore, the container storing those scores will be different from the input data, typically a `matrix` or a dense `DelayedMatrix` object. The latter will be particularly used when the total number of values exceeds 2\^31, which is the largest 32-bit standard integer value in R.
-   By default, when the input expression data is stored in a sparse data container, as it typically happens with single-cell data, then a slightly a slightly modified GSVA algorithm will run, if GSVA is the choice of algorithm, by which nonzero values are treated differently from zero values, leading to slightly different results than those obtained by applying the classical GSVA algorithm. If we set the parameter `sparse=FALSE` in the call to `gsvaParam()`, the classical GSVA algorithm will be used, which for a typical single-cell data set will result in longer running times and larger memory consumption than running it in the default sparse regime for this type of data.

In what follows, we will illustrate the use of GSVA on a publicly available single-cell transcriptomics data set of peripheral blood mononuclear cells (PBMCs) published by [Zheng et al. ([2017](#ref-zheng2017massively){role="doc-biblioref"})].

# [2] Import data

We import the PBMC data using the *[TENxPBMCData](https://bioconductor.org/packages/3.22/TENxPBMCData)* package, as a `SingleCellExperiment` object, defined in the *[SingleCellExperiment](https://bioconductor.org/packages/3.22/SingleCellExperiment)* package.

``` r
library(SingleCellExperiment)
library(TENxPBMCData)

sce <- TENxPBMCData(dataset="pbmc4k")
sce
class: SingleCellExperiment
dim: 33694 4340
metadata(0):
assays(1): counts
rownames(33694): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
  ENSG00000268674
rowData names(3): ENSEMBL_ID Symbol_TENx Symbol
colnames: NULL
colData names(11): Sample Barcode ... Individual Date_published
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):
```

# [3] Quality assessment and pre-processing

Here, we perform a quality assessment and pre-processing steps using the package *[scuttle](https://bioconductor.org/packages/3.22/scuttle)* [(McCarthy et al. [2017](#ref-mccarthy2017scater){role="doc-biblioref"})]. We start identifying mitochondrial genes.

``` r
library(scuttle)

is_mito <- grepl("^MT-", rowData(sce)$Symbol_TENx)
table(is_mito)
is_mito
FALSE  TRUE
33681    13
```

Calculate quality control (QC) metrics and filter out low-quality cells.

``` r
sce <- quickPerCellQC(sce, subsets=list(Mito=is_mito),
                           sub.fields="subsets_Mito_percent")
dim(sce)
[1] 33694  4147
```

Figure [1](#fig:cntxgene) below shows the empirical cumulative distribution of counts per gene in logarithmic scale.

``` r
cntxgene <- rowSums(assays(sce)$counts)+1
plot.ecdf(cntxgene, xaxt="n", panel.first=grid(), xlab="UMI counts per gene",
          log="x", main="", xlim=c(1, 1e5), las=1)
axis(1, at=10^(0:5), labels=10^(0:5))
abline(v=100, lwd=2, col="red")
```

Figure 1: [Empirical cumulative distribution of UMI counts per gene]\
The red vertical bar indicates a cutoff value of 100 UMI counts per gene across all cells, below which genes will be filtered out.

We filter out lowly-expressed genes, by selecting those with at least 100 UMI counts across all cells for downstream analysis.

``` r
sce <- sce[cntxgene >= 100, ]
dim(sce)
[1] 8823 4147
```

Calculate library size factors and normalized units of expression in logarithmic scale.

``` r
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)
assayNames(sce)
[1] "counts"    "logcounts"
```

# [4] Annotate cell types using GSVA

Here, we illustrate how to annotate cell types in the PBMC data using GSVA.

## [4.1] Read gene sets in GMT format

First, we fetch a collection of 22 leukocyte gene set signatures, containing a total 547 genes, which should help to distinguish among 22 mature human hematopoietic cell type populations isolated from peripheral blood or *in vitro* culture conditions, including seven T cell types: naïve and memory B cells, plasma cells, NK cell, and myeloid subsets. These gene sets have been used in the benchmarking publication by [Diaz-Mejia et al. ([2019](#ref-diaz2019evaluation){role="doc-biblioref"})], and were originally compiled by the [CIBERSORT](https://cibersortx.stanford.edu) developers, where they called it the LM22 signature [(Newman et al. [2015](#ref-newman2015robust){role="doc-biblioref"})]. The LM22 signature is stored in the *[GSVAdata](https://bioconductor.org/packages/3.22/GSVAdata)* experiment data package as a compressed text file in [GMT format](https://www.genepattern.org/file-formats-guide/#GMT), which can be read into R using the `readGMT()` function from the *[GSVA](https://bioconductor.org/packages/3.22/GSVA)* package, and will return the gene sets into a `GeneSetCollection` object, defined in the *[GSEABase](https://bioconductor.org/packages/3.22/GSEABase)* package.

``` r
library(GSEABase)
library(GSVA)

fname <- file.path(system.file("extdata", package="GSVAdata"),
                   "pbmc_cell_type_gene_set_signatures.gmt.gz")
gsets <- readGMT(fname)
gsets
GeneSetCollection
  names: B_CELLS_MEMORY, B_CELLS_NAIVE, ..., T_CELLS_REGULATORY_TREGS (22 total)
  unique identifiers: AIM2, BANK1, ..., SKAP1 (248 total)
  types in collection:
    geneIdType: SymbolIdentifier (1 total)
    collectionType: NullCollection (1 total)
```

## [4.2] Add gene identifier type metadata

Note that while gene identifers in the `sce` object correspond to [Ensembl stable identifiers](https://www.ensembl.org/info/genome/stable_ids/index.html) (`ENSG...`), the gene identifiers in the gene sets are [HGNC](https://www.genenames.org) gene symbols. This, in principle, precludes matching directly what gene in the single-cell data object `sce` corresponds to what gene set in the `GeneSetCollection` object `gsets`. However, the *[GSVA](https://bioconductor.org/packages/3.22/GSVA)* package can do that matching as long as the appropriate metadata is present in both objects.

In the case of a `GeneSetCollection` object, its `geneIdType` metadata slot stores the type of gene identifier. In the case of a `SingleCellExperiment` object, such as the previous `sce` object, such metadata is not present. However, using the function `gsvaAnnotation()` from the *[GSVA](https://bioconductor.org/packages/3.22/GSVA)* package, and the helper function `ENSEMBLIdentifier()` from the *[GSEABase](https://bioconductor.org/packages/3.22/GSEABase)* package, we add such metadata to the `sce` object as follows.

``` r
gsvaAnnotation(sce) <- ENSEMBLIdentifier("org.Hs.eg.db")
gsvaAnnotation(sce)
geneIdType: ENSEMBL (org.Hs.eg.db)
```

## [4.3] Build parameter object

We first build a parameter object using the function `gsvaParam()`. By default, the expression values in the `logocounts` assay will be selected for downstream analysis.

``` r
gsvapar <- gsvaParam(sce, gsets)
gsvapar
A GSVA::gsvaParam object
expression data:
  class: SingleCellExperiment
  dim: 8823 4147
  metadata(1): annotation
  assays(2): counts logcounts
  rownames(8823): ENSG00000279457 ENSG00000228463 ... ENSG00000198727
    ENSG00000273748
  rowData names(3): ENSEMBL_ID Symbol_TENx Symbol
  colnames: NULL
  colData names(22): Sample Barcode ... discard sizeFactor
  reducedDimNames(0):
  mainExpName: NULL
  altExpNames(0):
using assay: logcounts
using annotation:
  geneIdType: ENSEMBL (org.Hs.eg.db)
gene sets:
  GeneSetCollection
    names: B_CELLS_MEMORY, B_CELLS_NAIVE, ..., T_CELLS_REGULATORY_TREGS (22 total)
    unique identifiers: AIM2, BANK1, ..., SKAP1 (248 total)
    types in collection:
      geneIdType: SymbolIdentifier (1 total)
      collectionType: NullCollection (1 total)
gene set size: [1, Inf]
kcdf: auto
kcdfNoneMinSampleSize: 200
tau: 1
maxDiff: TRUE
absRanking: FALSE
sparse:  TRUE
checkNA: auto
missing data: didn't check
filterRows:  TRUE
ondisk:  auto
nonzero values: less than 2^31 (INT_MAX)
```

## [4.4] Calculate GSVA scores

While at this point, we could already run the entire GSVA algorithm with a call to the `gsva(gsvapar)` function. We show here how to do it in two steps. First we calculate GSVA rank values using the function `gsvaRanks()`.

``` r
gsvaranks <- gsvaRanks(gsvapar)
gsvaranks
A GSVA::gsvaRanksParam object
expression data:
  class: SingleCellExperiment
  dim: 8823 4147
  metadata(1): annotation
  assays(3): counts logcounts gsvaranks
  rownames(8823): ENSG00000279457 ENSG00000228463 ... ENSG00000198727
    ENSG00000273748
  rowData names(3): ENSEMBL_ID Symbol_TENx Symbol
  colnames: NULL
  colData names(22): Sample Barcode ... discard sizeFactor
  reducedDimNames(0):
  mainExpName: NULL
  altExpNames(0):
using assay: gsvaranks
using annotation:
  geneIdType: ENSEMBL (org.Hs.eg.db)
gene sets:
  GeneSetCollection
    names: B_CELLS_MEMORY, B_CELLS_NAIVE, ..., T_CELLS_REGULATORY_TREGS (22 total)
    unique identifiers: AIM2, BANK1, ..., SKAP1 (248 total)
    types in collection:
      geneIdType: SymbolIdentifier (1 total)
      collectionType: NullCollection (1 total)
gene set size: [1, Inf]
kcdf: auto
kcdfNoneMinSampleSize: 200
tau: 1
maxDiff: TRUE
absRanking: FALSE
sparse:  TRUE
checkNA: auto
missing data: didn't check
filterRows:  TRUE
ondisk:  auto
nonzero values: less than 2^31 (INT_MAX)
```

Second, we calculate the GSVA scores using the output of `gsvaRanks()` as input to the function `gsvaScores()`. By default, this function will calculate the scores for all gene sets specified in the input parameter object.

``` r
es <- gsvaScores(gsvaranks)
es
class: SingleCellExperiment
dim: 22 4147
metadata(0):
assays(1): es
rownames(22): B_CELLS_MEMORY B_CELLS_NAIVE ... T_CELLS_GAMMA_DELTA
  T_CELLS_REGULATORY_TREGS
rowData names(1): gs
colnames: NULL
colData names(22): Sample Barcode ... discard sizeFactor
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):
```

However, we could calculate the scores for another collection of gene sets by updating them in the `gsvaranks` object as follows.

``` r
geneSets(gsvaranks) <- geneSets(gsvapar)[1:2]
es2 <- gsvaScores(gsvaranks)
```

## [4.5] Using GSVA scores to assign cell types

Following [Amezquita et al. ([2020](#ref-amezquita2020orchestrating){role="doc-biblioref"})], and some of the steps described in "Chapter 5 Clustering" of the first version of the [OSCA book](https://bioconductor.org/books/3.16/OSCA.basic/clustering.html), we use GSVA scores to build a nearest-neighbor graph of the cells using the function `buildSNNGraph()` from the *[scran](https://bioconductor.org/packages/3.22/scran)* package [(Lun, McCarthy, and Marioni [2016](#ref-lun2016step){role="doc-biblioref"})]. The parameter `k` in the call to `buildSNNGraph()` specifies the number of nearest neighbors to consider during graph construction, and here we set `k=20` because it leads to a number of clusters close to the expected number of cell types.

``` r
library(scran)

g <- buildSNNGraph(es, k=20, assay.type="es")
```

Second, we use the function `cluster_walktrap()` from the *[igraph](https://CRAN.R-project.org/package=igraph)* package [(Csardi and Nepusz [2006](#ref-csardi2006igraph){role="doc-biblioref"})], to cluster cells by finding densely connected subgraphs. We store the resulting vector of cluster indicator values into the `sce` object using the function `colLabels()`.

``` r
library(igraph)

colLabels(es) <- factor(cluster_walktrap(g)$membership)
table(colLabels(es))

  1   2   3   4   5   6   7   8
495 601 502 525 972 191 345 516
```

Similarly to [Diaz-Mejia et al. ([2019](#ref-diaz2019evaluation){role="doc-biblioref"})], we apply a simple cell type assignment algorithm, which consists of selecting at each cell the gene set with highest GSVA score, tallying the selected gene sets per cluster, and assigning to the cluster the most frequent gene set, storing that assignment into the `sce` object with the function `colLabels()`.

``` r
## whmax <- apply(assay(es), 2, which.max)
whmax <- apply(assay(es), 2, function(x) which.max(as.vector(x)))
gsxlab <- split(rownames(es)[whmax], colLabels(es))
gsxlab <- names(sapply(sapply(gsxlab, table), which.max))
colLabels(es) <- factor(gsub("[0-9]\\.", "", gsxlab))[colLabels(es)]
table(colLabels(es))

     B_CELLS_NAIVE        EOSINOPHILS NK_CELLS_ACTIVATED   NK_CELLS_RESTING
               601               1027                495                191
 T_CELLS_CD4_NAIVE        T_CELLS_CD8
              1488                345
```

We can visualize the cell type assignments by projecting cells dissimilarity in two dimensions with a principal components analysis (PCA) on the GSVA scores, and coloring cells using the previously assigned clusters.

``` r
library(RColorBrewer)

res <- prcomp(assay(es))
varexp <- res$sdev^2 / sum(res$sdev^2)
nclusters <- nlevels(colLabels(es))
hmcol <- colorRampPalette(brewer.pal(nclusters, "Set1"))(nclusters)
par(mar=c(4, 5, 1, 1))
plot(res$rotation[, 1], res$rotation[, 2], col=hmcol[colLabels(es)], pch=19,
     xlab=sprintf("PCA 1 (%.0f%%)", varexp[1]*100),
     ylab=sprintf("PCA 2 (%.0f%%)", varexp[2]*100),
     las=1, cex.axis=1.2, cex.lab=1.5)
legend("topright", gsub("_", " ", levels(colLabels(es))), fill=hmcol, inset=0.01)
```

Figure 2: [Cell type assignments of PBMC scRNA-seq data, based on GSVA scores]\

Finally, if we want to better understand why a specific cell type is annotated to a given cell, we can use the `gsvaEnrichment()` function, which will show a GSEA enrichment plot. This function takes as input the output of `gsvaRanks()`, a given column (cell) in the input singl-cell data, and a given gene set. In Figure [3](#fig:gsvaenrichment) below, we show such a plot for the first cell annotated to the eosinophil cell type.

``` r
firsteosinophilcell <- which(colLabels(es) == "EOSINOPHILS")[1]
par(mar=c(4, 5, 1, 1))
gsvaEnrichment(gsvaranks, column=firsteosinophilcell, geneSet="EOSINOPHILS",
               cex.axis=1.2, cex.lab=1.5, plot="ggplot")
```

Figure 3: [GSVA enrichment plot of the EOSINOPHILS gene set in the expression profile of the first cell annotated to that cell type]\

In the previous call to `gsvaEnrichment()` we used the argument `plot="ggplot"` to produce a plot with the [ggplot2](https://cran.r-project.org/package=ggplot2) package. By default, if we call `gsvaEnrichment()` interactively, it will produce a plot using "base R", but either when we do it non-interactively, or when we set `plot="no"` it will return a `data.frame` object with the enrichment data.

# [5] Forthcoming features

These are features that we are working on and we expect to have them implemented in the near future (e.g., next release):

-   A specific implementation of the other methods ssGSEA, PLAGE and zscore to work on large datasets stored using a `DelayedArray` backend, such as HDF5, is not yet available.

We are still benchmarking and testing this version of GSVA for single-cell data. If you encounter problems or have suggestions, do not hesitate to contact us by opening an [issue](https://github.com/rcastelo/GSVA/issues) in the GSVA GitHub repo.

# Session information

Here is the output of `sessionInfo()` on the system on which this document was compiled running pandoc 2.7.3:

``` r
sessionInfo()
R version 4.5.2 (2025-10-31)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 24.04.3 LTS

Matrix products: default
BLAS:   /home/biocbuild/bbs-3.22-bioc/R/lib/libRblas.so
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_GB              LC_COLLATE=C
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

time zone: America/New_York
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods
[8] base

other attached packages:
 [1] igraph_2.2.2                scran_1.38.1
 [3] scuttle_1.20.0              TENxPBMCData_1.28.0
 [5] HDF5Array_1.38.0            h5mread_1.2.1
 [7] rhdf5_2.54.1                DelayedArray_0.36.0
 [9] SparseArray_1.10.9          S4Arrays_1.10.1
[11] abind_1.4-8                 Matrix_1.7-4
[13] SingleCellExperiment_1.32.0 sva_3.58.0
[15] BiocParallel_1.44.0         genefilter_1.92.0
[17] mgcv_1.9-4                  nlme_3.1-168
[19] RColorBrewer_1.1-3          edgeR_4.8.2
[21] limma_3.66.0                GSVAdata_1.46.0
[23] SummarizedExperiment_1.40.0 GenomicRanges_1.62.1
[25] Seqinfo_1.0.0               MatrixGenerics_1.22.0
[27] matrixStats_1.5.0           hgu95a.db_3.13.0
[29] GSEABase_1.72.0             graph_1.88.1
[31] annotate_1.88.0             XML_3.99-0.22
[33] org.Hs.eg.db_3.22.0         AnnotationDbi_1.72.0
[35] IRanges_2.44.0              S4Vectors_0.48.0
[37] Biobase_2.70.0              BiocGenerics_0.56.0
[39] generics_0.1.4              GSVA_2.4.8
[41] BiocStyle_2.38.0

loaded via a namespace (and not attached):
 [1] jsonlite_2.0.0           magrittr_2.0.4           magick_2.9.1
 [4] farver_2.1.2             rmarkdown_2.30           vctrs_0.7.1
 [7] memoise_2.0.1            tinytex_0.58             htmltools_0.5.9
[10] AnnotationHub_4.0.0      curl_7.0.0               BiocNeighbors_2.4.0
[13] Rhdf5lib_1.32.0          sass_0.4.10              bslib_0.10.0
[16] httr2_1.2.2              cachem_1.1.0             lifecycle_1.0.5
[19] pkgconfig_2.0.3          rsvd_1.0.5               R6_2.6.1
[22] fastmap_1.2.0            digest_0.6.39            dqrng_0.4.1
[25] irlba_2.3.7              ExperimentHub_3.0.0      RSQLite_2.4.6
[28] beachmat_2.26.0          filelock_1.0.3           labeling_0.4.3
[31] httr_1.4.8               compiler_4.5.2           bit64_4.6.0-1
[34] withr_3.0.2              S7_0.2.1                 DBI_1.3.0
[37] rappdirs_0.3.4           rjson_0.2.23             bluster_1.20.0
[40] tools_4.5.2              otel_0.2.0               glue_1.8.0
[43] rhdf5filters_1.22.0      grid_4.5.2               cluster_2.1.8.2
[46] memuse_4.2-3             gtable_0.3.6             BiocSingular_1.26.1
[49] ScaledMatrix_1.18.0      metapod_1.18.0           XVector_0.50.0
[52] BiocVersion_3.22.0       pillar_1.11.1            splines_4.5.2
[55] dplyr_1.2.0              BiocFileCache_3.0.0      lattice_0.22-9
[58] survival_3.8-6           bit_4.6.0                tidyselect_1.2.1
[61] locfit_1.5-9.12          Biostrings_2.78.0        knitr_1.51
[64] bookdown_0.46            xfun_0.56                statmod_1.5.1
[67] yaml_2.3.12              evaluate_1.0.5           codetools_0.2-20
[70] tibble_3.3.1             BiocManager_1.30.27      cli_3.6.5
[73] xtable_1.8-8             jquerylib_0.1.4          dichromat_2.0-0.1
[76] Rcpp_1.1.1               dbplyr_2.5.2             png_0.1-9
[79] parallel_4.5.2           ggplot2_4.0.2            blob_1.3.0
[82] sparseMatrixStats_1.22.0 SpatialExperiment_1.20.0 scales_1.4.0
[85] purrr_1.2.1              crayon_1.5.3             rlang_1.1.7
[88] KEGGREST_1.50.0
```

# References

Amezquita, Robert A, Aaron TL Lun, Etienne Becht, Vince J Carey, Lindsay N Carpp, Ludwig Geistlinger, Federico Marini, et al. 2020. "Orchestrating Single-Cell Analysis with Bioconductor." *Nature Methods* 17 (2): 137--45.

Csardi, Gabor, and Tamas Nepusz. 2006. "The Igraph Software Package for Complex Network Research." *InterJournal* Complex Systems: 1695. <https://igraph.org>.

Diaz-Mejia, J Javier, Elaine C Meng, Alexander R Pico, Sonya A MacParland, Troy Ketela, Trevor J Pugh, Gary D Bader, and John H Morris. 2019. "Evaluation of Methods to Assign Cell Type Labels to Cell Clusters from Single-Cell Rna-Sequencing Data." *F1000Research* 8: ISCB--Comm.

Hänzelmann, Sonja, Robert Castelo, and Justin Guinney. 2013. "GSVA: Gene Set Variation Analysis for Microarray and RNA-Seq Data." *BMC Bioinformatics* 14: 7. <https://doi.org/10.1186/1471-2105-14-7>.

Lun, Aaron TL, Davis J McCarthy, and John C Marioni. 2016. "A Step-by-Step Workflow for Low-Level Analysis of Single-Cell Rna-Seq Data with Bioconductor." *F1000Research* 5: 2122.

McCarthy, Davis J, Kieran R Campbell, Aaron TL Lun, and Quin F Wills. 2017. "Scater: Pre-Processing, Quality Control, Normalization and Visualization of Single-Cell Rna-Seq Data in R." *Bioinformatics* 33 (8): 1179--86.

Newman, Aaron M, Chih Long Liu, Michael R Green, Andrew J Gentles, Weiguo Feng, Yue Xu, Chuong D Hoang, Maximilian Diehn, and Ash A Alizadeh. 2015. "Robust Enumeration of Cell Subsets from Tissue Expression Profiles." *Nature Methods* 12 (5): 453--57.

Zheng, Grace XY, Jessica M Terry, Phillip Belgrader, Paul Ryvkin, Zachary W Bent, Ryan Wilson, Solongo B Ziraldo, et al. 2017. "Massively Parallel Digital Transcriptional Profiling of Single Cells." *Nature Communications* 8 (1): 14049.
