
# GSVA: gene set variation analysis

Robert Castelo[1\*], Axel Klenk[1\*\*] and Justin Guinney[2\*\*\*]

[1]Dept. of Medicine and Life Sciences, Universitat Pompeu Fabra, Barcelona, Spain\
[2]Tempus Labs, Inc.

[\*]<robert.castelo@upf.edu>\
[\*\*]<axel.klenk@gmail.com>\
[\*\*\*]<justin.guinney@tempus.com>

#### 18 March 2026

#### Abstract

Gene set variation analysis (GSVA) is a particular type of gene set enrichment method that works on single samples and enables pathway-centric analyses of molecular data by performing a conceptually simple but powerful change in the functional unit of analysis, from genes to gene sets. The GSVA package provides the implementation of four single-sample gene set enrichment methods, concretely *zscore*, *plage*, *ssGSEA* and its own called *GSVA*. While this methodology was initially developed for gene expression data, it can be applied to other types of molecular profiling data. In this vignette we illustrate how to use the GSVA package with bulk microarray and RNA-seq expression data.

#### Package

GSVA 2.4.8

**License**: Artistic-2.0

# [1] Quick start

*[GSVA](https://bioconductor.org/packages/3.22/GSVA)* is an R package distributed as part of the [Bioconductor](https://bioconductor.org) project. To install the package, start R and enter:

``` r
install.packages("BiocManager")
BiocManager::install("GSVA")
```

Once *[GSVA](https://bioconductor.org/packages/3.22/GSVA)* is installed, it can be loaded with the following command.

``` r
library(GSVA)
```

Given a gene expression data matrix, which we shall call `X`, with rows corresponding to genes and columns to samples, such as this one simulated from random Gaussian data:

``` r
p <- 10000 ## number of genes
n <- 30    ## number of samples
## simulate expression values from a standard Gaussian distribution
X <- matrix(rnorm(p*n), nrow=p,
            dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
X[1:5, 1:5]
            s1         s2        s3         s4          s5
g1  0.57980333  1.9460830 1.3822990  0.7803017 -1.61228280
g2  1.64977188 -0.4500512 0.3444646 -1.6849470  0.68413338
g3 -0.31084704  1.6284771 1.5367198  1.4805650  0.05926513
g4 -0.26310814  0.2276478 1.1820973 -0.8531259 -0.74497926
g5 -0.08792931 -0.6713186 0.5753177 -0.5327106  0.18114692
```

Given a collection of gene sets stored, for instance, in a `list` object, which we shall call `gs`, with genes sampled uniformly at random without replacement into 100 different gene sets:

``` r
## sample gene set sizes
gs <- as.list(sample(10:100, size=100, replace=TRUE))
## sample gene sets
gs <- lapply(gs, function(n, p)
                   paste0("g", sample(1:p, size=n, replace=FALSE)), p)
names(gs) <- paste0("gs", 1:length(gs))
```

We can calculate GSVA enrichment scores as follows. First we should build a parameter object for the desired methodology. Here we illustrate it with the GSVA algorithm of [Hänzelmann, Castelo, and Guinney ([2013](#ref-haenzelmann2013gsva){role="doc-biblioref"})] by calling the function `gsvaParam()`, but other parameter object constructor functions are available; see in the next section below.

``` r
gsvaPar <- gsvaParam(X, gs)
gsvaPar
A GSVA::gsvaParam object
expression data:
  matrix [10000, 30]
    rows: g1, g2, ..., g10000 (10000 total)
    cols: s1, s2, ..., s30 (30 total)
using assay: none
using annotation:
  geneIdType: Null
gene sets:
  list
    names: gs1, gs2, ..., gs100 (100 total)
    unique identifiers: g4578, g9822, ..., g427 (4040 total)
gene set size: [1, Inf]
kcdf: auto
kcdfNoneMinSampleSize: 200
tau: 1
maxDiff: TRUE
absRanking: FALSE
sparse:  FALSE
checkNA: auto
missing data: no
filterRows:  TRUE
ondisk:  auto
nonzero values: less than 2^31 (INT_MAX)
```

The first argument to the `gsvaParam()` function constructing this parameter object is the gene expression data matrix, and the second is the collection of gene sets. In this example, we provide expression data and gene sets into base R *matrix* and *list* objects, respectively, to the `gsvaParam()` function, but it can take also different specialized containers that facilitate the access and manipulation of molecular and phenotype data, as well as their associated metadata.

Second, we call the `gsva()` function with the parameter object as first argument. Other additional arguments to the `gsva()` function are `verbose` to control progress reporting and `BPPPARAM` to perform calculations in parallel through the package *[BiocParallel](https://bioconductor.org/packages/3.22/BiocParallel)*.

``` r
gsva.es <- gsva(gsvaPar, verbose=FALSE)
dim(gsva.es)
[1] 100  30
gsva.es[1:5, 1:5]
             s1          s2          s3          s4          s5
gs1 -0.10317431 -0.01519608  0.10542740 -0.19864230 -0.02928974
gs2 -0.06362930  0.18530341 -0.01913301 -0.02141608 -0.08000869
gs3 -0.08826469  0.04686503  0.31888132  0.05872338  0.32992926
gs4  0.08988159 -0.22954284  0.18200253  0.03204703  0.04312882
gs5 -0.17714013  0.22068397  0.07286054  0.12069209 -0.03722003
```

# [2] Introduction

Gene set variation analysis (GSVA) provides an estimate of pathway activity by transforming an input gene-by-sample expression data matrix into a corresponding gene-set-by-sample expression data matrix. This resulting expression data matrix can be then used with classical analytical methods such as differential expression, classification, survival analysis, clustering or correlation analysis in a pathway-centric manner. One can also perform sample-wise comparisons between pathways and other molecular data types such as microRNA expression or binding data, copy-number variation (CNV) data or single nucleotide polymorphisms (SNPs).

The GSVA package provides an implementation of this approach for the following methods:

-   *plage* [(Tomfohr, Lu, and Kepler [2005](#ref-tomfohr_pathway_2005){role="doc-biblioref"})]. Pathway level analysis of gene expression (PLAGE) standardizes expression profiles over the samples and then, for each gene set, it performs a singular value decomposition (SVD) over its genes. The coefficients of the first right-singular vector are returned as the estimates of pathway activity over the samples. Note that, because of how SVD is calculated, the sign of its singular vectors is arbitrary.

-   *zscore* [(Lee et al. [2008](#ref-lee_inferring_2008){role="doc-biblioref"})]. The z-score method standardizes expression profiles over the samples and then, for each gene set, combines the standardized values as follows. Given a gene set \\(\\gamma=\\{1,\\dots,k\\}\\) with standardized values \\(z_1,\\dots,z_k\\) for each gene in a specific sample, the combined z-score \\(Z\_\\gamma\\) for the gene set \\(\\gamma\\) is defined as: [\\[ Z\_\\gamma = \\frac{\\sum\_{i=1}\^k z_i}{\\sqrt{k}}\\,. \\]]

-   *ssgsea* [(Barbie et al. [2009](#ref-barbie_systematic_2009){role="doc-biblioref"})]. Single sample GSEA (ssGSEA) is a non-parametric method that calculates a gene set enrichment score per sample as the normalized difference in empirical cumulative distribution functions (CDFs) of gene expression ranks inside and outside the gene set. By default, the implementation in the GSVA package follows the last step described in [(Barbie et al. [2009](#ref-barbie_systematic_2009){role="doc-biblioref"}, online methods, pg. 2)] by which pathway scores are normalized, dividing them by the range of calculated values. This normalization step may be switched off using the argument `ssgsea.norm` in the call to the `gsva()` function; see below.

-   *gsva* [(Hänzelmann, Castelo, and Guinney [2013](#ref-haenzelmann2013gsva){role="doc-biblioref"})]. This is the default method of the package and similarly to ssGSEA, is a non-parametric method that uses the empirical CDFs of gene expression ranks inside and outside the gene set, but it starts by calculating an expression-level statistic that brings gene expression profiles with different dynamic ranges to a common scale, and combines the information of those empirical CDFs in a different way to provide an enrichment score.

The interested user may find full technical details about how these methods work in their corresponding articles cited above. If you use any of them in a publication, please cite them with the given bibliographic reference.

# [3] Overview of the GSVA functionality

The workhorse of the GSVA package is the function `gsva()`, which takes a *parameter object* as its main input. There are four classes of parameter objects corresponding to the methods listed above, and may have different additional parameters to tune, but all of them require at least the following two input arguments:

1.  A normalized gene expression dataset, which can be provided in one of the following containers:
    -   A `matrix` of expression values with genes corresponding to rows and samples corresponding to columns.
    -   An `ExpressionSet` object; see package *[Biobase](https://bioconductor.org/packages/3.22/Biobase)*.
    -   A `SummarizedExperiment` object, see package *[SummarizedExperiment](https://bioconductor.org/packages/3.22/SummarizedExperiment)*.
2.  A collection of gene sets; which can be provided in one of the following containers:
    -   A `list` object where each element corresponds to a gene set defined by a vector of gene identifiers, and the element names correspond to the names of the gene sets.
    -   A `GeneSetCollection` object; see package *[GSEABase](https://bioconductor.org/packages/3.22/GSEABase)*.

One advantage of providing the input data using specialized containers such as `ExpressionSet`, `SummarizedExperiment` and `GeneSetCollection` is that the `gsva()` function will automatically map the gene identifiers between the expression data and the gene sets (internally calling the function `mapIdentifiers()` from the package *[GSEABase](https://bioconductor.org/packages/3.22/GSEABase)*), when they come from different standard nomenclatures, i.e., *Ensembl* versus *Entrez*, provided the input objects contain the appropriate metadata; see next section.

If either the input gene expression data is provided as a `matrix` object or the gene sets are provided in a `list` object, or both, it is then the responsibility of the user to ensure that both objects contain gene identifiers following the same standard nomenclature.

Before the actual calculations take place, the `gsva()` function will apply the following filters:

1.  Discard genes in the input expression data matrix with constant expression.

2.  Discard genes in the input gene sets that do not map to a gene in the input gene expression data matrix.

3.  Discard gene sets that, after applying the previous filters, do not meet a minimum and maximum size, which by default is one for the minimum size and has no limit for the maximum size.

If, as a result of applying these three filters, either no genes or gene sets are left, the `gsva()` function will prompt an error. A common cause for such an error at this stage is that gene identifiers between the expression data matrix and the gene sets do not belong to the same standard nomenclature and could not be mapped. This may happen because either the input data were not provided using some of the specialized containers described above or the necessary metadata in those containers that allows the software to successfully map gene identifiers, is missing.

The method employed by the `gsva()` function is determined by the class of the parameter object that it receives as an input. An object constructed using the `gsvaParam()` function runs the method described by [Hänzelmann, Castelo, and Guinney ([2013](#ref-haenzelmann2013gsva){role="doc-biblioref"})], but this can be changed using the parameter constructor functions `plageParam()`, `zscoreParam()`, or `ssgseaParam()`, corresponding to the methods briefly described in the introduction; see also their corresponding help pages.

When using `gsvaParam()`, the user can additionally tune the following parameters, whose default values cover most of the use cases:

-   `kcdf`: The first step of the GSVA algorithm brings gene expression profiles to a common scale by calculating an expression statistic through the estimation of the CDF across samples. The way in which such an estimation is performed by GSVA is controlled by the `kcdf` parameter, which accepts the following four possible values: (1) `"auto"`, the default value, lets GSVA automatically decide the estimation method; (2) `"Gaussian"`, use a Gaussian kernel, suitable for continuous expression data, such as microarray fluorescent units in logarithmic scale and RNA-seq log-CPMs, log-RPKMs or log-TPMs units of expression; (2) `"Poisson"`, use a Poisson kernel, suitable for integer counts, such as those derived from RNA-seq alignments; (3) `"none"`, which will perform a direct estimation of the CDF without a kernel function.

-   `kcdfNoneMinSampleSize`: When `kcdf="auto"`, this parameter decides at what minimum sample size `kcdf="none"`, i.e., the estimation of the empirical cumulative distribution function (ECDF) of expression levels across samples is performed directly without using a kernel; see the previous `kcdf` parameter. By default `kcdfNoneMinSampleSize=200`.

-   `tau`: Exponent defining the weight of the tail in the random walk. By default `tau=1`.

-   `maxDiff`: The last step of the GSVA algorithm calculates the gene set enrichment score from two Kolmogorov-Smirnov random walk statistics. This parameter is a logical flag that allows the user to specify two possible ways to do such calculation: (1) `TRUE`, the default value, where the enrichment score is calculated as the magnitude difference between the largest positive and negative random walk deviations. This default value gives larger enrichment scores to gene sets whose genes are concordantly activated in one direction only; (2) `FALSE`, where the enrichment score is calculated as the maximum distance of the random walk from zero. This approach produces a distribution of enrichment scores that is bimodal, but it can be give large enrichment scores to gene sets whose genes are not concordantly activated in one direction only.

-   `absRanking`: Logical flag used only when `maxDiff=TRUE`. By default, `absRanking=FALSE` and it implies that a modified Kuiper statistic is used to calculate enrichment scores, taking the magnitude difference between the largest positive and negative random walk deviations. When `absRanking=TRUE` the original Kuiper statistic is used, by which the largest positive and negative random walk deviations are added together.

-   `sparse`: Logical flag used only when the input expression data is stored in a sparse matrix (e.g., a `dgCMatrix` or a `SingleCellExperiment` object storing the expression data in a `dgCMatrix`). In such as case, when `sparse=TRUE` (default), a sparse version of the GSVA algorithm will be applied. Otherwise, when `sparse=FALSE`, the classical version of the GSVA algorithm will be used.

In general, the default values for the previous parameters are suitable for most analysis settings, which usually consist of some kind of normalized continuous expression values.

# [4] Gene set definitions and containers

Gene sets constitute a simple, yet useful, way to define pathways because we use pathway membership definitions only, neglecting the information on molecular interactions. Gene set definitions are a crucial input to any gene set enrichment analysis because if our gene sets do not capture the biological processes we are studying, we will likely not find any relevant insights in our data from an analysis based on these gene sets.

There are multiple sources of gene sets, the most popular ones being [The Gene Ontology (GO) project](https://geneontology.org) and [The Molecular Signatures Database (MSigDB)](https://www.gsea-msigdb.org/gsea/msigdb). Sometimes gene set databases will not include the ones we need. In such a case we should either curate our own gene sets or use techniques to infer them from data.

The most basic data container for gene sets in R is the `list` class of objects, as illustrated before in the quick start section, where we defined a toy collection of three gene sets stored in a list object called `gs`:

``` r
class(gs)
[1] "list"
length(gs)
[1] 100
head(lapply(gs, head))
$gs1
[1] "g4578" "g9822" "g8896" "g8165" "g4646" "g6473"

$gs2
[1] "g555"  "g8418" "g7467" "g9278" "g1056" "g679"

$gs3
[1] "g5620" "g2821" "g6005" "g3346" "g6647" "g3735"

$gs4
[1] "g3696" "g2589" "g5303" "g8236" "g5183" "g21"

$gs5
[1] "g1156" "g9932" "g820"  "g6606" "g3764" "g2114"

$gs6
[1] "g1574" "g5192" "g1213" "g7009" "g4522" "g3442"
```

Using a Bioconductor organism-level package such as *[org.Hs.eg.db](https://bioconductor.org/packages/3.22/org.Hs.eg.db)* we can easily build a list object containing a collection of gene sets defined as GO terms with annotated Entrez gene identifiers, as follows:

``` r
library(org.Hs.eg.db)

goannot <- select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns="GO")
head(goannot)
  ENTREZID         GO EVIDENCE ONTOLOGY
1        1 GO:0002764      IBA       BP
2        1 GO:0005576      HDA       CC
3        1 GO:0005576      IDA       CC
4        1 GO:0005576      IEA       CC
5        1 GO:0005576      TAS       CC
6        1 GO:0005615      HDA       CC
genesbygo <- split(goannot$ENTREZID, goannot$GO)
length(genesbygo)
[1] 18864
head(genesbygo)
$`GO:0000009`
[1] "55650" "55650" "55650" "79087"

$`GO:0000012`
 [1] "1161"      "2074"      "3981"      "7141"      "7374"      "7515"
 [7] "7515"      "23411"     "54840"     "54840"     "54840"     "54840"
[13] "55247"     "55775"     "55775"     "55775"     "200558"    "100133315"

$`GO:0000014`
 [1] "2021"  "2067"  "2072"  "4361"  "4361"  "4361"  "5932"  "6419"  "6419"
[10] "6419"  "9941"  "10111" "28990" "64421"

$`GO:0000015`
[1] "2023"   "2023"   "2023"   "2026"   "2026"   "2027"   "2027"   "387712"
[9] "387712"

$`GO:0000016`
[1] "3938" "3938" "3938" "3938"

$`GO:0000017`
[1] "6523" "6523" "6523" "6524"
```

A more sophisticated container for gene sets is the `GeneSetCollection` object class defined in the *[GSEABase](https://bioconductor.org/packages/3.22/GSEABase)* package, which also provides the function `getGmt()` to import [gene matrix transposed (GMT)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) files such as those provided by [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) into a `GeneSetCollection` object. The experiment data package *[GSVAdata](https://bioconductor.org/packages/3.22/GSVAdata)* provides one such object with the old (3.0) version of the C2 collection of curated genesets from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb), which can be loaded as follows.

``` r
library(GSEABase)
library(GSVAdata)

data(c2BroadSets)
class(c2BroadSets)
[1] "GeneSetCollection"
attr(,"package")
[1] "GSEABase"
c2BroadSets
GeneSetCollection
  names: NAKAMURA_CANCER_MICROENVIRONMENT_UP, NAKAMURA_CANCER_MICROENVIRONMENT_DN, ..., ST_PHOSPHOINOSITIDE_3_KINASE_PATHWAY (3272 total)
  unique identifiers: 5167, 100288400, ..., 57191 (29340 total)
  types in collection:
    geneIdType: EntrezIdentifier (1 total)
    collectionType: BroadCollection (1 total)
```

The documentation of *[GSEABase](https://bioconductor.org/packages/3.22/GSEABase)* contains a description of the `GeneSetCollection` class and its associated methods.

# [5] Importing gene sets from GMT files

An important source of gene sets is the [Molecular Signatures Database (MSigDB)](https://www.gsea-msigdb.org/gsea/msigdb) [(Subramanian et al. [2005](#ref-subramanian_gene_2005){role="doc-biblioref"})], which stores them in plain text files following the so-called [*gene matrix transposed* (GMT) format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29). In the GMT format, each line stores a gene set with the following values separated by tabs:

-   A unique gene set identifier.
-   A gene set description.
-   One or more gene identifiers.

Because each different gene set may consist of a different number of genes, each line in a GMT file may contain a different number of tab-separated values. This means that the GMT format is not a tabular format, and therefore cannot be directly read with base R functions such as `read.table()` or `read.csv()`.

We need a specialized function to read GMT files. We can find such a function in the *[GSEABase](https://bioconductor.org/packages/3.22/GSEABase)* package with `getGmt()`, or in the *[qusage](https://bioconductor.org/packages/3.22/qusage)* package with `read.gmt()`.

GSVA also provides such a function called `readGMT()`, which takes as first argument the filename or URL of a, possibly compressed, GMT file. The call below illustrates how to read a GMT file from MSigDB providing its URL, concretely the one corresponding to the C7 collection of immunologic signature gene sets. Note that we also load the package *[GSEABase](https://bioconductor.org/packages/3.22/GSEABase)* because, by default, the value returned by `readGMT()` is a `GeneSetCollection` object defined in that package.

``` r
library(GSEABase)
library(GSVA)

URL <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c7.immunesigdb.v2024.1.Hs.symbols.gmt"
c7.genesets <- readGMT(URL)
```

    GeneSetCollection
      names: GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN, GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_UP, ..., KAECH_NAIVE_VS_MEMORY_CD8_TCELL_UP (4872 total)
      unique identifiers: ABCA2, ABCC5, ..., LINC00841 (20457 total)
      types in collection:
        geneIdType: SymbolIdentifier (1 total)
        collectionType: NullCollection (1 total)

By default, `readGMT()` returns a `GeneSetCollection` object, but this can be switched to `list` object by setting the argument `valueType="list"`. It will also attempt to figure out the type of identifier used in the gene set and set the corresponding metadata in the resulting object. However, this can be also manually set either through the parameter `geneIdType` or in a call to the setter method `gsvaAnnotation()`:

``` r
gsvaAnnotation(c7.genesets) <- SymbolIdentifier("org.Hs.eg.db")
```

This operation can actually also be done with `list` objects and it will add the metadata through an R attribute that later the `gsva()` function will be able to read. For understanding the different types of available gene identifier metadata constructor functions, please consult the help page of `GeneIdentifierType` in the *[GSEABase](https://bioconductor.org/packages/3.22/GSEABase)* package.

The specification of the [GMT format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) establishes that duplicated gene set names are not allowed. For this reason, the `getGmt()` function from the [GSEABase](https://bioconductor.org/packages/GSEABase) package prompts an error when duplicated gene names are found, while the `read.gmt()` function from the [qusage](https://bioconductor.org/packages/qusage) package silently accepts them in a list with duplicated element names.

The GSVA `readGMT()` function deals with duplicated gene set names as follows. By default, `readGMT()` warns the user about a duplicated gene set name and keeps only the first occurrence of the duplicated gene set in the returned object. We can illustrate this situation with an old GMT file from the MSigDB database that happens to have duplicated gene set names and which a small subset of it is stored in the *[GSVAdata](https://bioconductor.org/packages/3.22/GSVAdata)* package.

``` r
fname <- system.file("extdata", "c2.subsetdups.v7.5.symbols.gmt.gz",
                     package="GSVAdata")
c2.dupgenesets <- getGmt(fname, geneIdType=SymbolIdentifier())
Error in `validObject()`:
! invalid class "GeneSetCollection" object: each setName must be distinct
c2.dupgenesets
Error:
! object 'c2.dupgenesets' not found
```

We can see that `getGmt()` prompts an error. We can see below that this does not happen with `readGMT()` and that, by default, all but the first occurrence of the duplicated gene set have been removed.

``` r
c2.dupgenesets <- readGMT(fname, geneIdType=SymbolIdentifier())
Warning in deduplicateGmtLines(lines, deduplUse): GMT contains duplicated gene
set names; deduplicated using method: first
c2.dupgenesets
GeneSetCollection
  names: CORONEL_RFX7_DIRECT_TARGETS_UP, FOROUTAN_TGFB_EMT_UP, ..., CHANDRAN_METASTASIS_TOP50_DN (109 total)
  unique identifiers: ABAT, DIP2A, ..., MSRB2 (6929 total)
  types in collection:
    geneIdType: SymbolIdentifier (1 total)
    collectionType: NullCollection (1 total)
any(duplicated(names(c2.dupgenesets)))
[1] FALSE
```

The parameter `deduplUse` in the `readGMT()` function allow one to apply other policies to deal with duplicated gene set names, see the help page of `readGMT()` with `?readGMT` for full details on this parameter.

# [6] Quantification of pathway activity in bulk microarray and RNA-seq data

Here we illustrate how GSVA provides an analogous quantification of pathway activity in both microarray and RNA-seq data by using two such datasets that have been derived from the same biological samples. More concretely, we will use gene expression data of lymphoblastoid cell lines (LCL) from HapMap individuals that have been profiled using both technologies [(Huang et al. [2007](#ref-huang_genome-wide_2007){role="doc-biblioref"}, \@pickrell_understanding_2010)]. These data form part of the experimental package *[GSVAdata](https://bioconductor.org/packages/3.22/GSVAdata)* and the corresponding help pages contain details on how the data were processed. We start loading these data and verifying that they indeed contain expression data for the same genes and samples, as follows:

``` r
library(Biobase)

data(commonPickrellHuang)

stopifnot(identical(featureNames(huangArrayRMAnoBatchCommon_eset),
                    featureNames(pickrellCountsArgonneCQNcommon_eset)))
stopifnot(identical(sampleNames(huangArrayRMAnoBatchCommon_eset),
                    sampleNames(pickrellCountsArgonneCQNcommon_eset)))
```

Next, for the current analysis we use the subset of canonical pathways from the C2 collection of MSigDB Gene Sets v3.0. These correspond to the following pathways from KEGG, REACTOME and BIOCARTA:

``` r
canonicalC2BroadSets <- c2BroadSets[c(grep("^KEGG", names(c2BroadSets)),
                                      grep("^REACTOME", names(c2BroadSets)),
                                      grep("^BIOCARTA", names(c2BroadSets)))]
canonicalC2BroadSets
GeneSetCollection
  names: KEGG_GLYCOLYSIS_GLUCONEOGENESIS, KEGG_CITRATE_CYCLE_TCA_CYCLE, ..., BIOCARTA_ACTINY_PATHWAY (833 total)
  unique identifiers: 55902, 2645, ..., 8544 (6744 total)
  types in collection:
    geneIdType: EntrezIdentifier (1 total)
    collectionType: BroadCollection (1 total)
```

Additionally, we extend this collection of gene sets with two formed by genes with sex-specific expression, which also form part of the *[GSVAdata](https://bioconductor.org/packages/3.22/GSVAdata)* experiment data package. Here we use the constructor function `GeneSet` from the *[GSEABase](https://bioconductor.org/packages/3.22/GSEABase)* package to build the objects that we add to the `GeneSetCollection` object `canonicalC2BroadSets`.

``` r
data(genderGenesEntrez)

MSY <- GeneSet(msYgenesEntrez, geneIdType=EntrezIdentifier(),
               collectionType=BroadCollection(category="c2"),
               setName="MSY")
MSY
setName: MSY
geneIds: 266, 84663, ..., 353513 (total: 34)
geneIdType: EntrezId
collectionType: Broad
  bcCategory: c2 (Curated)mh (Mouse-Ortholog Hallmark)
  bcSubCategory: NA
details: use 'details(object)'
XiE <- GeneSet(XiEgenesEntrez, geneIdType=EntrezIdentifier(),
               collectionType=BroadCollection(category="c2"),
               setName="XiE")
XiE
setName: XiE
geneIds: 293, 8623, ..., 1121 (total: 66)
geneIdType: EntrezId
collectionType: Broad
  bcCategory: c2 (Curated)mh (Mouse-Ortholog Hallmark)
  bcSubCategory: NA
details: use 'details(object)'

canonicalC2BroadSets <- GeneSetCollection(c(canonicalC2BroadSets, MSY, XiE))
canonicalC2BroadSets
GeneSetCollection
  names: KEGG_GLYCOLYSIS_GLUCONEOGENESIS, KEGG_CITRATE_CYCLE_TCA_CYCLE, ..., XiE (835 total)
  unique identifiers: 55902, 2645, ..., 1121 (6810 total)
  types in collection:
    geneIdType: EntrezIdentifier (1 total)
    collectionType: BroadCollection (1 total)
```

We calculate now GSVA enrichment scores for these gene sets using first the normalized microarray data and then the normalized RNA-seq integer count data. Note that the only requirement to do the latter is to set the argument `kcdf="Poisson"`, which is `"Gaussian"` by default. Note, however, that if our RNA-seq normalized expression levels would be continuous, such as log-CPMs, log-RPKMs or log-TPMs, the default value of the `kcdf` argument should remain unchanged.

``` r
huangPar <- gsvaParam(huangArrayRMAnoBatchCommon_eset, canonicalC2BroadSets,
                      minSize=5, maxSize=500)
esmicro <- gsva(huangPar)
pickrellPar <- gsvaParam(pickrellCountsArgonneCQNcommon_eset,
                         canonicalC2BroadSets, minSize=5, maxSize=500,
                         kcdf="Poisson")
esrnaseq <- gsva(pickrellPar)
```

We are going to assess how gene expression profiles correlate between microarray and RNA-seq data and compare those correlations with the ones derived at pathway level. To compare gene expression values of both technologies, we will transform first the RNA-seq integer counts into log-CPM units of expression using the `cpm()` function from the *[edgeR](https://bioconductor.org/packages/3.22/edgeR)* package.

``` r
library(edgeR)

lcpms <- cpm(exprs(pickrellCountsArgonneCQNcommon_eset), log=TRUE)
```

We calculate Spearman correlations between gene expression profiles of the previous log-CPM values and the microarray RMA values.

``` r
genecorrs <- sapply(1:nrow(lcpms),
                    function(i, expmicro, exprnaseq)
                      cor(expmicro[i, ], exprnaseq[i, ], method="spearman"),
                    exprs(huangArrayRMAnoBatchCommon_eset), lcpms)
names(genecorrs) <- rownames(lcpms)
```

Now calculate Spearman correlations between GSVA enrichment scores derived from the microarray and the RNA-seq data.

``` r
pwycorrs <- sapply(1:nrow(esmicro),
                   function(i, esmicro, esrnaseq)
                     cor(esmicro[i, ], esrnaseq[i, ], method="spearman"),
                   exprs(esmicro), exprs(esrnaseq))
names(pwycorrs) <- rownames(esmicro)
```

Figure [1](#fig:compcorrs) below shows the two distributions of these correlations and we can see that GSVA enrichment scores provide an agreement between microarray and RNA-seq data comparable to the one observed between gene-level units of expression.

Figure 1: [Comparison of correlation values of gene and pathway expression profiles derived from microarray and RNA-seq data]\

Finally, in Figure [2](#fig:compsexgenesets) we compare the actual GSVA enrichment scores for two gene sets formed by genes with sex-specific expression. Concretely, one gene set (XIE) formed by genes that escape chromosome X-inactivation in females [(Carrel and Willard [2005](#ref-carrel_x-inactivation_2005){role="doc-biblioref"})] and another gene set (MSY) formed by genes located on the male-specific region of chromosome Y [(Skaletsky et al. [2003](#ref-skaletsky_male-specific_2003){role="doc-biblioref"})].

Figure 2: [Comparison of GSVA enrichment scores obtained from microarray and RNA-seq data for two gene sets formed by genes with sex-specific expression]\

We can see how microarray and RNA-seq single-sample GSVA enrichment scores correlate very well in these gene sets, with \\(\\rho=0.80\\) for the male-specific gene set and \\(\\rho=0.79\\) for the female-specific gene set. Male and female samples show higher GSVA enrichment scores in their corresponding gene set.

# [7] Example applications

## [7.1] Molecular signature identification

In [(Verhaak et al. [2010](#ref-verhaak_integrated_2010){role="doc-biblioref"})] four subtypes of glioblastoma multiforme (GBM) -proneural, classical, neural and mesenchymal- were identified by the characterization of distinct gene-level expression patterns. Using four gene set signatures specific to brain cell types (astrocytes, oligodendrocytes, neurons and cultured astroglial cells), derived from murine models by [Cahoy et al. ([2008](#ref-cahoy_transcriptome_2008){role="doc-biblioref"})], we replicate the analysis of [Verhaak et al. ([2010](#ref-verhaak_integrated_2010){role="doc-biblioref"})] by using GSVA to transform the gene expression measurements into enrichment scores for these four gene sets, without taking the sample subtype grouping into account. We start by having a quick glance to the data, which forms part of the *[GSVAdata](https://bioconductor.org/packages/3.22/GSVAdata)* package:

``` r
data(gbm_VerhaakEtAl)
gbm_eset
ExpressionSet (storageMode: lockedEnvironment)
assayData: 11861 features, 173 samples
  element names: exprs
protocolData: none
phenoData
  rowNames: TCGA.02.0003.01A.01 TCGA.02.0010.01A.01 ...
    TCGA.12.0620.01A.01 (173 total)
  varLabels: subtype
  varMetadata: labelDescription channel
featureData: none
experimentData: use 'experimentData(object)'
Annotation:
head(featureNames(gbm_eset))
[1] "AACS"    "FSTL1"   "ELMO2"   "CREB3L1" "RPS11"   "PNMA1"
table(gbm_eset$subtype)

  Classical Mesenchymal      Neural   Proneural
         38          56          26          53
data(brainTxDbSets)
lengths(brainTxDbSets)
      astrocytic_up        astroglia_up         neuronal_up oligodendrocytic_up
                 85                  88                  98                  70
lapply(brainTxDbSets, head)
$astrocytic_up
[1] "GRHL1"   "GPAM"    "PAPSS2"  "MERTK"   "BTG1"    "SLC46A1"

$astroglia_up
[1] "BST2"     "SERPING1" "ACTA2"    "C9orf167" "C1orf31"  "ANXA4"

$neuronal_up
[1] "STXBP1"  "JPH4"    "CACNG3"  "BRUNOL6" "CLSTN2"  "FAM123C"

$oligodendrocytic_up
[1] "DCT"    "ZNF536" "GNG8"   "ELOVL6" "NR2C1"  "RCBTB1"
```

GSVA enrichment scores for the gene sets contained in `brainTxDbSets` are calculated, in this case using `mx.diff=FALSE`, as follows:

``` r
gbmPar <- gsvaParam(gbm_eset, brainTxDbSets, maxDiff=FALSE)
gbm_es <- gsva(gbmPar)
```

Figure [3](#fig:gbmSignature) shows the GSVA enrichment scores obtained for the up-regulated gene sets across the samples of the four GBM subtypes. As expected, the *neural* class is associated with the neural gene set and the astrocytic gene sets. The *mesenchymal* subtype is characterized by the expression of mesenchymal and microglial markers, thus we expect it to correlate with the astroglial gene set. The *proneural* subtype shows high expression of oligodendrocytic development genes, thus it is not surprising that the oligodendrocytic gene set is highly enriched for ths group. Interestingly, the *classical* group correlates highly with the astrocytic gene set. In summary, the resulting GSVA enrichment scores recapitulate accurately the molecular signatures from [Verhaak et al. ([2010](#ref-verhaak_integrated_2010){role="doc-biblioref"})].

``` r
library(RColorBrewer)
subtypeOrder <- c("Proneural", "Neural", "Classical", "Mesenchymal")
sampleOrderBySubtype <- sort(match(gbm_es$subtype, subtypeOrder),
                             index.return=TRUE)$ix
subtypeXtable <- table(gbm_es$subtype)
subtypeColorLegend <- c(Proneural="red", Neural="green",
                        Classical="blue", Mesenchymal="orange")
geneSetOrder <- c("astroglia_up", "astrocytic_up", "neuronal_up",
                  "oligodendrocytic_up")
geneSetLabels <- gsub("_", " ", geneSetOrder)
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
hmcol <- hmcol[length(hmcol):1]

heatmap(exprs(gbm_es)[geneSetOrder, sampleOrderBySubtype], Rowv=NA,
        Colv=NA, scale="row", margins=c(3,5), col=hmcol,
        ColSideColors=rep(subtypeColorLegend[subtypeOrder],
                          times=subtypeXtable[subtypeOrder]),
        labCol="", gbm_es$subtype[sampleOrderBySubtype],
        labRow=paste(toupper(substring(geneSetLabels, 1,1)),
                     substring(geneSetLabels, 2), sep=""),
        cexRow=2, main=" \n ")
par(xpd=TRUE)
text(0.23,1.21, "Proneural", col="red", cex=1.2)
text(0.36,1.21, "Neural", col="green", cex=1.2)
text(0.47,1.21, "Classical", col="blue", cex=1.2)
text(0.62,1.21, "Mesenchymal", col="orange", cex=1.2)
mtext("Gene sets", side=4, line=0, cex=1.5)
mtext("Samples          ", side=1, line=4, cex=1.5)
```

Figure 3: [Heatmap of GSVA scores for cell-type brain signatures from murine models (y-axis) across GBM samples grouped by GBM subtype]\

## [7.2] Differential expression at pathway level

We illustrate here how to conduct a differential expression analysis at pathway level using the bulk RNA-seq data from [(Costa et al. [2021](#ref-costa2021genome){role="doc-biblioref"})]. This dataset consists of stranded 2x75nt paired-end reads sequenced from whole blood stored in [dried blood spots (DBS)](https://en.wikipedia.org/wiki/Dried_blood_spot). [Costa et al. ([2021](#ref-costa2021genome){role="doc-biblioref"})] generated these data from 21 DBS samples of extremely preterm newborns (neonates born before the 28th week of gestation), where 10 of them had been exposed to a fetal inflammatory response (FIR) before birth. A normalized matrix of logCPM units of expression of these data is stored in a [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment) object in the package *[GSVAdata](https://bioconductor.org/packages/3.22/GSVAdata)* and can be loaded as follows:

``` r
data(geneprotExpCostaEtAl2021)
se <- geneExpCostaEtAl2021
se
class: SummarizedExperiment
dim: 11279 21
metadata(0):
assays(1): logCPM
rownames(11279): 100 10000 ... 9994 9997
rowData names(1): Symbol
colnames(21): BS03 BS04 ... BS23 BS24
colData names(2): FIR Sex
```

To facilitate later on the automatic mapping of gene identifiers between gene sets and RNA-seq data, we should add annotation metadata to the `SummarizedExperiment` object as follows.

``` r
gsvaAnnotation(se) <- EntrezIdentifier("org.Hs.eg.db")
```

Here we have used the metadata constructor function `EntrezIdentifier()` because we can see that the gene identifiers in these expression data are entirely formed by numerical digits and these correspond to [NCBI Entrez gene identifiers](https://www.ncbi.nlm.nih.gov/gene). The sample (column) identifiers correspond to anonymized neonates, and the column (phenotype) metadata describes the exposure to FIR and the sex of the neonate. We can see that we have expression profiles for all four possible combinations of FIR exposure and sex.

``` r
colData(se)
DataFrame with 21 rows and 2 columns
          FIR      Sex
     <factor> <factor>
BS03      yes   female
BS04      yes   female
BS05      yes   male
BS06      no    female
BS07      no    female
...       ...      ...
BS20      no    female
BS21      yes   male
BS22      no    male
BS23      yes   male
BS24      no    male
table(colData(se))
     Sex
FIR   female male
  no       4    7
  yes      4    6
```

### [7.2.1] Data exploration at gene level

We do a brief data exploration at gene level, to have a sense of what we can expect in our analysis at pathway level. Figure [4](#fig:genelevelmds) below shows the projection in two dimensions of sample dissimilarity by means of a [multidimensional scaling (MDS)](https://en.wikipedia.org/wiki/Multidimensional_scaling) plot, produced with the `plotMDS()` function of the Bioconductor package *[limma](https://bioconductor.org/packages/3.22/limma)*. We can observe that sample dissimilarity in RNA expression from DBS samples is driven by the FIR and sex phenotypes, as shown in Fig. 1C of [Costa et al. ([2021](#ref-costa2021genome){role="doc-biblioref"})].

``` r
library(limma)

fircolor <- c(no="skyblue", yes="darkred")
sexpch <- c(female=19, male=15)
plotMDS(assay(se), col=fircolor[se$FIR], pch=sexpch[se$Sex])
```

Figure 4: [Gene-level exploration]\
Multidimensional scaling (MDS) plot at gene level. Red corresponds to `FIR=yes` and blue to `FIR=no`, while circles and squares correspond, respectively, to female and male neonates.

### [7.2.2] Filtering of immunologic gene sets

[Costa et al. ([2021](#ref-costa2021genome){role="doc-biblioref"})] report a postnatal activation of the innate immune system and an impairment of the adaptive immunity. For the purpose of exploring these results at pathway level, we will use the C7 collection of immunologic signature gene sets previously downloaded from the [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/human/collection_details.jsp#C7) database. We are going to futher filter this collection of gene sets to those formed by genes upregulated in innate leukocytes and adaptive mature lymphocytes, excluding those reported in studies on myeloid cells and the lupus autoimmune disease.

``` r
innatepat <- c("NKCELL_VS_.+_UP", "MAST_CELL_VS_.+_UP",
               "EOSINOPHIL_VS_.+_UP", "BASOPHIL_VS_.+_UP",
               "MACROPHAGE_VS_.+_UP", "NEUTROPHIL_VS_.+_UP")
innatepat <- paste(innatepat, collapse="|")
innategsets <- names(c7.genesets)[grep(innatepat, names(c7.genesets))]
length(innategsets)
[1] 53

adaptivepat <- c("CD4_TCELL_VS_.+_UP", "CD8_TCELL_VS_.+_UP", "BCELL_VS_.+_UP")
adaptivepat <- paste(adaptivepat, collapse="|")
adaptivegsets <- names(c7.genesets)[grep(adaptivepat, names(c7.genesets))]
excludepat <- c("NAIVE", "LUPUS", "MYELOID")
excludepat <- paste(excludepat, collapse="|")
adaptivegsets <- adaptivegsets[-grep(excludepat, adaptivegsets)]
length(adaptivegsets)
[1] 97

c7.genesets.filt <- c7.genesets[c(innategsets, adaptivegsets)]
length(c7.genesets.filt)
[1] 150
```

### [7.2.3] Running GSVA

To run GSVA on these data we build first the parameter object.

``` r
gsvapar <- gsvaParam(se, c7.genesets.filt, assay="logCPM", minSize=5,
                     maxSize=300)
```

Second, we run the GSVA algorithm by calling the `gsva()` function with the previoulsy built parameter object.

``` r
es <- gsva(gsvapar)
es
class: SummarizedExperiment
dim: 150 21
metadata(0):
assays(1): es
rownames(150):
  GSE18804_SPLEEN_MACROPHAGE_VS_BRAIN_TUMORAL_MACROPHAGE_UP
  GSE18804_SPLEEN_MACROPHAGE_VS_COLON_TUMORAL_MACROPHAGE_UP ...
  GSE7460_CD8_TCELL_VS_CD4_TCELL_ACT_UP
  GSE7460_CD8_TCELL_VS_TREG_ACT_UP
rowData names(1): gs
colnames(21): BS03 BS04 ... BS23 BS24
colData names(2): FIR Sex
```

Because the input expression data was provided in a `SummmarizedExperiment` object, the output of `gsva()` is again a `SummarizedExperiment` object, with two main differences with respect to the one given as input: (1) the one or more matrices of molecular data in the assay slot of the input object have been replaced by a single matrix of GSVA enrichment scores under the assay name `es`; and (2) the collection of mapped and filtered gene sets is included in the object and can be accessed using the methods `geneSets()` and `geneSetSizes()`.

``` r
assayNames(se)
[1] "logCPM"
assayNames(es)
[1] "es"
assay(es)[1:3, 1:3]
                                                                 BS03
GSE18804_SPLEEN_MACROPHAGE_VS_BRAIN_TUMORAL_MACROPHAGE_UP -0.06960698
GSE18804_SPLEEN_MACROPHAGE_VS_COLON_TUMORAL_MACROPHAGE_UP  0.28080100
GSE18804_SPLEEN_MACROPHAGE_VS_TUMORAL_MACROPHAGE_UP       -0.12573692
                                                                 BS04
GSE18804_SPLEEN_MACROPHAGE_VS_BRAIN_TUMORAL_MACROPHAGE_UP  0.05293804
GSE18804_SPLEEN_MACROPHAGE_VS_COLON_TUMORAL_MACROPHAGE_UP -0.13064317
GSE18804_SPLEEN_MACROPHAGE_VS_TUMORAL_MACROPHAGE_UP        0.16350068
                                                                 BS05
GSE18804_SPLEEN_MACROPHAGE_VS_BRAIN_TUMORAL_MACROPHAGE_UP -0.13353343
GSE18804_SPLEEN_MACROPHAGE_VS_COLON_TUMORAL_MACROPHAGE_UP -0.08696271
GSE18804_SPLEEN_MACROPHAGE_VS_TUMORAL_MACROPHAGE_UP       -0.16517679
```

``` r
head(lapply(geneSets(es), head))
$GSE18804_SPLEEN_MACROPHAGE_VS_BRAIN_TUMORAL_MACROPHAGE_UP
[1] "8540" "1645" "302"  "314"  "328"  "9582"

$GSE18804_SPLEEN_MACROPHAGE_VS_COLON_TUMORAL_MACROPHAGE_UP
[1] "51703" "136"   "29880" "240"   "64682" "313"

$GSE18804_SPLEEN_MACROPHAGE_VS_TUMORAL_MACROPHAGE_UP
[1] "22"     "9429"   "55860"  "1645"   "347902" "302"

$GSE22886_NEUTROPHIL_VS_DC_UP
[1] "101"    "84658"  "222487" "2334"   "57538"  "200315"

$GSE22886_NEUTROPHIL_VS_MONOCYTE_UP
[1] "8748"   "222487" "306"    "25825"  "699"    "762"

$GSE2585_THYMIC_MACROPHAGE_VS_MTEC_UP
[1] "10347" "10449" "2015"  "3268"  "10327" "334"
```

``` r
head(geneSetSizes(es))
GSE18804_SPLEEN_MACROPHAGE_VS_BRAIN_TUMORAL_MACROPHAGE_UP
                                                      122
GSE18804_SPLEEN_MACROPHAGE_VS_COLON_TUMORAL_MACROPHAGE_UP
                                                      141
      GSE18804_SPLEEN_MACROPHAGE_VS_TUMORAL_MACROPHAGE_UP
                                                      146
                             GSE22886_NEUTROPHIL_VS_DC_UP
                                                      105
                       GSE22886_NEUTROPHIL_VS_MONOCYTE_UP
                                                       82
                     GSE2585_THYMIC_MACROPHAGE_VS_MTEC_UP
                                                      142
```

### [7.2.4] Data exploration at pathway level

We do again a data exploration, this time at pathway level. Figure [5](#fig:pathwaylevelmds) below, shows an MDS plot of GSVA enrichment scores. We can see again that most variability is driven by the FIR phenotype, but this time the sex phenotype does not seem to affect sample dissimilarity at pathway level, probably because the collection of gene sets we have used does not include gene sets formed by genes with sex-specific expression.

``` r
plotMDS(assay(es), col=fircolor[es$FIR], pch=sexpch[es$Sex])
```

Figure 5: [Pathway-level exploration]\
Multidimensional scaling (MDS) plot at pathway level. Red corresponds to `FIR=yes` and blue to `FIR=no`, while circles and squares correspond, respectively, to female and male neonates.

### [7.2.5] Differential expression

We will perform now a differential expression analysis at pathway level using the Bioconductor packages *[limma](https://bioconductor.org/packages/3.22/limma)* [(Smyth [2004](#ref-Smyth_2004){role="doc-biblioref"})] and *[sva](https://bioconductor.org/packages/3.22/sva)*, the latter to adjust for sample heterogeneity using surrogate variable analysis [(Leek and Storey [2007](#ref-leek2007capturing){role="doc-biblioref"})]

``` r
library(sva)
library(limma)

## build design matrix of the model to which we fit the data
mod <- model.matrix(~ FIR, colData(es))
## build design matrix of the corresponding null model
mod0 <- model.matrix(~ 1, colData(es))
## estimate surrogate variables (SVs) with SVA
sv <- sva(assay(es), mod, mod0)
Number of significant surrogate variables is:  4
Iteration (out of 5 ):1  2  3  4  5
## add SVs to the design matrix of the model of interest
mod <- cbind(mod, sv$sv)
## fit linear models
fit <- lmFit(assay(es), mod)
## calculate moderated t-statistics using the robust regime
fit.eb <- eBayes(fit, robust=TRUE)
## summarize the extent of differential expression at 5% FDR
res <- decideTests(fit.eb)
summary(res)
       (Intercept) FIRyes
Down            25     14  38  41  33  53
NotSig         119     97  65  78  70  94
Up               6     39  47  31  47   3
```

As shown in Figure [6](#fig:esstdevxgssize) below, GSVA scores tend to have higher precision for larger gene sets11[[1] Thanks to Gordon Smyth for pointing this out to us.]{#sidenote-1 .sidenote}, albeit this trend breaks at the end of gene set sizes in this case. This trend is usually more clear when GSVA scores are derived from gene sets including smaller sizes (our smallest gene set here is about 100 genes), and from less heterogenous expression data. Here we use the getter method `geneSetSizes()` to obtain the vector of sizes of the gene sets filtered after the GSVA calculations on the output of the `gsva()` function.

``` r
gssizes <- geneSetSizes(es)
plot(sqrt(gssizes), sqrt(fit.eb$sigma), xlab="Sqrt(gene sets sizes)",
          ylab="Sqrt(standard deviation)", las=1, pch=".", cex=4)
lines(lowess(sqrt(gssizes), sqrt(fit.eb$sigma)), col="red", lwd=2)
```

Figure 6: [Pathway-level differential expression analysis]\
Residual standard deviation of GSVA scores as a function of gene set size. Larger gene sets tend to have higher precision.

When this trend is present, we may improve the statistical power to detect differentially expressed (DE) pathways by using the limma-trend pipeline [(Phipson et al. [2016](#ref-phipson2016robust){role="doc-biblioref"})]. More concretely, we should call the `eBayes()` function with the argument `trend=x`, where `x` is a vector of values corresponding to the square root of the sizes of the gene sets. As we have already seen, the values of these sizes can be easily obtained using GSVA's function `geneSetSizes()` on the output of the `gsva()` function. Here below, we call again `eBayes()` using the `trend` parameter. In this case, however, the change in the number of FIR DE pathways is negligible.

``` r
fit.eb.trend <- eBayes(fit, robust=TRUE, trend=sqrt(gssizes))
res <- decideTests(fit.eb.trend)
summary(res)
       (Intercept) FIRyes
Down            28     15  37  41  33  58
NotSig         116     94  66  75  70  88
Up               6     41  47  34  47   4
```

We can select DE pathways with FDR \< 5% as follows.

``` r
tt <- topTable(fit.eb.trend, coef=2, n=Inf)
DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.05]
length(DEpwys)
[1] 56
head(DEpwys)
[1] "GSE3039_CD4_TCELL_VS_ALPHABETA_CD8_TCELL_UP"
[2] "GSE3982_NEUTROPHIL_VS_CENT_MEMORY_CD4_TCELL_UP"
[3] "GSE3982_EOSINOPHIL_VS_CENT_MEMORY_CD4_TCELL_UP"
[4] "GSE36392_EOSINOPHIL_VS_NEUTROPHIL_IL25_TREATED_LUNG_UP"
[5] "GSE3982_NEUTROPHIL_VS_TH2_UP"
[6] "GSE22886_NEUTROPHIL_VS_DC_UP"
```

Figure [7](#fig:heatmapdepwys) below shows a heatmap of the GSVA enrichment scores of the subset of the 56 DE pathways, clustered by pathway and sample. We may observe that, consistently with the findings of [Costa et al. ([2021](#ref-costa2021genome){role="doc-biblioref"})], FIR-affected neonates display an enrichment of upregulated pathways associated with innate immunity, and an enrichment of downregulated pathways associated with adaptive immunity, with respect to FIR-unaffected neonates.

``` r
## get DE pathway GSVA enrichment scores, removing the covariates effect
DEpwys_es <- removeBatchEffect(assay(es[DEpwys, ]),
                               covariates=mod[, 2:ncol(mod)],
                               design=mod[, 1:2])
Coefficients not estimable: FIRyes
## cluster samples
sam_col_map <- fircolor[es$FIR]
names(sam_col_map) <- colnames(DEpwys_es)
sampleClust <- hclust(as.dist(1-cor(DEpwys_es, method="spearman")),
                      method="complete")

## cluster pathways
gsetClust <- hclust(as.dist(1-cor(t(DEpwys_es), method="pearson")),
                    method="complete")

## annotate pathways whether they are involved in the innate or in
## the adaptive immune response
labrow <- rownames(DEpwys_es)
mask <- rownames(DEpwys_es) %in% innategsets
labrow[mask] <- paste("(INNATE)", labrow[mask], sep="_")
mask <- rownames(DEpwys_es) %in% adaptivegsets
labrow[mask] <- paste("(ADAPTIVE)", labrow[mask], sep="_")
labrow <- gsub("_", " ", gsub("GSE[0-9]+_", "", labrow))

## pathway expression color scale from blue (low) to red (high)
library(RColorBrewer)
pwyexpcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
pwyexpcol <- pwyexpcol[length(pwyexpcol):1]

## generate heatmap
heatmap(DEpwys_es, ColSideColors=fircolor[es$FIR], xlab="Samples",
        ylab="Pathways", margins=c(2, 20), labCol="", labRow=labrow,
        col=pwyexpcol, scale="row", Colv=as.dendrogram(sampleClust),
        Rowv=as.dendrogram(gsetClust))
```

Figure 7: [Pathway-level signature of FIR]\
Heatmap of GSVA enrichment scores from pathways being called DE with 5% FDR between FIR-affected and unaffected neonates.

# [8] Interactive web app

The `gsva()` function can be also used through an interactive web app developed with *[shiny](https://CRAN.R-project.org/package=shiny)*. To start it just type on the R console:

``` r
res <- igsva()
```

It will open your browser with the web app shown here below. The button `SAVE & CLOSE` will close the app and return the resulting object on the R console. Hence, the need to call igsva() on the right-hand side of an assignment if you want to store the result in your workspace. Alternatively, you can use the `DOWNLOAD` button to download the result in a CSV file.

In the starting window of the web app, after running GSVA, a non-parametric kernel density estimation of sample profiles of GSVA scores will be shown. By clicking on one of the lines, the cumulative distribution of GSVA scores for the corresponding samples will be shown in the `GeneSets` tab, as illustrated in the image below.

# [9] Contributing

GSVA has benefited from contributions by multiple developers, see <https://github.com/rcastelo/GSVA/graphs/contributors> for a list of them. Contributions to the software codebase of GSVA are welcome as long as contributors abide to the terms of the [Bioconductor Contributor Code of Conduct](https://bioconductor.org/about/code-of-conduct). If you want to contribute to the development of GSVA please open an [issue](https://github.com/rcastelo/GSVA/issues) to start discussing your suggestion or, in case of a bugfix or a straightforward feature, directly a [pull request](https://github.com/rcastelo/GSVA/pulls).

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
 [1] sva_3.58.0                  BiocParallel_1.44.0
 [3] genefilter_1.92.0           mgcv_1.9-4
 [5] nlme_3.1-168                RColorBrewer_1.1-3
 [7] edgeR_4.8.2                 limma_3.66.0
 [9] GSVAdata_1.46.0             SummarizedExperiment_1.40.0
[11] GenomicRanges_1.62.1        Seqinfo_1.0.0
[13] MatrixGenerics_1.22.0       matrixStats_1.5.0
[15] hgu95a.db_3.13.0            GSEABase_1.72.0
[17] graph_1.88.1                annotate_1.88.0
[19] XML_3.99-0.22               org.Hs.eg.db_3.22.0
[21] AnnotationDbi_1.72.0        IRanges_2.44.0
[23] S4Vectors_0.48.0            Biobase_2.70.0
[25] BiocGenerics_0.56.0         generics_0.1.4
[27] GSVA_2.4.8                  BiocStyle_2.38.0

loaded via a namespace (and not attached):
 [1] blob_1.3.0                  Biostrings_2.78.0
 [3] fastmap_1.2.0               SingleCellExperiment_1.32.0
 [5] digest_0.6.39               rsvd_1.0.5
 [7] lifecycle_1.0.5             survival_3.8-6
 [9] statmod_1.5.1               KEGGREST_1.50.0
[11] RSQLite_2.4.6               magrittr_2.0.4
[13] compiler_4.5.2              rlang_1.1.7
[15] sass_0.4.10                 tools_4.5.2
[17] yaml_2.3.12                 knitr_1.51
[19] S4Arrays_1.10.1             bit_4.6.0
[21] DelayedArray_0.36.0         abind_1.4-8
[23] HDF5Array_1.38.0            memuse_4.2-3
[25] grid_4.5.2                  beachmat_2.26.0
[27] xtable_1.8-8                Rhdf5lib_1.32.0
[29] tinytex_0.58                cli_3.6.5
[31] rmarkdown_2.30              crayon_1.5.3
[33] otel_0.2.0                  httr_1.4.8
[35] rjson_0.2.23                DBI_1.3.0
[37] cachem_1.1.0                rhdf5_2.54.1
[39] splines_4.5.2               parallel_4.5.2
[41] BiocManager_1.30.27         XVector_0.50.0
[43] vctrs_0.7.1                 Matrix_1.7-4
[45] jsonlite_2.0.0              bookdown_0.46
[47] BiocSingular_1.26.1         bit64_4.6.0-1
[49] irlba_2.3.7                 magick_2.9.1
[51] locfit_1.5-9.12             h5mread_1.2.1
[53] jquerylib_0.1.4             glue_1.8.0
[55] codetools_0.2-20            ScaledMatrix_1.18.0
[57] pillar_1.11.1               htmltools_0.5.9
[59] rhdf5filters_1.22.0         R6_2.6.1
[61] sparseMatrixStats_1.22.0    evaluate_1.0.5
[63] lattice_0.22-9              png_0.1-9
[65] SpatialExperiment_1.20.0    memoise_2.0.1
[67] bslib_0.10.0                Rcpp_1.1.1
[69] SparseArray_1.10.9          xfun_0.56
[71] pkgconfig_2.0.3
```

# References

Barbie, David A., Pablo Tamayo, Jesse S. Boehm, So Young Kim, Susan E. Moody, Ian F. Dunn, Anna C. Schinzel, et al. 2009. "Systematic RNA Interference Reveals That Oncogenic KRAS-driven Cancers Require TBK1." *Nature* 462 (7269): 108--12. <https://doi.org/10.1038/nature08460>.

Cahoy, John D., Ben Emery, Amit Kaushal, Lynette C. Foo, Jennifer L. Zamanian, Karen S. Christopherson, Yi Xing, et al. 2008. "A Transcriptome Database for Astrocytes, Neurons, and Oligodendrocytes: A New Resource for Understanding Brain Development and Function." *J. Neurosci.* 28 (1): 264--78. <https://doi.org/10.1523/JNEUROSCI.4178-07.2008>.

Carrel, Laura, and Huntington F Willard. 2005. "X-Inactivation Profile Reveals Extensive Variability in X-Linked Gene Expression in Females." *Nature* 434 (7031): 400--404.

Costa, Daniel, Núria Bonet, Amanda Solé, José Manuel González de Aledo-Castillo, Eduard Sabidó, Ferran Casals, Carlota Rovira, et al. 2021. "Genome-Wide Postnatal Changes in Immunity Following Fetal Inflammatory Response." *The FEBS Journal* 288 (7): 2311--31. <https://doi.org/10.1111/febs.15578>.

Hänzelmann, Sonja, Robert Castelo, and Justin Guinney. 2013. "GSVA: Gene Set Variation Analysis for Microarray and RNA-Seq Data." *BMC Bioinformatics* 14: 7. <https://doi.org/10.1186/1471-2105-14-7>.

Huang, R Stephanie, Shiwei Duan, Wasim K Bleibel, Emily O Kistner, Wei Zhang, Tyson A Clark, Tina X Chen, et al. 2007. "A Genome-Wide Approach to Identify Genetic Variants That Contribute to Etoposide-Induced Cytotoxicity." *Proceedings of the National Academy of Sciences of the United States of America* 104 (23): 9758--63.

Lee, Eunjung, Han-Yu Chuang, Jong-Won Kim, Trey Ideker, and Doheon Lee. 2008. "Inferring Pathway Activity Toward Precise Disease Classification." *PLOS Comput Biol* 4 (11): e1000217. <https://doi.org/10.1371/journal.pcbi.1000217>.

Leek, Jeffrey T, and John D Storey. 2007. "Capturing Heterogeneity in Gene Expression Studies by Surrogate Variable Analysis." *PLoS Genetics* 3 (9): e161. <https://doi.org/10.1371/journal.pgen.0030161>.

Phipson, Belinda, Stanley Lee, Ian J Majewski, Warren S Alexander, and Gordon K Smyth. 2016. "Robust Hyperparameter Estimation Protects Against Hypervariable Genes and Improves Power to Detect Differential Expression." *The Annals of Applied Statistics* 10 (2): 946. <https://doi.org/10.1214/16-AOAS920>.

Pickrell, Joseph K., John C. Marioni, Athma A. Pai, Jacob F. Degner, Barbara E. Engelhardt, Everlyne Nkadori, Jean-Baptiste Veyrieras, Matthew Stephens, Yoav Gilad, and Jonathan K. Pritchard. 2010. "Understanding Mechanisms Underlying Human Gene Expression Variation with RNA Sequencing." *Nature* 464 (7289): 768--72.

Skaletsky, Helen, Tomoko Kuroda-Kawaguchi, Patrick J. Minx, Holland S. Cordum, LaDeana Hillier, Laura G. Brown, Sjoerd Repping, et al. 2003. "The Male-Specific Region of the Human Y Chromosome Is a Mosaic of Discrete Sequence Classes." *Nature* 423 (6942): 825--37.

Smyth, Gordon K. 2004. "Linear Models and Empirical Bayes Methods for Assessing Differential Expression in Microarray Experiments." *Stat Appl Genet Mol Biol* 3: Article3. <https://doi.org/10.2202/1544-6115.1027>.

Subramanian, Aravind, Pablo Tamayo, Vamsi K. Mootha, Sayan Mukherjee, Benjamin L. Ebert, Michael A. Gillette, Amanda Paulovich, et al. 2005. "Gene Set Enrichment Analysis: A Knowledge-Based Approach for Interpreting Genome-Wide Expression Profiles." *Proc. Natl. Acad. Sci. U.S.A.* 102 (43): 15545--50. <https://doi.org/10.1073/pnas.0506580102>.

Tomfohr, John, Jun Lu, and Thomas B Kepler. 2005. "Pathway Level Analysis of Gene Expression Using Singular Value Decomposition." *BMC Bioinformatics* 6 (September): 225. <https://doi.org/10.1186/1471-2105-6-225>.

Verhaak, Roel G W, Katherine A Hoadley, Elizabeth Purdom, Victoria Wang, Yuan Qi, Matthew D Wilkerson, C Ryan Miller, et al. 2010. "Integrated Genomic Analysis Identifies Clinically Relevant Subtypes of Glioblastoma Characterized by Abnormalities in PDGFRA, IDH1, EGFR, and NF1." *Cancer Cell* 17 (1): 98--110. <https://doi.org/10.1016/j.ccr.2009.12.020>.
