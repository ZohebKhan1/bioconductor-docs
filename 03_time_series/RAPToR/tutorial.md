
# RAPToR-showcase

Code []

-   [Show All Code](#)
-   [Hide All Code](#)

# `RAPToR` - Showcase

### RAPToR 1.2.0

Romain Bulteau

#### March 2023

# Contents

-   [[1] Improving Differential Expression analysis with age estimates](#improving-differential-expression-analysis-with-age-estimates)
    -   [[1.1] Data and strategy](#data-and-strategy)
    -   [[1.2] Estimating sample age](#estimating-sample-age)
    -   [[1.3] Impact of development on gene expression](#impact-of-development-on-gene-expression)
    -   [[1.4] Differential Expression analysis](#differential-expression-analysis)
    -   [[1.5] Further support for the increase in performance](#further-support-for-the-increase-in-performance)
    -   [[1.6] Functions and code](#functions-and-code)
        -   [[1.6.1] Code to generate objects](#code-to-generate-objects)
        -   [[1.6.2] DE functions](#de-functions)
-   [[2] Staging samples cross-species](#staging-samples-cross-species)
    -   [[2.1] Data and strategy](#data-and-strategy-1)
    -   [[2.2] Estimating the age of *C. elegans* embryos on a *D. melanogaster* reference](#estimating-the-age-of-c.-elegans-embryos-on-a-d.-melanogaster-reference)
    -   [[2.3] Adjusting the reference](#adjusting-the-reference)
    -   [[2.4] About the outliers](#about-the-outliers)
    -   [[2.5] Code to generate objects](#code-to-generate-objects-1)
-   [[3] Capturing tissue-specific development](#capturing-tissue-specific-development)
    -   [[3.1] Data and strategy](#data-and-strategy-2)
    -   [[3.2] Estimating sample global age](#estimating-sample-global-age)
    -   [[3.3] Caracterization of expression dynamics](#caracterization-of-expression-dynamics)
    -   [[3.4] Estimating tissue-specific age](#estimating-tissue-specific-age)
    -   [[3.5] Code to generate objects](#code-to-generate-objects-2)
-   [References](#references)
-   [SessionInfo](#sessioninfo)

\

This vignette showcases different uses for `RAPToR`.

\

# [1] Improving Differential Expression analysis with age estimates

Including estimated age as a covariate in Differential Expression (DE) analysis can substantially reduce previously unexplained variation between samples. In this example, we show that even when chronological age of different time points is known, using estimated age results in better model fits and consequently better power to detect DE genes.

## [1.1] Data and strategy

[Lehrbach et al. (2012)] profiled *C. elegans* wild-type (`wt`) and *pash-1(mj100)* mutants (`mut`) at 4 time points (0, 6, 12, and 24 hours past the L4 stage), each in triplicate (Accession: [E-MTAB-1333](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1333/), `dslehrbach2012`).

Given the time-series experimental design, searching for Differentially Expresses (DE) genes between strains requires including development in the model. While true in any time series context, it is particularly important here as late-larval development of *C. elegans* is known for drastic changes in gene expression within short time frames ([Snoek et al. (2014)]).

Using identical DE models with either chronological age or RAPToR age estimates as covariates, we can determine which of the two is the best predictor and thus gives more power to detect differential expression.

Code to generate the `dslehrbach2012` object can be found [at the end of this section](#code-to-generate-objects).

``` r
library(RAPToR)
library(wormRef)

library(stats)
library(parallel)
library(limma)
```

\

## [1.2] Estimating sample age

We start by applying a quantile-normalization and \\(log(X+1)\\) transformation to the data.

``` r
dslehrbach2012$g <- limma::normalizeBetweenArrays(dslehrbach2012$g,
                                                  method = "quantile")
dslehrbach2012$g <- log1p(dslehrbach2012$g)
```

Next, we select a reference and stage the samples. Since sample (chronological) age ranges from the fourth larval molt (L4) to 24h past L4, we select `Cel_YA_2` from `wormRef` that covers this range.

```
# load reference
r_ya <- prepare_refdata("Cel_YA_2", "wormRef", 400)

# estimate sample age
ae_lerhbach <- ae(dslehrbach2012$g, r_ya)
#>                 nb.genes
#> refdata            15710
#> samp               17441
#> intersect.genes    14952
#> Bootstrap set size is 4984
#> Performing age estimation...
#> Bootstrapping...
#>  Building gene subsets...
#>  Computing correlations...
#>  Performing age estimation...
#> Computing summary statistics...
```

We can see quite a lot of variation in the sample timings, especially at the 0- and 6-hour time points.

Another curious effect can be noted, easier to see when grouping the samples by batch. In the first two replicates, mutants are systematically (and statistically significantly) older than WT, while the opposite effect is true for the third replicate.

```
# compute age difference between wt & mut at each time point
isWT <- dslehrbach2012$p$strain=="wt"
ae_dif <- ae_lerhbach$age.estimates[isWT,1] - ae_lerhbach$age.estimates[!isWT,1]

# test for effect significance with a linear model
batch <- dslehrbach2012$p$rep[isWT]
summary(lm(ae_dif~batch))
#>
#> Call:
#> lm(formula = ae_dif ~ batch)
#>
#> Residuals:
#>     Min      1Q  Median      3Q     Max
#> -3.4889 -0.3440  0.0491  0.7862  2.4079
#>
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)  -2.1130     0.7867  -2.686   0.0250 *
#> batch2       -1.0811     1.1126  -0.972   0.3566
#> batch3        3.3169     1.1126   2.981   0.0154 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Residual standard error: 1.573 on 9 degrees of freedom
#> Multiple R-squared:  0.6535, Adjusted R-squared:  0.5765
#> F-statistic: 8.486 on 2 and 9 DF,  p-value: 0.008487
```

## [1.3] Impact of development on gene expression

Gene expression depends so strongly on development that correlation between samples is only predicted by their age difference, rather than by their genetic background.

The graph below shows correlation between all possible pairs of samples, with no clear effect of strain difference.

``` r
cc <- cor(dslehrbach2012$g, dslehrbach2012$g, method = "spearman")
```

## [1.4] Differential Expression analysis

We will use `limma` to perform a DE analysis with the following model:

[\\[ X \\sim strain \\times \\mathtt{ns(age, df = 2)} \\]] where \\(\\mathtt{age}\\) will either be chronological or estimated age, and the spline \\(\\mathtt{ns()}\\) will be used to handle non-linear expression dynamics along time.

To find differential expression of genes along development or between strains translates to comparing the following nested models: (2) vs. (1) for development, (3) vs. (1) for strain.

[\\[ \\begin{aligned} Y & = \\beta_0 + \\beta_1 I\_{strain} + (\\alpha_1 age\_{sp1} + \\alpha_2 age\_{sp2}) + (\\gamma_1 I\_{strain} age\_{sp1} + \\gamma_2 I\_{strain} age\_{sp2})\\\\ Y & = \\beta_0 + \\beta_1 I\_{strain}\\\\ Y & = \\beta_0 + (\\alpha_1 age\_{sp1} + \\alpha_2 age\_{sp2}) \\end{aligned} \\]] Genes will be considered DE with Benjamini-Holm adjusted p-values \\(\< 0.05\\) for the corresponding coefficients.

Note: for the purpose of this vignette, we don't consider whether genes are up- or down-regulated, but only if a significant effect is detected.

\

We now run the same DE analysis using either chronological age or RAPToR age estimates as predictor, and compare results to determine if there is an improvement. Code for `DGE()` can be found [at the end of this section](#de-functions).

```
# format gdata as log2 for limma input
X <- log2(exp(dslehrbach2012$g))

# with chron. age
dge.ca <- DGE(X = X, strain = dslehrbach2012$p$strain,
              age = dslehrbach2012$p$tpastL4,
              name = "dge.ca", return.model = T)
#> Loading required package: splines

# with RAPToR ae
dge.ae <- DGE(X = X, strain = dslehrbach2012$p$strain,
              age = ae_lerhbach$age.estimates[,1],
              name = "dge.ae", return.model = T)
```

Estimated age shows a clear increase in performance over chronological age when comparing adjusted p-values of each gene for detection of an effect. Red bars in the plots above correspond to the \\(0.05\\) threshold for significance and color-coded gene counts for each quadrant are indicated in the bottom-right.

``` r
# compute nb. of DE genes for mutation
mut_dif.ca <- sum(dge.ca$tmut$adj.P.Val < 0.05)
mut_dif.ae <- sum(dge.ae$tmut$adj.P.Val < 0.05)

# compute nb. of DE genes for development
age_dif.ca <- sum(dge.ca$tage$adj.P.Val < 0.05)
age_dif.ae <- sum(dge.ae$tage$adj.P.Val < 0.05)

# percentage increase of mutation DE genes
100 * (mut_dif.ae - mut_dif.ca) / (mut_dif.ca) # mut pct increase
#> [1] 96.02237

# percentage increase of development DE genes
100 * (age_dif.ae - age_dif.ca) / (age_dif.ca) # age pct increase
#> [1] 7.674787
```

Indeed, we detect nearly twice as many DE genes for strain using RAPToR age estimates compared to chronological age, and around \\(8\\%\\) more genes with development. Furthermore, the large proportion of genes changing with development (\\(57 - 62\\%\\)) compared to strain (\\(7-16\\%\\)) also shows how crucial it is to properly model developmental dynamics in DE analyses.

The increase in detected DE genes is partially explained by better model fits, shown below by the goodness-of-fit (GoF) distribution across genes. The goodness-of-fit (GoF) computed is an \\(R\^2 = 1 - \\frac{SS\_{res}}{SS\_{tot}}\\) per gene.

## [1.5] Further support for the increase in performance

If the asynchronicity that RAPToR detects between samples is erroneous, then adding random noise around the chronological age values should yield similar results to using our age estimates in the model. To test this, we simulate age sets with added noise of similar distribution to the age differences observed between chronological and estimated age.

```
set.seed(10) # for reproducibility

age_diffs <- (dslehrbach2012$p$tpastL4 + 50) - ae_lerhbach$age.estimates[,1]
# Note : we add 50 to shift tpastL4 to the age values and avoid negative
#        tpastl4 values. This has no impact on the DE analysis.

# estimate density function of age_diffs
d_ad <- density(age_diffs)
# generate 100 age sets of with random age_diffs-like noise
n <- 100
rd_ages <- lapply(seq_len(n), function(i){
  (dslehrbach2012$p$tpastL4 + 50) +
    sample(x = d_ad$x, size = nrow(dslehrbach2012$p),
           prob = d_ad$y, replace = T)
})
```

As done above, we run identical DE models and compare results with the model goodness-of-fit (GoF) per gene and look at the number of DE genes found for strain and development (BH-adjusted p-value \\(\< 0.05\\)).

```
# setup cluster for parallelization
cl <- parallel::makeCluster(6, "FORK")

# do DGE on all age sets
rd_dges <- parLapply(cl, seq_len(n), function(i){
  cat("\r", i,"/",n)
  DGE(X, dslehrbach2012$p$strain, rd_ages[[i]], name = paste0("rd.",i))
})

stopCluster(cl)
gc()
```

``` r
# get quantiles of GoF for plotting
qts <- seq(0,1, length.out = 100)
quants <- lapply(seq_len(n), function(i){
  quantile(rd_dges[[i]]$gof, probs = qts)
})
quants <- do.call(rbind, quants)
```

We can see that using estimated age systematically leads to better model fits and increased detection of DE genes both for strain and development.

Gene GoF using estimated age and binned by chronological age GoF quantiles (below, left) also clearly shows that for the same genes, we tend to have better fits with estimated age, while that is not the case when adding noise.

## [1.6] Functions and code

### [1.6.1] Code to generate objects

``` r
data_folder <- "../inst/extdata/"

requireNamespace("wormRef", quietly = T)
requireNamespace("utils", quietly = T)

requireNamespace("biomaRt", quietly = T) # bioconductor
requireNamespace("limma", quietly = T)   # bioconductor
requireNamespace("affy", quietly = T)    # bioconductor
requireNamespace("gcrma", quietly = T)   # bioconductor
```

*Note : set the `data_folder` variable to an existing path on your system where you want to store the objects.*

To generate `dslehrbach2012`:

``` r
# get probe ids from the arrayexpress
probe_ids <- read.table(
  "https://www.ebi.ac.uk/arrayexpress/files/A-AFFY-60/A-AFFY-60.adf.txt",
  h=T, comment.char = "", sep = "\t", skip = 17, as.is = T)

# pheno data
p_url <- paste0("https://www.ebi.ac.uk/arrayexpress/files/",
                "E-MTAB-1333/E-MTAB-1333.sdrf.txt")
p <- read.table(p_url, h = T, sep = "\t", comment.char = "", as.is = T)

# geno data
g_url <- unique(p$Comment..ArrayExpress.FTP.file.)
g_dir <- paste0(data_folder, "dslehrbach2012_raw/")
g_file <- paste0(g_dir, "raw.zip")

if(!file.exists(g_file)){
  dir.create(g_dir)
  utils::download.file(g_url, destfile = g_file)
}
utils::unzip(g_file, exdir = g_dir)

g <- affy::ReadAffy(filenames = p$Array.Data.File,
                    celfile.path = g_dir, phenoData = p)
g <- affy::expresso(g, bg.correct = F, normalize = F,
                    pmcorrect.method = "pmonly", summary.method = "median")
g <- 2^exprs(g) # expresso log2s the data
colnames(g) <- p$Source.Name

# format ids
g <- RAPToR::format_ids(g, probe_ids, from = 1, to = 7)
g <- RAPToR::format_ids(g, wormRef::Cel_genes, from = 3, to = 1)

# filter relevant fields
p <- p[, c(1,22,24)]
colnames(p) <- c("title", "tpastL4", "strain")
p$strain <- factor(p$strain,
                   levels = c('wild type', 'pash-1(mj100)'),
                   labels = c('wt', 'strain'))
p$rep <- factor(gsub("\\w+_h\\d+\\.(\\d)","\\1", p$title))

dslehrbach2012 <- list(g = g, p = p)
save(dslehrbach2012,
     file = paste0(data_folder, "dslerhbach2012.RData"), compress = "xz")

# cleanup
file.remove(g_file)
unlink(g_dir, recursive = T)

rm(g, g_url, g_dir, g_file, p, p_url, probe_ids)
```

### [1.6.2] DE functions

Functions to compute predictions (`pred_lmFit()`) and Goodness-of-fit (`GoF()`) from `limma` models.

``` r
pred_lmFit <- function(Fit){
# Predictions from a limma model
  tcrossprod(Fit$coefficients, Fit$design)
}

GoF <- function(Fit, X){
# Compute Goodness of Fit
  pred <- pred_lmFit(Fit)
  res <- (X - pred)
  ss <- apply(X, 1, function(ro) sum((ro - mean(ro))^2))
  Rsq <- sapply(seq_len(nrow(X)), function(i){
    1 - sum(res[i,]^2)/ss[i]
  })
  return(Rsq)
}
```

Function to run DE analysis with `limma`, as described in [*Differential Expression analysis*](#differential-expression-analysis) (`DGE()`).

``` r
DGE <- function(X, strain, age, df = 2, name = NULL, return.model = FALSE){
  require(splines)
  if(! length(strain) == ncol(X) | ! length(age) == ncol(X))
    stop("strain and age must be of length ncol(X).")

  # make pdat df
  pdat <- data.frame(strain = factor(strain),
                     age = as.numeric(age),
                     row.names = colnames(X))

  # build design matrix
  d <- model.matrix(~ 1 + ns(age, df = df) * strain, data = pdat)

  # fix colnames
  colnames(d) <- c("b0", paste0(rep("a", df), 1L:df),
                   "strainmut", paste0(rep("strainmut.a", df), 1L:df))

  # build contrast matrices for mut and age tests
  if(df == 2){
    cm.mut <- makeContrasts(mut = strainmut,
                            mut.i1 = strainmut.a1,
                            mut.i2 = strainmut.a2,
                            levels = d)
    cm.age <- makeContrasts(a1, a2,
                            a.i1 = strainmut.a1,
                            a.i2 = strainmut.a2,
                            levels = d)
  }

  if(df == 3){
    cm.mut <- makeContrasts(mut = strainmut,
                            mut.i1 = strainmut.a1,
                            mut.i2 = strainmut.a2,
                            mut.i3 = strainmut.a3,
                            levels = d)
    cm.age <- makeContrasts(a1, a2, a3,
                            a.i1 = strainmut.a1,
                            a.i2 = strainmut.a2,
                            a.i3 = strainmut.a3,
                            levels = d)
  }

  # fit model
  m.0 <- lmFit(object = X, design = d)
  # get GoF
  gof <- GoF(Fit = m.0, X = X)

  # find DE genes for mut
  m.m <- contrasts.fit(m.0, contrasts = cm.mut)
  m.m <- eBayes(m.m)

  # find DE genes for age
  m.a <- contrasts.fit(m.0, contrasts = cm.age)
  m.a <- eBayes(m.a)

  Tmut <- topTable(m.m, adjust.method = "BH",
                   number = Inf,
                   sort.by = "none")[, c("F", "P.Value", "adj.P.Val")]
  Tage <- topTable(m.a, adjust.method = "BH",
                   number = Inf,
                   sort.by = "none")[, c("F", "P.Value", "adj.P.Val")]

  res <- list(gof  = gof,
              tmut = Tmut,
              tage = Tage,
              name = name)
  if(return.model)
    res$model = m.0

  rm(m.0, m.m, m.a, Tmut, Tage, d, X, pdat, gof)
  gc(verbose = F)
  return(res)
}
```

\
\

# [2] Staging samples cross-species

Your organism of interest may not be well-studied or have an abundance of reference time-series data available. However, RAPToR still works when using a close organism as a reference, thanks to the conserved nature of developmental processes across species (particularly in early development).

Indeed, samples can be staged cross-species using ortholog genes.

## [2.1] Data and strategy

[Li et al. (2014)] defined a set of orthologs between *D. melanogaster* and *C. elegans* (hereafter, `glist`), which we will use to stage *C. elegans* single embryos on a *D. melanogaster* reference (of note, ensembl ortholog sets also work).

The *C. elegans* time-series of single-embryos was profiled and published by [Levin et al. (2016)] (Accession : [GSE60755](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60755), `dslevin2016cel`). The *Drosophila melanogaster* embryonic development time-series is part of the modENCODE project and published by [Graveley et al. (2011)] (data downloaded from [fruitfly.org](https://fruitfly.org/sequence/download.html), `dsgraveley2011`)

Code to generate `glist` (orthologs), `dslevin2016cel` (*C. elegans* data), and `dsgraveley2011` (*D. melanogaster* data) can be found [at the end of this section](#code-to-generate-objects-1).

``` r
library(RAPToR)
library(drosoRef)

library(limma)
library(stats)
```

## [2.2] Estimating the age of *C. elegans* embryos on a *D. melanogaster* reference

We start by applying a quantile-normalization and \\(log(X+1)\\) transformation to the data.

``` r
dsgraveley2011$g <- limma::normalizeBetweenArrays(dsgraveley2011$g,
                                                  method = "quantile")
dsgraveley2011$g <- log1p(dsgraveley2011$g)

dslevin2016cel$g <- limma::normalizeBetweenArrays(dslevin2016cel$g,
                                                  method = "quantile")
dslevin2016cel$g <- log1p(dslevin2016cel$g)
```

4 outliers in the *C. elegans* data are then removed from the analysis (see [*About the outliers*](about-the-outliers)).

``` r
# outlier samples (see below)
rem <- c(sample_0001 = 53L, sample_0002 = 54L,
         sample_0003 = 58L, sample_0004 = 59L)
```

To stage samples we must build a reference from the *Drosophila* data, which happens to be the `Dme_embryo` reference of the `drosoRef` package. It can thus be directly loaded with `prepare_refdata()`.

```
# reference built from dsgraveley2011
r_grav <- prepare_refdata("Dme_embryo", "drosoRef", 500)
```

We then convert the *C. elegans* gene IDs to their *D. melanogaster* orthologs, and stage the *C. elegans* embryos with the Drosophila reference. In this case, many-to-one orthologs are averaged and one-to-many are assigned to the first ID match.

```
dslevin2016cel$g_dmel <- format_ids(dslevin2016cel$g, glist,
                                    from = "wb_id", to = "fb_id")
#> Kept 5714 out of 20168 - aggregated into 4214

ae_cel_on_dmel <- ae(dslevin2016cel$g_dmel, r_grav)
#>                 nb.genes
#> refdata            12408
#> samp                4214
#> intersect.genes     3194
#> Bootstrap set size is 1065
#> Performing age estimation...
#> Bootstrapping...
#>  Building gene subsets...
#>  Computing correlations...
#>  Performing age estimation...
#> Computing summary statistics...
```

The gaps we see in the staging results are likely at timings where there are incompatible expression dynamics between the two species.

## [2.3] Adjusting the reference

By re-building the *Drosophila* reference on the first 2 components which are broad or monotonic, we can keep only these expression dynamics of development in the reference. This can improve staging because it filters out the incompatible dynamics.

```
# build same model, but restrict interpolation to 2 components instead of 8
m_grav2 <- ge_im(X = dsgraveley2011$g, p = dsgraveley2011$p,
                 formula = "X ~ s(age, bs = 'cr')", nc = 2)

# make adjusted reference object
r_grav2 <- make_ref(m_grav2, n.inter = 500,
                    t.unit = "h past egg-laying",
                    metadata = list("organism"="D. melanogaster"))
```

We then stage the *C. elegans* embryos on this adjusted reference.

```
ae_cel_on_dmel2 <- ae(dslevin2016cel$g_dmel, r_grav2)
#>                 nb.genes
#> refdata            12300
#> samp                4214
#> intersect.genes     3143
#> Bootstrap set size is 1048
#> Performing age estimation...
#> Bootstrapping...
#>  Building gene subsets...
#>  Computing correlations...
#>  Performing age estimation...
#> Computing summary statistics...
```

## [2.4] About the outliers

We removed 4 outlier samples from the analyses above because they have an erroneous noted age. This is evidenced by their position on principal components, where they appear to be around 2 hours "older" than their specified age.

```
pca_cel <- stats::prcomp(t(dslevin2016cel$g), rank = 10,
                         center = TRUE, scale = FALSE)
```

This is even further confirmed by their estimated age on the *Drosophila* reference, which places them as expected.

## [2.5] Code to generate objects

Required packages and variables:

``` r
data_folder <- "../inst/extdata/"

requireNamespace("wormRef", quietly = T)
requireNamespace("utils", quietly = T)
requireNamespace("GEOquery", quietly = T) # bioconductor
requireNamespace("Biobase", quietly = T)  # bioconductor
requireNamespace("biomaRt", quietly = T)  # bioconductor
```

*Note : set the `data_folder` variable to an existing path on your system where you want to store the objects.*

``` r
raw2tpm <- function(rawcounts, genelengths){
  if(nrow(rawcounts) != length(genelengths))
    stop("genelengths must match nrow(rawcounts).")
  x <- rawcounts/genelengths
  return(t( t(x) * 1e6 / colSums(x) ))
}

fpkm2tpm <- function(fpkm){
  return(exp(log(fpkm) - log(colSums(fpkm)) + log(1e6)))
}
```

\

To build `glist`, ortholog genes between *C. elegans* and *D. melanogaster* from [Li et al. (2014)] Supplementary Table 1:

``` r
tmp_file <- paste0(data_folder, "dmel_cel_orth.zip")
tmp_fold <- paste0(data_folder, "dmel_cel_orth/")
f_url <- paste0("https://genome.cshlp.org/content/suppl/2014/05/15/",
                "gr.170100.113.DC1/Supplemental_Files.zip")

utils::download.file(url = f_url, destfile = tmp_file)
utils::unzip(tmp_file, exdir = tmp_fold)

glist <- read.table(
  paste0(tmp_fold,
         "Supplementary\ files/TableS1\ fly-worm\ ortholog\ pairs.txt"),
  skip = 1, h=T, sep = "\t", as.is = T, quote = "\""
  )
colnames(glist) <- c("fb_id", "dmel_name", "cel_id", "cel_name")
glist$wb_id <- wormRef::Cel_genes[
  match(glist$cel_id, wormRef::Cel_genes$sequence_name), "wb_id"]

save(glist, file = paste0(data_folder, "sc2_glist.RData"), compress = "xz")

# cleanup
file.remove(tmp_file)
unlink(tmp_fold, recursive = T)
rm(tmp_file, tmp_fold, f_url)
```

To build `dsgraveley2011`, first get drosophila genes from ensembl:

``` r
mart <- biomaRt::useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
droso_genes <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                             "ensembl_transcript_id",
                                             "external_gene_name",
                                             "transcript_end",
                                             "transcript_start"),
                              mart = mart)
droso_genes$transcript_length <-
  droso_genes$transcript_end - droso_genes$transcript_start
droso_genes <- droso_genes[,c(1:3,6)]
colnames(droso_genes)[1:3] <- c("fb_id", "transcript_id", "gene_name")

rm(mart)
```

Then, download *D. melanogaster* data from [Graveley et al. (2011)] :

``` r
g_url_dsgraveley2011 <- paste0(
  "ftp://ftp.fruitfly.org/pub/download/modencode_expression_scores/",
  "Celniker_Drosophila_Annotation_20120616_1428_allsamps_",
  "MEAN_gene_expression.csv.gz")
g_file_dsgraveley2011 <- paste0(data_folder, "dsgraveley2011.csv.gz")
utils::download.file(g_url_dsgraveley2011, destfile = g_file_dsgraveley2011)

X_dsgraveley2011 <- read.table(gzfile(g_file_dsgraveley2011),
                               sep = ',', row.names = 1, h = T)

# convert gene ids to FBgn
X_dsgraveley2011 <- RAPToR::format_ids(X_dsgraveley2011, droso_genes,
                                       from = "gene_name", to = "fb_id")

# select embryo time series samples
X_dsgraveley2011 <- X_dsgraveley2011[,1:12]

P_dsgraveley2011 <- data.frame(
  sname = colnames(X_dsgraveley2011),
  age = as.numeric(gsub("em(\\d+)\\.\\d+hr", "\\1",
                        colnames(X_dsgraveley2011))),
  stringsAsFactors = FALSE)

dsgraveley2011 <- list(g = X_dsgraveley2011, p = P_dsgraveley2011)

save(dsgraveley2011,
     file = paste0(data_folder, "dsgraveley2011.RData"), compress = "xz")

# cleanup
file.remove(g_file_dsgraveley2011)
rm(g_url_dsgraveley2011, g_file_dsgraveley2011,
   X_dsgraveley2011, P_dsgraveley2011)
```

To build `dslevin2016cel`, *C. elegans* single-embryo data from [Levin et al. (2016)]

``` r
geo_dslevin2016cel <- "GSE60755"

g_url_dslevin2016cel <- GEOquery::getGEOSuppFiles(geo_dslevin2016cel,
                                                  makeDirectory = FALSE,
                                                  fetch_files = FALSE)
g_file_dslevin2016cel <- paste0(data_folder, "dslevin2016cel.txt.gz")
utils::download.file(url = as.character(g_url_dslevin2016cel$url[1]),
                     destfile = g_file_dslevin2016cel)

X_dslevin2016cel <- read.table(gzfile(g_file_dslevin2016cel),
                               h = T, sep = '\t', as.is = T, row.names = 1,
                               comment.char = "")

# filter poor quality samples
cm_dslevin2016cel <- RAPToR::cor.gene_expr(X_dslevin2016cel, X_dslevin2016cel)
f_dslevin2016cel <- which(0.67 > apply(cm_dslevin2016cel, 1,
                                       quantile, probs = .99))
X_dslevin2016cel <- X_dslevin2016cel[, -f_dslevin2016cel]

# convert to tpm & FBgn
X_dslevin2016cel <- X_dslevin2016cel[
  rownames(X_dslevin2016cel)%in%wormRef::Cel_genes$sequence_name,]
X_dslevin2016cel <- raw2tpm(
  rawcounts = X_dslevin2016cel,
  genelengths = wormRef::Cel_genes$transcript_length[
    match(rownames(X_dslevin2016cel), wormRef::Cel_genes$sequence_name)]
  )
X_dslevin2016cel <- RAPToR::format_ids(X_dslevin2016cel, wormRef::Cel_genes,
                                       from = "sequence_name", to = "wb_id")

# pheno data
P_dslevin2016cel <- Biobase::pData(GEOquery::getGEO(geo_dslevin2016cel,
                                                    getGPL = F)[[1]])

# filter relevant fields/samples
P_dslevin2016cel <- P_dslevin2016cel[
  , c("title", "geo_accession", "time point (minutes after 4-cell):ch1")]
colnames(P_dslevin2016cel)[3] <- "time"
P_dslevin2016cel$title <- as.character(P_dslevin2016cel$title)

P_dslevin2016cel <- P_dslevin2016cel[
  P_dslevin2016cel$title %in% colnames(X_dslevin2016cel),]

# formatting
P_dslevin2016cel$age <- as.numeric(P_dslevin2016cel$time) / 60
P_dslevin2016cel <- P_dslevin2016cel[order(P_dslevin2016cel$age),]
X_dslevin2016cel <- X_dslevin2016cel[, P_dslevin2016cel$title]

P_dslevin2016cel$title <- gsub('Metazome_CE_timecourse_', '',
                               P_dslevin2016cel$title)
colnames(X_dslevin2016cel) <- P_dslevin2016cel$title

# remove extra outlier
X_dslevin2016cel <- X_dslevin2016cel[,-127]
P_dslevin2016cel <- P_dslevin2016cel[-127,]

dslevin2016cel <- list(g = X_dslevin2016cel, p = P_dslevin2016cel)

save(dslevin2016cel,
     file = paste0(data_folder, "dslevin2016cel.RData"), compress = "xz")

# cleanup
file.remove(g_file_dslevin2016cel)
rm(geo_dslevin2016cel, g_url_dslevin2016cel, g_file_dslevin2016cel,
   f_dslevin2016cel, cm_dslevin2016cel, X_dslevin2016cel, P_dslevin2016cel)
```

\
\

# [3] Capturing tissue-specific development

In some species, tissues can develop at different rates, and these rates can vary between individuals. In *C. elegans*, [Perez et al. (2017)] have shown such developmental heterochrony between soma and germline tissues.

By restricting the genes RAPToR uses for staging to those associated with (or expressed in) specific tissues, it is possible to estimate development specific to these tissues.

## [3.1] Data and strategy

[Rockman, Skrovanek, and Kruglyak (2010)] published microarray profiling of 208 recombinant inbred lines of *C. elegans* N2 and Hawaii (CB4856) strains (Accession: [GSE23857](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23857), `dsrockman2010`). These 208 samples were described as "developmentally synchronized" in the original article. However, [Francesconi and Lehner (2014)] later showed there is significant developmental spread between samples, spanning around 20 hours of 20°C late-larval development, essentially making this dataset a very high-resolution time course.

\
We will start by staging samples with all available genes to get their "***global age***".

Then, we stage samples a second time, but restricting the input to a `germline` set of 2554 genes by combining the `germline_intrinsic`, `germline_oogenesis_enriched`, and `germline_sperm_enriched` categories defined in [Perez et al. (2017)] (based on previous work by [(**reinke2004genome?**)]). This will give us a ***germline age*** for the samples.

Staging the samples a third time, with a `soma` set of 2718 genes from the `osc` (molting) gene category defined in [Hendriks et al. (2014)], will give us the ***soma age***.

\

Finally, we evaluate the results of staging using principal (PCA) or independent (ICA) components, that summarize the expression dynamics of the data. For example, we can expect that components corresponding to particular tissues will have noticeably cleaner dynamics along their respective age estimates (*e.g.* oscillatory molting dynamics would be less noisy with soma age than germline age).

Code to generate `dsrockman2010` (RIL expression data) and `gsubset` (germline/soma gene sets) can be found [at the end of this section](#code-to-generate-objects-2)

``` r
library(RAPToR)
library(wormRef)

library(ica)
library(stats)
library(limma)

# for plotting
library(vioplot)
```

## [3.2] Estimating sample global age

We start by applying a quantile-normalization and \\(log(X+1)\\) transformation to the data.

``` r
dsrockman2010$g <- limma::normalizeBetweenArrays(dsrockman2010$g,
                                                 method = "quantile")
dsrockman2010$g <- log1p(dsrockman2010$g)
```

Then, we select an appropriate *C. elegans* reference (larval to young-adult), and stage the samples.

```
r_ya <- prepare_refdata("Cel_YA_1", "wormRef", n.inter = 400)

ae_dsrockman2010 <- ae(dsrockman2010$g, r_ya)
#>                 nb.genes
#> refdata            16327
#> samp               14440
#> intersect.genes    13288
#> Bootstrap set size is 4429
#> Performing age estimation...
#> Bootstrapping...
#>  Building gene subsets...
#>  Computing correlations...
#>  Performing age estimation...
#> Computing summary statistics...
```

Our estimates confirm the wide developmental spread and match the previously inferred ages by [Francesconi and Lehner (2014)].

## [3.3] Caracterization of expression dynamics

Expression dynamics in the data can be observed through ICA components, and enrichment in component gene loadings (or contributions) allows us to associate them to specific processes.

We find out how many independent components (IC) to extract from how many principal components are needed to explain \\(90\\%\\) of the variance in the data.

```
pca_rock <- summary(stats::prcomp(t(dsrockman2010$g),
                                  rank = 30, center = TRUE, scale = FALSE))
nc <- sum(pca_rock$importance["Cumulative Proportion",] < .90) + 1
nc
#> [1] 48
```

Then, we perform an ICA extracting the 48 components.

```
ica_rock <- ica::icafast(t(scale(t(dsrockman2010$g), center = T, scale = F)),
                         nc = nc)
```

We can plot the first few independent components along global estimated age of the samples :

It seems that past IC7, components don't have a definite link to developmental processes, which is expected given the (intended) genetic background heterogeneity in the samples. We can have a closer look at components clearly capturing development, and the enrichment of their loadings for the gene sets of interest.

```
dev_comps <- 1:7
```

From the gene loadings above, we can establish that

-   component IC2 is linked to oogenesis
-   component IC3 is clearly associated with spermatogenesis
-   components IC1, IC4, IC5, and IC7 have good contribution from the soma geneset, and show oscillatory dynamics.

## [3.4] Estimating tissue-specific age

Now, we stage the samples using only germline or soma gene subsets.

```
ae_soma <- ae(
  dsrockman2010$g[rownames(dsrockman2010$g) %in% gsubset$soma,],
  r_ya)
#>                 nb.genes
#> refdata            16327
#> samp                2005
#> intersect.genes     1815
#> Bootstrap set size is 605
#> Performing age estimation...
#> Bootstrapping...
#>  Building gene subsets...
#>  Computing correlations...
#>  Performing age estimation...
#> Computing summary statistics...

ae_germline <- ae(
  dsrockman2010$g[rownames(dsrockman2010$g) %in% gsubset$germline,],
  r_ya)
#>                 nb.genes
#> refdata            16327
#> samp                1856
#> intersect.genes     1666
#> Bootstrap set size is 555
#> Performing age estimation...
#> Bootstrapping...
#>  Building gene subsets...
#>  Computing correlations...
#>  Performing age estimation...
#> Computing summary statistics...
```

We notice the soma estimates of two samples (marked with \\(\\small{\\triangle}\\)) are quite off from their global or germline age. This caused by the combination of both the fact we used a small gene set for inferring age, and that very similar expression profiles can occur at different times in oscillatory profiles (which is the case for molting genes).

If we look at the global, soma, or germline correlation profiles of one of these samples, we can see 2 peaks of similar correlation for the soma, which is not the case for germline or global correlation.

This is the perfect situation to use a prior for staging. We can input the global age as a prior so that the correct peak of the soma correlation profile is favored in the erroneous samples.

```
ae_soma_prior <- ae(
  dsrockman2010$g[rownames(dsrockman2010$g)%in%gsubset$soma,],
  r_ya,
  prior = ae_dsrockman2010$age.estimates[,1], # gaussian prior values (mean)
  prior.params = 10                           # gaussian prior sd
  )
#>                 nb.genes
#> refdata            16327
#> samp                2005
#> intersect.genes     1815
#> Bootstrap set size is 605
#> Performing age estimation...
#> Bootstrapping...
#>  Building gene subsets...
#>  Computing correlations...
#>  Performing age estimation...
#> Computing summary statistics...
```

Indeed, the estimate is now shifted to the first (correct) peak. Note that the correlation profile itself (central black line) is unchanged (dotted lines can be different as they are derived from bootstrap estimates with random gene sets).

At the same time, all the other estimates are essentially unchanged.

```
# 80 & 141 are the offset samples
mean((ae_soma$age.estimates[,1] - ae_soma_prior$age.estimates[,1])[-c(80,141)])
#> [1] -0.014664
```

Now, we can observe the expression dynamics in the data along tissue-specific estimates.

As expected, components IC1, IC4, IC5, and IC7 that we previously associated with soma are much cleaner along the soma age, but very noisy along germline age.

At the same time, germline-associated components IC2 and IC3 appear quite noisy along soma age, but are very clean along germline age estimates.

This effect is the result of soma-germline heterochrony between the samples.

\
\

Note: In this vignette, we use `Cel_YA_1` to stage the samples, which is different from the RAPToR article ([Bulteau and Francesconi (2022)]) which uses `Cel_larv_YA`. The article shows a combination of soma-germline heterochrony *between* the reference and the samples, as well as among the samples. We shortened the analysis for the vignette to only show heterochrony among RILs (*e.g.* the ICA does not combine reference and samples, so components are different; however, enrichment of soma or germline genes to specific components/dynamics is still clear).

## [3.5] Code to generate objects

``` r
data_folder <- "../inst/extdata/"

requireNamespace("wormRef", quietly = T)
requireNamespace("utils", quietly = T)
requireNamespace("readxl", quietly = T)

requireNamespace("GEOquery", quietly = T) # bioconductor
requireNamespace("Biobase", quietly = T)  # bioconductor
requireNamespace("limma", quietly = T)    # bioconductor
```

*Note : set the `data_folder` variable to an existing path on your system where you want to store the objects.*

To generate `dsrockman2010`:

``` r
geo_dsrockman2010 <- "GSE23857"
geo_dsrockman2010 <- GEOquery::getGEO(geo_dsrockman2010, GSEMatrix = F)

# get pdata
P_dsrockman2010 <- do.call(
  rbind,
  lapply(GEOquery::GSMList(geo_dsrockman2010), function(go){
    unlist(GEOquery::Meta(go)[
      c("geo_accession", "characteristics_ch1", "characteristics_ch2",
        "label_ch1", "label_ch2")]
    )
  })
  )

P_dsrockman2010 <- as.data.frame(P_dsrockman2010, stringsAsFactors = F)

# get RG data from microarray
RG_dsrockman2010 <- list(
  R = do.call(cbind, lapply(GEOquery::GSMList(geo_dsrockman2010), function(go){
    GEOquery::Table(go)[, "R_MEAN_SIGNAL"]
  })),
  G = do.call(cbind, lapply(GEOquery::GSMList(geo_dsrockman2010), function(go){
    GEOquery::Table(go)[, "G_MEAN_SIGNAL"]
  }))
)

# normalize within channels
RG_dsrockman2010 <- limma::normalizeWithinArrays(RG_dsrockman2010,
                                                 method = "loess")
RG_dsrockman2010 <- limma::RG.MA(RG_dsrockman2010) # convert back from MA to RG

# get only the RIALs (not the mixed stage controls)
X_dsrockman2010 <- RG_dsrockman2010$R
X_dsrockman2010[, P_dsrockman2010$label_ch1 == "Cy5"] <-
  RG_dsrockman2010$G[, P_dsrockman2010$label_ch1 == "Cy5"]

# format probe/gene ids
gpl <- GEOquery::getGEO(GEOquery::Meta(
  GEOquery::GSMList(geo_dsrockman2010)[[1]]
  )$platform_id)
gpl <- GEOquery::Table(gpl)
gpl <- gpl[as.character(gpl$ID) %in% as.character(GEOquery::Table(
  GEOquery::GSMList(geo_dsrockman2010)[[1]]
  )$ID_REF), ]

sel <- gpl$ORF%in%wormRef::Cel_genes$sequence_name
gpl <- gpl[sel,]
X_dsrockman2010 <- X_dsrockman2010[sel,]

# filter bad quality samples
cm_dsrockman2010 <- cor(log1p(X_dsrockman2010), method = 'spearman')
f_dsrockman2010 <- which(0.95 > apply(cm_dsrockman2010, 1,
                                      quantile, probs = .95))

X_dsrockman2010 <- X_dsrockman2010[,-f_dsrockman2010]
P_dsrockman2010 <- P_dsrockman2010[-f_dsrockman2010,]

# format ids
rownames(X_dsrockman2010) <- gpl$ORF
X_dsrockman2010 <- RAPToR::format_ids(X_dsrockman2010, wormRef::Cel_genes,
                                      from = "sequence_name", to = "wb_id")

dsrockman2010 <- list(g = X_dsrockman2010, p = P_dsrockman2010)
save(dsrockman2010,
     file = paste0(data_folder, "dsrockman2010.RData"), compress = "xz")

rm(geo_dsrockman2010, gpl, sel,
   RG_dsrockman2010, X_dsrockman2010, P_dsrockman2010,
   cm_dsrockman2010, f_dsrockman2010)
```

To download the soma and germline gene sets (`gsubset`) :

``` r
library(readxl)
# germline sets from Perez et al. (2017)
germline_url <- paste0("https://static-content.springer.com/esm/",
                       "art%3A10.1038%2Fnature25012/MediaObjects/",
                       "41586_2017_BFnature25012_MOESM3_ESM.xlsx")
germline_file <- paste0(data_folder, "germline_gset.xlsx")
utils::download.file(url = germline_url, destfile = germline_file)

germline_set <- read_xlsx(germline_file, sheet = 3, na = "NA")[,c(1, 44:46)]
germline_set[is.na(germline_set)] <- FALSE
germline_set <- cbind(wb_id = germline_set[,1],
                      germline = apply(germline_set[, 2:4], 1,
                                       function(r) any(r)),
                      germline_set[, 2:4])

germline <- germline_set[germline_set$germline,1]
germline_intrinsic <- germline_set[germline_set$germline_intrinsic,1]
germline_oogenesis <- germline_set[germline_set$germline_oogenesis_enriched,1]
germline_sperm <- germline_set[germline_set$germline_sperm_enriched,1]

# soma set from Hendriks et al. (2014)
soma_url <- paste0("https://ars.els-cdn.com/content/image/",
                   "1-s2.0-S1097276513009039-mmc2.xlsx")
soma_file <- paste0(data_folder, "soma_gset.xlsx")
utils::download.file(url = soma_url, destfile = soma_file)

soma_set <- read_xlsx(soma_file, skip = 3, na = "NA")[,c(1, 4)]
soma_set$class <- factor(soma_set$class)

soma_set$soma <- soma_set$class == "osc"
soma_set <- soma_set[soma_set$soma, 1]

# save gene sets
gsubset <- list(germline = germline, soma = soma_set$`Gene WB ID`,
                germline_intrinsic = germline_intrinsic,
                germline_oogenesis = germline_oogenesis,
                germline_sperm = germline_sperm)

save(gsubset, file = paste0(data_folder, "sc3_gsubset.RData"), compress = "xz")

# cleanup
file.remove(germline_file)
file.remove(soma_file)
rm(germline_url, germline_file, germline_set,
   soma_url, soma_file, soma_set)
```

[Francesconi and Lehner (2014)] sample timings (`francesconi_time`):

``` r
# Copied from supp data of Francesconi & Lehner (2014)
francesconi_time <- data.frame(
  time = c(
  4.862660944, 4.957081545, 5.051502146, 5.145922747, 5.240343348, 5.334763948,
  5.429184549, 5.52360515, 5.618025751, 5.712446352, 5.806866953, 5.901287554,
  5.995708155, 6.090128755, 6.184549356, 6.278969957, 6.373390558, 6.467811159,
  6.56223176, 6.656652361, 6.751072961, 6.845493562, 6.939914163, 7.034334764,
  7.128755365, 7.223175966, 7.317596567, 7.412017167, 7.506437768, 7.600858369,
  7.69527897, 7.789699571, 7.884120172, 7.978540773, 8.072961373, 8.167381974,
  8.261802575, 8.356223176, 8.450643777, 8.545064378, 8.639484979, 8.733905579,
  8.82832618, 8.82832618, 8.82832618, 8.875536481, 8.875536481, 8.875536481,
  8.875536481, 8.875536481, 8.875536481, 8.875536481, 8.875536481, 8.875536481,
  8.969957082, 9.017167382, 9.017167382, 9.064377682, 9.064377682, 9.111587983,
  9.206008584, 9.206008584, 9.206008584, 9.300429185, 9.394849785, 9.489270386,
  9.489270386, 9.489270386, 9.489270386, 9.489270386, 9.583690987, 9.583690987,
  9.583690987, 9.583690987, 9.583690987, 9.583690987, 9.583690987, 9.583690987,
  9.630901288, 9.725321888, 9.819742489, 9.819742489, 9.819742489, 9.819742489,
  9.819742489, 9.819742489, 9.819742489, 9.91416309, 10.00858369, 10.05579399,
  10.05579399, 10.05579399, 10.05579399, 10.10300429, 10.19742489, 10.19742489,
  10.29184549, 10.29184549, 10.29184549, 10.38626609, 10.38626609, 10.38626609,
  10.43347639, 10.43347639, 10.43347639, 10.43347639, 10.43347639, 10.43347639,
  10.43347639, 10.43347639, 10.43347639, 10.43347639, 10.43347639, 10.43347639,
  10.527897, 10.6223176, 10.6223176, 10.6223176, 10.6223176, 10.6223176,
  10.6223176, 10.6223176, 10.6695279, 10.6695279, 10.6695279, 10.7639485,
  10.7639485, 10.7639485, 10.8583691, 10.8583691, 10.9527897, 10.9527897,
  10.9527897, 11.0472103, 11.1416309, 11.2360515, 11.2360515, 11.3304721,
  11.3304721, 11.3776824, 11.472103, 11.56652361, 11.66094421, 11.75536481,
  11.84978541, 11.94420601, 12.03862661, 12.13304721, 12.22746781, 12.32188841,
  12.41630901, 12.51072961, 12.60515021, 12.69957082, 12.79399142, 12.88841202,
  12.98283262, 13.07725322, 13.17167382, 13.26609442, 13.36051502, 13.45493562,
  13.54935622, 13.54935622, 13.54935622, 13.54935622, 13.54935622, 13.59656652,
  13.69098712, 13.78540773, 13.78540773, 13.78540773, 13.87982833, 13.97424893,
  14.06866953, 14.06866953, 14.06866953, 14.16309013, 14.25751073, 14.35193133,
  14.44635193, 14.54077253, 14.63519313, 14.72961373, 14.82403433, 14.82403433,
  14.82403433, 14.91845494, 14.96566524, 15.01287554, 15.10729614, 15.20171674,
  15.29613734, 15.39055794, 15.48497854, 15.57939914, 15.67381974, 15.76824034,
  15.86266094, 15.95708155, 16.05150215, 16.14592275, 16.24034335, 16.33476395,
  16.42918455, 16.52360515),
  geo_accession = c(
  "GSM588291", "GSM588174", "GSM588110", "GSM588097", "GSM588271", "GSM588203",
  "GSM588105", "GSM588200", "GSM588123", "GSM588122", "GSM588115", "GSM588100",
  "GSM588171", "GSM588190", "GSM588229", "GSM588206", "GSM588277", "GSM588129",
  "GSM588175", "GSM588151", "GSM588273", "GSM588216", "GSM588099", "GSM588117",
  "GSM588179", "GSM588164", "GSM588184", "GSM588092", "GSM588285", "GSM588272",
  "GSM588228", "GSM588121", "GSM588170", "GSM588194", "GSM588143", "GSM588149",
  "GSM588156", "GSM588220", "GSM588212", "GSM588089", "GSM588209", "GSM588253",
  "GSM588091", "GSM588113", "GSM588130", "GSM588202", "GSM588191", "GSM588244",
  "GSM588227", "GSM588197", "GSM588233", "GSM588292", "GSM588163", "GSM588196",
  "GSM588224", "GSM588283", "GSM588267", "GSM588257", "GSM588221", "GSM588274",
  "GSM588090", "GSM588114", "GSM588195", "GSM588265", "GSM588182", "GSM588093",
  "GSM588157", "GSM588251", "GSM588177", "GSM588188", "GSM588269", "GSM588145",
  "GSM588205", "GSM588162", "GSM588210", "GSM588166", "GSM588125", "GSM588252",
  "GSM588207", "GSM588173", "GSM588102", "GSM588286", "GSM588107", "GSM588238",
  "GSM588189", "GSM588106", "GSM588295", "GSM588192", "GSM588134", "GSM588183",
  "GSM588103", "GSM588198", "GSM588293", "GSM588218", "GSM588259", "GSM588234",
  "GSM588137", "GSM588152", "GSM588133", "GSM588250", "GSM588168", "GSM588235",
  "GSM588148", "GSM588279", "GSM588140", "GSM588241", "GSM588111", "GSM588231",
  "GSM588128", "GSM588131", "GSM588101", "GSM588088", "GSM588281", "GSM588159",
  "GSM588249", "GSM588290", "GSM588118", "GSM588154", "GSM588136", "GSM588268",
  "GSM588204", "GSM588160", "GSM588135", "GSM588098", "GSM588294", "GSM588225",
  "GSM588181", "GSM588248", "GSM588096", "GSM588217", "GSM588147", "GSM588176",
  "GSM588116", "GSM588146", "GSM588127", "GSM588104", "GSM588108", "GSM588262",
  "GSM588223", "GSM588161", "GSM588237", "GSM588172", "GSM588284", "GSM588256",
  "GSM588165", "GSM588211", "GSM588242", "GSM588169", "GSM588240", "GSM588264",
  "GSM588219", "GSM588287", "GSM588124", "GSM588178", "GSM588167", "GSM588258",
  "GSM588232", "GSM588141", "GSM588112", "GSM588208", "GSM588215", "GSM588132",
  "GSM588278", "GSM588275", "GSM588155", "GSM588153", "GSM588109", "GSM588185",
  "GSM588138", "GSM588094", "GSM588226", "GSM588236", "GSM588266", "GSM588119",
  "GSM588222", "GSM588246", "GSM588150", "GSM588261", "GSM588201", "GSM588247",
  "GSM588186", "GSM588280", "GSM588255", "GSM588288", "GSM588260", "GSM588139",
  "GSM588245", "GSM588270", "GSM588276", "GSM588199", "GSM588254", "GSM588120",
  "GSM588144", "GSM588243", "GSM588214", "GSM588180", "GSM588126", "GSM588282",
  "GSM588187", "GSM588158", "GSM588193", "GSM588213", "GSM588230", "GSM588263",
  "GSM588142", "GSM588095"),
  stringsAsFactors = F)
rownames(francesconi_time) <- francesconi_time$geo_accession
```

\
\
[Back to top](#top)

------------------------------------------------------------------------

# References

Bulteau, Romain, and Mirko Francesconi. 2022. "Real Age Prediction from the Transcriptome with RAPToR." *Nature Methods*, 1--7.

Francesconi, Mirko, and Ben Lehner. 2014. "The Effects of Genetic Variation on Gene Expression Dynamics During Development." *Nature* 505 (7482): 208.

Graveley, Brenton R, Angela N Brooks, Joseph W Carlson, Michael O Duff, Jane M Landolin, Li Yang, Carlo G Artieri, et al. 2011. "The Developmental Transcriptome of Drosophila Melanogaster." *Nature* 471 (7339): 473.

Hendriks, Gert-Jan, Dimos Gaidatzis, Florian Aeschimann, and Helge Großhans. 2014. "Extensive Oscillatory Gene Expression During c. Elegans Larval Development." *Molecular Cell* 53 (3): 380--92.

Lehrbach, Nicolas J, Cecilia Castro, Kenneth J Murfitt, Cei Abreu-Goodger, Julian L Griffin, and Eric A Miska. 2012. "Post-Developmental microRNA Expression Is Required for Normal Physiology, and Regulates Aging in Parallel to Insulin/IGF-1 Signaling in c. Elegans." *Rna* 18 (12): 2220--35.

Levin, Michal, Leon Anavy, Alison G Cole, Eitan Winter, Natalia Mostov, Sally Khair, Naftalie Senderovich, et al. 2016. "The Mid-Developmental Transition and the Evolution of Animal Body Plans." *Nature* 531 (7596): 637.

Li, Jingyi Jessica, Haiyan Huang, Peter J Bickel, and Steven E Brenner. 2014. "Comparison of d. Melanogaster and c. Elegans Developmental Stages, Tissues, and Cells by modENCODE RNA-Seq Data." *Genome Research* 24 (7): 1086--1101.

Perez, Marcos Francisco, Mirko Francesconi, Cristina Hidalgo-Carcedo, and Ben Lehner. 2017. "Maternal Age Generates Phenotypic Variation in Caenorhabditis Elegans." *Nature* 552 (7683): 106.

Rockman, Matthew V, Sonja S Skrovanek, and Leonid Kruglyak. 2010. "Selection at Linked Sites Shapes Heritable Phenotypic Variation in c. Elegans." *Science* 330 (6002): 372--76.

Snoek, L Basten, Mark G Sterken, Rita JM Volkers, Mirre Klatter, Kobus J Bosman, Roel PJ Bevers, Joost AG Riksen, Geert Smant, Andrew R Cossins, and Jan E Kammenga. 2014. "A Rapid and Massive Gene Expression Shift Marking Adolescent Transition in c. Elegans." *Scientific Reports* 4 (1): 1--5.

------------------------------------------------------------------------

# SessionInfo

```
sessionInfo()
```

    #> R version 4.1.2 (2021-11-01)
    #> Platform: x86_64-pc-linux-gnu (64-bit)
    #> Running under: Ubuntu 22.04.1 LTS
    #>
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so
    #>
    #> locale:
    #>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
    #>  [3] LC_TIME=C                  LC_COLLATE=en_US.UTF-8
    #>  [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=en_US.UTF-8
    #>  [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C
    #>  [9] LC_ADDRESS=C               LC_TELEPHONE=C
    #> [11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C
    #>
    #> attached base packages:
    #> [1] splines   parallel  stats     graphics  grDevices utils     datasets
    #> [8] methods   base
    #>
    #> other attached packages:
    #>  [1] vioplot_0.4.0      zoo_1.8-11         sm_2.2-5.7.1       ica_1.0-3
    #>  [5] drosoRef_0.2.0     limma_3.50.3       wormRef_0.5.0      RAPToR_1.2.0
    #>  [9] ggpubr_0.6.0       ggplot2_3.4.0.9000 BiocStyle_2.22.0
    #>
    #> loaded via a namespace (and not attached):
    #>  [1] Rcpp_1.0.10         lattice_0.20-45     tidyr_1.3.0
    #>  [4] digest_0.6.31       utf8_1.2.3          R6_2.5.1
    #>  [7] backports_1.4.1     evaluate_0.20       highr_0.10
    #> [10] pillar_1.8.1        Rdpack_2.4          rlang_1.0.6
    #> [13] rstudioapi_0.14     data.table_1.14.6   car_3.1-1
    #> [16] jquerylib_0.1.4     magick_2.7.3        Matrix_1.5-3
    #> [19] rmarkdown_2.20.1    labeling_0.4.2      stringr_1.5.0
    #> [22] tinytex_0.44        munsell_0.5.0       broom_1.0.3
    #> [25] compiler_4.1.2      xfun_0.37           pkgconfig_2.0.3
    #> [28] mgcv_1.8-41         htmltools_0.5.4     tidyselect_1.2.0
    #> [31] tibble_3.1.8        bookdown_0.32       codetools_0.2-18
    #> [34] fansi_1.0.4         dplyr_1.1.0         withr_2.5.0
    #> [37] rbibutils_2.2.13    grid_4.1.2          nlme_3.1-161
    #> [40] jsonlite_1.8.4      gtable_0.3.1        lifecycle_1.0.3
    #> [43] magrittr_2.0.3      scales_1.2.1        cli_3.6.0
    #> [46] stringi_1.7.12      cachem_1.0.6        carData_3.0-5
    #> [49] farver_2.1.1        ggsignif_0.6.4      pryr_0.1.6
    #> [52] bslib_0.4.2         generics_0.1.3      vctrs_0.5.2
    #> [55] tools_4.1.2         glue_1.6.2          beeswarm_0.4.0
    #> [58] purrr_1.0.1         abind_1.4-5         fastmap_1.1.0
    #> [61] yaml_2.3.7          colorspace_2.1-0    BiocManager_1.30.19
    #> [64] rstatix_0.7.2       knitr_1.42          sass_0.4.5

---

# RAPToR-DEcorrection

Code []

-   [Show All Code](#)
-   [Hide All Code](#)

# Correcting for age in differential expression analysis

### RAPToR 1.2.0

Romain Bulteau

#### March 2023

# Contents

-   [[1] Introduction](#introduction)
-   [[2] Quickstart](#quickstart)
    -   [[2.1] Quantifying age-driven expression changes with `ref_compare()`](#quantifying-age-driven-expression-changes-with-ref_compare)
    -   [[2.2] Correcting for age effects in differential expression analysis](#correcting-for-age-effects-in-differential-expression-analysis)
-   [[3] How it works](#how-it-works)
    -   [[3.1] Quantifying age-driven expression changes](#quantifying-age-driven-expression-changes)
    -   [[3.2] Correcting for age effects by integrating reference data in DE analysis](#correcting-for-age-effects-by-integrating-reference-data-in-de-analysis)
-   [[4] Demonstration of DE correction in a controlled case](#demonstration-of-de-correction-in-a-controlled-case)
    -   [[4.1] Data and strategy](#data-and-strategy)
    -   [[4.2] Defining age-matched and shifted sample sets](#defining-age-matched-and-shifted-sample-sets)
    -   [[4.3] Quantifying age-driven expression changes in the shifted sets](#quantifying-age-driven-expression-changes-in-the-shifted-sets)
    -   [[4.4] Effect of developmental shifts on DE analysis performance](#effect-of-developmental-shifts-on-de-analysis-performance)
        -   [[4.4.1] Gold-standard DE](#gold-standard-de)
        -   [[4.4.2] DE with developmental shifts](#de-with-developmental-shifts)
    -   [[4.5] DE analysis integrating reference data to correct age effect](#de-analysis-integrating-reference-data-to-correct-age-effect)
        -   [[4.5.1] Selecting reference windows](#selecting-reference-windows)
        -   [[4.5.2] DE models with reference data](#de-models-with-reference-data)
-   [[5] Functions and code](#functions-and-code)
    -   [[5.1] Code to generate `dsmiki2017`](#code-to-generate-dsmiki2017)
    -   [[5.2] Code to generate quickstart data](#code-to-generate-quickstart-data)
    -   [[5.3] Plotting functions](#plotting-functions)
    -   [[5.4] DE functions](#de-functions)
-   [References](#references)
-   [SessionInfo](#sessioninfo)

# [1] Introduction

The developmental speed of fast-growing organisms such as *C. elegans* is influenced by many different factors -- including experimental conditions themselves -- making it difficult or impossible to collect age-matched individuals between conditions. For example, if a mutant has delayed development but controls and mutants are collected at the same chronological (and therefore different physiological) time, the perturbation of interest will be completely confounded with development (**a**). Because of this, purely developmental gene expression changes can be mislabeled as the effect of a variable of interest (**b**).

When sample groups have age differences but still overlap, the developmental effect can simply be accounted for by including age as a covariate in the Differential Expression (DE) analysis. If there is no age overlap however, it becomes impossible to know whether an observed effect is due to the perturbation or age since both variables are completely confounded.

Using `RAPToR` reference data, we can bridge the gap between non-overlapping sample groups and rescue otherwise impossible DE analyses.

\

# [2] Quickstart

## [2.1] Quantifying age-driven expression changes with `ref_compare()`

```
# load libraries
library(RAPToR)
library(wormRef)
```

Given transcript per million (TPM) expression profiles from 3 control and 3 mutant samples.

``` r
head(qs$tpm)
#>                    ctr1     ctr2     ctr3     mut1     mut2     mut3
#> WBGene00000001 2.862195 3.118323 3.175077 3.340883 3.426089 3.546783
#> WBGene00000002 3.107584 2.579527 2.156750 1.783193 1.843059 1.780094
#> WBGene00000003 2.586109 3.101123 3.152630 3.335000 3.244302 2.833467
#> WBGene00000004 1.930230 1.801899 1.706351 1.938189 2.018981 2.191545
#> WBGene00000005 2.866435 2.483658 2.125052 2.360653 2.545615 2.647073
#> WBGene00000006 1.583126 1.693566 1.873643 2.296500 2.266956 2.017384
print(qs$pdat)
#>   sname strain
#> 1  ctr1    ctr
#> 2  ctr2    ctr
#> 3  ctr3    ctr
#> 4  mut1    mut
#> 5  mut2    mut
#> 6  mut3    mut
```

First, infer sample age with RAPToR.

```
# load reference
r_cel <- prepare_refdata("Cel_larv_YA", "wormRef", 600)

# estimate sample age
ae_qs <- ae(qs$tpm, r_cel)
```

```
plot(ae_qs, group = qs$pdat$strain, g.line=3, lmar=5)
```

Run `ref_compare()` to estimate developmental expression changes between groups.

```
qs_rc <- ref_compare(
    X      = qs$tpm,        # sample data, log(TPM+1)
    ref    = r_cel,         # ref object
    ae_obj = ae_qs,         # ae object
    group  = qs$pdat$strain # factor defining compared groups (wt/mut)
)
```

Plot sample log2 fold-changes (logFCs) vs. expected developmental logFCs inferred from the reference.

```
qs_lfc <- get_logFC(qs_rc)
#> Comparing ctr vs. mut

plot(qs_lfc, cex=.2, asp=1,
     xlab = "ctr vs. mut log2FCs", xlim = c(-10,10),
     ylab = "Developmental log2FCs", ylim = c(-10, 10))
```

Print the correlation between sample and reference logFCs, as well as the average age difference between groups.

```
print(qs_rc)
#> DE comparison with reference
#> ---
#>                 ref.logFC.r ref.logFC.r2 ae.avg.dif
#> ctr (ctrl, n=3)          NA           NA         NA
#> mut (n=3)             0.868        0.753      5.035
#> ---
```

In this example, ctr *vs.* mut logFCs have a correlation of \\(0.868\\) with the expected developmental logFCs computed from the reference data. Using \\(r\^2\\), this corresponds to approximately \\(75.3 \\%\\) of variance explained by development. Such a large developmental effect between the control and mutant groups is expected as they don't have any age overlap.

For an explanation of the analysis, see [*How it works*](#quantifying-age-driven-expression-changes).

A custom `ggplot` function to plot logFCs is also provided in [*Plotting functions*, below](#code-plotting_functions).

```
gg_logFC(qs_lfc, main = "Developmental effect",
         xlab = "ctr vs. mut log2FC", ylab = "Developmental log2FC")
#> Loading required package: viridisLite
```

## [2.2] Correcting for age effects in differential expression analysis

```
# load libraries
library(RAPToR)
library(wormRef)

library(DESeq2)
library(splines)
```

Given raw counts and transcript per million (TPM) expression profiles from 3 control and 3 mutant samples.

``` r
head(qs$count)
#>                ctr1 ctr2 ctr3 mut1 mut2 mut3
#> WBGene00000001 1000 1266 1843 1984 2310 2540
#> WBGene00000002 1387  765  658  386  442  398
#> WBGene00000003  620 1036 1500 1643 1594 1005
#> WBGene00000004  531  441  539  644  754  891
#> WBGene00000005 1361  872  803  947 1236 1339
#> WBGene00000006  423  469  799 1174 1211  886
head(qs$tpm)
#>                    ctr1     ctr2     ctr3     mut1     mut2     mut3
#> WBGene00000001 2.862195 3.118323 3.175077 3.340883 3.426089 3.546783
#> WBGene00000002 3.107584 2.579527 2.156750 1.783193 1.843059 1.780094
#> WBGene00000003 2.586109 3.101123 3.152630 3.335000 3.244302 2.833467
#> WBGene00000004 1.930230 1.801899 1.706351 1.938189 2.018981 2.191545
#> WBGene00000005 2.866435 2.483658 2.125052 2.360653 2.545615 2.647073
#> WBGene00000006 1.583126 1.693566 1.873643 2.296500 2.266956 2.017384
print(qs$pdat)
#>   sname strain
#> 1  ctr1    ctr
#> 2  ctr2    ctr
#> 3  ctr3    ctr
#> 4  mut1    mut
#> 5  mut2    mut
#> 6  mut3    mut
```

First, infer sample age with RAPToR.

```
# load reference
r_cel <- prepare_refdata("Cel_larv_YA", "wormRef", 600)

# estimate sample age
ae_qs <- ae(qs$tpm, r_cel)
```

```
plot(ae_qs, group = qs$pdat$strain, g.line=3, lmar=5)
```

Define a reference window to integrate.

```
r.window <- range(ae_qs$age.estimates[,1]) + c(1,-1) # +1h margin
r.idx <- get_refTP(r_cel, ae_values = r.window) # indices in the ref object
r.idx <- r.idx[1]:r.idx[2]

# extract time and expression values
r.time <- r_cel$time[r.idx]
r.tpm <- r_cel$interpGE[,r.idx]
```

Find the optimal spline degree-of-freedom (df) to model expression dynamics in the selected window.

```
dfs <- 1:8
ssq <- sapply(dfs, function(i){
  # compute SSQ of linear model on expression data
  sum(residuals(lm( t(r.tpm) ~ splines::ns(r.time, df=i) ))^2)
})

plot(dfs, ssq, type='b', lwd=2)
```

```
r.df <- 3 # select df=3
```

Transform reference TPMs into artificial counts using an arbitrary library size.

```
# get gene lengths
g.l <- wormRef::Cel_genes[match(rownames(r.tpm), wormRef::Cel_genes[,"wb_id"]),
                          "transcript_length"]
libsize <- 25e6

# convert tpm to artificial counts
r.count <- t( (t(exp(r.tpm)-1)/1e6) * (libsize/median(g.l))) * g.l
r.count[r.count < 0] <- 0
r.count <- round(r.count)
```

Combine sample and reference count matrices and pdata.

```
# filter low expressed genes (at least >5 in one sample)
qs$count <- qs$count[apply(qs$count, 1, max)>5, ]

# keep overlapping sample and ref genes, and merge in one matrix
comb.count <- do.call(cbind, format_to_ref(qs$count, r.count)[1:2])
#>                 nb.genes
#> refdata            19953
#> samp               17596
#> intersect.genes    17274

# combine pdata
comb.p <- data.frame(
  time = c(ae_qs$age.estimates[,1], r.time),
  strain = c(qs$pdat$strain, factor(rep('ctr', ncol(r.count)),
                                    levels = c('ctr', 'mut'))), # ref is ctr
  batch = factor(rep(c('samp', 'ref'), c(ncol(qs$count), ncol(r.count))),
                 levels = c("samp", "ref"))
)
```

Infer gene dispersions by fitting a model on sample data only.

```
# build DESeq object
dd0 <- DESeqDataSetFromMatrix(countData = comb.count[ ,1:6],
                              colData = comb.p[1:6, ],
                              design = ~strain) # no age yet
# compute gene dispersions
dd0 <- estimateSizeFactors(dd0)
dd0 <- estimateDispersions(dd0, fitType = "local")
d0 <- dispersions(dd0) # store dispersions
d0[is.na(d0)] <- 0 # remove NAs
```

Build a model with combined reference and sample data, injecting previously estimated gene dispersions.

```
dd1 <-  DESeqDataSetFromMatrix(
  countData = comb.count,
  colData = comb.p,
  design = ~ splines::ns(time, df=3)+batch+strain # use df found above
)
dd1 <- estimateSizeFactors(dd1)
# inject dispersions from sample-only model
dispersions(dd1) <- d0

dd1 <- nbinomWaldTest(dd1)
```

Extract DE results.

```
res <- results(dd1)
summary(res)
#>
#> out of 17274 with nonzero total read count
#> adjusted p-value < 0.1
#> LFC > 0 (up)       : 2961, 17%
#> LFC < 0 (down)     : 1894, 11%
#> outliers [1]       : 0, 0%
#> low counts [2]     : 0, 0%
#> (mean count < 0)
#> [1] see 'cooksCutoff' argument of ?results
#> [2] see 'independentFiltering' argument of ?results
```

```
DESeq2::plotMA(dd1)
```

Results are a [standard DESeq2 output](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) and can be manipulated as such.

For an explanation of the analysis, see [*How it works*](#correcting-for-age-effects-by-integrating-reference-data-in-de-analysis).

Custom functions to transform reference data into artificial counts and run DESeq2 as described above are given in [*DE functions*, below](#de-functions) (please note they may require tweaking). The following is equivalent to the code above.

```
dd1 <- run_DESeq2_ref(
  X = qs$count,         # sample counts
  p = qs$pdat,          # sample pdata
  formula = "~ strain", # formula (without age)
  ref = r_cel,          # ref. object
  ae_values = ae_qs$age.estimates[,1], # sample age estimates
  window.extend = 1,    # ref. window extension
  ns.df = 3             # spline df for ref. window
)
```

# [3] How it works

## [3.1] Quantifying age-driven expression changes

Between two conditions 'A' and 'B' where the sample groups have developmental differences, the expression changes (or log-fold changes, logFCs) observed between the groups will be a combination of perturbation and developmental effects (**a**).

We can determine the expression changes expected only from the difference in development using reference expression profiles matching the samples (**b**). Any correlation between observed (sample) and expected (reference) logFCs will then correspond to the developmental effect, with uncorrelated logFCs corresponding to the perturbation effect (**c**).

The `ref_compare()` function inputs:

-   sample data (expression matrix, expects \\(log(TPM+1)\\) as transcripts per million are more comparable across datasets)
-   the reference object used for age estimation,
-   the `ae` object (or age estimate values for the samples),
-   and a group variable (*e.g.* a factor defining WT and mutant).

```
qs_rc <- ref_compare(
    X      = qs$tpm,        # sample data, log(TPM+1)
    ref    = r_cel,         # ref object
    ae_obj = ae_qs,         # ae object
    group  = qs$pdat$strain # factor defining compared groups (wt/mut)
)
```

The resulting sample and reference (log2) logFCs can be extracted with `get_logFC()` and plotted. When there are more than 2 compared conditions, you can specify the factor levels to compare to the control (with `l`) when extracting the logFCs.

A custom `ggplot2` function for plotting logFCs can be found in the [*Plotting functions* section below](#code-plotting_functions).

```
gg_logFC(qs_lfc, main = "Developmental effect",
         xlab = "ctr vs. mut log2FC", ylab = "Developmental log2FC")
#> Loading required package: viridisLite
```

The \\(r\^2\\) between sample and reference logFCs is a rough estimate of the percentage of variance explained by development (`VarDev`, in the bottom right). In this case, over \\(75\\%\\) of expression changes can be explained by development, which is expected given the large age difference between the sample groups.

The average age difference between groups and correlation between sample and reference are also given directly in the output of `ref_compare()`.

```
print(qs_rc)
#> DE comparison with reference
#> ---
#>                 ref.logFC.r ref.logFC.r2 ae.avg.dif
#> ctr (ctrl, n=3)          NA           NA         NA
#> mut (n=3)             0.868        0.753      5.035
#> ---
```

## [3.2] Correcting for age effects by integrating reference data in DE analysis

When there is no developmental overlap between the sample groups to compare, integrating reference data in the DE analysis makes it possible to recover truly DE genes by bridging the gap, as shown in the cartoon below.

Most DE analysis tools input *raw counts* for their particular statistical properties, so the interpolated reference data must be converted from TPM to (artificial) counts assuming an arbitrary fixed library size (\\(25\\times 10\^6\\) counts). Then, because gene dispersions needed for statistical testing cannot be estimated from the noiseless artificial reference, they are inferred from a model without reference data and injected into the model with reference data.

Thus, resulting model coefficients (logFCs) between sample groups are corrected for development by the reference, and their respective statistical tests use dispersion values inferred only from samples. This approach improves upon what is presented in the RAPToR article ([Bulteau and Francesconi (2022)]), generating valid p-values.

**DISCLAIMER:** we provide the functions `find_df()`, `run_DESeq2_ref()` and sub-functions needed for age correction DE analysis only in this vignette and not as part of the package because (1) they go beyond the scope of RAPToR, and (2) they might need to be tweaked by the user to adapt to their needs and experimental designs (unlike `ref_compare()`). This approach should also work with tools other than `DEseq2` (e.g. [`edgeR`](https://bioconductor.org/packages/edgeR)).

# [4] Demonstration of DE correction in a controlled case

## [4.1] Data and strategy

[[Miki, Carl, and Großhans (2017)]](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5695088/) profiled time-series of *C. elegans* wild-type (WT) and *xrn-2* mutants ([GSE97775](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97775)). Code to download this data and generate the `dsmiki2017` object can be found [at the end of this vignette](#code-gen_dsmiki2017)

By selecting specific samples of matching development from both time-series, we define a gold-standard of truly DE genes. Then, by sliding the window of WT samples we can evaluate the effect of increasing developmental shifts between the groups on the DE analysis, as well as measure the advantages of integrating reference data in the model.

``` r
library(RAPToR)
library(wormRef)

library(DESeq2)
library(splines)
library(ROCR)

# for plotting
library(ggplot2)
library(ggpubr)
library(viridis)

col.palette1 <- c('grey20', 'firebrick', 'royalblue', 'forestgreen')
col.palette2 <- viridisLite::viridis(5)
```

## [4.2] Defining age-matched and shifted sample sets

We start by infering sample age.

```
# load reference
r_cel <- prepare_refdata("Cel_larv_YA", "wormRef", 600)

# estimate sample age
ae_miki <- ae(dsmiki2017$g, r_cel)
dsmiki2017$p$ae <- ae_miki$age.estimates[, 1]
```

``` r
plot(ae_miki, main = "Miki et al. (2017) samples",
     group = dsmiki2017$p$strain_long, show.boot_estimates = F,
     color = col.palette1[dsmiki2017$p$strain], lmar = 10, g.line= 1)
```

We then select a set of 3 WT and 3 mutant samples as the age-matched gold-standard, and define sets of 3 WT samples shifted by 1, 2, 3, 5, and 7 time points compared to the mutants.

``` r
GS_wt <- 8:10 # gold standard WT samples
GS_mut <- 19:21 # GS mutant samples

shifts <- -c(1,2,3,5,7) # number of timepoints to shift
subs <- c(list(gold.standard = c(GS_wt, GS_mut)),
          lapply(shifts, function(s) c(GS_wt + s, GS_mut)))
names(subs)[-1] <- paste0('s', abs(shifts))
```

\

## [4.3] Quantifying age-driven expression changes in the shifted sets

With `ref_compare()`, we quantify the effect of development in all the sample sets defined above, and plot the reference vs. sample logFCs.

``` r
rcs <- lapply(subs, function(s){
    RAPToR::ref_compare(
      X = dsmiki2017$g[, s],
      ref = r_cel,
      ae_values = ae_miki$age.estimates[s,1],
      group = dsmiki2017$p$strain[s]
      )
})
```

Unsurprisingly, age-matched reference logFCs (*i.e* development) account for increasing proportions of expression changes with larger age differences between groups.

\

## [4.4] Effect of developmental shifts on DE analysis performance

### [4.4.1] Gold-standard DE

Using the age-matched mutant and WT samples defined as gold-standard, we find a set truly DE genes using a standard [`DESeq2`](https://bioconductor.org/packages/DEseq2) workflow.

We filter genes to keep those overlapping with the reference, and with at least 5 counts in any sample.

``` r
g.filt <- dsmiki2017$g.raw[apply(dsmiki2017$g.raw, 1, function(r) max(r)>5),]
g.filt <- RAPToR::format_to_ref(g.filt, r_cel$interpGE)$samp
#>                 nb.genes
#> refdata            19953
#> samp               18073
#> intersect.genes    17659
```

We then build a DESeq model including age and strain and test for significant differences between mutant and WT groups. We define a `run_DESeq2_age()` function to do this (see [*DE functions* below for code](#code-de_functions)).

```
DE.GS <- run_DESeq2_age(g.filt[, subs$gold.standard],
                        dsmiki2017$p[subs$gold.standard, ])
res.GS <- get_DEres(DE.GS, coefname = "strain_xrn2_vs_wt") # get results table
```

We consider genes as DE with the thresholds below.

```
# Define DE (with p < thr AND |logFC| > thr)
thr.p <- 0.01
thr.logFC <- 1.0

# get true DE genes from gold-standard
truth <- (res.GS$padj < thr.p) & abs(res.GS$log2FoldChange) >= thr.logFC
table(truth)
#> truth
#> FALSE  TRUE
#> 16534  1125
```

With thresholds of \\(p.value \< 0.01\\) and \\(\|logFC\| \> 1\\), we find 1125 DE genes in the absence of developmental effects ("true" DE genes).\
\

### [4.4.2] DE with developmental shifts

We now apply the same model as above on the 5 sample subsets with shifted WT.

``` r
DE.shifts <- lapply(subs[-1], function(s){
  run_DESeq2_age(g.filt[, s], dsmiki2017$p[s, ])
  })
res.shifts <- lapply(DE.shifts, get_DEres)
```

[Precision and recall (P/R)](https://en.wikipedia.org/wiki/Precision_and_recall) metrics allow us to measure the performance the DE analysis (p-value and logFC) in discriminating truly DE genes (*i.e.* genes found in the gold-standard DE above) from non-DE genes. An area under the P/R curve of 1 is a perfect classifier.

``` r
# compute precision-recall curves
rocs_initial_age <- lapply(seq_along(subs[-1]), function(i){
  v <- res.shifts[[i]]$padj
  v[abs(res.shifts[[i]]$log2FoldChange)<=thr.logFC] <- 1
  ROCR::prediction(predictions = -v, labels = truth)
})
```

We see a sharp decline in the DE model performance in detecting true DE genes as the developmental shift increases, particularly once there is no more developmental overlap between the compared groups (starting at WT-3).

If we select DE genes with the same thresholds as the gold-standard, this drop in performance translates to both a decrease of true positives and an increase of false positives, as shown below.

\
\

## [4.5] DE analysis integrating reference data to correct age effect

### [4.5.1] Selecting reference windows

We select a reference window spanning the sample age range with an added margin of 1 hour on either side for each subset, and find the optimal spline degree-of-freedom (df) to model expression dynamics covered by the reference window.

We define a `find_df()` function which computes the residual sum of squares for a range of spline dfs and choose the optimal df values where a plateau is reached (see [*DE functions* below for code](#code-de_functions)).

```
dfs <- 2:8
dfs_ssq <- lapply(subs, function(s){
  find_df(r_cel, # ref object
          w = range(ae_miki$age.estimates[s,1]) + c(-1, 1), # time window
          dfs = dfs) # df range to test
})
# optimal df for each subset based on plot below
df_opti <- c(3,4,4,4,5,5)
```

The selected df increases with the developmental shift, which is expected since the reference window to include gets larger and may thus contain more complex dynamics.

### [4.5.2] DE models with reference data

Having defined the appropriate fit parameters for the reference, we can integrate in the DE model. We define `run_DESeq2_ref()` to do this (see [*DE functions* below for code](#code-de_functions)).

```
# run DESeq2 with reference data inclusion on shifted subsets
DEref.shifts <- lapply(seq_along(subs), function(i){
  s <-subs[[i]]  # select sample subset
  run_DESeq2_ref(X = g.filt[, s],                # sample counts
                 p = dsmiki2017$p[s, ],          # sample pdata
                 ae_values = dsmiki2017$p$ae[s], # age estimates
                 formula = ~strain,   # formula (age is added by the function)
                 ref = r_cel,         # ref object
                 ns.df = df_opti[i],  # df for spline fit
                 window.extend = 1)   # reference window extension
})
res_ref.shifts <- lapply(DEref.shifts, get_DEres) # get results table
```

As done above, let's assess the performance of these models in detecting truly DE genes.

``` r
# compute age-corrected precision-recall curves
rocs_corrected_age <- lapply(seq_along(subs), function(i){
  v <- res_ref.shifts[[i]]$padj
  v[abs(res_ref.shifts[[i]]$log2FoldChange)<=thr.logFC] <- 1
  ROCR::prediction(predictions = -v, labels = truth)
})
```

Compared to the non-corrected models we strongly rescue the performance of the analysis when there is no overlap between the compared sample groups (starting at `WT-3`).

In less shifted subsets or in the gold-standard, we see that the performance of the age-corrected DE analysis to find truly DE genes is comparable to applying no correction. This is expected since developmental dynamics can be properly modeled without reference data as long as there is overlap between the compared groups, nullifying the advantage of adding it to the model.

# [5] Functions and code

## [5.1] Code to generate `dsmiki2017`

[]

Required libraries and variables:

``` r
data_folder <- "../inst/extdata/"

library("GEOquery") # bioconductor package

requireNamespace("wormRef", quietly = T)
requireNamespace("utils", quietly = T)
```

*Note : set the `data_folder` variable to an existing path on your system where you want to store the objects.*

``` r
raw2tpm <- function(rawcounts, genelengths){
  if(nrow(rawcounts) != length(genelengths))
    stop("genelengths must match nrow(rawcounts).")
  x <- rawcounts/genelengths
  return(t( t(x) * 1e6 / colSums(x) ))
}

fpkm2tpm <- function(fpkm){
  return(exp(log(fpkm) - log(colSums(fpkm)) + log(1e6)))
}
```

\

``` r
geo_dsmiki2017 <- "GSE97775"
url_dsmiki2017 <- as.character(
  GEOquery::getGEOSuppFiles(geo_dsmiki2017,
                            makeDirectory = FALSE,
                            fetch_files = FALSE)[9,"url"]
  )
tmpf <- file.path(data_folder, "miki2017.txt.gz")
utils::download.file(url_dsmiki2017, destfile = tmpf)
g <- read.table(gzfile(tmpf), h=T, as.is = T, row.names = 1)
file.remove(tmpf)

# format genes to WB ids
g <- RAPToR::format_ids(g, wormRef::Cel_genes,
                        from = "wb_id", to = "wb_id",
                        aggr.fun = sum)

# store raw counts and convert to log(TPM+1)
g.raw <- g
g <- raw2tpm(g.raw,
             wormRef::Cel_genes$transcript_length[
               match(rownames(g.raw), wormRef::Cel_genes$wb_id)
               ])
g <- log1p(g)

# sample metadata
p <- Biobase::pData(
  GEOquery::getGEO(geo_dsmiki2017, getGPL = F)[[1]]
  )[14:35, c(1,2, 45)]
colnames(p)[3] <- "strain_long"
p[,3] <- factor(p[, 3], levels = unique(p[,3]),
                labels = c("N2 (wild-type)", "HW1660 (xrn-2(xe31))"))
p$strain <- p[,3]
levels(p$strain) <- c("wt", "xrn2")

# get chronological age from sample name
p$age <- as.numeric(gsub(".*_(\\d+)h", "\\1", p$title))
p$title <- colnames(g)

dsmiki2017 <- list(g = g, p = p, g.raw = g.raw)

save(dsmiki2017, file = file.path(data_folder, "dsmiki2017.RData"),
     compress = "xz")
rm(url_dsmiki2017, geo_dsmiki2017, tmpf, g, g.raw, p, raw2tpm)
```

## [5.2] Code to generate quickstart data

Given the `dsmiki2017` data generated above:

``` r
sel <- c(5:7, 19:21) # select 3 WT and 3 mut

qs <- list(
  tpm = dsmiki2017$g[, sel],
  count = dsmiki2017$g.raw[, sel],
  pdat = data.frame(
    sname = paste0(rep(c('ctr', 'mut'), e=3), rep(1:3,2)),
    strain = factor(rep(c('ctr', 'mut'), e=3)),
    stringsAsFactors = F)
)
colnames(qs$tpm) <- qs$pdat$sname
colnames(qs$count) <- qs$pdat$sname

save(qs, file = file.path(data_folder, 'sc4_qsdata.RData'), compress = "xz")
```

\

## [5.3] Plotting functions

[]

Plotting a logFC comparison as binned heatmap (`gg_logFC()`)

``` r
gg_logFC <- function(x, y=NULL, rg = range(xy, na.rm = T)*1.02,
                     l.breaks = log1p(c(0, 1, 5, 10, 50, 100, Inf)),
                     l.labels = c('1', '5', '10', '50', '100', '100+'),
                     nbins = 200,
                     xlab = "log2(FC) in x",
                     ylab = "log2(FC) in y",
                     main = "", get.r =T, add.vd = T,
                     DEgsel = NULL,
                     ...){
  # Make a 2d binned hexplot for showing logFC comparison
  require(ggplot2)
  require(viridisLite)

  if(is.null(y))
    xy <- x
  else
    xy <- cbind(x, y)
  xy <- as.data.frame(xy)

  colnames(xy) <- c("x", "y")

  g <- ggplot(data = xy, mapping = aes(x=x, y=y)) +
    stat_bin_hex(aes(fill = after_stat(cut(log1p(after_stat(count)), breaks = l.breaks,
                                 labels = F, right = T, include.lowest = T))),
             bins=nbins) +
    scale_fill_gradientn(colors = viridisLite::inferno(length(l.breaks)),
                         name = 'count', labels = l.labels) +
    theme_classic() + xlim(rg) + ylim(rg) +
    xlab(xlab) + ylab(ylab) + ggtitle(main) + coord_fixed()
    if(get.r){
    if(any(is.na(xy))){
      message("Warning: removed NA values to compute correlation coefs")
      rmsel <- apply(xy, 1, function(r) any(is.na(r)))
      xy <- xy[!rmsel,]
      if(!is.null(DEgsel))
        DEgsel <- DEgsel[!rmsel]
    }

    cc <- cor(xy)[1,2]
    cc2 <- cor(xy, method='spearman')[1,2]
    cctxt <- paste0("r = ", round(cc, 3), '\nrho = ', round(cc2, 3))
    if(add.vd)
      cctxt <- paste0(cctxt, "\n", round(100*(cc^2), 1), "% VarDev")
    if(!is.null(DEgsel)){
      ccDE <- cor(xy[DEgsel,])[1,2]
      cctxt <- paste0(cctxt, "\nr(DE) = ", round(ccDE, 3))
      if(add.vd)
        cctxt <- paste0(cctxt, "\n", round(100*(ccDE^2), 1), "% VarDev(DE)")
    }

    g <- g + annotate(geom="text", hjust=1, vjust=0,
                      x = rg[2],
                      y = rg[1],
                      label = cctxt)
  }

  return(g)
}
```

Make a color transparent (`transp()`) and add an axis with 2 tick marks `twoTicks()`

``` r
transp <- function(col, a=.5){
  # Make a color transparent
  colr <- col2rgb(col)
  return(rgb(colr[1,], colr[2,], colr[3,], a*255, maxColorValue = 255))
}

twoTicks <- function(side = c(1,2), col = NA,
                     col.ticks = 1, las = c(1,2), ...){
  # Generate 2 (edge) ticks on given axis of current plot
  # col = NA & col.ticks = 1 makes the axis line dissapear, but keeps ticks
  for(i in seq_along(side)){
    axt <- axTicks(side[i])
    axis(side[i], at = axt[c(1, length(axt))], col = col,
         col.ticks = col.ticks, las = las[i], ...)
  }
}
```

## [5.4] DE functions

[]

Running a standard DESeq2 workflow (`run_DESeq2_age()`) and extracting a results table (`get_DEres()`).

``` r
run_DESeq2_age <- function(X, p){
  # Run DEseq2 wt vs. mutant (with age covariate)
  require(DESeq2)
  rownames(p) <- colnames(X)
  p$aes <- scale(p$ae) # scale age value
  dds <- DESeqDataSetFromMatrix(countData = X,
                                colData = p,
                                design = ~aes+strain)
  dds <- DESeq(dds, test = "Wald", fitType = "local")
  return(dds)
}

get_DEres <- function(dds, coefname="strain_xrn2_vs_wt"){
  # Get results table from deseq output, managing NAs

  res <- results(dds, name=coefname)
  # manage NAs
  res$padj[is.na(res$padj)] <- 1
  res$log2FoldChange[is.na(res$log2FoldChange)] <- 0

  return(res)
}
```

\

Finding optimal df to fit a reference window (`find_df()`) and converting reference data to counts (`log1ptpm_2rawcounts()`, `ref_2counts()`).

``` r
find_df <- function(ref, w, dfs=2:8){
  # Compute ssq of spline fits to find optimal df

  w.idx <- seq(max(c(1L, which.min(abs(w[1]-ref$time))-1 )),
               min(c(length(ref$time), which.min(abs(w[2]-ref$time))+1)))
  # get time values of window
  ts <- ref$time[w.idx]
  # get PCs of reference window
  pc <- summary(stats::prcomp(t(ref$interpGE[,w.idx]), scale=F, center=T))
  # keep enough components for 99% var
  spc <- sum(pc$importance[3, ] < 0.99)+1
  pcfit <- pc$x[, 1L:spc]
  w <- pc$importance[2, 1L:spc]
  # compute ssq of fit residuals with different dfs
  # for each df, compute weighted ssq of spline fit on
  ssqs <- cbind(sapply(dfs, function(dfi){
    ssq <- sum(w * colSums(
      stats::residuals(stats::lm(pcfit~splines::ns(ts, df=dfi)))^2
    ))
  }))/(length(ts))
  return(ssqs)
}

log1ptpm_2rawcounts <- function(X, glengths, nreadbygl){
  # Transform log1p(tpm) to (artificial) raw counts
  # note : nreadbygl = colSums(rawcounts/genelengths)
  if(length(nreadbygl) != ncol(X))
    stop("nreadbygl != ncol(X)")
  if(length(glengths) != nrow(X))
    stop("glengths must be of length nrow(X)")
  X <- t( (t(exp(X) - 1)/1e6) * nreadbygl ) * glengths
  X[X<0] <- 0
  return(round(X))

}

ref_2counts <- function(ref, ae_values,
                        gltable = wormRef::Cel_genes[,c("wb_id",
                                                        "transcript_length")],
                        avg_librarysize = 25e6){
  # Get expression profiles of given age from a RAPToR reference as
  # (artificial) count data.
  # note : gltable must have WBids as col 1 and gene length as col 2.

  # ref expression profiles at given timepoints :
  rX <- RAPToR::get_refTP(ref, ae_values=ae_values, return.idx = F)

  # transform to counts
  gl <- gltable[match(rownames(rX), gltable[,1]), 2]
  rX <- log1ptpm_2rawcounts(rX, gl, nreadbygl = rep(avg_librarysize,
                                                    ncol(rX))/median(gl))

  return(rX)
}
```

\

Integrating reference data in a DESeq model (`run_DESeq2_ref()`).

``` r
run_DESeq2_ref <- function(X, p, formula, ref,
                           ae_values=NULL, window.extend=1, ns.df=3){
  # Run DEseq2 wt vs. mutant correcting for development with ref. data.

  # Do no not specify age in the formula, it is added directly by the function.
  # Age estimates should either be an 'ae' column of p or given as 'ae_values'.

  require(DESeq2)

  if(!any(colnames(p)=='ae') & is.null(ae_values)){
    stop("Age estimates should either be a column of p, or given to ae_values.")
  }
  if(!is.null(ae_values)){
    p$ae <- ae_values
  }

  ## Extract reference expression profiles in sample time window
  w.rg <- range(p$ae) + c(-window.extend, window.extend)
  w.idx <- seq(
    max(c(1, which.min(abs(w.rg[1]-ref$time))-1 )),
    min(c(length(ref$time), which.min(abs(w.rg[2]-ref$time))+1))
    )
  # ref window time values
  w.ts <- ref$time[w.idx]
  # ref window expression values
  w.GE <- ref_2counts(ref = ref, ae_values = w.ts)

  ## Join ref & sample data
  nX <- ncol(X)
  nR <- ncol(w.GE)

  # get overlapping genes & join expression data
  ovl <- RAPToR::format_to_ref(samp = X, refdata = w.GE, verbose = F)
  Xj <- as.matrix(cbind(ovl$refdata, ovl$samp))

  # get relevant fields from p
  f0 <- as.formula(formula)
  p2 <- p[, attr(terms(f0), "term.labels"), drop=F]
  # get 1st level of each variable
  lev0 <- lapply(p2, function(col) levels(col)[1])

  # join p data and add Time and ref/sample Batch
  pj <- cbind(
    Time = c(w.ts, p$ae),
    Batch = factor(rep(c('r', 's'), c(nR, nX))),
    rbind(do.call(cbind, lapply(lev0, rep, times=nR)), p2) # other terms
    )
  rownames(pj) <- colnames(Xj)

  # Estimate dispersions with samples only
  s0 <- pj$Batch=="s"
  dds0 <- DESeqDataSetFromMatrix(countData = Xj[,s0],
                                 colData = pj[s0,],
                                 design = f0)
  dds0 <- estimateSizeFactors(dds0)
  dds0 <- estimateDispersions(dds0, fitType = "local")
  dd <- dispersions(dds0) # store dispersions
  dd[is.na(dd)] <- 0 # remove NAs

  ## Build full DE model with reference
  # Add Time and ref/sample batch to model formula
  f1 <- update.formula(f0, substitute(
    ~ splines::ns(Time, df = ns.df) + Batch + .,
    list(ns.df=ns.df)
    ))
  dds <-  DESeqDataSetFromMatrix(countData = Xj,
                                 colData = pj,
                                 design = f1)
  dds <- estimateSizeFactors(dds)
  # inject dispersions from sample-only model
  dispersions(dds) <- dd

  dds <- nbinomWaldTest(dds)

  return(dds)
}
```

------------------------------------------------------------------------

# References

Bulteau, Romain, and Mirko Francesconi. 2022. "Real Age Prediction from the Transcriptome with RAPToR." *Nature Methods*, 1--7.

Miki, Takashi S, Sarah H Carl, and Helge Großhans. 2017. "Two Distinct Transcription Termination Modes Dictated by Promoters." *Genes & Development* 31 (18): 1870--79.

------------------------------------------------------------------------

# SessionInfo

```
sessionInfo()
```

    #> R version 4.1.2 (2021-11-01)
    #> Platform: x86_64-pc-linux-gnu (64-bit)
    #> Running under: Ubuntu 22.04.1 LTS
    #>
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so
    #>
    #> locale:
    #>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
    #>  [3] LC_TIME=C                  LC_COLLATE=en_US.UTF-8
    #>  [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=en_US.UTF-8
    #>  [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C
    #>  [9] LC_ADDRESS=C               LC_TELEPHONE=C
    #> [11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C
    #>
    #> attached base packages:
    #>  [1] stats4    splines   parallel  stats     graphics  grDevices utils
    #>  [8] datasets  methods   base
    #>
    #> other attached packages:
    #>  [1] viridis_0.6.2               ROCR_1.0-11
    #>  [3] DESeq2_1.34.0               SummarizedExperiment_1.24.0
    #>  [5] Biobase_2.54.0              MatrixGenerics_1.6.0
    #>  [7] matrixStats_0.63.0          GenomicRanges_1.46.1
    #>  [9] GenomeInfoDb_1.30.1         IRanges_2.28.0
    #> [11] S4Vectors_0.32.4            BiocGenerics_0.40.0
    #> [13] viridisLite_0.4.1           vioplot_0.4.0
    #> [15] zoo_1.8-11                  sm_2.2-5.7.1
    #> [17] ica_1.0-3                   drosoRef_0.2.0
    #> [19] limma_3.50.3                wormRef_0.5.0
    #> [21] RAPToR_1.2.0                ggpubr_0.6.0
    #> [23] ggplot2_3.4.0.9000          BiocStyle_2.22.0
    #>
    #> loaded via a namespace (and not attached):
    #>  [1] colorspace_2.1-0       ggsignif_0.6.4         pryr_0.1.6
    #>  [4] XVector_0.34.0         rstudioapi_0.14        farver_2.1.1
    #>  [7] hexbin_1.28.2          bit64_4.0.5            AnnotationDbi_1.56.2
    #> [10] fansi_1.0.4            codetools_0.2-18       cachem_1.0.6
    #> [13] geneplotter_1.72.0     knitr_1.42             jsonlite_1.8.4
    #> [16] broom_1.0.3            annotate_1.72.0        png_0.1-8
    #> [19] BiocManager_1.30.19    compiler_4.1.2         httr_1.4.4
    #> [22] backports_1.4.1        Matrix_1.5-3           fastmap_1.1.0
    #> [25] cli_3.6.0              htmltools_0.5.4        tools_4.1.2
    #> [28] gtable_0.3.1           glue_1.6.2             GenomeInfoDbData_1.2.7
    #> [31] dplyr_1.1.0            tinytex_0.44           Rcpp_1.0.10
    #> [34] carData_3.0-5          jquerylib_0.1.4        vctrs_0.5.2
    #> [37] Biostrings_2.62.0      nlme_3.1-161           xfun_0.37
    #> [40] stringr_1.5.0          rbibutils_2.2.13       lifecycle_1.0.3
    #> [43] rstatix_0.7.2          XML_3.99-0.13          zlibbioc_1.40.0
    #> [46] scales_1.2.1           RColorBrewer_1.1-3     yaml_2.3.7
    #> [49] gridExtra_2.3          memoise_2.0.1          sass_0.4.5
    #> [52] stringi_1.7.12         RSQLite_2.2.20         highr_0.10
    #> [55] genefilter_1.76.0      BiocParallel_1.28.3    Rdpack_2.4
    #> [58] rlang_1.0.6            pkgconfig_2.0.3        bitops_1.0-7
    #> [61] evaluate_0.20          lattice_0.20-45        purrr_1.0.1
    #> [64] labeling_0.4.2         cowplot_1.1.1          bit_4.0.5
    #> [67] tidyselect_1.2.0       magrittr_2.0.3         bookdown_0.32
    #> [70] R6_2.5.1               magick_2.7.3           generics_0.1.3
    #> [73] DelayedArray_0.20.0    DBI_1.1.3              pillar_1.8.1
    #> [76] withr_2.5.0            mgcv_1.8-41            survival_3.5-0
    #> [79] KEGGREST_1.34.0        abind_1.4-5            RCurl_1.98-1.10
    #> [82] tibble_3.1.8           crayon_1.5.2           car_3.1-1
    #> [85] utf8_1.2.3             rmarkdown_2.20.1       locfit_1.5-9.7
    #> [88] grid_4.1.2             data.table_1.14.6      blob_1.2.3
    #> [91] digest_0.6.31          xtable_1.8-4           tidyr_1.3.0
    #> [94] munsell_0.5.0          beeswarm_0.4.0         bslib_0.4.2

---

# RAPToR-refbuilding

Code []

-   [Show All Code](#)
-   [Hide All Code](#)

# `RAPToR` - Building References

### RAPToR 1.2.0

Romain Bulteau

#### March 2023

# Contents

-   [[1] Introduction](#introduction)
-   [[2] Selecting / Preparing a dataset](#selecting-preparing-a-dataset)
    -   [[2.1] Databases](#databases)
    -   [[2.2] What to look out for](#what-to-look-out-for)
    -   [[2.3] ID formatting](#id-formatting)
    -   [[2.4] Normalize and log expression](#normalize-and-log-expression)
    -   [[2.5] Exploring the data](#exploring-the-data)
        -   [[2.5.1] Correlation](#correlation)
        -   [[2.5.2] Plotting components](#plotting-components)
        -   [[2.5.3] Plotting random genes](#plotting-random-genes)
-   [[3] The gene expression interpolation model (GEIM)](#the-gene-expression-interpolation-model-geim)
    -   [[3.1] Fitting models in component space](#fitting-models-in-component-space)
    -   [[3.2] The GEIM interface](#the-geim-interface)
    -   [[3.3] Defining the appropriate model and parameters](#defining-the-appropriate-model-and-parameters)
        -   [[3.3.1] Model type](#model-type)
        -   [[3.3.2] Model performance](#model-performance)
        -   [[3.3.3] Number of components](#number-of-components)
        -   [[3.3.4] Comparing formulas](#comparing-formulas)
    -   [[3.4] Building a Reference object](#building-a-reference-object)
    -   [[3.5] Validating the reference](#validating-the-reference)
        -   [[3.5.1] Checking model predictions against components](#checking-model-predictions-against-components)
        -   [[3.5.2] Staging samples](#staging-samples)
    -   [[3.6] About aging references](#about-aging-references)
-   [[4] Reference-Building examples](#reference-building-examples)
    -   [[4.1] Example 1 - *C. elegans* larval development](#example-1---c.-elegans-larval-development)
        -   [[4.1.1] Data](#data)
        -   [[4.1.2] Normalization & exploration](#normalization-exploration)
        -   [[4.1.3] Model fitting](#model-fitting)
        -   [[4.1.4] Validation](#validation)
        -   [[4.1.5] Code to generate objects](#code-to-generate-objects)
    -   [[4.2] Example 2 - *D. melanogaster* embryonic development](#example-2---d.-melanogaster-embryonic-development)
        -   [[4.2.1] Data](#data-1)
        -   [[4.2.2] Normalization & exploration](#normalization-exploration-1)
        -   [[4.2.3] Model fitting](#model-fitting-1)
        -   [[4.2.4] Validation](#validation-1)
        -   [[4.2.5] Code to generate objects](#code-to-generate-objects-1)
    -   [[4.3] Example 3 - Uneven sampling of *D. rerio* embryos](#example-3---uneven-sampling-of-d.-rerio-embryos)
        -   [[4.3.1] Data](#data-2)
        -   [[4.3.2] Normalization & Quick look](#normalization-quick-look)
        -   [[4.3.3] Model fitting](#model-fitting-2)
        -   [[4.3.4] Validation](#validation-2)
        -   [[4.3.5] Model fitting without age transformation](#model-fitting-without-age-transformation)
        -   [[4.3.6] Code to generate objects](#code-to-generate-objects-2)
    -   [[4.4] Example 4 - *C. elegans* aging](#example-4---c.-elegans-aging)
        -   [[4.4.1] Data](#data-3)
        -   [[4.4.2] Normalization & Quick look](#normalization-quick-look-1)
        -   [[4.4.3] Filtering genes](#filtering-genes)
        -   [[4.4.4] Model fitting](#model-fitting-3)
        -   [[4.4.5] Validation](#validation-3)
        -   [[4.4.6] Code to generate objects](#code-to-generate-objects-3)
-   [References](#references)
-   [SessionInfo](#sessioninfo)

# [1] Introduction

This vignette is specifically focused on building `RAPToR` references, and assumes knowledge of general use of the package ([see the main `"RAPToR"` vignette](RAPToR.html)).

Building references is a key aspect of `RAPToR`, as samples need an appropriate reference to be staged. References are broadly usable once built, so they are worth the trouble to set up. This vignette, will go through the general workflow of building a reference from selecting an appropriate dataset, to validating a model for interpolation.

Along the explanations, examples will be given using *C. elegans* larval development time-series data published by [Aeschimann et al. (2017)] and [Hendriks et al. (2014)], stored in `dsaeschimann2017` and `dshendriks2014` respectively. [Code to create these objects can be found at the end of this vignette](#code-to-generate-objects).

Further examples are given at the end for broader coverage of reference-building needs.

# [2] Selecting / Preparing a dataset

We start from a gene expression profiling time-series spanning the developmental stages of interest for your organism. Thankfully, these time-series experiments are increasingly numerous in the literature and many models will already have some kind of useable data. You may also make (or already have) your own time-series on hand.

## [2.1] Databases

There are a few expression profiling databases to search in and download data from. The most well-known are the [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/) and the [Array Express](https://www.ebi.ac.uk/arrayexpress/).

Both of these databases have APIs to get their data directly from R (*e.g* [the `GEOquery` package](https://bioconductor.org/packages/release/bioc/html/GEOquery.html), used in multiple data loading scripts of the RAPToR vignettes).

## [2.2] What to look out for

Several points of the experimental design should be kept in mind when selecting data for a reference.

-   ***Are there replicates ?*** If so, good. This means you can confirm the dynamics in the data are not noise. A slightly sparser time-series with replicates could be a better candidate to build a reference than an experiment with one batch.
-   ***Is the time point sampling even ?*** Profiling is still expensive, so time-series sometimes account for dynamic ranges of development (*i.e.* early or fast-changing stages are more sampled). Even sampling is better for spline fitting, but data can still otherwise be used by transforming time values so they are more uniform (see [*Example 3*](#example-3---uneven-sampling-of-d.-rerio-embryo).
-   ***What's the developmental range ?*** Usually, the bigger, the better! Aging (as opposed to development) references need special handling (see [*About aging references*](#about-aging-references), and [*Example 4*](#example-4---c.-elegans-aging))
-   ***Is the profiling done on whole-organism, dissected parts, single-cells ?*** You should aim for profiling that matches your sample type. Whole-organism reference for whole-organism samples, dissected tissue/organ reference for similar samples. Dissected tissue (or single-cell) samples are often sparser than whole-organism for biological and technical reasons. To account for this, we recommend to filter out genes with median \\(log(TPM+1)\<0.5\\) across reference samples.
-   ***What's the profiling technology ?*** RNA-seq is much, *much* cleaner than microarray data, so it's preferable. There is however no trouble staging RNA-seq samples on microarray references and vice-versa. RNA-seq counts should have some within-sample normalization (*e.g.* TPM) to reflect gene expression accurately across samples.

## [2.3] ID formatting

Bioinformatics is a field plagued by diverse and fast-changing sets of IDs for genes, transcripts, etc. When you build a reference, you should always convert it to IDs that are **conventional and stable**. We like to use the organism-specific IDs (*e.g*, Wormbase for *C. elegans* : `WBGene00016153`, Flybase for *Drosophila* : `FBgn0010583`).

The [ensembl biomart](https://www.ensembl.org/info/data/biomart/index.html) or its associated R package [`biomaRt`](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) are a very useful resource to get gene, transcript or probe ID lists for conversion.

Below is a code snippet to get gene IDs for drosophila with `biomaRt`.

```
library(biomaRt)

# setup connection to ensembl
mart <- biomaRt::useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")

# get list of attributes
droso_genes <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                             "ensembl_transcript_id",
                                             "external_gene_name",
                                             "flybase_gene_id"),
                              mart = mart)

head(droso_genes)
#>   ensembl_gene_id ensembl_transcript_id external_gene_name flybase_gene_id
#> 1     FBgn0015380           FBtr0330209                drl     FBgn0015380
#> 2     FBgn0015380           FBtr0081195                drl     FBgn0015380
#> 3     FBgn0010356           FBtr0088240               Taf5     FBgn0010356
#> 4     FBgn0266023           FBtr0343232     lncRNA:CR44788     FBgn0266023
#> 5     FBgn0250732           FBtr0091512               gfzf     FBgn0250732
#> 6     FBgn0250732           FBtr0334671               gfzf     FBgn0250732
```

When multiple probe or transcript IDs match a single gene ID, we usually sum-aggregate counts and mean-aggregate other expression values (microarray or already-processed RNAseq as TPMs). This is taken care of with the `format_ids()` function.

## [2.4] Normalize and log expression

It's common practice to normalize expression data (*e.g.* to account for technical bias). You may deal with many different profiling technologies when building references, and may even join different datasets together for a single reference.

To stay as consistent as possible, we apply quantile-normalization on all data regardless of source or type. For this, we use the `normalizeBetweenArrays()` function from [`limma`](https://bioconductor.org/packages/release/bioc/html/limma.html).

We also log-transform the data with \\(log(X + 1)\\).

```
library(RAPToR)
library(limma)

dsaeschimann2017$g <- limma::normalizeBetweenArrays(dsaeschimann2017$g,
                                                    method = "quantile")
dsaeschimann2017$g <- log1p(dsaeschimann2017$g)

dshendriks2014$g <- limma::normalizeBetweenArrays(dshendriks2014$g,
                                                  method = "quantile")
dshendriks2014$g <- log1p(dshendriks2014$g)
```

## [2.5] Exploring the data

It's good practice to explore the data a little before moving on.

We start with quick look at the data: expression matrix and associated sample (or phenotypic) data.

```
dsaeschimann2017$g[1:5,1:3]
#>                let.7.n2853._18hr let.7.n2853._20hr let.7.n2853._22hr
#> WBGene00007063         2.1206604          2.469532          2.373273
#> WBGene00007064         2.1621558          2.260804          2.661102
#> WBGene00007065         2.7763061          2.847833          2.727037
#> WBGene00003525         0.9434159          2.466223          2.609585
#> WBGene00007067         1.0787531          1.081964          1.350796

head(dsaeschimann2017$p, n = 5)
#>                        title geo_accession           organism_ch1
#> GSM2113587 let.7.n2853._18hr    GSM2113587 Caenorhabditis elegans
#> GSM2113588 let.7.n2853._20hr    GSM2113588 Caenorhabditis elegans
#> GSM2113589 let.7.n2853._22hr    GSM2113589 Caenorhabditis elegans
#> GSM2113590 let.7.n2853._24hr    GSM2113590 Caenorhabditis elegans
#> GSM2113591 let.7.n2853._26hr    GSM2113591 Caenorhabditis elegans
#>                  strain time in development:ch1 age
#> GSM2113587 let-7(n2853)                18 hours  18
#> GSM2113588 let-7(n2853)                20 hours  20
#> GSM2113589 let-7(n2853)                22 hours  22
#> GSM2113590 let-7(n2853)                24 hours  24
#> GSM2113591 let-7(n2853)                26 hours  26
```

### [2.5.1] Correlation

Outliers in time series data can be revealed using sample-sample correlation heatmaps or boxplots. This also shows the clear correlation between samples of similar development.

```
cor_dsaeschimann2017 <- cor(dsaeschimann2017$g, method = "spearman")
```

``` r
ord <- order(dsaeschimann2017$p$age) # order by development
heatmap(cor_dsaeschimann2017[ord, ord], Colv = NA, Rowv = NA,
        scale = "none", keep.dendro = F, margins = c(1,1),
        # see Example 3 for transp() function.
        RowSideColors = transp(as.numeric(dsaeschimann2017$p$strain[ord])),
        labRow = "", labCol = "")
# add time labels
par(xpd = T)
mtext(text = paste0(unique(dsaeschimann2017$p$age), 'h'), side = 1, line = 4,
      at = seq(-.1,1.05, l = 11))
```

``` r
boxplot(cor_dsaeschimann2017~interaction(dsaeschimann2017$p$strain,
                                         dsaeschimann2017$p$age),
        col = transp(1:4), # see Example 3 for transp() function.
        xaxt = "n", at = seq(1,44, l = 55)[c(T,T,T,T,F)],
        ylab = "Spearman correlation", xlab = "Chronological age (h)")
#add time labels
axis(side = 1, at = seq(2,42, l = 11),
     labels = paste0(unique(dsaeschimann2017$p$age), 'h'))
legend(23,.87, fill = transp(1:4),
       legend = c("let-7", "lin-41", "let-7/lin-41", "N2"),
       bty = "n")
```

### [2.5.2] Plotting components

Dimension reductions such as PCA or ICA (Principal or Independent Component Analysis) yield components that summarize the gene expression landscape and dynamics ([Alter, Brown, and Botstein (2000)]). Plotting these components with respect to time is a good way to grasp general dynamics in the data.

For PCA, we want to perform a non-scaled, centered PCA. Centering is done gene-wise, not sample-wise (hence the matrix rotation below).

```
pca_dsaeschimann2017 <- stats::prcomp(t(dsaeschimann2017$g), rank = 25,
                                      center = TRUE, scale = FALSE)
```

``` r
par(mfrow = c(2,4), pty='s')

# for components 1:8
invisible(sapply(1:8, function(i){
  # plot points colored by strain
  plot(dsaeschimann2017$p$age, pca_dsaeschimann2017$x[,i],
       lwd = 2, col = dsaeschimann2017$p$strain,
       xlab = "Chronological age", ylab = "PC", main = paste0("PC", i))
  # connect the dots per strain
  sapply(levels(dsaeschimann2017$p$strain), function(l){
    s <- which(dsaeschimann2017$p$strain == l)
    points(dsaeschimann2017$p$age[s], pca_dsaeschimann2017$x[s,i],
           col = which(levels(dsaeschimann2017$p$strain) == l),
           type = 'l', lty = 2)
  })
  if(i == 1) # add legend on first panel
    legend("topleft", bty = 'n',
           legend = c("let-7", "lin-41", "let-7/lin-41", "N2"),
           pch = c(rep(1, 4)), lty = c(rep(NA, 4)),
           col = c(1:4), lwd = 3)
}))
```

In this *C. elegans* larval development data, we see very consistent dynamics across different strains. Also, PC2 and PC3 capture an oscillatory dynamic which is characteristic of molting in *C. elegans*.

### [2.5.3] Plotting random genes

We can also plot a few random genes, to get a grasp on the noise in the data.

``` r
# select random genes
set.seed(10)
gtp <- sample(nrow(dsaeschimann2017$g), size = 4)

par(mfrow = c(1,4), pty='s')
invisible(sapply(gtp, function(i){
  plot(dsaeschimann2017$p$age, dsaeschimann2017$g[i,],
       lwd = 2, col = dsaeschimann2017$p$strain,
       xlab = "Chronological age", ylab = "log(TPM+1)",
       main = rownames(dsaeschimann2017$g)[i])
  sapply(levels(dsaeschimann2017$p$strain), function(l){
    s <- which(dsaeschimann2017$p$strain == l)
    points(dsaeschimann2017$p$age[s], dsaeschimann2017$g[i,s],
           col = which(levels(dsaeschimann2017$p$strain) == l),
           type = 'l', lty = 2)
  })
  if(i == gtp[4])
    legend("topleft", bty = 'n',
           legend = c("let-7", "lin-41", "let-7/lin-41", "N2"),
           pch = c(rep(1, 4)), lty = c(rep(NA, 4)),
           col = c(1:4), lwd = 3)
}))
```

\
\

# [3] The gene expression interpolation model (GEIM)

Increasing the resolution of a profiling time series is a very unbalanced regression problem. We want to predict tens of thousands of dependent variables (genes) with a few independent variables (time, batch, ...).

## [3.1] Fitting models in component space

In order to fit a large number of output variable (genes), we project the data in a linear dimensionally-reduced space to interpolate at component-level, before re-projecting the data back to gene-level. We do this with Principal Components or Independant Components ( [Independant Component Analysis](https://en.wikipedia.org/wiki/Independent_component_analysis) ).

Both PCA and ICA perform the same type of linear transformation on the data (they just optimize different criteria). We get the following :

[\\[ X\_{(m\\times n)} = G\_{(m\\times c)}S\^{T}\_{(n\\times c)} \\]] with \\(X\\), the matrix of \\(m\\) genes by \\(n\\) samples, \\(G\\) the gene loadings (\\(m\\) genes by \\(c\\) components) and \\(S\^T\\) the sample scores (\\(n\\) samples by \\(c\\) components). When performing PCA (or ICA) on gene expression data, \\(S\\) is what's usually plotted (e.g. PC1 vs. PC2) to see how samples are grouped in the component space. It's what we plotted earlier in the section on observing data, for instance.

[Alter, Brown, and Botstein (2000)] demonstrated that singular value decomposition of gene expression data can be taken as "eigengenes", giving a global picture of the expression landscape and dynamics with a few components. GEIMs use this property. We fit a model on the columns of \\(S\^T\\) (eigengenes), predict in the component space, and reconstruct the gene expression data by a matrix product with the gene loadings.

We propose 2 model types : Generalized Additive Models (GAMs, the default) and Generalized Linear Models (GLMs). GAMs use `gam()` from the [`mgcv`](https://cran.r-project.org/web/packages/mgcv/index.html) package, and GLMs use the `glm()` function of the `stats` core package.

A standard R formula will be specified for the model. This formula can include the tools and functions generally used with `gam()` or `glm()`, most notably the variety of polynomial or smoothing splines implemented through the `s()` function of `mgcv` for GAMs. GLMs can also use splines from the `splines` core package, such as `ns()` for natural cubic splines.

## [3.2] The GEIM interface

Gene Expression Interpolation Models (GEIMs), are built with the `ge_im()` function, which outputs a `geim` object. This function inputs 3 key arguments :

-   `X` : the gene expression matrix of your time-series (genes as rows, samples as columns)
-   `p` : a dataframe of phenotypic data, samples as rows. This should include the age/time variable and any other covariates you want to include in the model (*e.g* batch, strain)
-   `formula` : the model formula. This should be a standard R formula, using terms found in `p`. **It must start with `X~`**.

Another important argument is the **n**umber of **c**omponents to interpolate on, `nc`.

For example, using the `dsaeschimann2017` dataset we could build the following model.

```
m_dsaeschimann2017 <- ge_im(X = dsaeschimann2017$g,
                            p = dsaeschimann2017$p,
                            formula = "X ~ s(age, bs = 'ts') + strain",
                            nc = 32)
```

Additional parameters are detailed in the documentation of the function `?ge_im()`.

Note that a single model formula is specified and applied to all the components, but models are fitted independently on each components.

Model predictions can be generated with the `predict()` function, as for any standard R model.

## [3.3] Defining the appropriate model and parameters

### [3.3.1] Model type

There are 5 types of GEIMs:

-   A GAM on PCA components (`method="gam",dim_red="pca"`) (default)
-   A GLM on PCA components (`method="glm",dim_red="pca"`)
-   A GAM on ICA components (`method="gam",dim_red="ica"`)
-   A GLM on ICA components (`method="glm",dim_red="ica"`)
-   A gene-by-gene linear model directly on the gene expression matrix (`method="limma"`)

Our default GEIM is fitting GAMs on PCA components, which is a robust choice when applying a smoothing spline to the data.

When using GAMs, there can be no interaction between terms (by definition). It is possible to include a `by` argument to the `s()` function, which essentially corresponds to separate model fits on each level of the specified group variable (and thus requires enough data to do so).

PCA and ICA interpolation usually yield near-identical results (ICA components tend to be more biologically meaningful or interpretable). ICA tends to outperform PCA when the data is very noisy (by design, since ICA essentially does signal extraction). It is however slower, especially if `nc` is large.

With the last option (`"limma"`), a model is fit per gene with the `lmFit()` function of the `limma` package. This approach makes no effort to reduce the dimensionality of the problem (`dim_red` and `nc` arguments are ignored) which means there is no information loss or bias introduced by dimension reduction. It is however very sensitive to noise.

### [3.3.2] Model performance

The `mperf()` function computes multiple indices to evaluate model performance.

In the formulas below, \\(X\\) corresponds to the input gene expression matrix (\\(m\\) genes as rows, \\(n\\) samples as columns), \\(\\hat{X}\\) to the model predictions. \\(x_i\\) corresponds to row \\(i\\) of matrix \\(X\\) and \\(x_i\^{(j)}\\) to sample \\(j\\) of that row. This notation is derived from the general regression problem, where \\(X\^T\\) corresponds to the set of \\(m\\) dependant variables to predict.

-   `aCC` : average Correlation Coefficient.

[\\[ aCC = \\frac{1}{m}\\sum\^{m}\_{i=1}{CC} = \\frac{1}{m}\\sum\^{m}\_{i=1}{\\cfrac{\\sum\^{n}\_{j=1}{(x_i\^{(j)}-\\bar{x}\_i)(\\hat{x}\_i\^{(j)}-\\bar{\\hat{x}}\_i)}}{\\sqrt{\\sum\^{n}\_{j=1}{(x_i\^{(j)}-\\bar{x}\_i)\^2(\\hat{x}\_i\^{(j)}-\\bar{\\hat{x}}\_i)\^2}}}} \\]]

-   `aRE` : average Relative Error.

[\\[ a\\delta = \\frac{1}{m}\\sum\^{m}\_{i=1}{\\delta} = \\frac{1}{m} \\sum\^{m}\_{i=1} \\frac{1}{n} \\sum\^{n}\_{j=1} \\cfrac{\| x_i\^{(j)} - \\hat{x}\_i\^{(j)} \| }{x_i\^{(j)}} \\]]

-   `MSE` : Mean Squared Error.

[\\[ MSE = \\frac{1}{m} \\sum\^{m}\_{i=1} \\frac{1}{n} \\sum\^{n}\_{j=1} (x_i\^{(j)} - \\hat{x}\_i\^{(j)} )\^2 \\]]

-   `aRMSE` : average Root MSE.

[\\[ aRMSE = \\frac{1}{m}\\sum\^{m}\_{i=1}{RMSE} = \\frac{1}{m} \\sum\^{m}\_{i=1} \\sqrt{\\cfrac{\\sum\^{n}\_{j=1} (x_i\^{(j)} - \\hat{x}\_i\^{(j)} )\^2}{n}} \\]]

Note that these indices are computed and averaged *with respect to variables (genes), not observations*. `mperf()` outputs either the overall (averaged) index, or values per-gene (`global` parameter).

```
# global values
g_mp <- mperf(dsaeschimann2017$g, predict(m_dsaeschimann2017), is.t = T)
g_mp
#> $aCC
#> [1] 0.7977299
#>
#> $aRE
#> [1] 0.1301014
#>
#> $MSE
#> [1] 0.01431891
#>
#> $aRMSE
#> [1] 0.1196617

# per gene values
ng_mp <- mperf(dsaeschimann2017$g, predict(m_dsaeschimann2017), is.t = T,
               global = F)
ng_mp <- lapply(ng_mp, na.omit) # remove NAs (eg. 0 variance genes)
ng_mp$aRE <- ng_mp$aRE[ng_mp$aRE < Inf] # remove Inf values (/0)
```

Model performance is reflected by the distributions of prediction indices over all genes.

``` r
par(mfrow = c(2,2))
invisible(sapply(names(ng_mp), function(idx){
  rg <- range(na.omit(ng_mp[[idx]]))

  # estimate density curve
  d <- density(na.omit(ng_mp[[idx]]), from = rg[1], to = rg[2])

  plot(d, main = paste0(gsub("a", "", idx, fixed = T),
                        " density (", length(ng_mp[[idx]]), " genes)"),
       xlab = idx, lwd = 2)
  # display global value
  abline(v = g_mp[[idx]], lty = 2, lwd = 2, col = "firebrick")
  text(g_mp[[idx]], .9*max(d$y), pos = 4, labels = idx,
       font = 2, col = "firebrick")
}))
```

This data is quite clean, so model fits are very good across the board.

### [3.3.3] Number of components

By default, the number of components to interpolate is set to the number of samples. However, we recommend setting a cutoff on explained variance of PCA (or ICA) components.

For example, on the `dsaeschimann2017` dataset, we set the threshold at \\(99\\%\\) :

```
nc <- sum(summary(pca_dsaeschimann2017)$importance[3,] < .99) + 1
nc
#> [1] 32
```

This threshold must be set in accordance with the noise in the data. For example, in very noisy data, would you consider that \\(99\\%\\) of the variance in the dataset corresponds to meaningful information ?

You can also define `nc` by plotting your components and stopping after the components stop capturing meaningful variation (dynamics) with respect to time/age. We define components with 'intelligible dynamics' with respect to time as those where a model fit explains \\(\>0.5\\) of the deviance (noisy components with no dynamics have poor fits).

In very noisy data, you may have to keep very few components (or even a single component) for the interpolation.

### [3.3.4] Comparing formulas

Choosing from different splines (and/or parameters) can be done with cross-validation (CV). The `ge_imCV()` function inputs the `X`, `p`, and a `formula_list` to test, as well as other potential CV parameters (*e.g.* training set size).

The default training/validation set ratio is `cv.s=0.8`, so \\(80\\%\\) of the data is used to build the model and the remainder for validation. When including (factor) covariates in the model, the training set is built such that all groups are proportionately represented in the training set (based on terms of the first formula in the list). The number of CV repeats is defined by `cv.n`.

The model type (GAM/GLM and PCA/ICA) is fixed for all formulas in one call of `ge_imCV()`.

`ge_imCV()` returns model performance indices from `mperf()` (excluding `aCC` to reduce computing time). Indices are computed on the validation set (CV Error, cve) *and* on the training set (Model PerFormance, mpf).

Below we choose between 4 available GAM smooth terms on the `dsaeschimann2017` GEIM.

```
smooth_methods <- c("tp", "ts", "cr", "ps")
flist <- as.list(paste0("X ~ s(age, bs = \'", smooth_methods, "\') + strain"))
cat(unlist(flist), sep = '\n')
#> X ~ s(age, bs = 'tp') + strain
#> X ~ s(age, bs = 'ts') + strain
#> X ~ s(age, bs = 'cr') + strain
#> X ~ s(age, bs = 'ps') + strain

set.seed(2)
cv_dsaeschimann2017 <- ge_imCV(X = dsaeschimann2017$g,
                               p = dsaeschimann2017$p,
                               formula_list = flist,
                               cv.n = 10, nc = nc,
                               nb.cores = 4)
#> CV on 4 models. cv.n = 10 | cv.s = 0.8
#>
#> ...Building training sets
#> ...Setting up cluster
#> ...Running CV
#> ...Cleanup and formatting
```

```
plot(cv_dsaeschimann2017, names = paste0("bs = ", smooth_methods), outline = F,
     swarmargs = list(cex = .8), boxwex=.5)
```

From the plots above, we can see the different splines perform similarly. All could work. We choose `ts` (a thin-plate regression spline), as it minimizes CV error without much impact on model performance.

Extra spline parameters can also be specified to the model. With `s()`, you can give the spline's basis dimension (\\(\\simeq\\) number of knots) to use with the `k` parameter. By default, the spline is a *penalized spline*, so it will not necessarily use `k` knots, but it will stay below that value. By setting `fx=TRUE`, a spline with `k` basis dimension is forced. Note this fits a spline of dimension `k` on *all* components, whereas the penalized spline will adjust. The `s()` or `choose.k` documentation gives further information.

In our experience, the parameter estimation done by `gam()` is usually sufficient and rarely requires tweaking. Below are examples of spline parameter tweaking with the `dsaeschimann2017` data.

```
ks <- c(4,6,8,10)
flistk <- as.list(c(
  "X ~ s(age, bs =  'cr') + strain",
  paste0("X ~ s(age, bs =  'cr', k = ", ks , ") + strain"),
  paste0("X ~ s(age, bs =  'cr', k = ", ks , ", fx=TRUE) + strain")
  ))
cat(unlist(flistk), sep = '\n')
#> X ~ s(age, bs =  'cr') + strain
#> X ~ s(age, bs =  'cr', k = 4) + strain
#> X ~ s(age, bs =  'cr', k = 6) + strain
#> X ~ s(age, bs =  'cr', k = 8) + strain
#> X ~ s(age, bs =  'cr', k = 10) + strain
#> X ~ s(age, bs =  'cr', k = 4, fx=TRUE) + strain
#> X ~ s(age, bs =  'cr', k = 6, fx=TRUE) + strain
#> X ~ s(age, bs =  'cr', k = 8, fx=TRUE) + strain
#> X ~ s(age, bs =  'cr', k = 10, fx=TRUE) + strain

cv_dsaeschimann2017k <- ge_imCV(X = dsaeschimann2017$g,
                                p = dsaeschimann2017$p,
                                formula_list = flistk,
                                cv.n = 10, nc = nc,
                                nb.cores = 4)
#> CV on 9 models. cv.n = 10 | cv.s = 0.8
#>
#> ...Building training sets
#> ...Setting up cluster
#> ...Running CV
#> ...Cleanup and formatting
```

``` r
par(mar = c(7,4,3,1))
plot(cv_dsaeschimann2017k,
     names = c("na", paste0("k=", ks, rep(c("", ", fx=T"), each = 4))),
     outline = F,
     col = transp(c("royalblue", rep(c(1, "firebrick"), each = 4))),
     tcol = c("royalblue", rep(c(1, "firebrick"), each = 4)),
     names.arrange = 5, swarmargs = list(cex = .8))
```

## [3.4] Building a Reference object

A `ref` object can be built from a GEIM using `make_ref()`, specifying interpolation resolution and relevant metadata (see `?plot.geim`).

On our `dsaeschimann2017` example :

```
r_dsaeschimann2017 <- make_ref(
  m = m_dsaeschimann2017,       # geim
  n.inter = 100,                # interpolation resolution
  t.unit = "h past egg-laying (25°C)", # time unit
  cov.levels = list("strain"= "N2"),  # covariate lvls to use for interpolation
  metadata = list("organism"= "C. elegans", # any metadata
                  "profiling"= "bulk RNAseq"))
```

As any R model, GEIMs have a `predict()` function (called internally by `make_ref()`) which can be used to manually get predictions either in component space or at the gene level. This can be useful for a deeper look at the model.

```
# first generate the new predictor data
n.inter <- 100 # nb of new timepoints
newdat <- data.frame(
  age = seq(min(dsaeschimann2017$p$age),
            max(dsaeschimann2017$p$age), l = n.inter),
  strain = rep("N2", n.inter) # we want to predict as N2
  )
head(newdat)
#>        age strain
#> 1 18.00000     N2
#> 2 18.20202     N2
#> 3 18.40404     N2
#> 4 18.60606     N2
#> 5 18.80808     N2
#> 6 19.01010     N2

# predict at gene level
pred_m_dsaeschimann2017 <- predict(m_dsaeschimann2017,
                                   newdata = newdat)
# predict at component level
pred_m_dsaeschimann2017_comp <- predict(m_dsaeschimann2017,
                                        newdata = newdat, as.c = TRUE)
```

## [3.5] Validating the reference

After building a reference, we check interpolation results by:

-   Checking model fits on components (plots)
-   Staging the samples on their own interpolated data, or better (if possible) stage another independent time-series on your reference for external validation.

### [3.5.1] Checking model predictions against components

Checking predictions against components allows you to immediately see if some dynamics get mishandled by the model, or if there is over fitting. It's acceptable to have slight offsets, no fit is perfect.

Plotting a model and a reference object (or equivalent metadata) shows component interpolation, with deviance explained (DE) and relative error (RE) for each component (this information is also returned by the plot function). DE can be used to define components with "intelligible" dynamics (*w.r.t.* time), when \\(DE\>0.5\\). In noisy data, this distinction can be useful to remove components which do not reflect meaningful developmental variation (but rather noise).

Predictions of the first few components from `dsaeschimann2017` are plotted below.

```
par(mfrow = c(2,4))
fit_vals <- plot(m_dsaeschimann2017, r_dsaeschimann2017, ncs=1:8,
                 col = dsaeschimann2017$p$strain, col.i = 'royalblue')
```

Of note, we are predicting model values as `N2` (lightblue). While all strains are shown on the plots, some model parameters depend on the selected `N2` strain.

The fit values are also returned by the plot function.

```
head(fit_vals)
#>     component.var.exp        r2 deviance.expl relative.err
#> PC1           0.68673 0.9972253     0.9972253    0.1658575
#> PC2           0.09546 0.9230361     0.9230361    0.3898380
#> PC3           0.08096 0.8783276     0.8783276    1.2464485
#> PC4           0.03242 0.9249808     0.9249808    5.3703222
#> PC5           0.01875 0.8510695     0.8510695    1.6352080
#> PC6           0.01274 0.9249475     0.9249475    2.9281560
```

You may notice some noisy components get "flattened", with a null model fitted. These components can be left in or removed as they generally have little to no impact on interpolation at the gene level (representing a minuscule part of total variance in the data). This can actually get rid of unwanted variation.

``` r
par(mfrow = c(1,3), pty='s')
fit_vals <- plot(m_dsaeschimann2017, r_dsaeschimann2017, ncs=11:13,
                 col = dsaeschimann2017$p$strain, col.i = 'royalblue',
                 l.pos = 'bottomleft')
```

### [3.5.2] Staging samples

Staging the samples used to build the reference on their interpolated version is a good first test. Then, staging another time-series from the literature on your reference is the best validation, if such data exists. This external validation confirms the interpolated dynamics indeed correspond to development processes and not a dataset-specific effect (which is unlikely, but not impossible).

We can use `dshendriks2014` for external validation.

```
ae_test_dsaeschimann2017 <- ae(dsaeschimann2017$g, r_dsaeschimann2017)
ae_test_dshendriks2014 <- ae(dshendriks2014$g, r_dsaeschimann2017)
```

## [3.6] About aging references

Unlike development, where gene expression changes are robust across individuals (and to an extent across species), expression changes along aging are much more subtle, more variable and noisier across experiments. Individuals can age at different rates, in different ways, and the process of aging itself is also still poorly understood. This makes it difficult for RAPToR to work with aging data "as is", and we recommend to strengthen the aging signal by restricting genes to an informative set.

Empirically, we have found that genes with monotonous trends along aging make for a robust choice. Genes whose expression correlates with age values above a given threshold can thus be selected to build a reference from a single component. We use this strategy in our article ([Bulteau and Francesconi (2022)]), and replicate it in [*Example 4*](#example-4---c.-elegans-aging).

Aging profiling data is particularly tricky to use for reference-building because on top of the many known factors that biologically influence aging, unknown differences in experimental procedures between labs can also influence aging, as evidenced by [Lithgow, Driscoll, and Phillips (2017)]. They further show that the t-zero can be different between labs (*e.g.* egg-laying, hatching, feeding), some start counting at "day 1", others at "day 0". It is therefore likely that inferred age units will differ from known chronological age of staged samples. Time-series should however show clear correlation between chronological and estimated age.

# [4] Reference-Building examples

Here are a few examples of reference building and validation on different organisms.

-   [*C. elegans* larval development](#example-1---c.-elegans-larval-development)
-   [*D. melanogaster* embryonic development](#example-2---d.-melanogaster-embryonic-development)
-   [*Danio rerio* embryonic development with uneven sampling](#example-3---uneven-sampling-of-d.-rerio-embryo)
-   [*C. elegans* aging](#example-4---c.-elegans-aging)

\
\

## [4.1] Example 1 - *C. elegans* larval development

### [4.1.1] Data

In this example, we use two *C. elegans* RNAseq time-series:

1.  A high-resolution time series of late larval development published by [Hendriks et al. (2014)], `dshendriks2014`, used to build the reference. (Accession : [GSE52861](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52861))
2.  A time series of larval development in 4 different strains published by [Aeschimann et al. (2017)], `dsaeschimann2017`, used for external validation. (Accession : [GSE80157](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80157))

The data is the same used in the vignette above, but we flip which dataset is the reference and which is used for external validation. Code to generate `dsaeschimann2017` and `dshendriks2014` can be found [at the end of this section](#code-to-generate-objects)

### [4.1.2] Normalization & exploration

We start by normalizing the data.

```
library(RAPToR)
library(limma)

dshendriks2014$g <- limma::normalizeBetweenArrays(dshendriks2014$g,
                                                  method = "quantile")
dshendriks2014$g <- log1p(dshendriks2014$g)

dsaeschimann2017$g <- limma::normalizeBetweenArrays(dsaeschimann2017$g,
                                                    method = "quantile")
dsaeschimann2017$g <- log1p(dsaeschimann2017$g)
```

Check the contents of the expression matrix and pheno data.

```
dshendriks2014$g[1:5,1:3]
#>                contDevA_N2_21h contDevA_N2_22h contDevA_N2_23h
#> WBGene00007063       1.4124620       1.5261100       1.4106924
#> WBGene00007064       1.8814221       2.0860505       2.2545866
#> WBGene00007065       2.7669727       2.6197855       2.5342036
#> WBGene00003525       0.3651526       0.9841291       1.7629401
#> WBGene00007067       1.0067865       0.9129110       0.7664744

head(dshendriks2014$p, n = 5)
#>                      title geo_accession time in development:ch1 age
#> GSM1277118 contDevA_N2_21h    GSM1277118                22 hours  22
#> GSM1277119 contDevA_N2_22h    GSM1277119                23 hours  23
#> GSM1277120 contDevA_N2_23h    GSM1277120                24 hours  24
#> GSM1277121 contDevA_N2_24h    GSM1277121                25 hours  25
#> GSM1277122 contDevA_N2_25h    GSM1277122                26 hours  26
```

Correlation between the samples.

``` r
cor_dshendriks <- cor(dshendriks2014$g, method = "spearman")
ord <- order(dshendriks2014$p$age)
heatmap(cor_dshendriks[ord, ord],
        Colv = NA, Rowv = NA, scale = "none",
        keep.dendro = F, margins = c(1,5),
        labRow = "", labCol = "")
# add time labels
par(xpd = T)
mtext(text = paste0(dshendriks2014$p$age),
      side = 1, line = 4.1, at = seq(-.17,.86, l = 16), cex = .8)

# add color legend
l.values <- seq(min(cor_dshendriks), max(cor_dshendriks), l = 10)
image(x = c(.95,1), y = seq(0.6,1, l = 10), useRaster = T,
      z = matrix(l.values, ncol = 10),
      col = hcl.colors(12, "YlOrRd", rev = TRUE), add = T)
text(.975, 1, pos = 3, labels = expression(rho), font = 2)
text(1, y = seq(0.6,1, l = 10), pos = 4,
     labels = round(l.values, 2), cex = .6)

mtext(at = 1.0, line = 4.1, side = 1, text = "(hours)", cex = .8)
```

``` r
boxplot(cor_dshendriks, names = paste0(dshendriks2014$p$age,'h'),
        boxwex=.5, las=1)
```

Now, plotting principal components.

```
pca_dshendriks <- summary(stats::prcomp(t(dshendriks2014$g), rank = 25,
                                        center = TRUE, scale = FALSE))
```

``` r
par(mfrow = c(2,3), pty='s')
invisible(sapply(1:6, function(i){
  plot(dshendriks2014$p$age, pca_dshendriks$x[,i], lwd = 2, type = 'b',
       xlab = "Chronological age", ylab = "PC", main = paste0("PC", i))
}))
```

### [4.1.3] Model fitting

We keep enough components to explain \\(99\\%\\) of the variance in the data, since it is very clean.

```
nc <- sum(pca_dshendriks$importance[3,] < .99) + 1
nc
#> [1] 10
```

We then fit a GAM model with a cubic spline on age on PCA components.

```
m_dshendriks2014 <- ge_im(X = dshendriks2014$g, p = dshendriks2014$p,
                          formula = "X ~ s(age, bs = 'cr')", nc = nc)
```

### [4.1.4] Validation

First, check global model performance.

```
# global model performance
mp_hendriks <- mperf(dshendriks2014$g, predict(m_dshendriks2014), is.t = T)
print(do.call(cbind, mp_hendriks))
#>            aCC       aRE         MSE      aRMSE
#> [1,] 0.8801989 0.1018511 0.005021943 0.07086567
```

And then per gene.

```
ng_mp_hendriks <-  mperf(dshendriks2014$g, predict(m_dshendriks2014),
                         is.t = T, global = F)
# remove NAs (eg. 0 variance genes) and Inf values (/0)
ng_mp_hendriks <- lapply(ng_mp_hendriks, na.omit)
ng_mp_hendriks$aRE <- ng_mp_hendriks$aRE[ng_mp_hendriks$aRE < Inf]
```

``` r
par(mfrow = c(2,2), mar=c(4,4,2,1))
invisible(sapply(names(ng_mp_hendriks), function(idx){
  rg <- range(na.omit(ng_mp_hendriks[[idx]]))
  # estimate density curve
  d <- density(na.omit(ng_mp_hendriks[[idx]]), from = rg[1], to = rg[2])
  plot(d, main = paste0(gsub("a", "", idx, fixed = T),
                        " density (", length(ng_mp_hendriks[[idx]]), " genes)"),
       xlab = idx, lwd = 2)
  # add global value
  abline(v = mp_hendriks[[idx]], lty = 2, lwd = 2, col = "firebrick")
  text(mp_hendriks[[idx]], .9*max(d$y), pos = 4, labels = idx,
       font = 2, col = "firebrick")
}))
```

Excellent fits.

Then, build the `ref` object.

```
r_dshendriks2014 <- make_ref(
  m = m_dshendriks2014,
  n.inter = 100,                  # interpolation resolution
  t.unit = "h past egg-laying, 25C",   # time unit
  metadata = list("organism" = "C. elegans", # any metadata
                  "profiling" = "bulk RNAseq"))
```

Plot model predictions on components.

``` r
par(mfrow = c(2,5), pty='s', mar = c(4,4,3,1))
plot(m_dshendriks2014, r_dshendriks2014, ncs = 1:10)
```

Finally, we stage the samples from the reference, as well as the external validation data (`dsaeschimann2017`). We exclude the earliest `dsaeschimann2017` samples (under 22h) from staging since they are outside the span covered by the reference

```
ae_dshendriks2014 <- ae(dshendriks2014$g, r_dshendriks2014)

too_young <- dsaeschimann2017$p$age <= 22
ae_dsaeschimann2017 <- ae(dsaeschimann2017$g[,!too_young], r_dshendriks2014)
```

``` r
par(mfrow = c(1,2), mar = c(4,5,3,1))
# dshendriks2014
dshendriks2014$p$age_est <- ae_dshendriks2014$age.estimates[,1]
rg <- range(c(dshendriks2014$p$age_est, dshendriks2014$p$age))
plot(age_est~age, data = dshendriks2014$p,
     xlab = "Chronological age", ylab = "Estimated age (r_dshendriks2014)",
     xlim = rg, ylim = rg,
     main = "Staging dshendriks2014
on r_dshendriks2014", lwd = 2)
points(age_est~age, data = dshendriks2014$p, type = 'l', lty = 2)
abline(a = 0, b = 1, lty = 3, lwd = 2) # x=y
legend("bottomright", legend = "x = y", lwd=3, col=1, lty = 3, bty='n')

# dsaeschimann2017
dsaeschimann2017$p[!too_young, "age_est"] <-
  ae_dsaeschimann2017$age.estimates[,1]
rg <- range(c(dsaeschimann2017$p$age_est[!too_young],
              dsaeschimann2017$p$age[!too_young]))
plot(age_est~age, data = dsaeschimann2017$p[!too_young,],
     xlab = "Chronological age", ylab = "Estimated age (r_dshendriks2014)",
     xlim = rg, ylim = rg,
     main = "Staging dsaeschimann2017
on r_dshendriks2014",
     lwd = 2, col = factor(dsaeschimann2017$p$strain))
invisible(sapply(levels(dsaeschimann2017$p$strain), function(l){
  s <- dsaeschimann2017$p$strain == l & !too_young
  points(age_est~age, data = dsaeschimann2017$p[s,],
         type = 'l', lty = 2,
         col = which(l==levels(factor(dsaeschimann2017$p$strain))))
}))
abline(a = 0, b = 1, lty = 3, lwd = 2)
legend("bottomright", bty='n',
       legend = c("let-7", "lin-41", "let-7/lin-41", "N2", "x = y"),
       lwd=3, col=c(1:4, 1), pch = c(1,1,1,1,NA), lty = c(rep(NA, 4), 3))
```

### [4.1.5] Code to generate objects

Required packages and variables:

``` r
data_folder <- "../inst/extdata/"

requireNamespace("wormRef", quietly = T)
requireNamespace("utils", quietly = T)
requireNamespace("GEOquery", quietly = T) # bioconductor
requireNamespace("Biobase", quietly = T)  # bioconductor
```

*Note : set `data_folder` to an existing path on your system where you want to store the objects.*

``` r
raw2tpm <- function(rawcounts, genelengths){
  if(nrow(rawcounts) != length(genelengths))
    stop("genelengths must match nrow(rawcounts).")
  x <- rawcounts/genelengths
  return(t( t(x) * 1e6 / colSums(x) ))
}

fpkm2tpm <- function(fpkm){
  return(exp(log(fpkm) - log(colSums(fpkm)) + log(1e6)))
}
```

To build `dshendriks2014`, *C. elegans* late larval development time series from [Hendriks et al. (2014)]

``` r
geo_dshendriks2014 <- "GSE52861"

g_url_dshendriks2014 <- GEOquery::getGEOSuppFiles(geo_dshendriks2014,
                                                  makeDirectory = FALSE,
                                                  fetch_files = FALSE)
g_file_dshendriks2014 <- paste0(data_folder, "dshendriks2014.txt.gz")
utils::download.file(url = as.character(g_url_dshendriks2014$url[2]),
                     destfile = g_file_dshendriks2014)

X_dshendriks2014 <- read.table(gzfile(g_file_dshendriks2014),
                               h=T, sep = '\t', stringsAsFactors = F,
                               row.names = 1)

# convert to tpm & wb_id
X_dshendriks2014 <- X_dshendriks2014[
  rownames(X_dshendriks2014)%in%wormRef::Cel_genes$wb_id,]
X_dshendriks2014 <- raw2tpm(
  rawcounts = X_dshendriks2014,
  genelengths = wormRef::Cel_genes$transcript_length[
    match(rownames(X_dshendriks2014), wormRef::Cel_genes$wb_id)])

# pheno data
P_dshendriks2014 <- Biobase::pData(GEOquery::getGEO(geo_dshendriks2014,
                                                    getGPL = F)[[1]])

# filter relevant fields/samples
P_dshendriks2014 <- P_dshendriks2014[
  (P_dshendriks2014$`strain:ch1` == 'N2') &
  (P_dshendriks2014$`growth protocol:ch1` == 'Continuous'), ]

P_dshendriks2014 <- P_dshendriks2014[, c("title", "geo_accession",
                                         "time in development:ch1")]

# get age from sample name
P_dshendriks2014$age <- as.numeric(
  sub('(\\d+)\\shours', '\\1', P_dshendriks2014$`time in development:ch1`))

# formatting
P_dshendriks2014$title <- gsub('RNASeq_polyA_', '',
                               gsub('hr', 'h',
                                    gsub('-', '.', fixed = T,
                                         as.character(P_dshendriks2014$title))))
colnames(X_dshendriks2014) <- gsub('RNASeq_polyA_','',
                                   colnames(X_dshendriks2014))
X_dshendriks2014 <- X_dshendriks2014[, P_dshendriks2014$title]

# save data
dshendriks2014 <- list(g = X_dshendriks2014, p = P_dshendriks2014)
save(dshendriks2014, file = paste0(data_folder, "dshendriks2014.RData"),
     compress = "xz")

# cleanup
file.remove(g_file_dshendriks2014)
rm(geo_dshendriks2014, g_url_dshendriks2014, g_file_dshendriks2014,
   X_dshendriks2014, P_dshendriks2014)
```

To build `dsaeschimann2017`, *C. elegans* larval development time series of 4 strains from [Aeschimann et al. (2017)]

``` r
geo_dsaeschimann2017 <- "GSE80157"

g_url_dsaeschimann2017 <- GEOquery::getGEOSuppFiles(geo_dsaeschimann2017,
                                                    makeDirectory = FALSE,
                                                    fetch_files = FALSE)
g_file_dsaeschimann2017 <- paste0(data_folder, "dsaeschimann2017.txt.gz")
utils::download.file(url = as.character(g_url_dsaeschimann2017$url[2]),
                     destfile = g_file_dsaeschimann2017)

X_dsaeschimann2017 <- read.table(gzfile(g_file_dsaeschimann2017),
                                 h=T, sep = '\t', stringsAsFactors = F,
                                 row.names = 1)

# convert to tpm & wb_id
X_dsaeschimann2017 <- X_dsaeschimann2017[
  rownames(X_dsaeschimann2017) %in% wormRef::Cel_genes$wb_id,]

X_dsaeschimann2017 <- raw2tpm(
  rawcounts = X_dsaeschimann2017,
  genelengths = wormRef::Cel_genes$transcript_length[
    match(rownames(X_dsaeschimann2017), wormRef::Cel_genes$wb_id)])

# pheno data
P_dsaeschimann2017 <- Biobase::pData(GEOquery::getGEO(geo_dsaeschimann2017,
                                                      getGPL = F)[[1]])
P_dsaeschimann2017[,10:34] <- NULL
P_dsaeschimann2017[, 3:8] <- NULL

colnames(P_dsaeschimann2017)[4] <- "strain"
P_dsaeschimann2017$strain <- factor(P_dsaeschimann2017$strain)
P_dsaeschimann2017$title <- make.names(P_dsaeschimann2017$title)

# get age from sample name
P_dsaeschimann2017$age <- as.numeric(
  sub('(\\d+)\\shours', '\\1', P_dsaeschimann2017$`time in development:ch1`))

# formatting
colnames(X_dsaeschimann2017) <- gsub('RNASeq_riboM_', '',
                                     colnames(X_dsaeschimann2017), fixed = T)
P_dsaeschimann2017$title <- gsub('RNASeq_riboM_', '',
                                 P_dsaeschimann2017$title, fixed = T)
X_dsaeschimann2017 <- X_dsaeschimann2017[, P_dsaeschimann2017$title]

# save data
dsaeschimann2017 <- list(g = X_dsaeschimann2017, p = P_dsaeschimann2017)

save(dsaeschimann2017, file = paste0(data_folder, "dsaeschimann2017.RData"),
     compress = "xz")

# cleanup
file.remove(g_file_dsaeschimann2017)
rm(geo_dsaeschimann2017, g_url_dsaeschimann2017, g_file_dsaeschimann2017,
   X_dsaeschimann2017, P_dsaeschimann2017)
```

\
\

## [4.2] Example 2 - *D. melanogaster* embryonic development

### [4.2.1] Data

In this example, we build a reference for *Drosophila melanogaster* embryo development. We use the following time-series:

1.  A bulk embryo development time-series, part of the modENCODE project and published by [Graveley et al. (2011)], `dsgraveley2011`, used to build the reference. (Data from [fruitfly.org](https://fruitfly.org/sequence/download.html))
2.  A high-resolution time-series of single embryos published by [Levin et al. (2016)], `dslevin2016dmel`, used for external validation. (Accession : [GSE60471](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60471))

Code to generate `dsgraveley2011` and `dslevin2016dmel` can be found [at the end of this section](#code-to-generate-objects-1)

### [4.2.2] Normalization & exploration

We start by normalizing the data.

```
library(RAPToR)
library(limma)

dsgraveley2011$g <- limma::normalizeBetweenArrays(dsgraveley2011$g,
                                                  method = "quantile")
dsgraveley2011$g <- log1p(dsgraveley2011$g)

dslevin2016dmel$g <- limma::normalizeBetweenArrays(dslevin2016dmel$g,
                                                   method = "quantile")
dslevin2016dmel$g <- log1p(dslevin2016dmel$g)
```

Check the contents of the expression matrix and pheno data.

```
dsgraveley2011$g[1:5, 1:5]
#>               em0.2hr   em2.4hr   em4.6hr   em6.8hr  em8.10hr
#> FBgn0000003 3.9688855 4.2758251 3.3179427 4.5643881 4.6995131
#> FBgn0000008 1.2963291 0.9221436 0.6952683 0.6472183 0.7443754
#> FBgn0000014 0.5099048 0.9511344 1.3960100 1.8617929 1.8433359
#> FBgn0000015 0.2435604 0.6427720 1.0519592 1.1089745 1.0195901
#> FBgn0000017 1.7970494 2.0885522 1.3390407 1.5331338 1.6763588

head(dsgraveley2011$p, n = 5)
#>      sname age
#> 1  em0.2hr   0
#> 2  em2.4hr   2
#> 3  em4.6hr   4
#> 4  em6.8hr   6
#> 5 em8.10hr   8
```

Correlation between the samples for outliers.

``` r
cor_dsgraveley2011 <- cor(dsgraveley2011$g, method = "spearman")
ord <- order(dsgraveley2011$p$age)

heatmap(cor_dsgraveley2011[ord, ord], Colv = NA, Rowv = NA,
        scale = "none", keep.dendro = F, margins = c(2,5),
        labRow = "", labCol = "")
# add time label
par(xpd = T)
mtext(text = dsgraveley2011$p$age, side = 1,
      line = 4, at = seq(-.16,.85, l = 12))
mtext(at = 1, line = 4, side = 1, text = "(hours)", cex = .8)

# add color legend
l.values <- seq(min(cor_dsgraveley2011), max(cor_dsgraveley2011), l = 9)
image(x = c(.95,1), y = seq(0.6,1, l = 9), useRaster = T,
      z = matrix(l.values, ncol = 9),
      col = hcl.colors(12, "YlOrRd", rev = TRUE), add = T)
text(.975, 1, pos = 3, labels = expression(rho), font = 2)
text(1, y = seq(0.6,1, l = 9)[c(T,F)], pos = 4,
     labels = round(l.values[c(T,F)], 2), cex = .6)
```

``` r
boxplot(cor_dsgraveley2011, names = paste0(dsgraveley2011$p$age,'h'),
        boxwex=.5, las=1)
```

Plotting principal components.

```
pca_dsgraveley2011 <- summary(stats::prcomp(t(dsgraveley2011$g), rank = 12,
                                            center = TRUE, scale = FALSE))
```

``` r
par(mfrow = c(2,4), pty = 's')
invisible(sapply(seq_len(8), function(i){
  plot(dsgraveley2011$p$age, pca_dsgraveley2011$x[,i], type = 'b', lwd = 2,
       xlab = "Chronological age", ylab = "PC", main = paste0("PC", i))
}))
```

### [4.2.3] Model fitting

We keep enough components to explain \\(99\\%\\) of the variance in the data.

```
nc <- sum(pca_dsgraveley2011$importance[3,] < .99) + 1
nc
#> [1] 8
```

We then fit a GAM on PCA components with a cubic spline on age.

```
m_dsgraveley2011 <- ge_im(X = dsgraveley2011$g,
                          p = dsgraveley2011$p,
                          formula = "X ~ s(age, bs = 'cr')", nc = nc)
```

### [4.2.4] Validation

First, check global model performance.

```
# global model performance
mp_grav <- mperf(dsgraveley2011$g, predict(m_dsgraveley2011), is.t = T)
print(do.call(cbind, mp_grav))
#>            aCC       aRE         MSE      aRMSE
#> [1,] 0.9400274 0.4040618 0.009908293 0.09954041
```

And then per gene.

```
ng_mp_grav <-  mperf(dsgraveley2011$g, predict(m_dsgraveley2011),
                     is.t = T, global = F)
# remove NAs (eg. 0 variance genes) and Inf values (/0)
ng_mp_grav <- lapply(ng_mp_grav, na.omit)
ng_mp_grav$aRE <- ng_mp_grav$aRE[ng_mp_grav$aRE < Inf]
```

``` r
par(mfrow = c(2,2), mar=c(4,4,2,1))
invisible(sapply(names(ng_mp_grav), function(idx){
  rg <- range(na.omit(ng_mp_grav[[idx]]))
  # estimate density curve
  d <- density(na.omit(ng_mp_grav[[idx]]), from = rg[1], to = rg[2])
  plot(d, main = paste0(gsub("a", "", idx, fixed = T),
                        " density (", length(ng_mp_grav[[idx]]), " genes)"),
       xlab = idx, lwd = 2)
  # add global value
  abline(v = mp_grav[[idx]], lty = 2, lwd = 2, col = "firebrick")
  text(mp_grav[[idx]], .9*max(d$y), pos = 4, labels = idx,
       font = 2, col = "firebrick")
}))
```

Then, build the `ref` object.

```
r_dsgraveley2011 <- make_ref(
  m = m_dsgraveley2011,
  n.inter = 100,                # interpolation resolution
  t.unit = "h past egg-laying", # time unit
  metadata = list("organism" = "D. melanogaster", # any metadata
                  "profiling" = "bulk RNAseq"))
```

Plot model predictions on components.

``` r
par(mfrow = c(2,4), pty='s', mar = c(4,4,3,1))
plot(m_dsgraveley2011, r_dsgraveley2011, col.i = 'royalblue')
```

No overfitting, good fits on the components explaining most of the variance. Some components flattened, but minimal associated variance.

Finally, we stage the samples from the reference, as well as the external validation data (`dslevin2016dmel`).

```
ae_dsgraveley2011 <- ae(dsgraveley2011$g, r_dsgraveley2011)
ae_dslevin2016dmel <- ae(dslevin2016dmel$g, r_dsgraveley2011)
```

``` r
par(mfrow = c(1,2), mar = c(4,5,3,1))
# dsgraveley2011
dsgraveley2011$p$age_est <- ae_dsgraveley2011$age.estimates[,1]
rg <- range(c(dsgraveley2011$p$age_est, dsgraveley2011$p$age))
plot(age_est~age, data = dsgraveley2011$p,
     xlab = "Chronological age", ylab = "Estimated age (r_dsgraveley2011)",
     xlim = rg, ylim = rg, lwd = 2,
     main = "Staging dsgraveley2011
on r_dsgraveley2011")
points(age_est~age, data = dsgraveley2011$p, type = 'l', lty = 2)
abline(a = 0, b = 1, lty = 3, lwd = 2) # x=y
legend("bottomright", legend = "x = y", lwd=3, col=1, lty = 3, bty='n')

# dslevin2016dmel
dslevin2016dmel$p$age_est <- ae_dslevin2016dmel$age.estimates[,1]
rg <- range(c(dslevin2016dmel$p$age_est , dslevin2016dmel$p$age))
plot(age_est~age, data = dslevin2016dmel$p,
     xlab = "Chronological age", ylab = "Estimated age (r_dsgraveley2011)",
     xlim = rg, ylim = rg, lwd = 2,
     main = "Staging dslevin2016dmel
on r_dsgraveley2011")
abline(a = 0, b = 1, lty = 3, lwd = 2) #x=y
legend("bottomright", legend = "x = y", lwd=3, col=1, lty = 3, bty='n')
```

We see that the validation data age estimates vary from their chronological age, particularly in later stages. However, this is not due to errors in the reference or staging, but instead to the fact this data is single-embryo profiling. Indeed, inter-individual variation in developmental speed results in striking differences with expected age, where it would otherwise be averaged out in bulk data.

The dynamics of the `dslevin2016dmel` data further confirm this.

```
pca_dslevin2016dmel <- stats::prcomp(t(dslevin2016dmel$g), rank = 20,
                                     center = TRUE, scale = FALSE)
```

``` r
par(mfrow = c(1,5), pty='s')
invisible(sapply(1:5, function(i){
    plot(dslevin2016dmel$p$age, pca_dslevin2016dmel$x[,i], lwd = 2,
         xlab = "Chronological age", ylab = "PC", main = paste0("PC", i))
}))
```

Inferred age however, clearly restores the dynamics cleanly.

``` r
par(mfrow = c(1,5), pty='s')
invisible(sapply(1:5, function(i){
    plot(dslevin2016dmel$p$age_est, pca_dslevin2016dmel$x[,i], lwd = 2,
         xlab = "Estimated age", ylab = "PC", main = paste0("PC", i),
         col = 'firebrick')
  box(col = 'firebrick', lwd=2)
}))
```

This example shows how inter-individual variation in developmental speed makes it difficult to profile high-temporal-resolution time series. Furthermore, the temporal resolution difference between the reference and validation data also demonstrates the effectiveness of gene expression interpolation.

### [4.2.5] Code to generate objects

Required packages and variables:

``` r
data_folder <- "../inst/extdata/"

requireNamespace("utils", quietly = T)
requireNamespace("GEOquery", quietly = T) # bioconductor
requireNamespace("Biobase", quietly = T)  # bioconductor
requireNamespace("biomaRt", quietly = T)  # bioconductor
```

*Note : set `data_folder` to an existing path on your system where you want to store the objects.*

``` r
raw2tpm <- function(rawcounts, genelengths){
  if(nrow(rawcounts) != length(genelengths))
    stop("genelengths must match nrow(rawcounts).")
  x <- rawcounts/genelengths
  return(t( t(x) * 1e6 / colSums(x) ))
}

fpkm2tpm <- function(fpkm){
  return(exp(log(fpkm) - log(colSums(fpkm)) + log(1e6)))
}
```

Get *Drosophila* gene ID table from ensembl.

``` r
mart <- biomaRt::useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
droso_genes <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                             "ensembl_transcript_id",
                                             "external_gene_name",
                                             "transcript_end",
                                             "transcript_start"),
                              mart = mart)
droso_genes$transcript_length <-
  droso_genes$transcript_end - droso_genes$transcript_start
droso_genes <- droso_genes[,c(1:3,6)]
colnames(droso_genes)[1:3] <- c("fb_id", "transcript_id", "gene_name")

rm(mart)
```

To build `dsgraveley2011`, bulk *D. melanogaster* embryo time series from [Graveley et al. (2011)]

``` r
g_url_dsgraveley2011 <- paste0(
  "ftp://ftp.fruitfly.org/pub/download/modencode_expression_scores/",
  "Celniker_Drosophila_Annotation_20120616_1428_allsamps_",
  "MEAN_gene_expression.csv.gz")
g_file_dsgraveley2011 <- paste0(data_folder, "dsgraveley2011.csv.gz")
utils::download.file(g_url_dsgraveley2011, destfile = g_file_dsgraveley2011)

X_dsgraveley2011 <- read.table(gzfile(g_file_dsgraveley2011),
                               sep = ',', row.names = 1, h = T)

# convert gene ids to FBgn
X_dsgraveley2011 <- RAPToR::format_ids(X_dsgraveley2011, droso_genes,
                                       from = "gene_name", to = "fb_id")

# select embryo time series samples
X_dsgraveley2011 <- X_dsgraveley2011[,1:12]

P_dsgraveley2011 <- data.frame(
  sname = colnames(X_dsgraveley2011),
  age = as.numeric(gsub("em(\\d+)\\.\\d+hr", "\\1",
                        colnames(X_dsgraveley2011))),
  stringsAsFactors = FALSE)

dsgraveley2011 <- list(g = X_dsgraveley2011, p = P_dsgraveley2011)

save(dsgraveley2011,
     file = paste0(data_folder, "dsgraveley2011.RData"), compress = "xz")

# cleanup
file.remove(g_file_dsgraveley2011)
rm(g_url_dsgraveley2011, g_file_dsgraveley2011,
   X_dsgraveley2011, P_dsgraveley2011)
```

To build `dslevin2016dmel`, single-embryo *D. melanogaster* embryo time series from [Levin et al. (2016)]

``` r
geo_dslevin2016dmel <- "GSE60471"
g_url_dslevin2016dmel <- GEOquery::getGEOSuppFiles(geo_dslevin2016dmel,
                                                   makeDirectory = FALSE,
                                                   fetch_files = FALSE)
g_file_dslevin2016dmel <- paste0(data_folder, "dslevin2016dmel.txt.gz")
utils::download.file(url = as.character(g_url_dslevin2016dmel$url[3]),
                     destfile = g_file_dslevin2016dmel)

X_dslevin2016dmel <- read.table(gzfile(g_file_dslevin2016dmel), h = T,
                                sep = '\t', as.is = T, row.names = 1,
                                comment.char = "")

# filter poor quality samples
cm_dslevin2016dmel <- RAPToR::cor.gene_expr(X_dslevin2016dmel,
                                            X_dslevin2016dmel)
f_dslevin2016dmel <- which(0.6 > apply(cm_dslevin2016dmel, 1,
                                       quantile, probs = .99))
X_dslevin2016dmel <- X_dslevin2016dmel[, -f_dslevin2016dmel]

# convert to tpm & FBgn
X_dslevin2016dmel <- X_dslevin2016dmel[
  rownames(X_dslevin2016dmel)%in%droso_genes$fb_id,]
X_dslevin2016dmel <- raw2tpm(
  rawcounts = X_dslevin2016dmel,
  genelengths = droso_genes$transcript_length[
    match(rownames(X_dslevin2016dmel), droso_genes$fb_id)])

# pheno data
P_dslevin2016dmel <- Biobase::pData(GEOquery::getGEO(geo_dslevin2016dmel,
                                                     getGPL = F)[[1]])

# filter relevant fields/samples
P_dslevin2016dmel <- P_dslevin2016dmel[,
  c("title", "geo_accession", "time (minutes cellularization stage):ch1")]
colnames(P_dslevin2016dmel)[3] <- "time"
P_dslevin2016dmel$title <- as.character(P_dslevin2016dmel$title)

P_dslevin2016dmel <- P_dslevin2016dmel[
  P_dslevin2016dmel$title %in% colnames(X_dslevin2016dmel),]
X_dslevin2016dmel <- X_dslevin2016dmel[, P_dslevin2016dmel$title]

# formatting
P_dslevin2016dmel$title <- gsub('Metazome_Drosophila_timecourse_', '',
                                P_dslevin2016dmel$title)
colnames(X_dslevin2016dmel) <- P_dslevin2016dmel$title
P_dslevin2016dmel$age <- as.numeric(P_dslevin2016dmel$time) / 60

# save data
dslevin2016dmel <- list(g = X_dslevin2016dmel, p = P_dslevin2016dmel)
save(dslevin2016dmel, file = paste0(data_folder, "dslevin2016dmel.RData"), compress = "xz")

# cleanup
file.remove(g_file_dslevin2016dmel)
rm(geo_dslevin2016dmel, g_url_dslevin2016dmel, g_file_dslevin2016dmel,
   X_dslevin2016dmel, P_dslevin2016dmel,
   cm_dslevin2016dmel, f_dslevin2016dmel)
```

\
\

## [4.3] Example 3 - Uneven sampling of *D. rerio* embryos

### [4.3.1] Data

This example uses two *Danio rerio* (zebrafish) embryo development time-series:

1.  A time-series of zebrafish embryonic development (with uneven time sampling) published by [White et al. (2017)], `dswhite2017`, used to build the reference. ([Data accessible from the publication](https://elifesciences.org/articles/30860))
2.  A high-resolution time-series of embryonic development published by [Levin et al. (2016)], `dslevin2016zeb`, used for external validation. (Accession : [GSE60619](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60619))

Code to generate `dswhite2017` and `dslevin2016zeb` can be found [at the end of this section](#code-to-generate-objects-2)

The reference data (`dswhite2017`) has uneven time sampling, as can often be the case to account for differences in dynamic ranges of expression: later time points are more sparse. We can apply a transformation on the time values (using `asinh()` in this case) so they become more uniform and thus avoid interpolation bias.

### [4.3.2] Normalization & Quick look

We start by normalizing the data.

```
library(RAPToR)
library(limma)

dswhite2017$g <- limma::normalizeBetweenArrays(dswhite2017$g,
                                               method = "quantile")
dswhite2017$g <- log1p(dswhite2017$g)

dslevin2016zeb$g <- limma::normalizeBetweenArrays(dslevin2016zeb$g,
                                                  method = "quantile")
dslevin2016zeb$g <- log1p(dslevin2016zeb$g)
```

Check the contents of the expression matrix and pheno data.

```
dswhite2017$g[1:5, 1:4]
#>                    zmp_ph133_B zmp_ph133_D zmp_ph133_E zmp_ph133_F
#> ENSDARG00000000001    2.192007    2.019082    1.929426    2.031762
#> ENSDARG00000000002    1.149510    1.188959    0.900076    1.185358
#> ENSDARG00000000018    2.456661    2.534134    2.224970    2.364784
#> ENSDARG00000000019    4.432509    4.529970    4.608232    4.533400
#> ENSDARG00000000068    4.406696    4.460862    4.267657    4.294028

head(dswhite2017$p, n = 5)
#>        sample accession_number         stage stageName sampleName age
#> 1 zmp_ph133_B       ERS1079239 Zygote:1-cell    1-cell   1-cell-1   0
#> 2 zmp_ph133_D       ERS1079240 Zygote:1-cell    1-cell   1-cell-2   0
#> 3 zmp_ph133_E       ERS1079241 Zygote:1-cell    1-cell   1-cell-3   0
#> 4 zmp_ph133_F       ERS1079243 Zygote:1-cell    1-cell   1-cell-4   0
#> 5 zmp_ph133_G       ERS1079244 Zygote:1-cell    1-cell   1-cell-5   0
#>   batch
#> 1     1
#> 2     2
#> 3     3
#> 4     4
#> 5     5
```

Correlation between samples for outliers.

``` r
cor_dswhite2017 <- cor(dswhite2017$g, method = "spearman")
ord <- order(dswhite2017$p$age)
heatmap(cor_dswhite2017[ord, ord], Colv = NA, Rowv = NA,
        scale = "none", keep.dendro = F, margins = c(1,5),
        RowSideColors = transp(as.numeric(dswhite2017$p$batch[ord])),
        labRow = "", labCol = "")
par(xpd = T)
mtext(text = unique(dswhite2017$p$age), side = 1, line = c(3.8, 4),
      at = seq(-.12, .915, l = length(unique(dswhite2017$p$age))), cex = .6)

# color key
l.values <- seq(min(cor_dswhite2017), max(cor_dswhite2017), l = 10)
image(x = c(.95,1), y = seq(0.6,1, l = 10), useRaster = T,
      z = matrix(l.values, ncol = 10),
      col = hcl.colors(12, "YlOrRd", rev = TRUE), add = T)

text(.975, 1, pos = 3, labels = expression(rho), font = 2)
text(1, y = seq(0.6,1, l = 10), pos = 4,
     labels = round(l.values, 2), cex = .6)

xlp <- 1.025
batch_legend <- as.character(1:5)
text(xlp, .5, labels = "batch", font = 2, cex = .8, adj = .5)
text(xlp, seq(.3,.48, l = 5), labels = batch_legend, adj = 1, pos = 1,
     col = levels(dswhite2017$p$batch), font = 2, cex = .7)

mtext(at = xlp, line = 4, side = 1, text = "(hours)", cex = .8)
```

``` r
boxplot(cor_dswhite2017~interaction(dswhite2017$p$batch,
                                    dswhite2017$p$age),
        col = transp(1:5), # see 'Code to generate objects' for transp()
        xaxt = "n", at = seq(1,90, l = 90*(6/5))[c(T,T,T,T,T,F)], boxwex=.5,
        ylab = "Spearman correlation", xlab = "Chronological age (h)", las=1)
#add time labels
axis(side = 1, at = seq(2.5,87.5, l = 90/5),
     labels = paste0(unique(dswhite2017$p$age), 'h'))
legend('bottom', fill = transp(1:5), title = "Batch",
       legend = c(1:5), horiz = T, bty = "n")
```

No outliers.

Plotting principal components.

```
pca_dswhite2017 <- stats::prcomp(t(dswhite2017$g), rank = 25,
                                 center = TRUE, scale = FALSE)
```

``` r
par(mfrow = c(2,4), pty='s')
invisible(sapply(1:8, function(i){
  plot(dswhite2017$p$age, pca_dswhite2017$x[,i],
       lwd = 2, col = dswhite2017$p$batch,
       xlab = "Chronological age", ylab = "PC", main = paste0("PC", i))
  # add dotted lines
  sapply(seq_along(levels(dswhite2017$p$batch)), function(l){
    s <- which(dswhite2017$p$batch == levels(dswhite2017$p$batch)[l])
    points(dswhite2017$p$age[s], pca_dswhite2017$x[s,i], col = l,
           type = 'l', lty = 2)
  })
  # legend
  if(i == 1)
    legend("bottomright", bty = 'n', legend = batch_legend, title = "Batch",
           pch = c(rep(1, 5)), lty = c(rep(NA, 5)), col = c(1:5), lwd = 3)
}))
```

Sampling is sparser towards the end of the time series, and expression dynamics are also "wider". Because of this, fitting splines along age on these components will poorly fit the earlier time points.

To bypass this, we can use `asinh(age)`, which has a similar effect to a logarithm to stretch the age values closer to a uniform scale. We also add an intercept to avoid overstretching the first few time points. The relationship between `age` and `asinh(1+age)` is shown below.

```
plot(asinh(1+age)~age, data=dswhite2017$p)
# add a grid
abline(v=dswhite2017$p$age, h=asinh(1+dswhite2017$p$age),
       col = 'grey80', lty=3)
```

The transformation will be specified in the model formula so that we can predict new points from time values on the `age` scale (otherwise, we would need to predict with using transformed values such that they correspond to a uniform time scale when transformed back).

### [4.3.3] Model fitting

We keep enough components to explain \\(99\\%\\) of the variance in the data.

```
nc <- sum(summary(pca_dswhite2017)$importance[3,] < .99) + 1
nc
#> [1] 67
```

We will fit a GAM on PCA components. The type of spline to fit as well as the value of the intercept in `asinh()` will be determined using cross-validation. Since sample batch is clearly absent from the components, we exclude it for a more parsimonious model.

```
set.seed(3)
# intercept values to try
intercepts <- c(0,1,2,5)

# list of formulas to test
flist <- as.list(paste0(
  "X ~ s(",
  c("age", # age without transformation
    paste0("asinh(", intercepts, "+age)")), # asinh with intercept
  ", bs = '",
  rep(c("tp", "gp"), e=1+length(intercepts)), # 2 different splines
  "', k=9)"))
# print formula list
cat(unlist(flist), sep = '\n')
#> X ~ s(age, bs = 'tp', k=9)
#> X ~ s(asinh(0+age), bs = 'tp', k=9)
#> X ~ s(asinh(1+age), bs = 'tp', k=9)
#> X ~ s(asinh(2+age), bs = 'tp', k=9)
#> X ~ s(asinh(5+age), bs = 'tp', k=9)
#> X ~ s(age, bs = 'gp', k=9)
#> X ~ s(asinh(0+age), bs = 'gp', k=9)
#> X ~ s(asinh(1+age), bs = 'gp', k=9)
#> X ~ s(asinh(2+age), bs = 'gp', k=9)
#> X ~ s(asinh(5+age), bs = 'gp', k=9)

cv_dswhite2017 <- ge_imCV(X = dswhite2017$g, p = dswhite2017$p,
                          formula_list = flist, nc = nc,
                          cv.n = 10, nb.cores = 6)
#> CV on 10 models. cv.n = 10 | cv.s = 0.8
#>
#> ...Building training sets
#> ...Setting up cluster
#> ...Running CV
#> ...Cleanup and formatting
```

``` r
par(mar = c(7,4,3,1))
cpal <- c("grey30", "#420A68FF",
          "#932667FF", "#DD513AFF", "#FCA50AFF") # color palette
plot(cv_dswhite2017, to_plot = c("aRE", "MSE"),
     names = paste0("s:", rep(c("tp", "gp"), e=1+length(intercepts)),
                    "/", c("nT", paste0("i=", intercepts))),
     col = NA, lwd=2, names.arrange = 5, boxwex = .5, outline=F,
     border = cpal, tcol = cpal, swarmargs = list(cex = .5, col = cpal))
```

We first see that with no `asinh()` transformation on time (grey, `nT` boxes), model performance and CV error are worse than on transformed data. Next, we note that than thin plate splines (`tp`) perform slightly better than Gaussian process splines (`gp`) from MSE. Finally, it seems that an intercept of 1 minimizes the overall CV error and model MSE.

We select `s:tp/i=1`, a thin plate spline with transformed age and intercept of 1.

```
m_dswhite2017 <- ge_im(X = dswhite2017$g,
                       p = dswhite2017$p,
                       formula = "X ~ s(asinh(1+age), bs = 'gp', k=9)",
                       nc = nc)
```

### [4.3.4] Validation

First, check global model performance.

```
# global model performance
mp_white <- mperf(dswhite2017$g, predict(m_dswhite2017), is.t = T)
print(do.call(cbind, mp_white))
#>            aCC      aRE       MSE     aRMSE
#> [1,] 0.8299206 1.090387 0.0515087 0.2269553
```

And then per gene.

```
ng_mp_white <-  mperf(dswhite2017$g, predict(m_dswhite2017),
                     is.t = T, global = F)
# remove NAs (eg. 0 variance genes) and Inf values (/0)
ng_mp_white <- lapply(ng_mp_white, na.omit)
ng_mp_white$aRE <- ng_mp_white$aRE[ng_mp_white$aRE < Inf]
```

``` r
par(mfrow = c(2,2), mar=c(4,4,2,1))
invisible(sapply(names(ng_mp_white), function(idx){
  rg <- range(na.omit(ng_mp_white[[idx]]))
  # estimate density curve
  d <- density(na.omit(ng_mp_white[[idx]]), from = rg[1], to = rg[2])
  plot(d, main = paste0(gsub("a", "", idx, fixed = T),
                        " density (", length(ng_mp_white[[idx]]), " genes)"),
       xlab = idx, lwd = 2)
  # add global value
  abline(v = mp_white[[idx]], lty = 2, lwd = 2, col = "firebrick")
  text(mp_white[[idx]], .9*max(d$y), pos = 4, labels = idx,
       font = 2, col = "firebrick")
}))
```

Then, build the `ref` object.

```
# make a 'ref' object'
r_dswhite2017 <- make_ref(
  m_dswhite2017,
  n.inter = 500,
  t.unit = "h post-fertilization",     # time unit
  metadata = list("organism" = "D. rerio", # any metadata
                  "profiling" = "single-embryo RNAseq",
                  "note" = "asinh(1+age) interpolation"))
```

Plot component interpolation.

``` r
par(mfrow = c(1,4), pty='s')
plot(m_dswhite2017, r_dswhite2017, ncs=1:4, l.pos = 'bottomright')
```

We can zoom in on the early time points to ensure they are indeed properly fit.

``` r
par(mfrow = c(1,4), pty='s')
plot(m_dswhite2017, r_dswhite2017, ncs=1:4, l.pos = 'bottomright', xlim = c(0,40))
```

```
# stage samples
ae_dswhite2017 <- ae(dswhite2017$g, r_dswhite2017)
ae_dslevin2016zeb <- ae(dslevin2016zeb$g, r_dswhite2017)
#>                 nb.genes
#> refdata            27754
#> samp               21429
#> intersect.genes    21429
```

``` r
par(mfrow = c(1,2), mar = c(4,5,3,1), pty='s')
rg <- range(c(ae_dswhite2017$age.estimates[,1], dswhite2017$p$age))
plot(ae_dswhite2017$age.estimates[,1]~dswhite2017$p$age,
     xlab = "Chronological age", ylab = "Estimated age (dswhite2017)",
     xlim = rg, ylim = rg,
     main = "Chron. vs Estimated ages for dswhite2017
(on r_dswhite2017)",
     lwd = 2)
abline(a = 0, b = 1, lty = 3, lwd = 1)

plot(ae_dslevin2016zeb$age.estimates[,1]~dslevin2016zeb$p$age,
     xlab = "Chronological age", ylab = "Estimated age (dswhite2017)",
     xlim = rg, ylim = rg,
     main = "Chron. vs Estimated ages for dslevin2016zeb
(on r_dswhite2017)", lwd = 2)
abline(a = 0, b = 1, lty = 3, lwd = 2)

legend("bottomright", legend = "x = y", lwd=3, col=1, lty = 3, bty='n')
```

### [4.3.5] Model fitting without age transformation

We'll build the same model without considering uneven sampling, for comparison. This is purposely not the optimal model and will show interpolation issues to look out for.

```
m_dswhite2017_nT <- ge_im(X = dswhite2017$g,
                          p = dswhite2017$p,
                          formula = "X ~ s(age, bs = 'tp', k=9)",
                          nc = nc)
```

    #>            aCC       aRE      MSE    aRMSE
    #> [1,] 0.7776173 0.5979607 2.672777 1.634863

Build a `ref` object

```
# make a 'ref' object'
r_dswhite2017_nT <- make_ref(
  m_dswhite2017_nT,
  n.inter = 500,
  t.unit = "h post-fertilization",     # time unit
  metadata = list("organism" = "D. rerio", # any metadata
                  "profiling" = "single-embryo RNAseq",
                  "note" = "Not the best model!"))
```

Plot component interpolation.

``` r
par(mfrow = c(1,4), pty='s')
plot(m_dswhite2017_nT, r_dswhite2017_nT, ncs=1:4, l.pos = 'bottomright')
```

``` r
par(mfrow = c(1,2), pty='s')
plot(m_dswhite2017_nT, r_dswhite2017_nT, ncs=3:4, show.legend = F,
     xlim = c(0,50))
```

Consequences can be seen on when staging external data (*i.e.* the validation data), *but not necessarily when staging the reference data* as we'll see below. Validating references with independent/external data is thus highly recommended.

```
ae_dswhite2017_nT <- ae(dswhite2017$g, r_dswhite2017_nT)
ae_dslevin2016zeb_nT <- ae(dslevin2016zeb$g, r_dswhite2017_nT)
```

``` r
par(mfrow = c(1,2), mar = c(4,5,3,1), pty='s')
rg <- range(c(ae_dswhite2017_nT$age.estimates[,1], dswhite2017$p$age))
plot(ae_dswhite2017_nT$age.estimates[,1]~dswhite2017$p$age,
     xlab = "Chronological age", ylab = "Estimated age (r_dswhite2017_nT)",
     xlim = rg, ylim = rg,
     main = "Staging dswhite2017
(on r_dswhite2017_nT)",
     lwd = 2)
abline(a = 0, b = 1, lty = 3, lwd = 1)

plot(ae_dslevin2016zeb_nT$age.estimates[,1]~dslevin2016zeb$p$age,
     xlab = "Chronological age", ylab = "Estimated age (r_dswhite2017_nT)",
     xlim = rg, ylim = rg,
     main = "Staging dslevin2016zeb
(on r_dswhite2017_nT)", lwd = 2)
abline(a = 0, b = 1, lty = 3, lwd = 2)

legend("bottomright", legend = "x = y", lwd=3, col=1, lty = 3, bty='n')
```

The "gaps" or "steps" in the validation data estimates are caused by interpolation bias (overfitting and poorly fit dynamics mentioned above). These model errors create local "unlikely/unrealistic" gene expression areas in the interpolation, which do not correlate with the samples of corresponding age.

Gaps will most often be in-between original time points of the reference data, meaning the effect can't be seen by staging reference samples. Validation data of sufficient temporal resolution will however have clear 'blank' ranges of age estimates, around which are clustered the samples of corresponding age.

The sub-optimal model we used had clear red flags on component plots, but interpolation bias can be much more subtle. In such cases, only staging external data will reveal the problem.

### [4.3.6] Code to generate objects

Required packages and variables:

``` r
data_folder <- "../inst/extdata/"

requireNamespace("utils", quietly = T)
requireNamespace("GEOquery", quietly = T) # bioconductor
requireNamespace("Biobase", quietly = T)  # bioconductor
requireNamespace("biomaRt", quietly = T)  # bioconductor
```

*Note : set `data_folder` to an existing path on your system where you want to store the objects.*

``` r
raw2tpm <- function(rawcounts, genelengths){
  if(nrow(rawcounts) != length(genelengths))
    stop("genelengths must match nrow(rawcounts).")
  x <- rawcounts/genelengths
  return(t( t(x) * 1e6 / colSums(x) ))
}

fpkm2tpm <- function(fpkm){
  return(exp(log(fpkm) - log(colSums(fpkm)) + log(1e6)))
}
```

``` r
mart <- biomaRt::useMart("ensembl", dataset = "drerio_gene_ensembl")
zeb_genes <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                           "transcript_length"),
                            mart = mart)
rm(mart)
```

To build `dswhite2017`, *D. rerio* embryo time series from [White et al. (2017)]

``` r
p_url_dswhite2017 <- paste0("http://europepmc.org/articles/PMC5690287/",
                            "bin/elife-30860-supp1.tsv")
g_url_dswhite2017 <- apste0("http://europepmc.org/articles/PMC5690287/",
                            "bin/elife-30860-supp2.tsv")
g_file_dswhite2017 <- paste0(data_folder, "dswhite2017.tsv")
utils::download.file(g_url_dswhite2017, destfile = g_file_dswhite2017)

X_dswhite2017 <- read.table(g_file_dswhite2017, h = T, sep  ="\t",
                            as.is = T, quote = "\"")
rownames(X_dswhite2017) <- X_dswhite2017$Gene.ID
X_dswhite2017 <- X_dswhite2017[,-(1:8)]

# convert to tpm & ensembl_id
X_dswhite2017 <- X_dswhite2017[
  rownames(X_dswhite2017)%in%zeb_genes$ensembl_gene_id,]
X_dswhite2017 <- raw2tpm(
  rawcounts = X_dswhite2017,
  genelengths = zeb_genes$transcript_length[
    match(rownames(X_dswhite2017), zeb_genes$ensembl_gene_id)])

# pheno data
P_dswhite2017 <- read.table(p_url_dswhite2017, h = T, sep = "\t", as.is = T)
P_dswhite2017 <- P_dswhite2017[P_dswhite2017$sequencing == "RNASeq",
                               c("sample", "accession_number", "stage",
                                 "stageName", "sampleName")]

# timings of stages (from White et al. eLife (2017)).
# in hours post-fertilization
timepoints <- data.frame(stage = unique(P_dswhite2017$stageName),
                         hours_pf = c(0, .75, 2.25, 3, 4.3, 5.25, 6, 8, 10.3,
                                      16, 19, 24, 30, 36, 48, 72, 96, 120),
                         stringsAsFactors = F, row.names = "stage")
P_dswhite2017$age <- timepoints[P_dswhite2017$stageName, "hours_pf"]

# formatting
P_dswhite2017$batch <- factor(gsub(".*-(\\d)$", "\\1",
                                   P_dswhite2017$sampleName))
X_dswhite2017 <- X_dswhite2017[, P_dswhite2017$sample]

# save data
dswhite2017 <- list(g = X_dswhite2017, p = P_dswhite2017)
save(dswhite2017, file = paste0(data_folder, "dswhite2017.RData"),
     compress = "xz")

# cleanup
file.remove(g_file_dswhite2017)
rm(p_url_dswhite2017, g_url_dswhite2017, g_file_dswhite2017,
   X_dswhite2017, P_dswhite2017, timepoints)
```

To build `dslevin2016zeb`, *D. rerio* embryo time series from [Levin et al. (2016)]

``` r
geo_dslevin2016zeb <- "GSE60619"

g_url_dslevin2016zeb <- GEOquery::getGEOSuppFiles(geo_dslevin2016zeb,
                                                  makeDirectory = FALSE,
                                                  fetch_files = FALSE)
g_file_dslevin2016zeb <- paste0(data_folder, "dslevin2016zeb.txt.gz")
utils::download.file(url = as.character(g_url_dslevin2016zeb$url[2]),
                     destfile = g_file_dslevin2016zeb)

X_dslevin2016zeb <- read.table(gzfile(g_file_dslevin2016zeb), h = T,
                               sep = '\t', as.is = T, row.names = 1,
                               comment.char = "")

# convert to tpm & ensembl_id
X_dslevin2016zeb <- X_dslevin2016zeb[
  rownames(X_dslevin2016zeb)%in%zeb_genes$ensembl_gene_id,]
X_dslevin2016zeb <- raw2tpm(
  rawcounts = X_dslevin2016zeb,
  genelengths = zeb_genes$transcript_length[
    match(rownames(X_dslevin2016zeb), zeb_genes$ensembl_gene_id)])

# pheno data
P_dslevin2016zeb <- Biobase::pData(GEOquery::getGEO(geo_dslevin2016zeb,
                                                    getGPL = F)[[1]])

# filter relevant fields/samples
P_dslevin2016zeb <- P_dslevin2016zeb[, c("title", "geo_accession",
                                         "time (min after fertilization):ch1")]
colnames(P_dslevin2016zeb)[3] <- "time"
P_dslevin2016zeb$title <- as.character(P_dslevin2016zeb$title)

P_dslevin2016zeb <- P_dslevin2016zeb[
  P_dslevin2016zeb$title %in% colnames(X_dslevin2016zeb),]
X_dslevin2016zeb <- X_dslevin2016zeb[, P_dslevin2016zeb$title]

# formatting
P_dslevin2016zeb$title <- gsub('Metazome_ZF_timecourse_', '',
                               P_dslevin2016zeb$title)
colnames(X_dslevin2016zeb) <- P_dslevin2016zeb$title

P_dslevin2016zeb$age <- as.numeric(P_dslevin2016zeb$time) / 60

# save data
dslevin2016zeb <- list(g = X_dslevin2016zeb, p = P_dslevin2016zeb)
save(dslevin2016zeb, file = paste0(data_folder, "dslevin2016zeb.RData"),
     compress = "xz")

# cleanup
file.remove(g_file_dslevin2016zeb)
rm(geo_dslevin2016zeb, g_url_dslevin2016zeb, g_file_dslevin2016zeb,
   X_dslevin2016zeb, P_dslevin2016zeb)
```

Function to make a color transparent (`transp()`)

``` r
transp <- function(col, a=.5){
  colr <- col2rgb(col)
  return(rgb(colr[1,], colr[2,], colr[3,], a*255, maxColorValue = 255))
}
```

\
\

## [4.4] Example 4 - *C. elegans* aging

### [4.4.1] Data

Two *C. elegans* aging time-series are used in this example:

1.  An unpublished time-series of an RNAi-hypersensitive strain produced by Byrne *et al.*, `dsbyrne2020`, used to build the reference. (Accession : [GSE93826](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93826))
2.  A time-series of aging in 3 diet conditions published [Hou et al. (2016)], `dshou2016`, used for validation. (Accession : [GSE77110](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77110))

Code to generate `dsbyrne2020` and `dshou2016` can be found [at the end of this section](#code-to-generate-objects-3)

### [4.4.2] Normalization & Quick look

We start by normalizing the data.

```
library(RAPToR)
library(limma)

dsbyrne2020$g <- limma::normalizeBetweenArrays(dsbyrne2020$g,
                                               method = "quantile")
dsbyrne2020$g <- log1p(dsbyrne2020$g)

dshou2016$g <- limma::normalizeBetweenArrays(dshou2016$g,
                                             method = "quantile")
dshou2016$g <- log1p(dshou2016$g)
```

Check the contents of the expression matrix and pheno data.

```
dsbyrne2020$g[1:5, 1:4]
#>                GSM2463097 GSM2463098 GSM2463099 GSM2463100
#> WBGene00000001  3.8512729  3.9854008  3.8985759  4.0972464
#> WBGene00000002  1.2460447  1.4028785  1.2457411  0.9150779
#> WBGene00000003  1.5149659  1.5658889  1.6733056  1.7358959
#> WBGene00000004  1.8531309  1.7911185  1.9029840  2.2112104
#> WBGene00000005  0.8301227  0.7806835  0.7145102  0.3913988

head(dsbyrne2020$p, n = 5)
#>               title geo_accession age
#> GSM2463097 d3_ev_R1    GSM2463097   3
#> GSM2463098 d3_ev_R3    GSM2463098   3
#> GSM2463099 d3_ev_R4    GSM2463099   3
#> GSM2463100 d6_ev_R1    GSM2463100   6
#> GSM2463101 d6_ev_R3    GSM2463101   6
```

Correlation between samples for outliers.

``` r
cor_dsbyrne2020 <- cor(dsbyrne2020$g, method = "spearman")
diag(cor_dsbyrne2020) <- NA
ord <- order(dsbyrne2020$p$age)
cols <- as.character(1:3)
heatmap(cor_dsbyrne2020[ord, ord], Colv = NA, Rowv = NA,
        scale = "none", keep.dendro = F, margins = c(1,5),
        labRow = "", labCol = "")
par(xpd = T)
agetab <- table(dsbyrne2020$p$age)
mtext(text = names(agetab), side = 1, line = c(4),
      at = seq(-.17, .86, l = length(dsbyrne2020$p$age))[
        round(cumsum(agetab) - agetab/2 +0.1)
      ], cex = .6)

# color key and labels
l.values <- seq(min(cor_dsbyrne2020, na.rm = T),
                max(cor_dsbyrne2020, na.rm = T), l = 10)
image(x = c(.95,1), y = seq(0.6,1, l = 10), useRaster = T,
      z = matrix(l.values, ncol = 10),
      col = hcl.colors(12, "YlOrRd", rev = TRUE), add = T)
text(.975, 1, pos = 3, labels = expression(rho), font = 2)
text(1, y = seq(0.6,1, l = 10), pos = 4,
     labels = round(l.values, 2), cex = .6)
mtext(at = 1.025, line = 4, side = 1, text = "(days)", cex = .8)
```

``` r
boxplot(cor_dsbyrne2020, xaxt="n",
        lwd=2, boxwex=.5, las=1, col = 0,
        at = 1:sum(agetab)+rep(1:length(agetab), agetab),
        ylab = "Spearman correlation",
        xlab = "Chronological age (days of adulthood)")
#add time labels
axis(side = 1, at = (1:sum(agetab)+rep(1:length(agetab), agetab))[
  round(cumsum(agetab) - agetab/2 +0.1)],
  labels = paste0(unique(dsbyrne2020$p$age), 'd'))
```

No outliers. We note the overall sample-sample correlation is lower than with development time-series experiments.

### [4.4.3] Filtering genes

We will keep genes with a strong monotonic aging signal, which we define as those with absolute spearman correlation with age \\(\> \\sqrt(1/3)\\). This roughly corresponds to monotonous genes for which aging explains over \\(1/3\\) of the variance.

```
# compute correlation
gn_cor <- apply(dsbyrne2020$g, 1, cor, y=dsbyrne2020$p$age, method = "spearman")

selg <- (abs(gn_cor)>sqrt(1/3))
table(selg)
#> selg
#> FALSE  TRUE
#> 22245  8704
```

Plotting principal components.

```
pca_dsbyrne2020 <- stats::prcomp(t(dsbyrne2020$g[selg,]),
                                 center = TRUE, scale = FALSE)
```

``` r
par(mfrow = c(1,4), pty='s')
invisible(sapply(1:4, function(i){
  plot(dsbyrne2020$p$age, pca_dsbyrne2020$x[,i],
       lwd = 2, xlab = "Chronological age (days)", ylab = "PC",
       main = paste0("PC", i))
}))
```

### [4.4.4] Model fitting

We only keep the first component, which explains \\(72.71 \\%\\) of the variance in the data. Given the previous selection of genes, this component is monotonous.

We fit a GAM on the component.

```
m_dsbyrne2020 <- ge_im(X = dsbyrne2020$g[selg,],
                       p = dsbyrne2020$p,
                       formula = "X ~ s(age, bs = 'cr', k=4)",
                       nc = 1)
```

### [4.4.5] Validation

Check global model performance.

    #>            aCC       aRE        MSE    aRMSE
    #> [1,] 0.7675197 0.1954797 0.02296315 0.151536

And then per gene.

```
ng_mp_hou <-  mperf(dsbyrne2020$g[selg, ], predict(m_dsbyrne2020),
                     is.t = T, global = F)
# remove NAs (eg. 0 variance genes) and Inf values (/0)
ng_mp_hou <- lapply(ng_mp_hou, na.omit)
ng_mp_hou$aRE <- ng_mp_hou$aRE[ng_mp_hou$aRE < Inf]
```

``` r
par(mfrow = c(2,2), mar=c(4,4,2,1))
invisible(sapply(names(ng_mp_hou), function(idx){
  rg <- range(na.omit(ng_mp_hou[[idx]]))
  # estimate density curve
  d <- density(na.omit(ng_mp_hou[[idx]]), from = rg[1], to = rg[2])
  plot(d, main = paste0(gsub("a", "", idx, fixed = T),
                        " density (", length(ng_mp_hou[[idx]]), " genes)"),
       xlab = idx, lwd = 2)
  # add global value
  abline(v = mp_hou[[idx]], lty = 2, lwd = 2, col = "firebrick")
  text(mp_hou[[idx]], .9*max(d$y), pos = 4, labels = idx,
       font = 2, col = "firebrick")
}))
```

Then, build the `ref` object.

```
# make a 'ref' object'
r_dsbyrne2020 <- make_ref(
  m_dsbyrne2020,
  by.inter = .01,
  t.unit = "days of adulthood, 20C",     # time unit
  metadata = list("organism" ="C. elegans", # any metadata
                  "profiling"="bulk RNAseq",
                  "condition"="rrf-3 mutant, liquid culture"))
```

Plot component interpolation.

``` r
par(mfrow = c(1,1), pty='s')
plot(m_dsbyrne2020, r_dsbyrne2020, l.pos = 'bottomright')
```

Looks good.

```
# stage samples
ae_dsbyrne2020 <- ae(dsbyrne2020$g, r_dsbyrne2020)
ae_dshou2016 <- ae(dshou2016$g, r_dsbyrne2020)
```

``` r
par(mfrow = c(1,2), mar = c(4,5,3,1), pty='s')
rg <- range(c(ae_dsbyrne2020$age.estimates[,1], dsbyrne2020$p$age))
plot(ae_dsbyrne2020$age.estimates[,1]~dsbyrne2020$p$age,
     xlab = "Chronological age (d of adulthood)",
     ylab = "Estimated age (r_dsbyrne2020)",
     xlim = rg, ylim = rg,
     main = "Chron. vs Estimated ages for dsbyrne2020
(on r_dsbyrne2020)",
     lwd = 2)
abline(a = 0, b = 1, lty = 3, lwd = 1)
legend("topleft", legend = c("x = y"), lwd=3, col=1, lty = 3, bty='n')

rg <- range(c(ae_dshou2016$age.estimates[,1], dshou2016$p$age))
plot(ae_dshou2016$age.estimates[,1]~dshou2016$p$age,
     col = dshou2016$p$diet, xlim = rg, ylim = rg,
     xlab = "Chronological age (d of adulthood)",
     ylab = "Estimated age (r_dsbyrne2020)",
     main = "Chron. vs Estimated ages for dshou2016
(on r_dsbyrne2020)", lwd = 2)
abline(a = 0, b = 1, lty = 3, lwd = 2)
legend("topleft", lty=NA, lwd=3, pch=1, col=1:3, bty='n',
       legend = c("ad ilibitum", "carloric restriction",
                  "intermittent fasting"))
```

We expected the youngest 2-day-old samples of `dshou2016` to be outside the scope of the reference, which starts at 3 days of adulthood. However, beyond this we also note a consistent difference of around 3 days between chronological and estimated age.

This offset could be explained by a combination of many factors, such as the few listed below.

  --------------------------------------------------------------------------
                            dsbyrne2020                   dshou2016
  ---------------- ------------------------------ --------------------------
  Strain            rrf-3 (RNAi hyper-sensitive)              N2

  Embryo signal               fertile              sterile (FUDR treatment)

  Culture medium           liquid culture                 NGM plates

  Food source            EV HT115 (E. coli)        UV-killed OP50 (E. coli)
  --------------------------------------------------------------------------

The difference could also come from unknown differences in experimental procedures between labs, or even from a difference in their definition of *"days of adulthood"*, as [Lithgow, Driscoll, and Phillips (2017)] have previously discussed.

### [4.4.6] Code to generate objects

Required packages and variables:

``` r
data_folder <- "../inst/extdata/"

requireNamespace("RAPToR", quietly = T)
requireNamespace("wormRef", quietly = T)
requireNamespace("utils", quietly = T)
requireNamespace("GEOquery", quietly = T) # bioconductor
requireNamespace("Biobase", quietly = T)  # bioconductor
requireNamespace("biomaRt", quietly = T)  # bioconductor
requireNamespace("affy", quietly = T)  # bioconductor
```

*Note : set `data_folder` to an existing path on your system where you want to store the objects.*

``` r
raw2tpm <- function(rawcounts, genelengths){
  if(nrow(rawcounts) != length(genelengths))
    stop("genelengths must match nrow(rawcounts).")
  x <- rawcounts/genelengths
  return(t( t(x) * 1e6 / colSums(x) ))
}

fpkm2tpm <- function(fpkm){
  return(exp(log(fpkm) - log(colSums(fpkm)) + log(1e6)))
}
```

To build `dsbyrne2020`, *C. elegans* unpublished aging time-series by Byrne et al., [GSE93826](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93826)

``` r
geo_id <- "GSE93826"
geo_obj <- GEOquery::getGEO(geo_id)[[1]]

# get expression data
sfs <- GEOquery::getGEOSuppFiles(geo_id, fetch_files = F, makeDirectory = F)

tmpfolder <- file.path(data_folder, "byrne2020tmp")
dir.create(tmpfolder)
destfile <- file.path(tmpfolder, sfs$fname[1])
download.file(sfs$url[1], destfile = destfile)

untar(destfile, exdir = tmpfolder)

flist <- list.files(tmpfolder)
g <- lapply(seq_along(flist)[-1], function(i){
  gi <- read.table(gzfile(file.path(tmpfolder,flist[i])),
                   h=F, row.names = 1)
  colnames(gi) <- gsub("(GSM246\\d+)_.*", "\\1", flist[i])
  return(gi)
})

g <- do.call(cbind, g)

# format and convert to tpm
g <- RAPToR::format_ids(g, wormRef::Cel_genes, from = "wb_id",
                        to = "wb_id", aggr.fun = sum)

g <- raw2tpm(g, wormRef::Cel_genes$transcript_length[
  match(rownames(g), wormRef::Cel_genes$wb_id)])
g <- g[apply(g, 1, sum) > 0, ] # keep expressed genes only
g <- g[-3123, ] # remove 0-variance artefact gene

# get pheno data
p <-  Biobase::pData(geo_obj)
p <- p[, c("title", "geo_accession", "age:ch1")]
colnames(p)[3] <- "age"
p$age <- as.numeric(as.character(
  factor(p$age, levels = paste0('day ', c(3,6,9,12,15,18)),
         labels = c(3,6,9,12,15,18))
  ))

p <- p[order(p$age), ]
g <- g[, p$geo_accession]

# save data
dsbyrne2020 <- list(g=g, p=p)
save(dsbyrne2020, file = file.path(data_folder, "dsbyrne2020.RData"),
     compress = "xz")

# cleanup
file.remove(file.path(tmpfolder, flist))
file.remove(tmpfolder, destfile)
rm(g, p, tmpfolder, flist, geo_id, geo_obj, sfs, destfile)
```

To build `dshou2016`, *C. elegans* aging time-series from [Hou et al. (2016)]

``` r
geo_id <- "GSE77110"
geo_obj <- GEOquery::getGEO(geo_id)[[1]]

# get pheno data
p <- Biobase::pData(geo_obj)
p <- p[, c("title", "geo_accession", "age:ch1",
           "diet:ch1", "source_name_ch1")]
p$age <- as.numeric(gsub("adult day ", "", as.character(p$`age:ch1`)))
p <- p[, -3]
colnames(p) <- c("title", "geo_accession", "diet", "source_name", "age")
p$diet <- factor(p$diet, levels = c("ad libitum",
                                    "calorie restriction",
                                    "intermittent fasting"))

# get microarray probe IDs
mart <- biomaRt::useMart(biomart = "ensembl",
                         dataset = "celegans_gene_ensembl")
probe_ids <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                           "affy_c_elegans"),
                            mart = mart)

# download data
sfile <- GEOquery::getGEOSuppFiles(geo_id, makeDirectory = F, fetch_files = F)
# You may need to increase the allowed connection size to download the dataset
# Sys.setenv("VROOM_CONNECTION_SIZE") <- 131072*4

tarfolder <- paste0(data_folder,"raw_hou")
dir.create(tarfolder, showWarnings = F)
tarfile <- paste0(tarfolder,'/', as.character(sfile$fname[1]))
utils::download.file(url = as.character(sfile$url[1]), destfile = tarfile)
untar(tarfile = tarfile, exdir = tarfolder)

flist <- paste0(tarfolder,'/',list.files(tarfolder))
for (f in flist[grepl(".gz", flist)]){
  GEOquery::gunzip(filename = f, destname = gsub("(.*).gz", "\\1", f),
         overwrite = T, remove = T)
}

# load and format gdata
flist <- list.files(tarfolder)
flist <- flist[sapply(p$geo_accession, function(gg) which(grepl(gg, flist)))]
g <- affy::ReadAffy(filenames = flist, celfile.path = tarfolder, phenoData = p)
g <- affy::expresso(g, bg.correct = F, normalize = F,
                    pmcorrect.method = "pmonly", summary.method = "median")
g <- 2^Biobase::exprs(g) # expresso log2s the data
g <- RAPToR::format_ids(g, probe_ids, from = 2, to = 1)

# save data
dshou2016 <- list(g=g, p=p)
save(dshou2016, file = file.path(data_folder, "dshou2016.RData"),
     compress = "xz")

# cleanup
file.remove(file.path(tarfolder, list.files(tarfolder)), tarfolder)
file.remove(tarfile)

rm(geo_id, geo_obj, g, p, probe_ids, mart,
   tarfolder, tarfile, flist, sfile, f)
```

[Back to top](#top)

------------------------------------------------------------------------

# References

Aeschimann, Florian, Pooja Kumari, Hrishikesh Bartake, Dimos Gaidatzis, Lan Xu, Rafal Ciosk, and Helge Großhans. 2017. "LIN41 Post-Transcriptionally Silences mRNAs by Two Distinct and Position-Dependent Mechanisms." *Molecular Cell* 65 (3): 476--89.

Alter, Orly, Patrick O Brown, and David Botstein. 2000. "Singular Value Decomposition for Genome-Wide Expression Data Processing and Modeling." *Proceedings of the National Academy of Sciences* 97 (18): 10101--6.

Bulteau, Romain, and Mirko Francesconi. 2022. "Real Age Prediction from the Transcriptome with RAPToR." *Nature Methods*, 1--7.

Graveley, Brenton R, Angela N Brooks, Joseph W Carlson, Michael O Duff, Jane M Landolin, Li Yang, Carlo G Artieri, et al. 2011. "The Developmental Transcriptome of Drosophila Melanogaster." *Nature* 471 (7339): 473.

Hendriks, Gert-Jan, Dimos Gaidatzis, Florian Aeschimann, and Helge Großhans. 2014. "Extensive Oscillatory Gene Expression During c. Elegans Larval Development." *Molecular Cell* 53 (3): 380--92.

Hou, Lei, Dan Wang, Di Chen, Yi Liu, Yue Zhang, Hao Cheng, Chi Xu, et al. 2016. "A Systems Approach to Reverse Engineer Lifespan Extension by Dietary Restriction." *Cell Metabolism* 23 (3): 529--40.

Levin, Michal, Leon Anavy, Alison G Cole, Eitan Winter, Natalia Mostov, Sally Khair, Naftalie Senderovich, et al. 2016. "The Mid-Developmental Transition and the Evolution of Animal Body Plans." *Nature* 531 (7596): 637.

Lithgow, Gordon J, Monica Driscoll, and Patrick Phillips. 2017. "A Long Journey to Reproducible Results." *Nature* 548 (7668): 387--88.

White, Richard J, John E Collins, Ian M Sealy, Neha Wali, Christopher M Dooley, Zsofia Digby, Derek L Stemple, et al. 2017. "A High-Resolution mRNA Expression Time Course of Embryonic Development in Zebrafish." *Elife* 6: e30860.

------------------------------------------------------------------------

# SessionInfo

```
sessionInfo()
```

    #> R version 4.1.2 (2021-11-01)
    #> Platform: x86_64-pc-linux-gnu (64-bit)
    #> Running under: Ubuntu 22.04.1 LTS
    #>
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so
    #>
    #> locale:
    #>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
    #>  [3] LC_TIME=C                  LC_COLLATE=en_US.UTF-8
    #>  [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=en_US.UTF-8
    #>  [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C
    #>  [9] LC_ADDRESS=C               LC_TELEPHONE=C
    #> [11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C
    #>
    #> attached base packages:
    #>  [1] stats4    splines   parallel  stats     graphics  grDevices utils
    #>  [8] datasets  methods   base
    #>
    #> other attached packages:
    #>  [1] viridis_0.6.2               ROCR_1.0-11
    #>  [3] DESeq2_1.34.0               SummarizedExperiment_1.24.0
    #>  [5] Biobase_2.54.0              MatrixGenerics_1.6.0
    #>  [7] matrixStats_0.63.0          GenomicRanges_1.46.1
    #>  [9] GenomeInfoDb_1.30.1         IRanges_2.28.0
    #> [11] S4Vectors_0.32.4            BiocGenerics_0.40.0
    #> [13] viridisLite_0.4.1           vioplot_0.4.0
    #> [15] zoo_1.8-11                  sm_2.2-5.7.1
    #> [17] ica_1.0-3                   drosoRef_0.2.0
    #> [19] limma_3.50.3                wormRef_0.5.0
    #> [21] RAPToR_1.2.0                ggpubr_0.6.0
    #> [23] ggplot2_3.4.0.9000          BiocStyle_2.22.0
    #>
    #> loaded via a namespace (and not attached):
    #>  [1] colorspace_2.1-0       ggsignif_0.6.4         pryr_0.1.6
    #>  [4] XVector_0.34.0         rstudioapi_0.14        farver_2.1.1
    #>  [7] hexbin_1.28.2          bit64_4.0.5            AnnotationDbi_1.56.2
    #> [10] fansi_1.0.4            codetools_0.2-18       cachem_1.0.6
    #> [13] geneplotter_1.72.0     knitr_1.42             jsonlite_1.8.4
    #> [16] broom_1.0.3            annotate_1.72.0        png_0.1-8
    #> [19] BiocManager_1.30.19    compiler_4.1.2         httr_1.4.4
    #> [22] backports_1.4.1        Matrix_1.5-3           fastmap_1.1.0
    #> [25] cli_3.6.0              htmltools_0.5.4        tools_4.1.2
    #> [28] gtable_0.3.1           glue_1.6.2             GenomeInfoDbData_1.2.7
    #> [31] dplyr_1.1.0            tinytex_0.44           Rcpp_1.0.10
    #> [34] carData_3.0-5          jquerylib_0.1.4        vctrs_0.5.2
    #> [37] Biostrings_2.62.0      nlme_3.1-161           xfun_0.37
    #> [40] stringr_1.5.0          rbibutils_2.2.13       lifecycle_1.0.3
    #> [43] rstatix_0.7.2          XML_3.99-0.13          zlibbioc_1.40.0
    #> [46] scales_1.2.1           RColorBrewer_1.1-3     yaml_2.3.7
    #> [49] gridExtra_2.3          memoise_2.0.1          sass_0.4.5
    #> [52] stringi_1.7.12         RSQLite_2.2.20         highr_0.10
    #> [55] genefilter_1.76.0      BiocParallel_1.28.3    Rdpack_2.4
    #> [58] rlang_1.0.6            pkgconfig_2.0.3        bitops_1.0-7
    #> [61] evaluate_0.20          lattice_0.20-45        purrr_1.0.1
    #> [64] labeling_0.4.2         cowplot_1.1.1          bit_4.0.5
    #> [67] tidyselect_1.2.0       magrittr_2.0.3         bookdown_0.32
    #> [70] R6_2.5.1               magick_2.7.3           generics_0.1.3
    #> [73] DelayedArray_0.20.0    DBI_1.1.3              pillar_1.8.1
    #> [76] withr_2.5.0            mgcv_1.8-41            survival_3.5-0
    #> [79] KEGGREST_1.34.0        abind_1.4-5            RCurl_1.98-1.10
    #> [82] tibble_3.1.8           crayon_1.5.2           car_3.1-1
    #> [85] utf8_1.2.3             rmarkdown_2.20.1       locfit_1.5-9.7
    #> [88] grid_4.1.2             data.table_1.14.6      blob_1.2.3
    #> [91] digest_0.6.31          xtable_1.8-4           tidyr_1.3.0
    #> [94] munsell_0.5.0          beeswarm_0.4.0         bslib_0.4.2

---

# RAPToR-datapkgs

Code []

-   [Show All Code](#)
-   [Hide All Code](#)

# `RAPToR` - Data-packages

### RAPToR 1.2.0

Romain Bulteau

#### March 2023

# Contents

-   [[1] Introduction](#introduction)
-   [[2] What's a data-package ?](#whats-a-data-package)
-   [[3] Reference data](#reference-data)
-   [[4] Data-package interface with RAPToR](#data-package-interface-with-raptor)
    -   [[4.1] `.prepref_` functions](#prepref_-functions)
    -   [[4.2] `ref_list` object](#ref_list-object)
    -   [[4.3] `.plot_refs()` function](#plot_refs-function)
    -   [[4.4] Other objects](#other-objects)
-   [SessionInfo](#sessioninfo)

# [1] Introduction

This vignette is aimed at those who've already familiarized themselves with [`RAPToR` reference building](RAPToR-refbuilding.html). It's also for us to keep track of guidelines to continue improving `RAPToR` with new references.

If you've already built some references and want to make them available to the world (or use them more easily yourself), you're in the right place. You will need some basic knowledge of R package development. Data-packages are, after all, *packages*. This document only details how to set up your data-package to interact properly with `RAPToR`.

# [2] What's a data-package ?

By definition, a "data-package" is an R package in which one stores large datasets (over a few Mo). This is a good practice for several reasons.

1.  If your data rarely or never changes, updates to the data-package (and thus, download of the data) will be minimal. If included in a standard package, large data can be a burden during install.
2.  The data may never be used. Why have users download data they won't need ?
3.  CRAN standards limit package size to 5MB (documentation included). A large dataset is better off separated from methods and functions that may need it.

Hadley Wickam gives thorough advice on organizing data in packages in his [*R packages* book](http://r-pkgs.had.co.nz/data.html).

# [3] Reference data

`RAPToR` uses references which can sometimes be tedious to build, so we want to give users access to pre-built references as much as possible. References are stored as `.RData` objects, which must include everything needed to make the interpolated reference:

-   gene expression data,
-   optimal interpolation parameters,
-   some metadata (*e.g.* time units).

There can be as many references as needed in a data-package; we group packages by organism.

Being consistent, clear and concise with naming can help users (or yourself !) find their way around references. For example, `wormRef` references are named with an organism code `Cel` (*C. elegans*) followed by the developmental period covered by the reference *e.g.* `larval`.

Document the references. It is standard practice to document data. What is the data? Where is the data from? Publication? etc.

The structure of the `Cel_larval` reference object is detailed below as an example.

``` r
str(wormRef::Cel_larval, vec.len = 2)
#> List of 6
#>  $ g          : num [1:18718, 1:62] 2.72 3.06 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:18718] "WBGene00007063" "WBGene00007064" ...
#>   .. ..$ : chr [1:62] "DH2_N2_0" "DH2_N2_2" ...
#>  $ p          :'data.frame': 62 obs. of  5 variables:
#>   ..$ sname    : chr [1:62] "DH2_N2_0" "DH2_N2_2" ...
#>   ..$ age      : num [1:62] 0 2 4 6 8 ...
#>   ..$ cov      : Factor w/ 3 levels "O.20","O.25",..: 1 1 1 1 1 ...
#>   ..$ age_ini  : num [1:62] 0 2 4 6 8 ...
#>   ..$ accession: chr [1:62] "GSM1192801" "GSM1192802" ...
#>  $ geim_params:List of 4
#>   ..$ formula: chr "X ~ s(age, bs = 'ds') + cov"
#>   ..$ method : chr "gam"
#>   ..$ dim_red: chr "pca"
#>   ..$ nc     : num 40
#>  $ t.unit     : chr "h past egg-laying (20C)"
#>  $ cov.levels :List of 1
#>   ..$ cov: chr "O.20"
#>  $ metadata   :List of 3
#>   ..$ organism  : chr "C. elegans"
#>   ..$ profiling : chr "whole-organism, bulk"
#>   ..$ technology: chr "RNAseq"
```

With

-   `g` The gene expression matrix (genes as rows, samples as columns).
-   `p` A dataframe of phenotypic data on the samples :
    -   `sname` sample names,
    -   `age` *developmental* age of the samples (scaled),
    -   `cov` covariate, factor indicating which of 3 time series,
    -   `age_ini` *chronological* age of the samples,
    -   `accession` sample accession ID for GEO.
-   `geim_params` A list with necessary parameters for interpolation
-   `t.unit` the time unit.
-   `cov.levels` A named list with covariate levels to interpolate as.
-   `metadata` A named list with any extra metadata

The documentation for `Cel_larval` can be accessed with `?Cel_larval`.

# [4] Data-package interface with RAPToR

For a user to access data-package information and references directly from `RAPToR`, we've set up a standard system using reference names.

A few objects are *necessary* for this interface to work.

## [4.1] `.prepref_` functions

`.prepref_` functions (note the "dot") are the key of the interface : they **prep**are the **ref**erence for the user. They must respect the naming convention `.prepref_ref_name()` (*e.g.* `.prepref_Cel_larval()`) and take `n.inter` and/or `by.inter` as arguments.

These functions are the backbone called by `prepare_refdata()` when fetching a reference, and thus should output the reference `ref` object with the specified parameters. This means building the GEIM model, and calling `make_ref()` with the appropriate parameters and metadata.

We have made a [function factory](https://adv-r.hadley.nz/function-factories.html) to generate these functions. It inputs the reference data object described above, and returns the corresponding prepref function.

```
.prepref_skel <- function(data, from=NULL, to=NULL){
  # .prepref function factory
  f <- function(n.inter=NULL, by.inter=NULL){
    m <- RAPToR::ge_im(
      X = data$g,
      p = data$p,
      formula = data$geim_params$formula,
      method = data$geim_params$method,
      dim_red = data$geim_params$dim_red,
      nc = data$geim_params$nc
    )
    return(RAPToR::make_ref(m,
                            n.inter = n.inter,
                            by.inter = by.inter,
                            from = from,
                            to = to,
                            t.unit = data$t.unit,
                            cov.levels = data$cov.levels,
                            metadata = data$metadata)
    )
  }
  return(f)
}
```

To make `.prepref_Cel_larval()`, we simply include the following code in the data-package, along with the function factory.

```
.prepref_Cel_larval <- .prepref_skel(wormRef::Cel_larval)
```

## [4.2] `ref_list` object

`RAPToR` expects a `ref_list` object in the data-package. This is what's displayed when calling the `list_refs(datapkg)` function.

```
library(RAPToR)
library(wormRef)

list_refs(datapkg = "wormRef")
#>          name               organism                 description
#> 1  Cel_embryo Caenorhabditis elegans       Embryonic development
#> 2  Cel_larval Caenorhabditis elegans          Larval development
#> 3 Cel_larv_YA Caenorhabditis elegans Larval to adult development
#> 4    Cel_YA_1 Caenorhabditis elegans              YA development
#> 5    Cel_YA_2 Caenorhabditis elegans              YA development
#>                                range
#> 1 1-cell(-50) to 840 min past 4-cell
#> 2      0 to 55 h post-hatching (20C)
#> 3      5 to 76 h post-hatching (20C)
#> 4   41.5 to 72 h post-hatching (20C)
#> 5     48 to 87 h post-hatching (20C)
```

The form/layout of this object is free, but reference names should be included somewhere, as the user needs them to access the reference through `prepare_refdata()`.

## [4.3] `.plot_refs()` function

This function is optional, but very useful to guide users to the correct reference for their samples. Including a `.plot_refs()` function (note the "dot") in a data-package will allow it to be called by the `plot_refs(datapkg)` function in `RAPToR`.

```
plot_refs(datapkg = "wormRef")
```

## [4.4] Other objects

You're free to include any extra objects in your data-package that may be useful. For example, the `wormRef` package has a `Cel_devstages` object with information on key developmental stages of *C. elegans* (which is used for building the plot in `.plot_refs()` above).

# SessionInfo

```
sessionInfo()
```

    #> R version 4.1.2 (2021-11-01)
    #> Platform: x86_64-pc-linux-gnu (64-bit)
    #> Running under: Ubuntu 22.04.1 LTS
    #>
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so
    #>
    #> locale:
    #>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
    #>  [3] LC_TIME=C                  LC_COLLATE=en_US.UTF-8
    #>  [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=en_US.UTF-8
    #>  [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C
    #>  [9] LC_ADDRESS=C               LC_TELEPHONE=C
    #> [11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C
    #>
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods
    #> [7] base
    #>
    #> other attached packages:
    #>  [1] vioplot_0.4.0      zoo_1.8-11         sm_2.2-5.7.1
    #>  [4] ica_1.0-3          drosoRef_0.2.0     limma_3.50.3
    #>  [7] wormRef_0.5.0      ggpubr_0.6.0       ggplot2_3.4.0.9000
    #> [10] BiocStyle_2.22.0   RAPToR_1.2.0
    #>
    #> loaded via a namespace (and not attached):
    #>  [1] Rcpp_1.0.10         lattice_0.20-45     tidyr_1.3.0
    #>  [4] digest_0.6.31       utf8_1.2.3          R6_2.5.1
    #>  [7] backports_1.4.1     evaluate_0.20       highr_0.10
    #> [10] pillar_1.8.1        Rdpack_2.4          rlang_1.0.6
    #> [13] rstudioapi_0.14     data.table_1.14.6   car_3.1-1
    #> [16] jquerylib_0.1.4     magick_2.7.3        Matrix_1.5-3
    #> [19] rmarkdown_2.20.1    labeling_0.4.2      splines_4.1.2
    #> [22] stringr_1.5.0       munsell_0.5.0       tinytex_0.44
    #> [25] broom_1.0.3         compiler_4.1.2      xfun_0.37
    #> [28] pkgconfig_2.0.3     mgcv_1.8-41         htmltools_0.5.4
    #> [31] tidyselect_1.2.0    tibble_3.1.8        bookdown_0.32
    #> [34] codetools_0.2-18    fansi_1.0.4         dplyr_1.1.0
    #> [37] withr_2.5.0         rbibutils_2.2.13    grid_4.1.2
    #> [40] nlme_3.1-161        jsonlite_1.8.4      gtable_0.3.1
    #> [43] lifecycle_1.0.3     magrittr_2.0.3      scales_1.2.1
    #> [46] cli_3.6.0           stringi_1.7.12      cachem_1.0.6
    #> [49] carData_3.0-5       farver_2.1.1        ggsignif_0.6.4
    #> [52] pryr_0.1.6          bslib_0.4.2         generics_0.1.3
    #> [55] vctrs_0.5.2         tools_4.1.2         glue_1.6.2
    #> [58] beeswarm_0.4.0      purrr_1.0.1         abind_1.4-5
    #> [61] parallel_4.1.2      fastmap_1.1.0       yaml_2.3.7
    #> [64] colorspace_2.1-0    BiocManager_1.30.19 rstatix_0.7.2
    #> [67] knitr_1.42          sass_0.4.5

---
