
# Effect size estimation with apeglm

#### Anqi Zhu, Joseph G. Ibrahim, and Michael I. Love

#### 10/29/2025

Abstract

*apeglm* provides empirical Bayes shrinkage estimators for effect sizes for a variety of GLM models; *apeglm* stands for "Approximate Posterior Estimation for GLM". apeglm package version: 1.32.0

# Typical RNA-seq call from DESeq2

**Note:** the typical RNA-seq workflow for users would be to call *apeglm* estimation from within the `lfcShrink` function from the *DESeq2* package. The unevaluated code chunk shows how to obtain *apeglm* shrinkage estimates after running `DESeq`. See the DESeq2 vignette for more details. The `lfcShrink` wrapper function takes care of many details below, and unifies the interface for multiple shrinkage estimators. The coefficient to shrink can be specified either by name or by number (following the order in `resultsNames(dds)`). Be aware that *DESeq2*'s `lfcShrink` interface provides LFCs on the log2 scale, while *apeglm* provides coefficients on the natural log scale.

```
res <- lfcShrink(dds, coef=2, type="apeglm")
```

# Acknowledgments

Joshua Zitovsky contributed fast C++ code for the beta-binomial likelihood, demonstrated in a later section of the vignette.

We have benefited in the development of `apeglm` from feedback or contributions from the following individuals:

Wolfgang Huber, Cecile Le Sueur, Charlotte Soneson

# Example RNA-seq analysis

Here we show example code which mimics what will happen inside the `lfcShrink` function when using the *apeglm* method [(Zhu, Ibrahim, and Love 2018)].

Load a prepared `SummarizedExperiment`:

```
library(airway)
data(airway)
head(assay(airway))
```

    ##                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
    ## ENSG00000000003        679        448        873        408       1138
    ## ENSG00000000005          0          0          0          0          0
    ## ENSG00000000419        467        515        621        365        587
    ## ENSG00000000457        260        211        263        164        245
    ## ENSG00000000460         60         55         40         35         78
    ## ENSG00000000938          0          0          2          0          1
    ##                 SRR1039517 SRR1039520 SRR1039521
    ## ENSG00000000003       1047        770        572
    ## ENSG00000000005          0          0          0
    ## ENSG00000000419        799        417        508
    ## ENSG00000000457        331        233        229
    ## ENSG00000000460         63         76         60
    ## ENSG00000000938          0          0          0

For demonstration, we will use 3000 genes of the `airway` dataset, the first from those genes with at least 10 counts across all samples.

```
keep <- head(which(rowSums(assay(airway)) >= 10), 3000)
airway <- airway[keep,]
```

First run a *DESeq2* differential expression analysis [(Love, Huber, and Anders 2014)] (size factors, and dispersion estimates could similarly be estimated using *edgeR*):

```
library(DESeq2)
dds <- DESeqDataSet(airway, ~cell + dex)
dds$dex <- relevel(dds$dex, "untrt")
dds <- DESeq(dds)
res <- results(dds)
```

Defining data and parameter objects necessary for `apeglm`. We must multiply the coefficients from *DESeq2* by a factor, because *apeglm* provides natural log coefficients. Again, this would be handled inside of `lfcShrink` in *DESeq2* for a typical RNA-seq analysis.

```
x <- model.matrix(design(dds), colData(dds))
param <- dispersions(dds)
mle <- log(2) * cbind(res$log2FoldChange, res$lfcSE)
offset <- matrix(log(sizeFactors(dds)),
                 ncol=ncol(dds),
                 nrow=nrow(dds),byrow=TRUE)
```

# Running apeglm

Here `apeglm` on 3000 genes takes a few seconds on a laptop. It scales with number of genes, the number of samples and the number of variables in the design formula, where here we have 5 coefficients (one for the four cell cultures and one for the difference due to dexamethasone treatment).

We provide `apeglm` with the *SummarizedExperiment* although the function can also run on a *matrix* of counts or other observed data. We specify a `coef` as well as a `threshold` which we discuss below. Note that we multiple the `threshold` by `log(2)` to convert from log2 scale to natural log scale.

The original interface to apeglm looked as follows, but was slow, we therefore recommend using the faster interfaces shown below.

```
# original interface, also may give Lapack error
fit <- apeglm(Y=airway, x=x, log.lik=logLikNB, param=param, coef=ncol(x),
              threshold=log(2) * 1, mle=mle, offset=offset)
```

There are better, faster implementations of apeglm specifically for negative binomial likelihoods. The version `nbinomR` is \~5 times faster than the default `method="general"`.

We will use this run for downstream analysis. Here we pass the *SummarizedExperiment*, which will pass back the ranges. For faster `apeglm` runs, provide only the counts matrix (as seen in subsequent chunks).

```
library(apeglm)
system.time({
  fitR <- apeglm(Y=airway, x=x, log.lik=NULL, param=param, coef=ncol(x),
                 threshold=log(2) * 1, mle=mle, offset=offset, method="nbinomR")
})
```

    ##    user  system elapsed
    ##   3.078   0.115   3.194

```
fit <- fitR
names(fit)
```

    ## [1] "map"           "sd"            "prior.control" "fsr"
    ## [5] "svalue"        "interval"      "thresh"        "diag"
    ## [9] "ranges"

```
str(fit$prior.control)
```

    ## List of 7
    ##  $ no.shrink            : int [1:4] 1 2 3 4
    ##  $ prior.mean           : num 0
    ##  $ prior.scale          : num 0.304
    ##  $ prior.df             : num 1
    ##  $ prior.no.shrink.mean : num 0
    ##  $ prior.no.shrink.scale: num 15
    ##  $ prior.var            : num 0.0921

The version `nbinomCR` is \~10 times faster than the default `general`.

```
system.time({
  fitCR <- apeglm(Y=assay(airway), x=x, log.lik=NULL, param=param, coef=ncol(x),
                 threshold=log(2) * 1, mle=mle, offset=offset, method="nbinomCR")
})
```

    ##    user  system elapsed
    ##   1.674   0.056   1.730

The version `nbinomC` returns only the MAP coefficients and can be \~50-100 times faster than the default `general`. The MAP coefficients are the same as returned by `nbinomCR` above, we just skip the calculation of posterior SD. A variant of `nbinomC` is `nbinomC*` which includes random starts.

```
system.time({
  fitC <- apeglm(Y=assay(airway), x=x, log.lik=NULL, param=param, coef=ncol(x),
                 threshold=log(2) * 1, mle=mle, offset=offset, method="nbinomC")
})
```

    ##    user  system elapsed
    ##   0.121   0.000   0.121

Among other output, we have the estimated coefficients attached to the ranges of the *SummarizedExperiment* used as input:

```
class(fit$ranges)
```

    ## [1] "CompressedGRangesList"
    ## attr(,"package")
    ## [1] "GenomicRanges"

```
mcols(fit$ranges, use.names=TRUE)
```

    ## DataFrame with 3000 rows and 5 columns
    ##                 X.Intercept. cellN061011 cellN080611 cellN61311     dextrt
    ##                    <numeric>   <numeric>   <numeric>  <numeric>  <numeric>
    ## ENSG00000000003      6.64788   0.0546137   0.2242529 -0.1570180 -0.2602199
    ## ENSG00000000419      6.23638  -0.0912014  -0.0206726 -0.0526312  0.1203700
    ## ENSG00000000457      5.45524   0.0540325  -0.0611456  0.0435641  0.0133132
    ## ENSG00000000460      3.77417   0.5394114   0.2564594  0.3450027 -0.0420874
    ## ENSG00000000971      8.53432   0.1733404   0.1364808 -0.4750196  0.2732967
    ## ...                      ...         ...         ...        ...        ...
    ## ENSG00000107731      6.63160  0.00129427  -0.3588497  0.2327120 -0.6946087
    ## ENSG00000107736      2.06827 -0.47986419  -0.5995966 -0.4171680 -0.0224223
    ## ENSG00000107738      8.62310 -0.14041395  -0.8643571 -0.2631779  0.2137555
    ## ENSG00000107742      2.06485  0.10029416  -0.4432328  0.0207201 -0.0116752
    ## ENSG00000107745      7.45809  0.01255731   0.0417601 -0.1318218 -0.0729778

We can compare the coefficients from *apeglm* with the `"normal"` shrinkage type from the original *DESeq2* paper (2014). This method, which makes use of a Normal-based prior, is no longer the default shrinkage estimator for `lfcShrink`. *apeglm* provides coefficients on the natural log scale, so we must convert to log2 scale by multiplying by `log2(exp(1))`. Note that *DESeq2*'s `lfcShrink` function converts *apeglm* coefficients to the log2 scale internally.

```
system.time({
  res.shr <- lfcShrink(dds, coef=5, type="normal")
})
```

    ##    user  system elapsed
    ##   1.469   0.025   1.494

```
DESeq2.lfc <- res.shr$log2FoldChange
apeglm.lfc <- log2(exp(1)) * fit$map[,5]
```

Here we plot *apeglm* estimators against *DESeq2*:

```
plot(DESeq2.lfc, apeglm.lfc)
abline(0,1)
```

Here we plot *MLE*, *DESeq2* and *apeglm* estimators against the mean of normalized counts:

```
par(mfrow=c(1,3))
lims <- c(-8,8)
hline <- function() abline(h=c(-4:4 * 2),col=rgb(0,0,0,.2))
xlab <- "mean of normalized counts"
plot(res$baseMean, res$log2FoldChange, log="x",
     ylim=lims, main="MLE", xlab=xlab)
hline()
plot(res$baseMean, DESeq2.lfc, log="x",
     ylim=lims, main="DESeq2", xlab=xlab)
hline()
plot(res$baseMean, apeglm.lfc, log="x",
     ylim=lims, main="apeglm", xlab=xlab)
hline()
```

# Specific coefficients

Note that p-values and FSR define different events, and are not on the same scale. An FSR of 0.5 means that the estimated sign is as bad as random guess.

```
par(mfrow=c(1,2),mar=c(5,5,1,1))
plot(res$pvalue, fit$fsr, col="blue",
     xlab="DESeq2 pvalue", ylab="apeglm local FSR",
     xlim=c(0,1), ylim=c(0,.5))
abline(0,1)
plot(-log10(res$pvalue), -log10(fit$fsr),
     xlab="-log10 DESeq2 pvalue", ylab="-log10 apeglm local FSR",
     col="blue")
abline(0,1)
```

The s-value was proposed by [Stephens (2016)], as a statistic giving the aggregate false sign rate for tests with equal or lower s-value than the one considered. We recommend using a lower threshold on s-values than typically used for adjusted p-values, for example one might be interested in sets with 0.01 or 0.005 aggregate FSR.

```
plot(res$padj, fit$svalue, col="blue",
     xlab="DESeq2 padj", ylab="apeglm svalue",
     xlim=c(0,.2), ylim=c(0,.02))
```

More scrutiny can be applied by using an LFC threshold greater than zero, and asking for the probability of a "false-sign-or-small" (FSOS) event: that the effect size is not further from zero in distance than the threshold amount. We can run the `svalue` function on these per-gene probabilities to produce s-values that bound the FSOS rate for sets of genes. By specifying `threshold=log(2) * 1` above, `apeglm` will then output a vector `thresh` in the results list that gives the per-gene probabilities of false-sign-or-small events.

```
s.val <- svalue(fit$thresh)
cols <- ifelse(s.val < .01, "red", "black")
keep <- res$baseMean > 10 # filter for plot
plot(res$baseMean[keep],
     log2(exp(1)) * fit$map[keep,5],
     log="x", col=cols[keep],
     xlab=xlab, ylab="LFC")
abline(h=c(-1,0,1), col=rgb(1,0,0,.5), lwd=2)
```

# Modeling zero-inflated counts

We have created a separate GitHub repository giving an example of how the *apeglm* estimator can be used for Zero-Inflated Negative Binomial data. This approach uses the *zinbwave* method and Bioconductor package to estimate the probability of each zero belonging to the excess zero component. We compare using a Negative Binomial likelihood with the excess zeros down-weighted and using a Zero-Inflated Negative Binomial likelihood. These two approaches with *apeglm* perform similarly but we note that the first approach involves less additional code and is faster to compute.

<https://github.com/mikelove/zinbwave-apeglm>

# Modeling ratios of counts

We also show an short example using an alternative likelihood to the negative binomial. Suppose we have allele-specific counts for n=20 vs 20 samples across 5000 genes. We can define a binomial model and test for the allelic balance across groups of samples.

Here we will *simulate* allele counts from our existing dataset for demonstration. We spike in 10 genes with strong allelic imbalance (instead of an allelic ratio close to 0.5, these will have a ratio of 0.75).

```
library(emdbook)
n <- 20
f <- factor(rep(1:2,each=n))
mu <- ifelse(res$baseMean > 50, res$baseMean, 50)
set.seed(1)
cts <- matrix(rnbinom(nrow(dds)*2*n,
                      mu=mu,
                      size=1/dispersions(dds)),
              ncol=2*n)
theta <- runif(nrow(cts),1,1000)
prob <- rnorm(nrow(cts),.5,.05) # close to 0.5
ase.cts <- matrix(rbetabinom(prod(dim(cts)), prob=prob,
                             size=cts, theta=rep(theta,ncol(cts))),
                  nrow=nrow(cts))
idx <- 1:10
idx2 <- which(f == 2)
theta[idx] <- 1000
prob[idx] <- 0.75
# the spiked in genes have an allelic ratio of 0.75
ase.cts[idx,idx2] <- matrix(rbetabinom(length(idx)*length(idx2), prob=prob[idx],
                                       size=cts[idx,idx2], theta=theta[idx]),
                            nrow=length(idx))
```

We first need to estimate MLE coefficients and standard errors.

One option to run `apeglm` would be to define a beta-binomial likelihood function which uses the total counts as a parameter, and the logit function as a link. And then this function could be provided to the `log.lik` argument.

```
betabinom.log.lik <- function(y, x, beta, param, offset) {
  xbeta <- x %*% beta
  p.hat <- (1+exp(-xbeta))^-1
  dbetabinom(y, prob=p.hat, size=param[-1], theta=param[1], log=TRUE)
}
```

However, `apeglm` has faster C++ implementations for the beta-binomial, which were implemented and tested by Joshua Zitovsky. Here, using `method="betabinCR"` is 3 times faster than using the R-defined `log.lik`, and the `"betabinCR"` method also scales significantly better with more samples and more coefficients. As with the negative binomial, `method="betabinC"` can be used if the standard errors are not needed, and this will be 5 times faster than the R-defined `log.lik` approach on this dataset.

The following code performs two iterations of estimating the MLE coefficients, and then estimating the beta-binomial dispersion.

```
theta.hat <- 100 # rough initial estimate of dispersion
x <- model.matrix(~f)
niter <- 3
system.time({
  for (i in 1:niter) {
    param <- cbind(theta.hat, cts)
    fit.mle <- apeglm(Y=ase.cts, x=x, log.lik=NULL, param=param,
                      no.shrink=TRUE, log.link=FALSE, method="betabinCR")
    theta.hat <- bbEstDisp(success=ase.cts, size=cts,
                           x=x, beta=fit.mle$map,
                           minDisp=.01, maxDisp=5000)
  }
})
```

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ##    user  system elapsed
    ##  10.783   0.060  10.847

We can then plot the MLE estimates over the mean:

```
coef <- 2
xlab <- "mean of normalized counts"
plot(res$baseMean, fit.mle$map[,coef], log="x", xlab=xlab, ylab="log odds")
points(res$baseMean[idx], fit.mle$map[idx,coef], col="dodgerblue", cex=3)
```

Now we run the posterior estimation, including a prior on the second coefficient:

```
mle <- cbind(fit.mle$map[,coef], fit.mle$sd[,coef])
param <- cbind(theta.hat, cts)
system.time({
  fit2 <- apeglm(Y=ase.cts, x=x, log.lik=NULL, param=param,
                 coef=coef, mle=mle, threshold=0.5,
                 log.link=FALSE, method="betabinCR")
})
```

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, possibly due to insufficient
    ## numeric precision

    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value
    ## Warning in betabinGLM(x = x, Y = YNZ, sizes = sizeNZ, thetas = theta, weights =
    ## weightsNZ, : the line search routine failed, unable to sufficiently decrease
    ## the function value

    ##    user  system elapsed
    ##   2.167   0.016   2.183

In the `apeglm` plot, we color in red the genes with a low aggregate probability of false-sign-or-small (FSOS) events (s-value \< .01), where we've again defined "small" on the log odds scale using the `threshold` argument above.

```
par(mfrow=c(1,2))
ylim <- c(-1,1.5)
s.val <- svalue(fit2$thresh) # small-or-false-sign value
plot(res$baseMean, fit.mle$map[,coef], main="MLE",
     log="x", xlab=xlab, ylab="log odds", ylim=ylim)
points(res$baseMean[idx], fit.mle$map[idx,coef], col="dodgerblue", cex=3)
abline(h=0,col=rgb(1,0,0,.5))
cols <- ifelse(s.val < .01, "red", "black")
plot(res$baseMean, fit2$map[,coef], main="apeglm",
     log="x", xlab=xlab, ylab="log odds", col=cols, ylim=ylim)
points(res$baseMean[idx], fit2$map[idx,coef], col="dodgerblue", cex=3)
abline(h=0,col=rgb(1,0,0,.5))
```

```
logit <- function(x) log(x/(1-x))
logit(.75)
```

    ## [1] 1.098612

```
table(abs(logit(prob[s.val < .01])) > .5)
```

    ##
    ## TRUE
    ##   10

# Session Info

```
sessionInfo()
```

    ## R version 4.5.1 Patched (2025-08-23 r88802)
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
    ##  [1] emdbook_1.3.14              apeglm_1.32.0
    ##  [3] DESeq2_1.50.0               airway_1.29.0
    ##  [5] SummarizedExperiment_1.40.0 Biobase_2.70.0
    ##  [7] GenomicRanges_1.62.0        Seqinfo_1.0.0
    ##  [9] IRanges_2.44.0              S4Vectors_0.48.0
    ## [11] BiocGenerics_0.56.0         generics_0.1.4
    ## [13] MatrixGenerics_1.22.0       matrixStats_1.5.0
    ##
    ## loaded via a namespace (and not attached):
    ##  [1] sass_0.4.10         SparseArray_1.10.0  lattice_0.22-7
    ##  [4] digest_0.6.37       magrittr_2.0.4      evaluate_1.0.5
    ##  [7] grid_4.5.1          RColorBrewer_1.1-3  mvtnorm_1.3-3
    ## [10] fastmap_1.2.0       plyr_1.8.9          jsonlite_2.0.0
    ## [13] Matrix_1.7-4        scales_1.4.0        numDeriv_2016.8-1.1
    ## [16] codetools_0.2-20    jquerylib_0.1.4     abind_1.4-8
    ## [19] cli_3.6.5           bbmle_1.0.25.1      rlang_1.1.6
    ## [22] XVector_0.50.0      cachem_1.1.0        DelayedArray_0.36.0
    ## [25] yaml_2.3.10         S4Arrays_1.10.0     tools_4.5.1
    ## [28] parallel_4.5.1      BiocParallel_1.44.0 coda_0.19-4.1
    ## [31] bdsmatrix_1.3-7     dplyr_1.1.4         ggplot2_4.0.0
    ## [34] locfit_1.5-9.12     vctrs_0.6.5         R6_2.6.1
    ## [37] lifecycle_1.0.4     MASS_7.3-65         pkgconfig_2.0.3
    ## [40] bslib_0.9.0         pillar_1.11.1       gtable_0.3.6
    ## [43] glue_1.8.0          Rcpp_1.1.0          tidyselect_1.2.1
    ## [46] tibble_3.3.0        xfun_0.53           knitr_1.50
    ## [49] dichromat_2.0-0.1   farver_2.1.2        htmltools_0.5.8.1
    ## [52] rmarkdown_2.30      compiler_4.5.1      S7_0.2.0

# References

Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." *Genome Biology* 15 (12): 550. <https://doi.org/10.1186/s13059-014-0550-8>.

Stephens, Matthew. 2016. "False Discovery Rates: A New Deal." *Biostatistics* 18 (2). <https://doi.org/10.1093/biostatistics/kxw041>.

Zhu, Anqi, Joseph G. Ibrahim, and Michael I. Love. 2018. "Heavy-Tailed Prior Distributions for Sequence Count Data: Removing the Noise and Preserving Large Differences." *Bioinformatics*. <https://doi.org/10.1093/bioinformatics/bty895>.
