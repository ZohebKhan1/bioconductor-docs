# Using a BioMart other than Ensembl

Steffen Durinck, Wolfgang Huber, Mike Smith

#### 13 March 2026

#### Package

biomaRt 2.66.2

# Contents

-   [[1] Introduction](#introduction)
-   [[2] Using a BioMart other than Ensembl](#using-a-biomart-other-than-ensembl)
    -   [[2.1] Wormbase](#wormbase)
    -   [[2.2] Phytozome](#phytozome)
        -   [[2.2.1] Version 12](#version-12)
        -   [[2.2.2] Version 13](#version-13)
-   [[3] Session Info](#session-info)

# [1] Introduction

In recent years a wealth of biological data has become available in public data repositories. Easy access to these valuable data resources and firm integration with data analysis is needed for comprehensive bioinformatics data analysis. The *[biomaRt](https://bioconductor.org/packages/3.22/biomaRt)* package, provides an interface to a growing collection of databases implementing the [BioMart software suite](https://www.ensembl.org/info/data/biomart/index.html). The package enables retrieval of large amounts of data in a uniform way without the need to know the underlying database schemas or write complex SQL queries. Examples of BioMart databases are Ensembl, Uniprot and HapMap. These major databases give *[biomaRt](https://bioconductor.org/packages/3.22/biomaRt)* users direct access to a diverse set of data and enable a wide range of powerful online queries from R.

# [2] Using a BioMart other than Ensembl

There are a small number of non-Ensembl databases that offer a BioMart interface to their data. The *[biomaRt](https://bioconductor.org/packages/3.22/biomaRt)* package can be used to access these in a very similar fashion to Ensembl. The majority of *[biomaRt](https://bioconductor.org/packages/3.22/biomaRt)* functions will work in the same manner, but the construction of the initial Mart object requires slightly more setup. In this section we demonstrate the setting requires to query [Wormbase ParaSite](https://parasite.wormbase.org/index.html) and [Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html). First we need to load *[biomaRt](https://bioconductor.org/packages/3.22/biomaRt)*.

``` r
library(biomaRt)
```

## [2.1] Wormbase

To demonstrate the use of the *[biomaRt](https://bioconductor.org/packages/3.22/biomaRt)* package with non-Ensembl databases the next query is performed using the Wormbase ParaSite BioMart. In this example, we use the `listMarts()` function to find the name of the available marts, given the URL of Wormbase. We use this to connect to Wormbase BioMart using the `useMart()` function.11[[1] Note that we use the `https` address and must provide the port as `443`. Queries to WormBase will fail without these options.]

``` r
listMarts(host = "parasite.wormbase.org")
```

    ## Error in `req_perform()`:
    ## ! Failed to perform HTTP request.
    ## Caused by error in `curl::curl_fetch_memory()`:
    ## ! Unsupported protocol [parasite.wormbase.org]:
    ## Received HTTP/0.9 when not allowed

``` r
wormbase <- useMart(
  biomart = "parasite_mart",
  host = "https://parasite.wormbase.org",
  port = 443
)
```

We can then use functions described earlier in this vignette to find and select the gene dataset, and print the first 6 available attributes and filters. Then we use a list of gene names as filter and retrieve associated transcript IDs and the transcript biotype.

``` r
listDatasets(wormbase)
```

    ##     dataset          description version
    ## 1 wbps_gene All Species (WBPS19)      19

``` r
wormbase <- useDataset(mart = wormbase, dataset = "wbps_gene")
head(listFilters(wormbase))
```

    ##                  name     description
    ## 1     species_id_1010          Genome
    ## 2 nematode_clade_1010  Nematode Clade
    ## 3     chromosome_name Chromosome name
    ## 4               start           Start
    ## 5                 end             End
    ## 6              strand          Strand

``` r
head(listAttributes(wormbase))
```

    ##                      name        description         page
    ## 1          species_id_key      Internal Name feature_page
    ## 2    production_name_1010     Genome project feature_page
    ## 3       display_name_1010        Genome name feature_page
    ## 4        taxonomy_id_1010        Taxonomy ID feature_page
    ## 5 assembly_accession_1010 Assembly accession feature_page
    ## 6     nematode_clade_1010     Nematode clade feature_page

``` r
getBM(
  attributes = c(
    "external_gene_id",
    "wbps_transcript_id",
    "transcript_biotype"
  ),
  filters = "gene_name",
  values = c("unc-26", "his-33"),
  mart = wormbase
)
```

    ##   external_gene_id wbps_transcript_id transcript_biotype
    ## 1           his-33         F17E9.13.1     protein_coding
    ## 2           unc-26          JC8.10a.1     protein_coding
    ## 3           unc-26          JC8.10b.1     protein_coding
    ## 4           unc-26          JC8.10c.1     protein_coding
    ## 5           unc-26          JC8.10d.1     protein_coding

## [2.2] Phytozome

### [2.2.1] Version 12

The Phytozome 12 BioMart was [retired](https://jgi.doe.gov/more-intuitive-phytozome-interface/) in August 2021 and can not longer be accessed.

### [2.2.2] Version 13

Version 13 of Phytozome can be found at https://phytozome-next.jgi.doe.gov/ and if you wish to query that version the URL used to create the Mart object must reflect that.

``` r
phytozome_v13 <- useMart(
  biomart = "phytozome_mart",
  dataset = "phytozome",
  host = "https://phytozome-next.jgi.doe.gov"
)
```

Once this is set up the usual *[biomaRt](https://bioconductor.org/packages/3.22/biomaRt)* functions can be used to interrogate the database options and run queries.

``` r
getBM(
  attributes = c("organism_name", "gene_name1"),
  filters = "gene_name_filter",
  values = "82092",
  mart = phytozome_v13
)
```

    ##          organism_name gene_name1
    ## 1 Smoellendorffii_v1.0      82092

# [3] Session Info

``` r
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
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB
    ##  [4] LC_COLLATE=C               LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C
    ## [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
    ##
    ## time zone: America/New_York
    ## tzcode source: system (glibc)
    ##
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base
    ##
    ## other attached packages:
    ## [1] biomaRt_2.66.2   BiocStyle_2.38.0
    ##
    ## loaded via a namespace (and not attached):
    ##  [1] rappdirs_0.3.4       sass_0.4.10          generics_0.1.4       xml2_1.5.2
    ##  [5] RSQLite_2.4.6        stringi_1.8.7        hms_1.1.4            digest_0.6.39
    ##  [9] magrittr_2.0.4       evaluate_1.0.5       bookdown_0.46        fastmap_1.2.0
    ## [13] blob_1.3.0           jsonlite_2.0.0       progress_1.2.3       AnnotationDbi_1.72.0
    ## [17] DBI_1.3.0            BiocManager_1.30.27  httr_1.4.8           purrr_1.2.1
    ## [21] Biostrings_2.78.0    httr2_1.2.2          jquerylib_0.1.4      cli_3.6.5
    ## [25] rlang_1.1.7          crayon_1.5.3         XVector_0.50.0       dbplyr_2.5.2
    ## [29] Biobase_2.70.0       bit64_4.6.0-1        withr_3.0.2          cachem_1.1.0
    ## [33] yaml_2.3.12          otel_0.2.0           tools_4.5.2          memoise_2.0.1
    ## [37] dplyr_1.2.0          filelock_1.0.3       BiocGenerics_0.56.0  curl_7.0.0
    ## [41] png_0.1-8            vctrs_0.7.1          R6_2.6.1             stats4_4.5.2
    ## [45] BiocFileCache_3.0.0  lifecycle_1.0.5      Seqinfo_1.0.0        KEGGREST_1.50.0
    ## [49] stringr_1.6.0        S4Vectors_0.48.0     IRanges_2.44.0       bit_4.6.0
    ## [53] pkgconfig_2.0.3      pillar_1.11.1        bslib_0.10.0         glue_1.8.0
    ## [57] xfun_0.56            tibble_3.3.1         tidyselect_1.2.1     knitr_1.51
    ## [61] htmltools_0.5.9      rmarkdown_2.30       compiler_4.5.2       prettyunits_1.2.0

``` r
warnings()
```
