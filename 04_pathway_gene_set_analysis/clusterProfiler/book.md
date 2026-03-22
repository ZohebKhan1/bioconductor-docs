# Welcome

#  Introduction

##  Motivation

The book is meant as a guide for mining biological knowledge to elucidate or interpret molecular mechanisms using a suite of R packages, including , , , , ,  and . Hence, if you are starting to read this book, we assume you have a working knowledge of how to use R.

##  Citation

If you use the software suite in published research, please cite the most appropriate paper(s) from this list:

1.  **G Yu**`\*`. [Background bias in functional enrichment analysis: Insights from clusterProfiler](https://doi.org/10.59717/j.xinn-life.2025.100181). ***The Innovation Life***. 2026, 4(1):100181.
2.  **G Yu**`\*`. [Thirteen years of clusterProfiler](https://doi.org/10.1016/j.xinn.2024.100722). ***The Innovation***. 2024, 5(6):100722.
3.  S Xu`\#`, E Hu`\#`, Y Cai`\#`, Z Xie`\#`, X Luo`\#`, L Zhan, W Tang, Q Wang, B Liu, R Wang, W Xie, T Wu, L Xie, **G Yu**`\*`. [Using clusterProfiler to characterise Multi-Omics Data](https://www.nature.com/articles/s41596-024-01020-z). ***Nature Protocols***. 2024, 19(11):3292-3320.
4.  Q Wang`\#`, M Li`\#`, T Wu, L Zhan, L Li, M Chen, W Xie, Z Xie, E Hu, S Xu, **G Yu**`\*`. [Exploring epigenomic datasets by ChIPseeker](https://doi.org/10.1002/cpz1.585). ***Current Protocols***, 2022, 2(10): e585.
5.  T Wu`\#`, E Hu`\#`, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo`\*`, **G Yu**`\*`. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. ***The Innovation***. 2021, 2(3):100141. doi: [10.1016/j.xinn.2021.100141](https://doi.org/10.1016/j.xinn.2021.100141)
6.  **G Yu**^\*^. Gene Ontology Semantic Similarity Analysis Using GOSemSim. In: Kidder B. (eds) Stem Cell Transcriptional Networks. ***Methods in Molecular Biology***. 2020, 2117:207-215. Humana, New York, NY. doi: [10.1007/978-1-0716-0301-7_11](https://doi.org/10.1007/978-1-0716-0301-7_11)
7.  **G Yu**^\*^. Using meshes for MeSH term enrichment and semantic analyses. ***Bioinformatics***. 2018, 34(21):3766--3767. doi: [10.1093/bioinformatics/bty410](https://doi.org/10.1093/bioinformatics/bty410)
8.  **G Yu**, QY He^\*^. ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. ***Molecular BioSystems***. 2016, 12(2):477-479. doi: [10.1039/C5MB00663E](https://doi.org/10.1039/C5MB00663E)
9.  **G Yu**^\*^, LG Wang, and QY He^\*^. ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. ***Bioinformatics***. 2015, 31(14):2382-2383. doi: [10.1093/bioinformatics/btv145](https://doi.org/10.1093/bioinformatics/btv145)
10. **G Yu**^\*^, LG Wang, GR Yan, QY He^\*^. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. ***Bioinformatics***. 2015, 31(4):608-609. doi: [10.1093/bioinformatics/btu684](https://doi.org/10.1093/bioinformatics/btu684)
11. **G Yu**, LG Wang, Y Han and QY He^\*^. clusterProfiler: an R package for comparing biological themes among gene clusters. ***OMICS: A Journal of Integrative Biology***. 2012, 16(5):284-287. doi: [10.1089/omi.2011.0118](https://doi.org/10.1089/omi.2011.0118)
12. **G Yu**, F Li, Y Qin, X Bo^\*^, Y Wu, S Wang^\*^. GOSemSim: an R package for measuring semantic similarity among GO terms and gene products. ***Bioinformatics***. 2010, 26(7):976-978. doi: [10.1093/bioinformatics/btq064](https://doi.org/10.1093/bioinformatics/btq064)

##  Book structure

-   Part 1 (Semantic similarity analysis) describes ,  and  packages for measuring semantic similarity of genes or gene products based on Gene Ontology, Disease Ontology and Medical Subject Headings.
-   Part 2 (Enrichment analysis) introduces over-reprensentation analysis and gene set enrichment analysis using  (supports GO, KEGG, MSigDb, WikiPathway, and many others via universal interface),  (DO, Disease-Gene Network, Network of Cancer Genes),  (MeSH), and  (Reactome pathway). Functional enrichment analysis of Genomic coordination is supported via  and comparison among multiple conditions is also supported by . We implemented a number of visualization methods in the  package to help users to interpret their results.
-   Part 3 (Miscellaneous topics) describes useful utilities including translating gene IDs and manipulating enrichment results.

##  Want to help?

The book's source code is hosted on GitHub, at <https://github.com/YuLab-SMU/biomedical-knowledge-mining-book>. Any feedback on the book is very welcome. Feel free to [open an issue](https://github.com/YuLab-SMU/biomedical-knowledge-mining-book/issues/new) on GitHub or send me a pull request if you notice typos or other issues (I'm not a native English speaker ;) ).

---

# Semantic Similarity Analysis

Functional similarity of gene products can be estimated by controlled biological vocabularies, such as Gene Ontology (GO), Disease Ontology (DO) and Medical Subject Headings (MeSH).

Four methods including Resnik [@philip_semantic_1999], Jiang [@jiang_semantic_1997], Lin [@lin_information-theoretic_1998] and Schlicker [@schlicker_new_2006] have been presented to determine the semantic similarity of two GO terms based on the annotation statistics of their common ancestor terms. Wang [@wang_new_2007] proposed a method to measure the similarity based on the graph structure of GO. Each of these methods has its own advantages and weaknesses and can be applied to other ontologies that have similar structure (i.e. directed acyclic graph).

## Information content-based methods

Four methods proposed by Resnik [@philip_semantic_1999], Jiang [@jiang_semantic_1997], Lin [@lin_information-theoretic_1998] and Schlicker [@schlicker_new_2006] are information content (IC) based, which depend on the frequencies of two GO terms involved and that of their closest common ancestor term in a specific corpus of GO annotations. The information content of a GO term is computed by the negative log probability of the term occurring in GO corpus. A rarely used term contains a greater amount of information.

The frequency of a term t is defined as:

$$p(t) = \frac{n_{t'}}{N} | t' \in \left\{t, \; children\: of\: t \right\}$$

where $n_{t'}$ is the number of term $t'$, and $N$ is the total number of terms in GO corpus.

Thus the information content is defined as:

$$IC(t) = -\log(p(t))$$

As GO allows multiple parents for each concept, two terms can share parents through multiple paths. IC-based methods calculate the similarity of two GO terms based on the information content of their closest common ancestor term, which is also called the most informative common ancestor (MICA).

### Resnik method

The Resnik method is defined as:

$$sim_{Resnik}(t_1,t_2) = IC(MICA)$$

### Lin method

The Lin method is defined as:

$$sim_{Lin}(t_1,t_2) = \frac{2IC(MICA)}{IC(t_1)+IC(t_2)}$$

### Rel method

The Relevance method, which was proposed by Schlicker, combines Resnik's and Lin's methods and is defined as:

$$sim_{Rel}(t_1,t_2) = \frac{2IC(MICA)(1-p(MICA))}{IC(t_1)+IC(t_2)}$$

### Jiang method

The Jiang and Conrath's method is defined as:

$$sim_{Jiang}(t_1,t_2) = 1-\min(1, IC(t_1) + IC(t_2) - 2IC(MICA))$$

## Graph-based methods

Graph-based methods use the topology of the GO graph structure to compute semantic similarity. Formally, a GO term A can be represented as $DAG_{A}=(A,T_{A},E_{A})$ where $T_{A}$ is the set of GO terms in $DAG_{A}$, including term A and all of its ancestor terms in the GO graph, and $E_{A}$ is the set of edges connecting the GO terms in $DAG_{A}$.

### Wang method

To encode the semantics of a GO term in a measurable format to enable quantitative comparison, Wang[@wang_new_2007] first defined the semantic value of term A as the aggregate contribution of all terms in $DAG_{A}$ to the semantics of term A, where terms closer to term A in $DAG_{A}$ contribute more to its semantics. Thus, the contribution of a GO term $t$ to the semantics of GO term $A$ is defined as the S-value of GO term $t$ related to term $A$.

For any of term $t$ in $DAG_{A}$, its S-value related to term $A$, $S_{A}(\textit{t})$ is defined as:

$$\left\{\begin{array}{l} S_{A}(A)=1 \\ S_{A}(\textit{t})=\max\{w_{e} \times S_{A}(\textit{t}') | \textit{t}' \in children \: of(\textit{t}) \} \; if \: \textit{t} \ne A \end{array} \right.$$

where $w_{e}$ is the semantic contribution factor for edge $e \in E_{A}$ linking term $t$ with its child term $t'$. Term $A$'s contribution to itself is defined as 1. After obtaining the S-values for all terms in $DAG_{A}$, the semantic value of GO term A, $SV(A)$, is calculated as:

$$SV(A)=\displaystyle\sum_{t \in T_{A}} S_{A}(t)$$

Thus given two GO terms A and B, the semantic similarity between these two terms is defined as:

$$sim_{Wang}(A, B) = \frac \cap T_{B}}{S_{A}(t) + S_{B}(t)}}{SV(A) + SV(B)}$$

where $S_{A}(\textit{t})$ is the S-value of GO term $t$ related to term $A$ and $S_{B}(\textit{t})$ is the S-value of GO term $t$ related to term $B$.

This method proposed by Wang [@wang_new_2007] determines the semantic similarity of two GO terms based on both the locations of these terms in the GO graph and their relations with their ancestor terms.

## Combine methods

Since a gene product can be annotated by multiple GO terms, semantic similarity among gene products needs to be aggregated from different semantic similarity scores of multiple GO terms associated with genes, including `max`, `avg`, `rcmax` and `BMA`.

### max

The `max` method calculates the maximum semantic similarity score over all pairs of GO terms between these two GO term sets.

$$sim_{max}(g_1, g_2) = \displaystyle\max_{1 \le i \le m, 1 \le j \le n} sim(go_{1i}, go_{2j})$$

### avg

The `avg` calculates the average semantic similarity score over all pairs of GO terms.

$$sim_{avg}(g_1, g_2) = \frac^m\sum_{j=1}^nsim(go_{1i}, go_{2j})}{m \times n}$$

### rcmax

Similarities among two sets of GO terms form a matrix, the `rcmax` method uses the maximum of `RowScore` and `ColumnScore`, where `RowScore` (or `ColumnScore`) is the average of maximum similarity on each row (or column).

$$sim_{rcmax}(g_1, g_2) = \max(\frac^m \max_{1 \le j \le n} sim(go_{1i}, go_{2j})}{m},\frac^n \max_{1 \le i \le m} sim(go_{1i},go_{2j})}{n})$$

### BMA

The `BMA` method, used the **B**est-**M**atch **A**verage strategy, calculates the average of all maximum similarities on each row and column, and is defined as:

$$sim_{BMA}(g_1, g_2) = \frac^m \max_{1 \le j \le n}sim(go_{1i}, go_{2j}) + \displaystyle\sum_{1=j}^n \max_{1 \le i \le m}sim(go_{1i}, go_{2j})} {m+n}$$

## Summary

The idea behind semantic similarity measurement is the notion that genes with similar functions should have similar annotation vocabularies and be closely related in the ontology structure. Measuring similarity is critical for expanding knowledge, since similar objects tend to exhibit similar behaviors, which supports many bioinformatics applications to infer gene/protein functions, miRNA functions, genetic interactions, protein-protein interactions, miRNA-mRNA interactions, and cellular localization.

We developed several Bioconductor packages, including  [@yu2010; @yu_gosemsim_2020] for computing semantic similarity among GO terms, sets of GO terms, gene products and gene clusters (see also [Chapter 2](#GOSemSim)),  [@yu_dose_2015] for Disease Ontology (DO) (see also [Chapter 3](#DOSE-semantic-similarity)) and  [@yu_meshes_2018] that based on Medical Subject Headings (MeSH) (see also [Chapter 4](#meshes-semantic-similarity)).

---

# Enrichment Analysis

    ## Terminology

    ### Gene sets and pathway

    A gene set is an unordered collection of genes that are functionally related. A pathway can be interpreted as a gene set by ignoring functional relationships among genes.

    ### Gene Ontology (GO)

    [Gene Ontology](http://www.geneontology.org/) defines concepts/classes used to describe gene function, and relationships between these concepts. It classifies functions along three aspects:

    + MF: Molecular Function
      - molecular activities of gene products
    + CC: Cellular Component
      - where gene products are active
    + BP: Biological Process
      - pathways and larger processes made up of the activities of multiple gene products

    GO terms are organized in a directed acyclic graph, where edges between terms represent parent-child relationship.

    ### Kyoto Encyclopedia of Genes and Genomes (KEGG)

    [KEGG](https://www.genome.jp/kegg/) is a collection of manually drawn pathway maps representing molecular interaction and reaction networks. These pathways cover a wide range of biochemical processes that can be divided into [7 broad categories](https://www.genome.jp/kegg/pathway.html):

    1. Metabolism
    2. Genetic information processing
    3. Environmental information processing
    4. Cellular processes
    5. Organismal systems
    6. Human diseases
    7. Drug development.

    <!-- https://pathview.uncc.edu/data/khier.tsv -->

    ### Other gene sets

    GO and KEGG are the most frequently used for functional analysis. They are typically the first choice because of their long-standing curation and availability for a wide range of species.

    Other gene sets include but are not limited to Disease Ontology ([DO](http://disease-ontology.org/)), Disease Gene Network ([DisGeNET](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4397996/)), [wikiPathways](https://www.wikipathways.org), Molecular Signatures Database ([MSigDb](http://software.broadinstitute.org/gsea/msigdb)).

    ## Over Representation Analysis

    Over Representation Analysis (ORA) [@boyle2004] is a widely used approach to determine whether known biological functions or processes are over-represented (= enriched) in an experimentally-derived gene list, *e.g.* a list of differentially expressed genes (DEGs).

    The _p_-value can be calculated by hypergeometric distribution.

    $p = 1 - \displaystyle\sum_{i = 0}^{k-1}\frac{{M \choose i}{{N-M} \choose {n-i}}} {{N \choose n}}$

    In this equation, `N` is the total number of genes in the background distribution,
    `M` is the number of genes within that distribution that are annotated (either directly or indirectly) to the gene set of interest,
    `n` is the size of the list of genes of interest and `k` is the number of genes within that list which are annotated to the gene set. The background distribution by default is all the genes that have annotation. _P_-values should be adjusted for [multiple comparison](https://en.wikipedia.org/wiki/Multiple_comparisons_problem).

    **Example:** Suppose we have 17,980 genes detected in a Microarray study and 57 genes were differentially expressed. From these, 2641 are annotated to gene set of interest. Among the differentially expressed genes, 28 belong to the gene set^[example adopted from <https://guangchuangyu.github.io/2012/04/enrichment-analysis/>].

    d <- data.frame(non_DE_gene=c(2613, 15310), DE_gene=c(28, 29))
    row.names(d) <- c("In_category", "not_in_category")
    d

Whether the gene set of interest is significantly over-represented in the differentially expressed genes can be assessed using a hypergeometric distribution. This corresponds to a one-sided version of Fisher's exact test.

`{r} #| label: fisher-test fisher.test(d, alternative = "greater")`

## Gene Set Enrichment Analysis

A common approach to analyzing gene expression profiles is identifying differentially expressed genes that are deemed interesting. The [ORA enrichment analysis](#ora-algorithm) is based on these differentially expressed genes. This approach will find genes where the difference is large and will fail where the difference is small, but evidenced in coordinated way in a set of related genes. Gene Set Enrichment Analysis (GSEA)[@subramanian_gene_2005] directly addresses this limitation. All genes can be used in GSEA; GSEA aggregates the per gene statistics across genes within a gene set, therefore making it possible to detect situations where all genes in a predefined set change in a small but coordinated way. This is important since it is likely that many relevant phenotypic differences are manifested by small but consistent changes in a set of genes.

Genes are ranked based on their phenotypes. Given a priori defined set of gene *S* (e.g., genes sharing the same *DO* category), the goal of GSEA is to determine whether the members of *S* are randomly distributed throughout the ranked gene list (*L*) or primarily found at the top or bottom.

There are three key elements of the GSEA method:

-   Calculation of an Enrichment Score.
    -   The enrichment score (*ES*) represents the degree to which a set *S* is over-represented at the top or bottom of the ranked list *L*. The score is calculated by walking down the list *L*, increasing a running-sum statistic when we encounter a gene in *S* and decreasing when it is not encountered. The magnitude of the increment depends on the gene statistics (e.g., correlation of the gene with phenotype). The *ES* is the maximum deviation from zero encountered in the random walk; it corresponds to a weighted Kolmogorov-Smirnov(KS)-like statistic [@subramanian_gene_2005].
-   Estimation of Significance Level of *ES*.
    -   The *p*-value of the *ES* is calculated using a permutation test. Specifically, we permute the gene labels of the gene list *L* and recompute the *ES* of the gene set for the permutated data, which generates a null distribution for the *ES*. The *p*-value of the observed *ES* is then calculated relative to this null distribution.
-   Adjustment for Multiple Hypothesis Testing.
    -   When the entire gene sets are evaluated, the estimated significance level is adjusted to account for multiple hypothesis testing and also *q*-values are calculated for FDR control.

We implemented the GSEA algorithm proposed by Subramanian [@subramanian_gene_2005]. Alexey Sergushichev implemented an algorithm for fast GSEA calculation in the  [@korotkevich_fast_2019] package. In our packages (, ,  and ), users can use the GSEA algorithm implemented in `DOSE` or `fgsea` by specifying the parameter `by="DOSE"` or `by="fgsea"`. By default, the `fgsea` method will be used since it is much faster.

## Leading edge analysis and core enriched genes

Leading edge analysis reports `Tags` to indicate the percentage of genes contributing to the enrichment score, `List` to indicate where in the list the enrichment score is attained and `Signal` for enrichment signal strength.

It would also be very interesting to get the core enriched genes that contribute to the enrichment. Our packages (, ,  and ) support leading edge analysis and report core enriched genes in GSEA analysis.

## References

---

# GO enrichment analysis

    GO comprises three orthogonal ontologies, i.e. molecular function (MF), biological
    process (BP), and cellular component (CC).

    opts_chunk$set(message=FALSE, warning=FALSE, eval=TRUE, echo=TRUE, cache=TRUE)

## Supported organisms

GO analyses (`groupGO()`, `enrichGO()` and `gseGO()`) support organisms that have an `OrgDb` object available (see also [session 2.2](#gosemsim-supported-organisms)).

If a user has GO annotation data (in a `data.frame` format with the first column as gene ID and the second column as GO ID), they can use the `enricher()` and `GSEA()` functions to perform an over-representation test and gene set enrichment analysis. Users can also read GO annotation from GMT files using the `read.gmt()` function, which parses the file into a `data.frame` suitable for these functions.

If the genes are annotated by direct annotation, they should also be annotated by their ancestor GO nodes (indirect annotation). If a user only has direct annotation, they can pass their annotation to the `buildGOmap` function, which will infer indirect annotation and generate a `data.frame` that is suitable for both `enricher()` and `GSEA()`.

## Direct vs indirect GO annotation

One common issue users encounter is the discrepancy between direct gene-to-GO mapping and the results from enrichment analysis. For instance, a user might find that a specific GO term is enriched, but when they check the input gene list against the direct annotations, the term doesn't appear to be associated with many genes.

This happens because GO structure is a Directed Acyclic Graph (DAG). If a gene is annotated to a specific term (e.g., "glucose metabolic process"), it is implicitly annotated to all its ancestor terms (e.g., "metabolic process"). `enrichGO` and `gseGO` account for this **indirect annotation** by propagating annotations up the graph.

If you are providing your own annotation file (e.g., using `enricher` or `GSEA`), you must ensure that it includes these indirect annotations. If your annotation file only contains direct annotations, you can use the `buildGOmap()` function to generate the full annotation set (direct + indirect) before performing enrichment analysis.

# Use the result for enrichment

enricher(gene, TERM2GENE = gomap_with_ancestors, ...)

    This step is crucial for accurate GO enrichment analysis when not using `OrgDb` objects.

    ## GO classification

    In , the `groupGO()` function is designed for gene classification based on GO distribution at a specific level. Here we use the dataset `geneList` provided by .

    data(geneList, package="DOSE")
    gene <- names(geneList)[abs(geneList) > 2]

    # Entrez gene ID
    head(gene)

    ggo <- groupGO(gene     = gene,
                   OrgDb    = org.Hs.eg.db,
                   ont      = "CC",
                   level    = 3,
                   readable = TRUE)

    head(ggo)

The `gene` parameter is a vector of gene IDs (can be any ID type that is supported by the corresponding `OrgDb`, see also @sec-id-convert). If `readable` is set to `TRUE`, the input gene IDs will be converted to gene symbols.

## GO over-representation analysis

The  package implements `enrichGO()` for gene ontology over-representation test.

`{r} #| label: enrichGO-example ego <- enrichGO(gene          = gene,                 universe      = names(geneList),                 OrgDb         = org.Hs.eg.db,                 ont           = "CC",                 pAdjustMethod = "BH",                 pvalueCutoff  = 0.01,                 qvalueCutoff  = 0.05,         readable      = TRUE) head(ego)`

Any gene ID type that is supported in `OrgDb` can be directly used in GO analyses. Users need to specify the `keyType` parameter to specify the input gene ID type.

However, the GO annotation in `OrgDb` is indexed by ENTREZID. If the input ID type is not ENTREZID, `clusterProfiler` will first map the input IDs to ENTREZID using the `bitr()` function. This mapping step may lead to a loss of information, as some IDs might not have a corresponding ENTREZID or multiple IDs might map to the same ENTREZID (many-to-one mapping).

To address this issue, `enrichGO` and `gseGO` perform the mapping and use the mapped ENTREZIDs for enrichment analysis. The results are then mapped back to the original input ID type. This ensures that the gene counts and gene ratios in the final result accurately reflect the input IDs, avoiding discrepancies caused by the intermediate mapping step.

For example, if we use ENSEMBL IDs as input:

ego2 \<- enrichGO(gene = gene.df\$ENSEMBL, OrgDb = org.Hs.eg.db, keyType = 'ENSEMBL', ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05) head(ego2, 3)

    Gene IDs can be mapped to gene Symbols by using the parameter `readable=TRUE` or [`setReadable()` function](#sec-setReadable).

    ## GO Gene Set Enrichment Analysis

    The  package provides the `gseGO()` function for [gene set enrichment analysis](#gsea-algorithm) using gene ontology.

    ego3 <- gseGO(geneList     = geneList,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "CC",
                  minGSSize    = 100,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)

The format of input data, `geneList`, is documented in the @sec-id-convert. Beware that only gene sets with size in `[minGSSize, maxGSSize]` will be tested.

## GO analysis for non-model organisms

Both the `enrichGO()` and `gseGO()` functions require an `OrgDb` object as the background annotation. For organisms that don't have `OrgDb` provided by `Bioconductor`, users can query one (if available) online via . If there is no `OrgDb` available, users can obtain GO annotation from other sources, e.g. from , or annotate the genes using [Blast2GO](https://www.blast2go.com/) or the [Trinotate](https://rnabio.org/module-07-trinotate/0007/02/01/Trinotate/) pipeline. Then the `enricher()` or `GSEA()` functions can be used to perform GO analysis for these organisms, similar to the examples using WikiPathways and MSigDB. Another solution is to create an `OrgDb` on your own using the  package.

Here is an example of querying GO annotation from Ensembl using .

Alternatively, you can use `AnnotationHub` to query and retrieve `OrgDb` objects for a wide range of organisms, including non-model species.

**Example: GO enrichment analysis for Maize (Zea mays)**

# Query for Zea mays OrgDb

query(hub, "zea") \# Retrieve the OrgDb (e.g., AH55736, assume it's the correct record) maize \<- hub\[\['AH55736'\]\]

# Check keys

length(keys(maize)) columns(maize)

# Perform enrichment analysis

# Assuming 'sample_genes' is a vector of Entrez IDs

res \<- enrichGO(sample_genes, OrgDb=maize, pvalueCutoff=0.05, qvalueCutoff=0.05)

    This approach allows you to use the standard `enrichGO` workflow with organisms not included in the default Bioconductor `OrgDb` packages.

    ## Visualization of GO enrichment results

    The `clusterProfiler` package provides several visualization methods to help interpret enrichment results, including `barplot`, `dotplot`, `cnetplot`, `emapplot`, `goplot` and `plotGOgraph`.

    Please refer to the @sec-enrichplot chapter for details.

    ## Filtering and simplifying GO terms

    GO enrichment analysis often results in many redundant terms due to the hierarchical structure of GO. Parent terms and child terms often share a large proportion of genes, leading to multiple significant results that represent similar biological processes.

    ### Removing specific terms or levels

    The `dropGO()` function can be used to remove specific GO terms or GO levels from the results. The `enrichGO()` function tests the whole GO corpus and enriched result may contains very general terms. If users want to restrict the result to a specific GO level (e.g., level 3 or 4), they can use the `gofilter()` function. Both functions work with results obtained from `enrichGO()`, `gseGO()` and `compareCluster()`.

    ### Simplifying enriched GO terms

    To reduce redundancy, `clusterProfiler` provides a `simplify` method. It uses semantic similarity (via `GOSemSim`) to calculate the similarity between enriched terms and removes highly similar terms, keeping the most significant one.

    The `simplify` method applies `select_fun` (which can be a user-defined function) to the feature specified by `by` (e.g., p.adjust) to select one representative term from redundant terms (which have similarity higher than `cutoff`). This results in a cleaner and more interpretable visualization.

    In @fig-simplify-emapplot(A), we can found that there are many redundant terms form a highly condense network. After removing redundant terms using the `simplify()` method, the result is more clear to view the whole story.

    data(geneList, package="DOSE")
    de <- names(geneList)[abs(geneList) > 2]
    bp <- enrichGO(de, ont="BP", OrgDb = 'org.Hs.eg.db')
    bp <- pairwise_termsim(bp)
    bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
    p1 <- emapplot(bp)
    p2 <- emapplot(bp2)
    plot_list(p1, p2, ncol=2, tag_levels = 'A')

Alternatively, users can use slim version of GO and use the `enricher()` or `gseGO()` functions to analyze.

**Note on `enricher` results**: The `simplify()` method is designed for `enrichGO()` results where the ontology (BP, CC, or MF) is known. If you use `enricher()` to perform GO analysis (e.g., with custom annotation), the result object does not strictly define the ontology. To use `simplify()`, you must explicitly set the ontology slot in the result object:

## Troubleshooting

### No results found

It is common to receive a message saying "No gene set have size \> 10 ... return NULL" or finding no significant terms.

1.  **Gene Set Size**: The default `minGSSize` is 10 and `maxGSSize` is 500. If your background annotations for a specific term have fewer than 10 genes or more than 500, they are excluded. You can adjust these parameters.
2.  **Universe**: If you provide a custom `universe` (background gene list), ensure it is large enough. A small universe can make it statistically difficult to find enrichment.
3.  **P-value Adjustment**: `clusterProfiler` uses multiple testing correction (e.g., BH) by default. Some web tools might use raw p-values or different thresholds, leading to more (but potentially false positive) results. `clusterProfiler` tends to be more conservative and reliable.

### Annotation quality

Sometimes, the `OrgDb` provided by Bioconductor might be outdated compared to the latest online databases (e.g., TAIR for Arabidopsis).

-   **Check Dates**: Always check the metadata of the `OrgDb` object to see the source and date of the annotation.
-   **Custom Annotation**: If the `OrgDb` is too old, consider downloading the latest annotation (e.g., GAF or GOSLIM files) from the organism's primary database and using `enricher()`/`GSEA()` with `buildGOmap()` as described in the [Direct vs indirect GO annotation](#direct-vs-indirect-go-annotation) section.

## Summary

GO semantic similarity can be calculated by  [@yu2010]. We can use it to cluster genes/proteins into different clusters based on their functional similarity and can also use it to measure the similarities among GO terms to reduce the redundancy of GO enrichment results.

## References

---

# KEGG enrichment analysis

    The KEGG FTP service has not been freely available for academic use since 2012, and there are many software packages using outdated KEGG annotation data. The  package supports downloading the latest online version of KEGG data using the [KEGG website](https://www.kegg.jp), which is freely available for academic users. Both KEGG pathways and modules are supported in .

    ## Supported organisms

    The  package supports all organisms that have KEGG annotation data available in the KEGG database. Users should pass an abbreviation of academic name to the `organism` parameter. The full list of KEGG supported organisms can be accessed via <http://www.genome.jp/kegg/catalog/org_list.html>. [KEGG Orthology](https://www.genome.jp/kegg/ko.html) (KO) Database is also supported by specifying `organism = "ko"`.

    The  package provides the `search_kegg_organism()` function to help search for supported organisms.

    search_kegg_organism('ece', by='kegg_code')
    ecoli <- search_kegg_organism('Escherichia coli', by='scientific_name')
    dim(ecoli)
    head(ecoli)

## Finding correct parameters: A case study with Rice

The `enrichKEGG()` function requires three key parameters: `gene` (gene IDs), `organism` (species code), and `keyType` (type of gene ID). Many users struggle with determining the correct values for `organism` and `keyType`, especially for non-model organisms. Here, we use rice (*Oryza sativa*) as an example to demonstrate how to find these parameters.

### Determine the `organism` code

The `organism` parameter accepts a KEGG organism code (e.g., 'hsa' for human, 'mmu' for mouse). To find the code for your species, you can search the [KEGG Organism Catalog](http://www.genome.jp/kegg/catalog/org_list.html) or use the `search_kegg_organism()` function as mentioned above.

For rice, searching for "rice" or "Oryza sativa" in the catalog reveals two main subspecies: \* *Oryza sativa japonica* (Japanese rice) -\> code: `osa` \* *Oryza sativa indica* (Indian rice) -\> code: `dosa`

Choose the one that matches your data. For example, if we choose `dosa`.

### Determine the `keyType`

Once the `organism` is determined, we need to know what gene ID types (`keyType`) are supported for that organism. A practical trick is to use `enrichKEGG()` as a query tool. By inputting a dummy gene ID and a likely `keyType`, the function will return an error message listing the expected ID format if the input is invalid.

# Try 'ncbi-geneid'

enrichKEGG(gene = "abc", organism = "dosa", keyType = "ncbi-geneid") \# Output: \# Error in KEGG_convert("kegg", keyType, species) : \# ncbi-geneid is not supported for dosa ...

# Try 'uniprot'

enrichKEGG(gene = "abc", organism = "dosa", keyType = "uniprot") \# Output: \# Error in KEGG_convert("kegg", keyType, species) : \# uniprot is not supported for dosa ...

    From the output, we learn that for `dosa`, the supported `keyType` is `kegg`, and the expected gene ID format looks like `Os06t0664200-01` (RAP-ID). For `osa`, the gene IDs are typically NCBI Gene IDs (e.g., `3131385`).

    ### ID Conversion

    If your gene list uses a different ID type (e.g., `LOC4334374` or `LOC_Os01g01010.1`), you need to convert them to the supported KEGG ID.

    For `LOC4334374`, you can use `AnnotationHub` to retrieve the rice `OrgDb` and convert symbols to Entrez IDs (which correspond to `osa` KEGG IDs).

    ah <- AnnotationHub()
    # Query for Oryza sativa
    query(ah, 'oryza')
    oryza <- ah[['AH55775']] # Example AH ID, check for latest
    # Check keys
    head(keys(oryza, keytype = "SYMBOL"))
    # Convert SYMBOL to ENTREZID
    # ...

For RAP-ID conversion (e.g., `LOC_Os01g01010.1` to `Os01t0100100-01`), you might need to download mapping files from the [Rice Annotation Project Database (RAP-DB)](https://rapdb.dna.affrc.go.jp/) or [Oryzabase](https://shigen.nig.ac.jp/rice/oryzabase/) and perform the conversion manually or using scripts.

### Perform Enrichment Analysis

With the correct `organism`, `keyType`, and converted gene IDs, you can now run the analysis:

# For dosa (using RAP-IDs)

kk \<- enrichKEGG(gene = gene_list_rap, organism = "dosa", keyType = "kegg", pvalueCutoff = 0.05)

    This approach of identifying `organism` code and verifying `keyType` via error messages applies to any species available in the KEGG database.

    ## KEGG Data Localization

    The  package retrieves KEGG data online via the KEGG API, which ensures that users analyze their data with the latest knowledge. However, this dependency on internet access can be a limitation if the network is unstable or unavailable (e.g., in some cloud environments). Furthermore, for the sake of reproducibility, it is often desirable to freeze the KEGG data version used in an analysis.

    To address these issues, we provide two methods to perform KEGG analysis locally.

    ### Using `gson`

    The  package provides the `gson_KEGG()` function to download the latest KEGG pathway and module data for a specific organism and store it in a `GSON` object. This object contains all necessary information for enrichment analysis and can be saved to a file for offline use.

    # download KEGG data for human
    kk_gson <- gson_KEGG(species = "hsa")

    # save to a file
    write.gson(kk_gson, file = "kegg_hsa.gson")

    # read from file
    kk_gson <- read.gson("kegg_hsa.gson")

The `GSON` object or the file can be used in `enricher()` and `GSEA()` functions (via `gson` parameter or by parsing it). For more detailed information on using the GSON format, including how to perform generalized enrichment analysis with multiple knowledge bases using `gsonList`, please refer to @sec-gson-enrichment.

### Using `createKEGGdb`

Another approach is to create a local `KEGG.db` package containing the KEGG data. While the official `KEGG.db` has not been updated since 2011, we can generate a new one using the `createKEGGdb` package.

First, install `createKEGGdb` from GitHub:

Then, use `create_kegg_db()` to download data and package it. You can specify one or more species, or even "all" to download data for all supported species.

This command will generate a `KEGG.db_1.0.tar.gz` file in the current directory. Install it as a standard R package:

Once installed, you can perform enrichment analysis using the local data by setting `use_internal_data = TRUE` in `enrichKEGG()`:

This ensures that the analysis runs offline and is fully reproducible using the installed `KEGG.db` version.

## KEGG ID Conversion

The  package provides the `bitr_kegg()` function to support ID conversion via the KEGG API. This is particularly useful for species that do not have an `OrgDb` object but are supported by KEGG.

### Basic Usage

Here is an example of converting KEGG IDs to NCBI Protein IDs and then to UniProt IDs for human genes.

# Convert KEGG ID to NCBI Protein ID

eg2np \<- bitr_kegg(hg, fromType='kegg', toType='ncbi-proteinid', organism='hsa') head(eg2np)

# Convert NCBI Protein ID to UniProt ID

np2up \<- bitr_kegg(eg2np\[,2\], fromType='ncbi-proteinid', toType='uniprot', organism='hsa') head(np2up)

    The ID type (both `fromType` and `toType`) should be one of 'kegg', 'ncbi-geneid', 'ncbi-proteinid', or 'uniprot'. The 'kegg' ID is the primary ID used in the KEGG database. A rule of thumb is that 'kegg' ID corresponds to `entrezgene` ID for eukaryote species and `Locus` ID for prokaryotes.

    For many prokaryote species, `entrezgene` IDs are not available. For example, `ece:Z5100` (E. coli O157:H7 EDL933) has `NCBI-ProteinID` and `UniProt` links, but not `NCBI-GeneID`. Attempting to convert `Z5100` to `ncbi-geneid` will result in an error stating that `ncbi-geneid` is not supported. However, conversion to `ncbi-proteinid` or `uniprot` is possible.

    # Example for prokaryote ID conversion
    bitr_kegg("Z5100", fromType="kegg", toType='ncbi-proteinid', organism='ece')
    bitr_kegg("Z5100", fromType="kegg", toType='uniprot', organism='ece')

### Converting KO numbers to Pathways

A common confusion exists between `K` numbers (KEGG Orthology) and `ko` numbers (KEGG Pathway maps). `K` numbers represent orthologous groups, while `ko` numbers represent pathway maps. Users often want to map `K` numbers to pathways.

`{r} #| label: bitr-kegg-path-ko bitr_kegg("K00844", "kegg", "Path", "ko")`

To retrieve the names of the pathways (e.g., "Glutathione metabolism" instead of "ko00480"), we can define a simple helper function `ko2name` to query the KEGG API.

`{r} #| label: ko2name-example x <- bitr_kegg("K00799", "kegg", "Path", "ko") y <- ko2name(x$Path) merge(x, y, by.x='Path', by.y='ko')`

## KEGG pathway over-representation analysis

kk \<- enrichKEGG(gene = gene, organism = 'hsa', pvalueCutoff = 0.05) head(kk)

    Input ID type can be `kegg`, `ncbi-geneid`, `ncbi-proteinid` or `uniprot` (see also @sec-kegg-id-conversion). Unlike `enrichGO()`, there is no `readable` parameter for `enrichKEGG()`. However, users can use the `setReadable()` function (see also @sec-setReadable) if there is an `OrgDb` available for the species.

    ## Translating Gene IDs to Names

    For GO analysis, the `readable` parameter controls whether to translate the IDs to human-readable gene names. This parameter is not available for KEGG analysis. However, the `setReadable()` function can translate input gene IDs to gene names if the corresponding `OrgDb` object is available.

    # 'kk' is the enrichKEGG result from the previous section
    kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    head(kk)

## KEGG pathway gene set enrichment analysis

`{r} kk2 <- gseKEGG(geneList     = geneList,                organism     = 'hsa',                minGSSize    = 120,                pvalueCutoff = 0.05,                verbose      = FALSE) head(kk2)`

## KEGG Pathway Classification

KEGG pathways are organized into a hierarchical structure with top-level categories such as "Metabolism", "Genetic Information Processing", "Environmental Information Processing", "Cellular Processes", "Organismal Systems", "Human Diseases", and "Drug Development". This classification helps users interpret enrichment results by providing higher-level biological context. For instance, users might be interested specifically in signaling pathways ("Environmental Information Processing") or metabolic pathways ("Metabolism") while excluding disease-related pathways.

The  package integrates this classification information directly into the enrichment results. The output `data.frame` contains `category` and `subcategory` columns, which can be used to filter or group the results.

### Data Update Strategy

While  queries the KEGG API for species-specific gene sets to ensure the latest data is used, the pathway classification structure is relatively stable and universal across species. To optimize performance and avoid unnecessary dependencies (e.g., web scraping packages),  caches this classification data within the package.

To prevent this cached data from becoming outdated, the package utilizes **GitHub Actions** for automated maintenance. A workflow is triggered weekly to crawl the latest pathway classification from the KEGG website. If any discrepancies are found between the online data and the cached data, the workflow automatically creates a Pull Request to update the package. This automation ensures that the classification data in  remains synchronized with KEGG updates efficiently and reliably.

## KEGG module over-representation analysis

[KEGG Module](http://www.genome.jp/kegg/module.html) is a collection of manually defined functional units. In some situations, KEGG Modules have a more straightforward interpretation.

`{r} #| label: kegg-module-example mkk <- enrichMKEGG(gene = gene,                    organism = 'hsa',                    pvalueCutoff = 1,                    qvalueCutoff = 1) head(mkk)`

## KEGG module gene set enrichment analysis

`{r} mkk2 <- gseMKEGG(geneList = geneList,                  organism = 'hsa',                  pvalueCutoff = 1) head(mkk2)`

## Visualize enriched KEGG pathways

The  package implements [several methods](#sec-enrichplot) to visualize enriched terms. Most of them are general methods that can be used on GO, KEGG, MSigDb, and other gene set annotations. Here, we introduce the `clusterProfiler::browseKEGG()` and `pathview::pathview()` functions to help users explore enriched KEGG pathways with genes of interest.

To view the KEGG pathway, users can use the `browseKEGG()` function, which will open a web browser and highlight enriched genes. See @fig-browseKEGG for a screenshot example.

`{r} #| label: browseKEGG-example #| eval: false browseKEGG(kk, 'hsa04110')`

Users can also use the `pathview()` function from the  [@luo_pathview] to visualize enriched KEGG pathways identified by the  package [@yu2012]. See @fig-pathview for a rendered pathway image.

The following example illustrates how to visualize the "hsa04110" pathway, which was enriched in our previous analysis.

## References

---

# Analysis with Other Knowledge Bases

    ## WikiPathways analysis

    [WikiPathways](https://www.wikipathways.org) is a continuously updated pathway database curated by a community of researchers and pathway enthusiasts. WikiPathways produces monthly releases of GMT files for supported organisms at [data.wikipathways.org](http://data.wikipathways.org/current/gmt/). The  package [@yu2012] supports enrichment analysis (either ORA or GSEA) for WikiPathways using the `enrichWP()` and `gseWP()` functions. These functions will automatically download and parse the latest WikiPathways GMT file for the selected organism.

    Supported organisms can be listed by:

    get_wp_organisms()

As an alternative to manually downloading gmt files, install the  to gain scripting access to the latest gmt files using the `downloadPathwayArchive()` function.

```r
## supported organisms can be accessed via the following command:
## rWikiPathways::listOrganisms()

wpgmtfile <- rWikiPathways::downloadPathwayArchive(organism="Homo sapiens", format = "gmt")
```

Once the GMT file is downloaded, we can use `read.gmt.wp()` to parse it. Note that the `read.gmt.wp()` function is designed for WikiPathways GMT files. For ordinary GMT files, please use the `read.gmt()` function.

```r
wp2gene <- read.gmt.wp(wpgmtfile)
#TERM2GENE
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene)
#TERM2NAME
wpid2name <- wp2gene %>% dplyr::select(wpid, name)

ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewp)

ewp2 <- GSEA(geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE)
head(ewp2)
```
-->
```
### WikiPathways over-representation analysis

enrichWP(gene, organism = "Homo sapiens")

    ### WikiPathways gene set enrichment analysis

    gseWP(geneList, organism = "Homo sapiens")

If your input gene ID type is not Entrez gene ID, you can use the [`bitr()`](#bitr) function to convert gene ID. If you want to convert the gene IDs in output result to gene symbols, you can use the `setReadable()` function, see also @sec-setReadable.

## DAVID functional analysis

[DAVID](https://david.ncifcrf.gov/) is a popular bioinformatics resource for functional annotation and enrichment analysis. Although we recommend using `enrichGO` and `enrichKEGG` which use up-to-date data maintained by Bioconductor,  provides `enrichDAVID` to support DAVID analysis for users who prefer it.

Users need to register an account on DAVID website and use the email address to access the web service.

`{r} #| label: enrichdavid-example #| eval: false require(clusterProfiler) data(geneList, package="DOSE") gene = names(geneList)[abs(geneList) > 2] david = enrichDAVID(gene = gene,                      idType="ENTREZ_GENE_ID",                      listType="Gene",                      annotation="KEGG_PATHWAY",                     david.user = "clusterProfiler@hku.hk")`

The result can be visualized using the same functions as other enrichment results.

`{r} #| label: visualize-david #| eval: false barplot(david) cnetplot(david, foldChange=geneList)`

### DAVID Web Service Limitations

DAVID Web Service has the following limitations:

-   A job with more than 3000 genes to generate gene or term cluster report will not be handled by DAVID due to resource limit.
-   No more than 200 jobs in a day from one user or computer.
-   DAVID Team reserves right to suspend any improper uses of the web service without notice.

For more details, please refer to <http://david.abcc.ncifcrf.gov/content.jsp?file=WS.html>.

### Custom Background

The `enrichDAVID` function also supports a custom background (universe) using the `universe` parameter.

`{r} #| label: enrichdavid-universe #| eval: false enrichDAVID(gene = gene,             universe = names(geneList),             idType = "ENTREZ_GENE_ID",             listType = "Gene",             annotation = "KEGG_PATHWAY")`

### Comparing DAVID Functional Profiles

Comparison of functional profiles among different gene clusters is also supported via `compareCluster`.

`{r} #| label: compare-david #| eval: false data(gcSample) x=compareCluster(gcSample, fun="enrichDAVID", annotation="KEGG_PATHWAY") plot(x)`

### Note on Data Updates

As highlighted in a Nature Methods study (@wadi_impact_2016), DAVID has a very low update frequency, often going more than 5 years between updates, leading to a continuous decline in knowledge completeness. This means analysis results may only reflect outdated biological knowledge, which could introduce biases in current research. As shown in previous comparisons, using `clusterProfiler` with Bioconductor annotation data often yields more comprehensive results (more annotated genes and enriched pathways) compared to DAVID.

## Generalized Enrichment Analysis with GSON format

The `gson` (Gene Set Object Notation) format provides a standardized way to store gene set information along with metadata, making it easier to share and use custom gene sets for enrichment analysis. The  package defines the `GSON` class to store this information.

### Helper functions for GSON

The  and related packages provide a series of helper functions to download and convert data from various sources into `GSON` objects. These functions typically start with the `gson_` prefix. Users can use tab completion in R (e.g., typing `gson_` and pressing Tab) to explore available functions. Some common examples include:

-   `gson_KEGG()`: Download KEGG pathway/module data.
-   `gson_WP()`: Download WikiPathways data.
-   `gson_GO()`: Download Gene Ontology data.

Additionally, users can use the `gson::gson()` function to manually create `GSON` objects from their own data tables.

### Analyzing multiple knowledge bases with `gsonList`

One of the powerful features of the GSON format is the ability to perform enrichment analysis on multiple knowledge bases simultaneously. The `gsonList()` function allows users to combine multiple `GSON` objects into a list, which can then be passed to the `enricher()` function.

# Example: Prepare GSON objects (assuming you have downloaded or created them)

# kegg_gson \<- gson_KEGG(species = "hsa")

# wp_gson \<- gson_WP(species = "Homo sapiens")

# go_gson \<- gson_GO(OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP")

# Combine them into a list

# gsons \<- gsonList(kegg = kegg_gson, wp = wp_gson, go = go_gson)

# Perform enrichment analysis

# res \<- enricher(gene = geneList, gson = gsons)

    The result of `enricher()` with a `gsonList` input is effectively a list of enrichment results, where each component corresponds to one of the input knowledge bases. You can access individual results by name or index.

    ### Visualization with `autofacet`

    To visualize the results from multiple knowledge bases together,  provides the `autofacet()` function. This function can be added to a `dotplot` to automatically facet the results by knowledge base.

    # dotplot(res) + autofacet()

This approach allows for a comprehensive view of enrichment results across different databases, facilitating comparison and interpretation. For example, you might observe that GO terms cover a broader range of biological processes due to higher annotation coverage, while KEGG pathways might provide more specific mechanistic insights.

## References

---

# Reactome enrichment analysis

     is designed for reactome pathway based analysis [@yu_reactomepa_2016]. Reactome is an open-source, open access, manually curated and peer-reviewed pathway database.

    ## Supported organisms

    Currently  supports several model organisms, including 'celegans', 'fly', 'human', 'mouse', 'rat', 'yeast' and 'zebrafish'. The input gene ID should be Entrez gene ID. We recommend using [`clusterProfiler::bitr()`](#bitr) to convert biological IDs.

    ## Reactome pathway over-representation analysis

    Enrichment analysis is a widely used approach to identify biological
    themes.  implemented `enrichPathway()` that uses [hypergeometric model](#ora-algorithm) to assess whether the number of selected genes associated with a reactome pathway is larger than expected.

    data(geneList, package="DOSE"
    )
    de <- names(geneList)[abs(geneList) > 1.5]
    head(de)
    x <- enrichPathway(gene=de, pvalueCutoff = 0.05, readable=TRUE)
    head(x)

## Reactome pathway gene set enrichment analysis

## Pathway Visualization

 implemented the `viewPathway()` to visualize selected reactome pathways. More general purpose visualization methods for ORA and GSEA results are provided in the  package and are documented in @sec-enrichplot. See @fig-viewpathway for an example.

`{r} #| label: fig-viewpathway #| fig-height: 6 #| fig-width: 8 #| fig-cap: "**Visualize reactome pathway.**" #| fig-scap: "Visualize reactome pathway." viewPathway("E2F mediated regulation of DNA replication",              readable = TRUE,              foldChange = geneList)`

```r fig.height=6, fig.width=12}
dotplot(x, showCategory=15)
```

Enrichment map can be visualized by the `enrichMap` function:
```r fig.height=10, fig.width=10}
emapplot(x)
```

In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories, we developed __*cnetplot*__ function to extract the complex association between genes and diseases.
```r fig.height=8, fig.width=8}
cnetplot(x, categorySize="pvalue", foldChange=geneList)
```

## Visualize GSEA result

```r fig.height=8, fig.width=8}
emapplot(y, color="pvalue")
```

```r fig.height=7, fig.width=10}
gseaplot(y, geneSetID = "R-HSA-69242")
```

## Comparing enriched reactome pathways among gene clusters with clusterProfiler

We have developed an `R` package [@yu_clusterprofiler:_2012] for comparing biological themes among gene clusters.  works fine with  and can compare biological themes at reactome pathway perspective.

```r fig.height=8, fig.width=13, eval=FALSE}
require(clusterProfiler)
data(gcSample)
res <- compareCluster(gcSample, fun="enrichPathway")
dotplot(res)
```

-->
```
## References

---

# Disease enrichment analysis

    We developed  [@yu_dose_2015] package to promote the investigation of diseases.  provides five methods for [measuring semantic similarities among DO terms and gene products](#DOSE-semantic-similarity), hypergeometric model and gene set enrichment analysis (GSEA) for associating disease with gene list and extracting disease association insight from genome wide expression profiles.

    ## Disease over-representation analysis

     supports enrichment analysis of Disease Ontology (DO) [@schriml_disease_2011], [Network of Cancer Gene](http://ncg.kcl.ac.uk/) [@omer_ncg] and [Disease Gene Network](http://disgenet.org/) (DisGeNET) [@janet_disgenet]. In addition, several visualization methods were provided by  to help interpreting semantic and enrichment results.

    ### Over-representation analysis for disease ontology

    In the following example, we selected genes with fold change above 1.5 as the differential genes and analyzed their disease association.

    data(geneList)
    gene <- names(geneList)[abs(geneList) > 1.5]
    head(gene)
    x <- enrichDO(gene          = gene,
                  ont           = "HDO",
                  pvalueCutoff  = 0.05,
                  pAdjustMethod = "BH",
                  universe      = names(geneList),
                  minGSSize     = 5,
                  maxGSSize     = 500,
                  qvalueCutoff  = 0.05,
                  readable      = FALSE)
    head(x)

The `enrichDO()` function requires an entrezgene ID vector as input, which is mostly the differential gene list from gene expression profile studies. Please refer to @sec-id-convert if you need to convert other gene ID types to entrezgene ID.

The `ont` parameter can be "HDO" (Human Disease Ontology), "HPO" (Human Phenotype Ontology) or "MPO" (Mouse Phenotype Ontology). `pvalueCutoff` setting the cutoff value of *p* value and adjusted *p* value; `pAdjustMethod` setting the *p* value correction methods, include the Bonferroni correction ("bonferroni"), Holm ("holm"), Hochberg ("hochberg"), Hommel ("hommel"), Benjamini & Hochberg ("BH") and Benjamini & Yekutieli ("BY") while `qvalueCutoff` is used to control *q*-values.

The `universe` sets the background gene universe for testing. If users do not explicitly set this parameter, `enrichDO()` will set the universe to all human genes that have DO annotation.

The `minGSSize` (and `maxGSSize`) indicates that only those DO terms that have more than `minGSSize` (and less than `maxGSSize`) annotated genes will be tested.

The `readable` is a logical parameter that indicates whether the entrezgene IDs will be mapped to gene symbols or not, see also @sec-setReadable.

### Over-representation analysis for the network of cancer gene

[Network of Cancer Gene](http://ncg.kcl.ac.uk/) (NCG) [@omer_ncg] is a manually curated repository of cancer genes. NCG release 5.0 (Aug. 2015) collects 1,571 cancer genes from 175 published studies.  supports analyzing gene list and determine whether they are enriched in genes known to be mutated in a given cancer type.

`{r} #| label: enrichncg-example gene2 <- names(geneList)[abs(geneList) > 3] ncg <- enrichNCG(gene2)  head(ncg)`

### Over-representation analysis for the disease gene network

[DisGeNET](http://disgenet.org/)[@janet_disgenet] is an integrative and comprehensive resource of gene-disease associations from several public data sources and the literature. It contains gene-disease associations and snp-gene-disease associations.

The enrichment analysis of disease-gene associations is supported by the `enrichDGN` function and analysis of snp-gene-disease associations is supported by the `enrichDGNv` function.

snp \<- c("rs1401296", "rs9315050", "rs5498", "rs1524668", "rs147377392", "rs841", "rs909253", "rs7193343", "rs3918232", "rs3760396", "rs2231137", "rs10947803", "rs17222919", "rs386602276", "rs11053646", "rs1805192", "rs139564723", "rs2230806", "rs20417", "rs966221") dgnv \<- enrichDGNv(snp) head(dgnv)

    ## Disease gene set enrichment analysis

    ### `gseDO` function

    In the following example, in order to speed up the compilation of this document, only gene sets with size above 120 were tested and only 100 permutations were performed.

    data(geneList)
    y <- gseDO(geneList,
               minGSSize     = 120,
               pvalueCutoff  = 0.2,
               pAdjustMethod = "BH",
               verbose       = FALSE)
    head(y, 3)

### `gseNCG` function

`{r} #| label: gsencg-example ncg <- gseNCG(geneList,               pvalueCutoff  = 0.5,               pAdjustMethod = "BH",               verbose       = FALSE) ncg <- setReadable(ncg, 'org.Hs.eg.db') head(ncg, 3)`

### `gseDGN` function

`{r} #| label: gsedgn-example dgn <- gseDGN(geneList,               pvalueCutoff  = 0.2,               pAdjustMethod = "BH",               verbose       = FALSE) dgn <- setReadable(dgn, 'org.Hs.eg.db') head(dgn, 3)`

## References

---

# Universal enrichment analysis

    The  package [@yu2012] supports both hypergeometric test and gene set enrichment analyses of many ontologies/pathways, but it's still not enough as users may want to analyze their data with unsupported organisms, slim versions of GO, novel functional annotations (e.g. GO via BlastGO or KEGG via KAAS), unsupported ontologies/pathways, or customized annotations.

    The  package provides `enricher()` function for hypergeometric test and `GSEA()` function for gene set enrichment analysis that are designed to accept user-defined annotations. They accept two additional parameters `TERM2GENE` and `TERM2NAME`. As indicated in the parameter names, `TERM2GENE` is a data.frame with the first column of term ID and the second column of corresponding mapped genes, and `TERM2NAME` is a `data.frame` with the first column of term ID and the second column of corresponding term names. `TERM2NAME` is optional.

    ## Input data

    For over representation analysis, all we need is a gene vector, that is a vector of gene IDs. These gene IDs can be obtained by differential expression analysis (*e.g.* with the  package).

    For gene set enrichment analysis, we need a ranked list of genes.  provides an example dataset `geneList` which was derived from `R` package  that contained 200 samples, including 29 samples in grade I, 136 samples in grade II and 35 samples in grade III. We computed the ratios of geometric means of grade III samples versus geometric means of grade I samples. Logarithm of these ratios (base 2) were stored in `geneList` dataset. If you want to prepare your own `geneList`, please refer to the [FAQ](#sec-genelist).

    We can load the sample data into R via:

    data(geneList, package="DOSE")
    head(geneList)

Suppose we define fold change greater than 2 as DEGs:

`{r} #| label: define-de-genes gene <- names(geneList)[abs(geneList) > 2] head(gene)`

## Cell Marker

## cell_marker_data \<- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt')

url \<- "http://yikedaxue.slwshop.cn/Cell_marker_Human.xlsx" f \<- tempfile(fileext = ".xlsx") download.file(url, f)

cell_marker_data \<- readxl::read_excel(f, 1)

## instead of `cellName`, users can use other features (e.g. `cancerType`)

cells \<- cell_marker_data %\>% dplyr::select('Cell name', GeneID) %\>% dplyr::mutate(GeneID = strsplit(GeneID, ',')) %\>% tidyr::unnest(cols = c(GeneID))

    ### Cell Marker over-representation analysis

    x <- enricher(gene, TERM2GENE = cells)
    head(x)

### Cell Marker gene set enrichment analysis

## MSigDb analysis

[Molecular Signatures Database](http://software.broadinstitute.org/gsea/msigdb) is a collection of annotated gene sets. It contains 8 major collections:

-   H: hallmark gene sets
-   C1: positional gene sets
-   C2: curated gene sets
-   C3: motif gene sets
-   C4: computational gene sets
-   C5: GO gene sets
-   C6: oncogenic signatures
-   C7: immunologic signatures

Users can download [GMT files](www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) from [Broad Institute](http://software.broadinstitute.org/gsea/msigdb) and use the `read.gmt()` function to parse the file to be used in `enricher()` and `GSEA()`.

There is an R package, [msigdbr](https://cran.r-project.org/package=msigdbr), that already packed the MSigDB gene sets in tidy data format that can be used directly with  [@yu2012].

It supports several species:

`{r} #| label: msigdbr-species library(msigdbr) msigdbr_species()`

We can retrieve all human gene sets:

`{r} #| label: msigdbr-human m_df <- msigdbr(species = "Homo sapiens") head(m_df, 2) |> as.data.frame()`

Or specific collection. Here we use C6, oncogenic gene sets as an example:

`{r} #| label: msigdbr-c6 m_t2g <- m_df |>   dplyr::filter(gs_collection == "C6") |>   dplyr::select(gs_name, ncbi_gene) head(m_t2g)`

### MSigDb over-representation analysis

`{r} #| label: msigdb-ora-example em <- enricher(gene, TERM2GENE=m_t2g) head(em)`

### MSigDb gene set enrichment analysis

In over-representation analysis, we use oncogenic gene sets (i.e. C6) to test whether the DE genes are involved in the process that leads to cancer. In this example, we will use the C3 category to test whether genes are up/down-regulated by sharing specific motifs using the GSEA approach.

em2 \<- GSEA(geneList, TERM2GENE = C3_t2g) head(em2) \`\`\`

## References

---

# Biological theme comparison

    The  package was developed for biological theme comparison [@yu2012; @wu_clusterprofiler_2021], and it provides a function, `compareCluster`, to automatically calculate enriched functional profiles of each gene clusters and aggregate the results into a single object. Comparing functional profiles can reveal functional consensus and differences among different experiments and helps in identifying differential functional modules in omics datasets.

    ## Comparing multiple gene lists

    The `compareCluster()` function applies selected function (via the `fun` parameter) to perform enrichment analysis for each gene list.

    data(gcSample)
    str(gcSample)

Users can use a named list of gene IDs as the input that passed to the `geneCluster` parameter.

`{r} #| label: compare-multiple-gene-lists-enrichGO ck <- compareCluster(geneCluster = gcSample, fun = enrichGO,                      OrgDb = org.Hs.eg.db, ont="BP") ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID") head(ck)`

## Formula interface of compareCluster

As an alternative to using named list, the `compareCluster()` function also supports passing a formula to describe more complicated experimental designs (*e.g.*, $Gene \sim time + treatment$).

formula_res \<- compareCluster(Entrez\~group+othergroup, data=mydf, fun=enrichGO, OrgDb = org.Hs.eg.db, ont="BP")

head(formula_res)

    ## Functional analysis of single-cell marker genes

    The `compareCluster()` function is highly versatile and can be directly integrated into single-cell analysis workflows. For example, after identifying marker genes for each cluster using `Seurat`, the results can be directly used for functional enrichment analysis.

    ### Using Seurat results

    Typically, `FindAllMarkers()` from `Seurat` returns a data frame where rows are genes and columns include cluster information and statistical metrics.

    # Load example data
    data("pbmc3k")
    sce <- pbmc3k.final

    # Identify marker genes
    sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE,
                                  min.pct = 0.25,
                                  thresh.use = 0.25)

    # Filter markers
    markers <- sce.markers |>
        group_by(cluster) |>
        filter(p_val_adj < 0.001) |>
        ungroup()

    # ID conversion (if needed)
    gid <- bitr(unique(markers$gene), 'SYMBOL', 'ENTREZID', OrgDb = 'org.Hs.eg.db')
    markers <- full_join(markers, gid, by = c('gene' = 'SYMBOL'))

    # Perform comparison using formula interface
    x <- compareCluster(ENTREZID ~ cluster, data = markers, fun = 'enrichKEGG')

    # Visualization
    dotplot(x, label_format = 40) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Using COSG results

Methods like `COSG` return marker genes as a list or data frame where columns represent clusters. Since a data frame is essentially a list of equal-length vectors, it can be passed directly to `compareCluster()`.

# Identify markers using COSG

marker_cosg \<- cosg( sce, groups = 'all', assay = 'RNA', slot = 'data', mu = 1, n_genes_user = 100 )

# The first element is a data frame of gene symbols

# Columns correspond to clusters

head(marker_cosg\[\[1\]\])

# Directly use the data frame for enrichment

y \<- compareCluster(marker_cosg\[\[1\]\], fun = 'enrichGO', OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL', ont = "MF")

# Visualization

dotplot(y, label_format = 60) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

    This demonstrates that `compareCluster()` can seamlessly handle various data structures commonly produced by single-cell analysis tools, simplifying the downstream functional interpretation.

    ## Visualization of functional profile comparison

    ### Dot plot

    We can visualize the result using the `dotplot()` method.

    dotplot(ck)
    dotplot(formula_res)

#\| echo: false #\| fig-keep: last #\| dev: 'png' #\| fig-width: 20 #\| fig-height: 14 #\| out-width: '100%' #\| fig-scap: \| #\| Comparing enrichment results of multiple gene lists. #\| fig-cap: \| #\| **Comparing enrichment results of multiple gene lists.** #\| (A) Using a named list of gene clusters, the results were #\| displayed as multiple columns with each one represents #\| an enrichment result of a gene cluster. #\| (B) Using formula interface, the columns represent gene clusters defined by the formula.

library(ggplot2) library(aplot)

p1 \<- dotplot(ck) + scale_y_discrete(labels=function(x) yulab.utils::str_wrap(x, width=50)) + scale_size(range = c(3, 10)) + theme_dose(14)

p2 \<- dotplot(formula_res) + scale_y_discrete(labels=function(x) yulab.utils::str_wrap(x, width=30)) + scale_size(range = c(5, 15)) + theme_dose(14) + theme(axis.text.x=element_text(angle=30, hjust=1))

plot_list(p1, p2, ncol=2, tag_levels='A')

    The formula interface allows more complicated gene cluster definition. In @fig-ccn(B), the gene clusters were defined by two variables (i.e. `group` that divides genes into `upregulated` and `downregulated` and `othergroup` that divides the genes into two categories of `A` and `B`.). The `dotplot()` function allows us to use one variable to divide the result into different facet and plot the result with other variables in each facet panel (@fig-ccf).

    dotplot(formula_res, x="group") + facet_grid(~othergroup)

By default, only top 5 (most significant) categories of each cluster are plotted. Users can change the parameter `showCategory` to specify how many categories of each cluster to be plotted, and if `showCategory` is set to `NULL`, the whole result will be plotted. The `showCategory` parameter also allows passing a vector of selected categories to [plot pathways of interest](#showing-specific-pathways)

The `dotplot()` function tries to make the comparison among different clusters more informative and reasonable. After extracting *e.g.* 10 categories for each clusters,  tries to collect overlap of these categories among clusters. For example, `term A` is enriched in all the gene clusters (*e.g.*, `g1` and `g2`) and is in the 10 most significant categories of `g1` but not `g2`.  will capture this information and include `term A` in `g2` cluster to make the comparison in `dotplot` more reasonable.

This feature ensures that if a term is significant and selected for one cluster, its statistics in other clusters are also displayed, allowing for a valid comparison. Consequently, the number of categories shown for some clusters might exceed the `showCategory` limit (e.g., 15 categories shown for *g2* while `showCategory` is 10). Disabling this by setting `includeAll = FALSE` in `dotplot()` may result in a plot that suggests zero overlap between clusters, which can be misleading (see @fig-includeAll A vs B) and is not recommended.

`` {r} #| label: fig-includeAll #| fig-width: 14 #| fig-height: 10 #| fig-cap: | #|   **Comparison of dotplot with and without includeAll parameter.**  #|   (A) `includeAll=FALSE` produces a plot that strictly follows `showCategory` but may hide overlapping terms.  #|   (B) The default behavior (`includeAll=TRUE`) recovers the missing overlap information, making the comparison more reasonable. p1 <- dotplot(ck, showCategory=5, includeAll=FALSE) + ggtitle("includeAll=FALSE") p2 <- dotplot(ck, showCategory=5) + ggtitle("default") plot_list(p1, p2, ncol=2, tag_levels='A') ``

The `dotplot()` function accepts a parameter `size` for setting the scale of dot sizes. The default parameter `size` is set to `geneRatio`, which corresponds to the `GeneRatio` column of the output. If it is set to `count`, the comparison will be based on gene counts, while if set to `rowPercentage`, the dot sizes will be normalized by `count/(sum of each row)`. Users can also map the dot size to other variables or derived variables (see [Chapter 16](#clusterProfiler-dplyr)).

To provide the full information, we also provide number of identified genes in each category (numbers in parentheses) when `by` is set to `rowPercentage` and number of gene clusters in each cluster label (numbers in parentheses) when `by` is set to `geneRatio`, as shown in @fig-ccn.

The p-values indicate which categories are more likely to have biological meanings. The dots in the plot are color-coded based on their corresponding adjusted p-values. Color gradient ranging from red to blue corresponds to the order of increasing adjusted p-values. That is, red indicates low p-values (high enrichment), and blue indicates high p-values (low enrichment). Adjusted p-values were filtered out by the threshold given by the parameter `pvalueCutoff`, and FDR can be estimated by `qvalue`.

### Gene-Concept Network

The [cnetplot](#cnetplot) also works with `compareCluster()` result.

`` {r} #| label: fig-mcnetplot   #| fig-width: 12 #| fig-height: 8 #| fig-scap: | #|   `cnetplot()` for comparing functional profiles of multiple gene clusters. #| fig-cap: | #|   **`cnetplot()` for comparing functional profiles of multiple gene clusters.**  #|   Genes and functional categories (*i.e.*, pathways) are encoded as pies to distinguish different gene clusters. cnetplot(ck) ``

## Summary

The comparison function was designed as a framework for comparing gene clusters of any kind of ontology associations, not only [groupGO](#go-classification), [enrichGO](#clusterprofiler-go-ora), [enrichKEGG](#clusterprofiler-kegg-pathway-ora), [enrichMKEGG](#clusterprofiler-kegg-module-ora), [enrichWP](#clusterprofiler-wikipathway-ora) and [enricher](#universal-api) that were provided in this package, but also other biological and biomedical ontologies, including but not limited to [enrichPathway](#reactomepa-ora), [enrichDO](#dose-do-ora), [enrichNCG](#dose-ncg-ora), [enrichDGN](#dose-dgn-ora) and [enrichMeSH](#meshes-ora).

In [@yu2012], we analyzed the publicly available expression dataset of breast tumor tissues from 200 patients (GSE11121, Gene Expression Omnibus) [@schmidt2008]. We identified 8 gene clusters from differentially expressed genes, and used the `compareCluster()` function to compare these gene clusters by their enriched biological process. In [@wu_clusterprofiler_2021], we analyzed the GSE8057 dataset which contains expression data from ovarian cancer cells at multiple time points and under two treatment conditions. Eight groups of DEG lists were analyzed simultaneously using `compareCluster()` with WikiPathways. The result indicates that the two drugs have distinct effects at the beginning but consistent effects in the later stages (Fig. 4 of [@wu_clusterprofiler_2021]).

## References

---

# Visualization of functional enrichment result

    The  package implements several visualization methods to help interpret enrichment results. It supports visualizing enrichment results obtained from  [@yu_dose_2015],  [@yu2012; @wu_clusterprofiler_2021],
     [@yu_reactomepa_2016] and  [@yu_meshes_2018]. Both over representation analysis (ORA) and gene set enrichment analysis (GSEA) are supported.

    Note: Several visualization methods were first implemented in  and rewrote from scratch using . If you want to use the [old methods](https://www.biostars.org/p/375555), you can use the [doseplot](https://github.com/GuangchuangYu/doseplot) package.

    ## Bar Plot

    Bar plot is the most widely used method to visualize enriched terms. It depicts the enrichment scores (*e.g.* p values) and gene count or ratio as bar height
    and color @fig-Barplot. Users can specify the number of terms (most significant) or selected terms (see also the [FAQ](#showing-specific-pathways)) to display via the `showCategory` parameter.

    data(geneList)
    de <- names(geneList)[abs(geneList) > 2]

    edo <- enrichDGN(de)

Other variables that derived using `mutate` can also be used as bar height or color as demonstrated in @sec-mutate and @fig-Barplot.

mutate(edo, qscore = -log(p.adjust, base=10)) \|\> barplot(x="qscore")

    p1 <- barplot(edo, showCategory=20)
    p2 <- mutate(edo, qscore = -log(p.adjust, base=10)) |>
        barplot(x="qscore")

    plot_list(p1, p2, ncol=2, tag_levels='A')

### Split top categories by ontology

When visualizing GO enrichment across all ontologies (`ont = "ALL"`), it is often more informative to show the top terms within each ontology (BP/CC/MF) rather than selecting the top terms globally. `barplot()` forwards additional arguments via `...`, and one useful option is `split`. Setting `split = "ONTOLOGY"` will perform the "top N per ontology" selection and makes it straightforward to color (and later facet) the results by the ontology.

library(clusterProfiler) library(org.Hs.eg.db) ego_all \<- enrichGO( gene = de, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", readable = TRUE )

    p_nofacet <- barplot(
      ego_all,
      x = "Count",
      showCategory = 15,
      split = "ONTOLOGY"
    ) +
      aes(fill = ONTOLOGY) +
      scale_fill_brewer(palette = "Set2")

    p_nofacet +
      geom_text(aes(label = Count), nudge_x = 2) +
      scale_y_discrete()

Note that `barplot()` wraps long term descriptions onto multiple lines by default (controlled by `label_format`, default is 30). If you prefer to keep the y-axis labels on a single line, you can override the default scale by adding `+ scale_y_discrete()` as shown above.

\`\``{r} #| label: fig-barplot-ontology-facet #| fig-width: 12 #| fig-height: 9 #| fig-scap: | #|   Bar plot of GO terms with`autofacet()`splitting by ontology. #| fig-cap: | #|   **Bar plot of GO terms with`autofacet()\` splitting by ontology.\*\* p_facet \<- barplot( ego_all, x = "Count", showCategory = 15, split = "ONTOLOGY" ) + aes(fill = ONTOLOGY) + scale_fill_brewer(palette = "Set2") + enrichplot::autofacet(by = "row", scales = "free")

p_facet + scale_y_discrete()

    If the fortified data contain an `ONTOLOGY` column (as in the output generated by `split = "ONTOLOGY"`), you can facet the plot directly with `+ enrichplot::autofacet()`, and it will automatically split panels by ontology.

    ## Dot plot

    Dot plot is similar to bar plot with the capability to encode another score as dot size.

    edo2 <- gseDO(geneList)
    dotplot(edo, showCategory=30, label_format=NULL) + ggtitle("dotplot for ORA")
    dotplot(edo2, showCategory=30, label_format=NULL) + ggtitle("dotplot for GSEA")

#\| echo: false #\| fig-width: 12 #\| fig-height: 10 #\| fig-scap: \| #\| Dot plot of enriched terms. #\| fig-cap: \| #\| **Dot plot of enriched terms.** library(ggplot2)

edo2 \<- gseDO(geneList) p1 \<- dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA") p2 \<- dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA") plot_list(p1, p2, ncol=2, tag_levels = 'A')

    Note: The `dotplot()` function also works with [`compareCluster()` output](#compare-dotplot).

    ### Highlighting specific pathways

    The `label_format` parameter in `dotplot()` accepts a function to format y-axis labels, which can be used to highlight specific pathways. This is useful for emphasizing pathways of interest in publication figures.

    First, perform enrichment analysis:

    data(geneList)
    de <- names(geneList)[1:100]
    x <- enrichDO(de)

Define a function to highlight pathways containing "cancer":

`{r} #| label: fig-label-format-func f <- function(ids) `

To render the HTML formatting, use the `ggtext` package with `element_markdown()`:

`{r} #| label: fig-dotplot-highlight library(ggplot2) library(ggtext) library(enrichplot) dotplot(x, label_format=f) + theme(axis.text.y = element_markdown())`

You can also highlight specific pathways by index:

`</i>") ids\[2\] \<- paste0("`<i style='color:#0072B2'>**", ids\[2\], "**`</i>") return(ids) }

dotplot(x, label_format=f2) + theme(axis.text.y = element_markdown())

    Or apply different colors to all pathways:

    f3 <- function(ids) {
        # cols <- rcartocolor::carto_pal(length(ids), "Vivid")
        # cols <- colorspace::rainbow_hcl(length(ids))
        cols <- rainbow(length(ids))
        ids <- paste0("<i style='color:", cols, "'>**", ids, "**</i>")
        return(ids)
    }

    dotplot(x, label_format=f3) + theme(axis.text.y = element_markdown())

### Formula interface of dotplot

The `x` variable of `dotplot()` supports a formula interface, allowing users to use derived variables. For example, we can calculate `GeneRatio/BgRatio` or `-log(p.adjust)` and use them as the x-axis variable.

`` {r} #| label: fig-dotplot-formula #| fig-width: 12 #| fig-height: 8 #| fig-cap: | #|   **Dotplot with formula interface.** #|   (A) `x = ~GeneRatio/BgRatio`. #|   (B) `x = ~ -log(p.adjust)`. p1 <- dotplot(edo, x = ~GeneRatio/BgRatio) p2 <- dotplot(edo, x = ~ -log(p.adjust)) plot_list(p1, p2, ncol=2, tag_levels='A') ``

This feature provides great flexibility in visualization. For instance, the **Rich Factor** is a common metric defined as the ratio of the number of differentially expressed genes annotated in a pathway to the total number of genes annotated in that pathway.

$$Rich Factor = \frac{Count}{BgRatio \times N}$$

where $N$ is the total number of genes in the background distribution.

We can easily visualize the **Rich Factor** using the formula interface.

Alternatively, you can precompute derived variables with `mutate()` (e.g., add `neglog10p = -log10(p.adjust)` or `richFactor = Count / (BgRatio * N)`) and then pass the new column name to `x`. This approach is useful when you want to reuse the computed variables across multiple plots. See the dedicated section on using dplyr verbs with enrichment results in @sec-mutate.

### Adjusting dot sizes

To increase the size of the dots in the dot plot, we can use the `scale_size` function from the `ggplot2` package. This allows us to specify the range of the dot sizes.

`` {r} #| label: fig-dotplot-size #| fig-width: 12 #| fig-height: 8 #| fig-cap: | #|   **Dotplot with adjusted dot sizes.** #|   (A) Default dotplot. #|   (B) Dotplot with `scale_size(range=c(2, 20))`. p1 <- dotplot(edo, showCategory=10) + ggtitle("Default") p2 <- dotplot(edo, showCategory=10) + scale_size(range=c(2, 20)) + ggtitle("Adjusted size") plot_list(p1, p2, ncol=2, tag_levels='A') ``

## Dotplot2: Comparing two clusters

The `dotplot2()` function is a specialized version of dotplot designed for comparing two selected clusters from `compareCluster` results. This function provides focused visualization when users want to directly contrast specific biological conditions or experimental groups.

# Example using compareCluster result

library(clusterProfiler) data(gcSample) xx \<- compareCluster(gcSample, fun="enrichKEGG", organism="hsa", pvalueCutoff=0.05)

# Compare cluster A and cluster B

dotplot2(xx, vars = c("A", "B"))

    The `dotplot2()` function accepts the same parameters as `dotplot()` but requires `cluster1` and `cluster2` arguments to specify which clusters to compare. This focused visualization helps in identifying differential enrichment patterns between specific experimental conditions.

    ## Gene-Concept Network

    data(geneList, package="DOSE")
    de <- names(geneList)[abs(geneList) > 2]
    edo <- enrichDGN(de)

Both the `barplot()` and `dotplot()` only displayed most significant or selected enriched terms, while users may want to know which genes are involved in these significant terms. To consider the potential biological complexities in which a gene may belong to multiple annotation categories and provide information about numeric changes when available, we developed the `cnetplot()` function to extract the complex associations. The `cnetplot()` depicts the linkages of genes and biological concepts (*e.g.* GO terms or KEGG pathways) as a network. GSEA result is also supported with only core enriched genes displayed. For GSEA results, `cnetplot` specifically visualizes the core enriched genes identified through leading edge analysis, which represent the subset of genes most responsible for driving the enrichment signal.

Users can use the 'fc_threshold' parameter in the `cnetplot` function to filter genes to only include those with `|foldChange| > fc_threshold`. This allows users to plot only high fold-change genes (e.g., `fc_threshold = 1`). In addition, users can use the `node_label_gene` parameter to label genes that are shared between groups of terms (i.e. `node_label_gene = "share"`).

`{r} #| label: fig-Networkplot   #| fig-width: 16 #| fig-height: 18 #| fig-scap: | #|   Network plot of enriched terms. #| fig-cap: | #|   **Network plot of enriched terms.** #|  ## convert gene ID to Symbol edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID') p1 <- cnetplot(edox, foldChange=geneList) p2 <- cnetplot(edox, categorySizeBy=~ -log10(pvalue), foldChange=geneList) ## only plot high fold-change genes p3 <- cnetplot(edox, foldChange=geneList, fc_threshold = 2) p4 <- cnetplot(edox, foldChange=geneList, node_label = "share") plot_list(p1, p2, p3, p4, ncol=2, tag_levels = 'A')`

If you would like label subset of the nodes, you can use the `node_label` parameter, which supports 4 possible selections (i.e. "category", "gene", "all" and "none"), as demonstrated in @fig-cnetNodeLabel.

The `node_label` parameter also supports enhanced functionality:

1.  **Vector selection**: Specify specific genes to label using a character vector, e.g., `node_label = c("gene1", "gene2")`
2.  **'exclusive' labeling**: Label only genes that belong exclusively to one category using `node_label = "exclusive"`
3.  **'share' labeling**: Label genes that are shared between multiple categories using `node_label = "share"` (as shown in the example above)
4.  **Conditional filtering**: Use comparison operators to filter genes based on fold change, e.g., `node_label = "> 1"` or `node_label = "< -1"` to label genes with absolute fold change greater than 1

`{r} #| label: fig-cnetNodeLabel   #| fig-width: 12 #| fig-height: 16 #| fig-scap: | #|   Labelling nodes by selected subset. #| fig-cap: | #|   **Labelling nodes by selected subset.**  #|   gene category (A), gene name (B),  #|   both gene category and gene name (C, default)  #|   and not to label at all (D). p1 <- cnetplot(edox, node_label="category")  p2 <- cnetplot(edox, node_label="gene")  p3 <- cnetplot(edox, node_label="all")  p4 <- cnetplot(edox, node_label="none",          color_category='firebrick',          color_item='steelblue')  plot_list(p1, p2, p3, p4, ncol=2, tag_levels = 'A')`

The `cnetplot` function can be used as a general method to visualize data relationships in a network diagram. Please refer to the vignette of .

### Customizing gene colors in cnetplot

Users can color specific genes by providing a named vector of fold changes containing only those genes.

`{r} #| label: fig-cnetplot-custom-color foldChange <- c(rep(1, 6), rep(-1, 4)) names(foldChange) <- c("MARCO", "GZMB", "CXCL11", "CXCL10",                        "LAG3", "CCL8", "PDK1", "GABRP", "MELK", "CENPE") p <- cnetplot(edox, foldChange=foldChange) p`

If you want to remove the color legend:

`{r} #| label: fig-cnetplot-no-legend p <- p + guides(color='none') p`

To create a custom legend, you can use a dummy data frame and `geom_point`:

g + guides(alpha=guide_legend( override.aes=list(color=c("blue", "red"), alpha=1, size=3), title = "VIP genes", reverse = TRUE ))

    Note: The `cnetplot()` function also works with [`compareCluster()` output](#compare-cnetplot).

    ## Heatmap-like functional classification

    The `heatplot` is similar to `cnetplot`, but displays the relationships as a
    heatmap. The gene-concept network may become too complicated if users want to
    show a large number of significant terms. The `heatplot` can simplify the result
    and make it easier to identify expression patterns.

    p1 <- heatplot(edox, showCategory=5)
    p2 <- heatplot(edox, foldChange=geneList, showCategory=5)
    plot_list(p1, p2, ncol=1, tag_levels = 'A')

The `showTop` parameter can be used to limit the number of genes displayed in the heatmap. This is particularly useful when dealing with large gene sets where only the top genes based on fold change or significance need to be visualized. For example, `showTop = 20` will display only the top 20 genes in the heatmap.

## Tree plot

The `treeplot()` function performs hierarchical clustering of enriched terms. It relies on the pairwise similarities of the enriched terms calculated by the `pairwise_termsim()` function, which by default uses Jaccard's similarity index (JC). Users can also use semantic similarity values when supported (*e.g.*, [GO](#GOSemSim), [DO](#DOSE-semantic-similarity), and [MeSH](#meshes-semantic-similarity)).

The default agglomeration method in `treeplot()` is `ward.D`, and users can specify other methods via the `cluster_method` parameter (*e.g.*, 'average', 'complete', 'median', 'centroid', *etc.*; see also the documentation of the `hclust()` function). The `treeplot()` function will cut the tree into several subtrees (specified by the `nCluster` parameter, default is 5) and label subtrees using high-frequency words. This reduces the complexity of the enriched result and improves user interpretation ability.

For fine-grained control over text appearance, the `leave_fontsize` and `clade_fontsize` parameters allow users to adjust the font size of leaf labels and clade labels respectively. For example, `leave_fontsize = 8, clade_fontsize = 10` will set leaf labels to 8pt and clade labels to 10pt.

\`\``{r} #| label: fig-treeplot   #| fig-width: 20 #| fig-height: 10 #| fig-scap: | #|   Tree plot of enriched terms. #| fig-cap: | #|   **Tree plot of enriched terms.**  #|   default (A),`hclust_method = "average"\` (B) library(ggtree)

edox2 \<- pairwise_termsim(edox) p1 \<- treeplot(edox2, cladelab_offset=8, tiplab_offset=.3, fontsize_cladelab =5) + hexpand(.2) p2 \<- treeplot(edox2, cluster_method = "average", cladelab_offset=14, tiplab_offset=.3, fontsize_cladelab =5) + hexpand(.3) aplot::plot_list(p1, p2, tag_levels='A', ncol=2)

    ## Semantic Space Plot

    While `treeplot()` visualizes semantic similarity using a hierarchical structure, the `ssplot()` (Semantic Space Plot) projects enriched terms into a low-dimensional space (e.g., using Multidimensional Scaling, MDS). This provides a complementary spatial view where terms with high semantic similarity are clustered together.

    Like `treeplot()`, `ssplot()` requires `pairwise_termsim()` to be run first. It automatically groups terms into clusters (default `nCluster` is determined automatically) and labels them with representative words.

    ssplot(edox2, nCluster=5) + ggtitle("ssplot (nCluster=5)")

`ssplot()` is particularly useful for visualizing the overall semantic landscape of enrichment results and identifying distinct functional modules in a continuous space.

## Enrichment Map

Enrichment map organizes enriched terms into a network with edges connecting overlapping gene sets. In this way, mutually overlapping gene sets are tend to cluster together, making it easy to identify functional module.

### Handling GO term redundancy

GO annotations often contain redundant terms that can dominate enrichment results, potentially obscuring other biological stories. The `simplify()` function from the `clusterProfiler` package uses semantic similarity (via the GOSemSim package) to remove redundant GO terms, providing a clearer view of distinct functional modules.

# Prepare example data

data(geneList, package="DOSE") de \<- names(geneList)\[abs(geneList) \> 2\]

# Perform GO enrichment analysis

ego \<- enrichGO(de, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)

# Remove redundant GO terms using simplify()

ego_simplified \<- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)

# Visualize both original and simplified results

ego \<- pairwise_termsim(ego) ego_simplified \<- pairwise_termsim(ego_simplified)

p1 \<- emapplot(ego, node_label_size=.8, size_edge=.5) + scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust", guide = guide_colorbar(reverse = TRUE, order=1), trans='log10') + ggtitle("Original GO terms")

p2 \<- emapplot(ego_simplified, node_label_size=.8, size_edge=.5) + scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust", guide = guide_colorbar(reverse = TRUE, order=1), trans='log10') + ggtitle("After removing redundant terms")

# Combine plots

library(patchwork) p1 + p2 + plot_layout(ncol = 2)

    The `simplify()` function removes redundant GO terms based on semantic similarity (default cutoff = 0.7). This reveals distinct functional modules that might be obscured by redundant terms in the original enrichment results. The `pairwise_termsim()` function calculates pairwise similarities between terms, which is required for `emapplot()` visualization.

    The `emapplot` function supports results obtained from hypergeometric test and gene set enrichment analysis.
    The `size_category` parameter can be used to resize nodes and the `layout` parameter can adjust the layout, as demonstrated in @fig-Enrichment.

    edo <- pairwise_termsim(edo)
    p1 <- emapplot(edo) # node_label = "category" (default)
    p2 <- emapplot(edo, node_label = "none")
    p3 <- emapplot(edo, node_label = "none", size_category=1.5)
    p4 <- emapplot(edo, node_label = "none", layout="with_fr")

    plot_list(p1, p2, p3, p4,
            ncol=2, tag_levels = 'A',
            design="AAAAAA\nBBCCDD",
            heights = c(1, .3))

The `node_label` parameter controls how the labels were displayed. The enriched terms will be displayed by default with `node_label="category"` and it can be disabled by setting `node_label="none"`.

The `node_label_size` parameter allows users to adjust the font size of node labels in the enrichment map. This is particularly useful when dealing with many overlapping terms or when labels need to be more readable. For example, `node_label_size = 3` will increase the label size compared to the default.

If `node_label="group"`, the `emapplot` function will cluster the enriched terms into different groups and only group names (determined by wordcloud) will be displayed. If `node_label="all"`, then the enriched terms and group names will be displayed simultaneously.

\`\``{r} #| label: fig-Enrichment-node-label   #| fig-width: 16 #| fig-height: 21 #| fig-scap: | #|`emapplot`with enriched terms clustering. #| fig-cap: | #|   **`emapplot`with enriched terms clustering.** #|`node_label="group"`(A) and  #|`node_label="all"\` (B). p5 \<- emapplot(edo, node_label = "group") p6 \<- emapplot(edo, node_label = "all")

plot_list(p5, p6, ncol=1, tag_levels = 'A')

    ## Biological theme comparison

    The `emapplot` function also supports results obtained from `compareCluster` function of `clusterProfiler` package. In addition to `size_category` and `layout` parameters, the number of circles in the bottom left corner can be adjusted  using the `legend_n` parameteras, and proportion of clusters in the pie chart can be adjusted using the `pie` parameter, when `pie="count"`, the proportion of  clusters in the pie chart is determined by the number of genes, as demonstrated in @fig-Enrichment2.

    data(gcSample)
    xx <- compareCluster(gcSample, fun="enrichKEGG",
                         organism="hsa", pvalueCutoff=0.05)
    xx <- pairwise_termsim(xx)
    p1 <- emapplot(xx)
    p2 <- emapplot(xx)
    p3 <- emapplot(xx, pie="count")
    p4 <- emapplot(xx, pie="count", size_category=1.5, layout="kk")
    plot_list(p1, p2, p3, p4, ncol=2, tag_levels = 'A')

## UpSet Plot

The `upsetplot` is an alternative to `cnetplot` for visualizing the complex association between genes and gene sets. It emphasizes the gene overlapping among different gene sets.

`{r} #| label: fig-upsetORA   #| fig-width: 12 #| fig-height: 5 #| fig-scap: | #|   Upsetplot for over-representation analysis. #| fig-cap: | #|   **Upsetplot for over-representation analysis.**  upsetplot(edo)`

For over-representation analysis, `upsetplot` will calculate the overlaps among different gene sets as demonstrated in @fig-upsetORA. For GSEA result, it will plot the fold change distributions of different categories (e.g. unique to pathway, overlaps among different pathways).

`{r} #| label: fig-upsetGSEA   #| fig-width: 12 #| fig-height: 5 #| fig-scap: | #|   Upsetplot for gene set enrichment analysis. #| fig-cap: | #|   **Upsetplot for gene set enrichment analysis.**  kk2 <- gseKEGG(geneList     = geneList,                organism     = 'hsa',                minGSSize    = 120,                pvalueCutoff = 0.05,                verbose      = FALSE) upsetplot(kk2)`

## ridgeline plot for expression distribution of GSEA result

The `ridgeplot` will visualize expression distributions of core enriched genes for GSEA enriched categories. It helps users to interpret up/down-regulated pathways.

`{r} #| label: fig-ridgeplot   #| fig-width: 12 #| fig-height: 12 #| fig-scap: | #|   Ridgeplot for gene set enrichment analysis. #| fig-cap: | #|   **Ridgeplot for gene set enrichment analysis.**  ridgeplot(edo2)`

## running score and preranked list of GSEA result

Running score and preranked list are traditional methods for visualizing GSEA result. The  package supports both of them to visualize the distribution of the gene set and the enrichment score.

`` {r} #| label: fig-gseaplot   #| fig-width: 12 #| fig-height: 16 #| fig-scap: | #|   gseaplot for GSEA result(`by = "runningScore"`). #| fig-cap: | #|   **gseaplot for GSEA result(`by = "runningScore"`).**  #|   `by = "runningScore"` (A), `by = "preranked"` (B), default (C) p1 <- gseaplot(edo2, geneSetID = 1, by = "runningScore", title = edo2$Description[1]) p2 <- gseaplot(edo2, geneSetID = 1, by = "preranked", title = edo2$Description[1]) p3 <- gseaplot(edo2, geneSetID = 1, title = edo2$Description[1]) plot_list(p1, p2, p3, ncol=1, tag_levels='A') ``

The `gseaplot` function also allows users to customize the colors of the running score line and the rank line.

`` {r} #| label: fig-gseaplot-custom-color #| fig-width: 12 #| fig-height: 8 #| fig-cap: | #|   **gseaplot with custom colors.** #|   (A) Default color scheme. #|   (B) Customized colors using `color`, `color.line`, and `color.vline`. p_default <- gseaplot(edo2, geneSetID = 1, title = edo2$Description[1]) p_custom <- gseaplot(edo2, geneSetID = 1, color="#DAB546", color.line='firebrick',                      color.vline="steelblue", title = edo2$Description[1]) plot_list(p_default, p_custom, ncol=2, tag_levels='A') ``

### Unified color setting with `set_enrichplot_color()`

For consistent color settings across different plot types in `enrichplot`, the `set_enrichplot_color()` function provides a unified approach to customize color mappings. This avoids repetitive color code specifications in different plotting functions.

#### Function introduction

The `set_enrichplot_color()` function can be added to any `enrichplot` visualization using the `+` operator (ggplot2 style). It offers flexible control over color transformations and palettes.

#### Parameters

-   `type`: Type of color aesthetic to modify (default: `"fill"`, can also be `"color"` for line colors)
-   `transform`: Transformation to apply to the values used for color mapping
    -   `"log10"`: Apply log10 transformation (default for recent versions)
    -   `"identity"`: Use original values without transformation
-   `colors`: Custom color palette as a vector of colors (e.g., `c("red", "blue")`)

#### Usage examples

\`\``{r} #| label: fig-set-enrichplot-color #| fig-width: 12 #| fig-height: 8 #| fig-cap: | #|   **Color customization using`set_enrichplot_color()\`.\*\* #\| (A) Default log10 transformation. #\| (B) Custom red-blue palette with log10 transformation. #\| (C) Identity transformation (no transformation).

# Example data

library(DOSE) data(geneList) de \<- names(geneList)\[abs(geneList) \> 2\] edo \<- enrichDGN(de)

# A: Default log10 transformation

p1 \<- dotplot(edo, showCategory=15) + set_enrichplot_color(type='fill', transform='log10') + ggtitle("Default log10 transform")

# B: Custom colors with log10 transformation

p2 \<- dotplot(edo, showCategory=15) + set_enrichplot_color(type='fill', transform='log10', colors=c("red", "blue")) + ggtitle("Red-blue palette")

# C: Identity transformation (no transformation)

p3 \<- dotplot(edo, showCategory=15) + set_enrichplot_color(type='fill', transform='identity') + ggtitle("Identity transform")

library(aplot) plot_list(p1, p2, p3, ncol=3, tag_levels='A')

    **Note**: Starting from enrichplot v1.29.2, log10 transformation of p-values is applied by default in functions like `dotplot()`. The `set_enrichplot_color()` function allows users to override this default behavior or customize color palettes consistently across all visualization types.

    Another method to plot GSEA result is the `gseaplot2` function:

    gseaplot2(edo2, geneSetID = 1, title = edo2$Description[1])

The `gseaplot2` also supports multile gene sets to be displayed on the same figure:

`{r} #| label: fig-gseaplot22   #| fig-width: 12 #| fig-height: 8 #| fig-scap: | #|   Gseaplot2 for GSEA result of multile gene sets. #| fig-cap: | #|   **Gseaplot2 for GSEA result of multile gene sets.**  gseaplot2(edo2, geneSetID = 1:3)`

User can also displaying the pvalue table on the plot via `pvalue_table` parameter:

`{r} #| label: fig-gseaplot23   #| fig-width: 12 #| fig-height: 8 #| fig-scap: | #|   Gseaplot2 for GSEA result of multile gene sets(add pvalue_table). #| fig-cap: | #|   **Gseaplot2 for GSEA result of multile gene sets(add pvalue_table).**  gseaplot2(edo2, geneSetID = 1:3, pvalue_table = TRUE,           color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")`

The `pvalue_table` can be customized using the following parameters:

-   `pvalue_table_rownames` to specify row names of the table (if NULL, no row names will be displayed)
-   `pvalue_table_columns` to specify column names of the table

`{r} #| label: fig-gseaplot-pvalue   #| fig-width: 12 #| fig-height: 8 #| fig-scap: | #|   Gseaplot2 with customized pvalue table. #| fig-cap: | #|   **Gseaplot2 with customized pvalue table.**  #|  gseaplot2(edo2, geneSetID = 1, pvalue_table = TRUE,           pvalue_table_rownames = NULL,           pvalue_table_columns = c("ID", "NES", "p.adjust"))`

User can specify `subplots` to only display a subset of plots:

`` {r} #| label: fig-gseaplot24 #| fig-width: 12 #| fig-height: 8 #| fig-scap: | #|   Gseaplot2 for GSEA result of multile gene sets(add subplots). #| fig-cap: | #|   **Gseaplot2 for GSEA result of multile gene sets(add subplots).**  #|   `subplots = 1` (A),`subplots = 1:2` (B) #|  p1 <- gseaplot2(edo2, geneSetID = 1:3, subplots = 1) p2 <- gseaplot2(edo2, geneSetID = 1:3, subplots = 1:2) plot_list(p1, p2, ncol=1, tag_levels = 'A') ``

### Labeling genes in GSEA plot

Users can use `geom_gsea_gene()` to label specific genes in the GSEA plot.

# Get gene set ID

id \<- edo2\$ID\[1\]

# Randomly select genes to label

set.seed(123) genes \<- sample(edo2\[\[id\]\], 5)

# Label genes on gseaplot2

p \<- gseaplot2(edo2, geneSetID = 1, title = edo2\$Description\[1\])

# Add geom_gsea_gene layer to the first subplot (running score)

p\[\[1\]\] \<- p\[\[1\]\] + geom_gsea_gene(genes, geom=geom_label) p

    If users prefer to label with gene symbols, they can convert the gene IDs to symbols (e.g., using `setReadable`) before plotting.

    # Assuming org.Hs.eg.db is available
    if (require("org.Hs.eg.db")) {
        edo2_symbol <- setReadable(edo2, 'org.Hs.eg.db', 'ENTREZID')
        id <- edo2_symbol$ID[1]
        genes_symbol <- sample(edo2_symbol[[id]], 5)

        p_symbol <- gseaplot2(edo2_symbol, geneSetID = 1, title = edo2_symbol$Description[1])
        p_symbol[[1]] <- p_symbol[[1]] + geom_gsea_gene(genes_symbol, geom=geom_text_repel)
        p_symbol
    }

The `gsearank` function plot the ranked list of genes belong to the specific gene set.

`{r} #| label: fig-gsearank #| fig-width: 12 #| fig-height: 8 #| fig-scap: | #|   Ranked list of genes belong to the specific gene set. #| fig-cap: | #|   **Ranked list of genes belong to the specific gene set.** gsearank(edo2, 1, title = edo2[1, "Description"])`

Multiple gene sets can be aligned using `cowplot`:

pp \<- lapply(1:3, function(i) { anno \<- edo2\[i, c("NES", "pvalue", "p.adjust")\] lab \<- paste0(names(anno), "=", round(anno, 3), collapse="`\n`{=tex}")

    gsearank(edo2, i, edo2[i, 2]) + xlab(NULL) +ylab(NULL) +
        annotate("text", 10000, edo2[i, "enrichmentScore"] * .75, label = lab, hjust=0, vjust=0)

}) plot_list(gglist=pp, ncol=1)

    ### Extracting data from gsearank plot

    The `gsearank` function can also output the data used for plotting by setting `output = "table"`. This allows users to inspect the running score and other metrics for genes in the gene set.

    gsearank(edo2, 1, output = "table") |> head()

Users can also merge this table with gene information (e.g. Symbol) using `bitr` or other methods.

rank_table \<- gsearank(edo2, 1, output = "table") \# Assuming 'gene' column contains Entrez IDs if (require("org.Hs.eg.db") && require("clusterProfiler")) { gene_info \<- bitr(rank_table\$gene, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db") rank_table_symbol \<- merge(rank_table, gene_info, by.x="gene", by.y="ENTREZID") head(rank_table_symbol) }

    ## pubmed trend of enriched terms

    One of the problem of enrichment analysis is to find pathways for further
    investigation. Here, we provide `pmcplot` function to plot the number/proportion
    of publications trend based on the query result from PubMed Central. Of course,
    users can use `pmcplot` in other scenarios. All text that can be queried on PMC
    is valid as input of `pmcplot`.

    terms <- edo$Description[1:5]
    p <- pmcplot(terms, 2010:2020)
    p2 <- pmcplot(terms, 2010:2020, proportion=FALSE)
    plot_list(p, p2, ncol=2)

## Volcano plot for ORA results

The `volplot()` function provides volcano plot visualization for over-representation analysis (ORA) results, allowing users to visualize the relationship between fold change (or other effect size metrics) and statistical significance.

# Example using enrichment result with fold change information

# Note: volplot requires results with fold change values

volplot(edo)

    The `volplot()` function accepts standard ggplot2 aesthetics and can be customized with additional parameters for point size, color, and labeling thresholds.

    ## Horizontal GSEA plot

    The `hplot()` function creates horizontal versions of GSEA plots, which can be useful for comparing multiple gene sets side-by-side or for presentations where horizontal layout is preferred.

    # Example using GSEA result
    edo2 <- gseDO(geneList)
    hplot(edo2, geneSetID = 1:3)

The `hplot()` function supports the same parameters as `gseaplot2()` but arranges the output horizontally. This is particularly useful for comparing multiple pathways or when vertical space is limited.

## GO Graph

The `goplot()` function can accept the output of `enrichGO` and visualize the enriched GO induced graph. See @fig-goplot for an example.

`{r} library(clusterProfiler) library(org.Hs.eg.db) data(geneList, package = "DOSE") de <- names(geneList)[abs(geneList) > 2] ego <- enrichGO(de, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)`

`{r} #| label: fig-goplot #| fig-height: 12 #| fig-width: 8 #| fig-cap: "**Goplot of enrichment analysis.**" #| fig-scap: "Goplot of enrichment analysis." goplot(ego)`

The `plotGOgraph()` function is another option to visualize the GO topology.

## References

---

# dplyr verbs for manipulating enrichment result

library(ggplot2) theme_set(theme_grey())

select \<- clusterProfiler:::select.enrichResult slice \<- clusterProfiler:::slice.enrichResult

    data(geneList)
    de = names(geneList)[1:100]
    x = enrichDO(de)

## filter

`{r} #| label: filter-example filter(x, p.adjust < .05, qvalue < 0.2)`

## arrange

`{r} #| label: arrange-example mutate(x, geneRatio = parse_ratio(GeneRatio)) %>%   arrange(desc(geneRatio))`

## select

`{r} #| label: select-example select(x, -geneID) %>% head`

## mutate

#\| fig-width: 7 #\| fig-height: 7 #\| fig-scap: \| #\| Visualizing rich factor of enriched terms using lolliplot. #\| fig-cap: \| #\| **Visualizing rich factor of enriched terms using lolliplot.**

# k/M

y \<- mutate(x, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) y

library(ggplot2) library(forcats) library(enrichplot)

ggplot(y, showCategory = 20, aes(richFactor, fct_reorder(Description, richFactor))) + geom_segment(aes(xend=0, yend = Description)) + geom_point(aes(color=p.adjust, size = Count)) + scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) + scale_size_continuous(range=c(2, 10)) + theme_minimal() + xlab("rich factor") + ylab(NULL) + ggtitle("Enriched Disease Ontology")

    A very similar concept is Fold Enrichment, which is defined as the ratio of two proportions, (k/n) / (M/N). Using `mutate` to add the fold enrichment variable is also easy:

    ```r
    mutate(x, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))

Here, the calculation of rich factor and fold enrichment is only for demonstration purposes. The `enrichplot` package provides the `dotplot` function that can directly visualize these two values without adding them to the enrichment result.

## slice

We can use `slice` to choose rows by their ordinal position in the enrichment result. Grouped result use the ordinal position with the group.

In the following example, a GSEA result of Reactome pathway was sorted by the absolute values of NES and the result was grouped by the sign of NES. We then extracted first 5 rows of each groups. The result was displayed in @fig-dot-sliceNES.

#\| fig-width: 9 #\| fig-height: 5 #\| fig-scap: \| #\| Choose pathways by ordinal positions. #\| fig-cap: \| #\| **Choose pathways by ordinal positions.**

library(ReactomePA) x \<- gsePathway(geneList)

y \<- arrange(x, abs(NES)) %\>% group_by(sign(NES)) %\>% slice(1:5)

library(forcats) library(ggplot2) library(enrichplot)

ggplot(y, aes(NES, fct_reorder(Description, NES), fill=qvalue), showCategory=10) + geom_col(orientation='y') + scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + theme_minimal() + ylab(NULL)

    ## summarise

    pbar <- function(x) {
      pi=seq(0, 1, length.out=11)

      mutate(x, pp = cut(p.adjust, pi)) |>
        group_by(pp) |>
        summarise(cnt = n()) |>
        ggplot(aes(pp, cnt)) + geom_col() +
        theme_minimal() +
        xlab("p value intervals") +
        ylab("Frequency") +
        ggtitle("p value distribution")
    }

    x <- enrichDO(de, pvalueCutoff=1, qvalueCutoff=1)
    set.seed(2020-09-10)
    random_genes <- sample(names(geneList), 100)
    y <- enrichDO(random_genes, pvalueCutoff=1, qvalueCutoff=1)
    p1 <- pbar(x)
    p2 <- pbar(y)
    aplot::plot_list(p1, p2, ncol=1, tag_levels = 'A')

## Alternative filtering approaches

While this chapter focuses on using `dplyr` verbs for manipulating enrichment results, users can also employ base R operations for filtering. For detailed information on using `[`, `[[`, and `$` operators with the `asis` parameter to preserve enrichment result objects, please refer to the [Data frame interface](#sec-data-frame-interface) section in the [Useful utilities](#useful-utilities) chapter.

The `dplyr` approach provides a more readable and chainable syntax for complex filtering operations, while base R operators offer familiarity for users accustomed to traditional data frame manipulation.

---

# AI-Assisted Biological Interpretation

Traditional enrichment analysis typically results in a list of significant pathways or GO terms. While statistically sound, these lists often leave researchers asking, "**So what?**" What is the underlying biological mechanism? Who are the key drivers? Is this a pro-survival or pro-death signal?

To bridge the gap between statistical results and biological insights, `clusterProfiler` introduces an AI-powered interpretation module. By leveraging Large Language Models (LLMs) and a multi-agent system, `clusterProfiler` can now act as a virtual bioinformatician, converting dry enrichment lists into coherent, evidence-based biological narratives.

## The `interpret` Function

The core function for this feature is `interpret()`. It accepts enrichment results (e.g., from `enrichKEGG`, `enrichGO`, or `compareCluster`) and uses an LLM to generate a structured report.

To use this feature, you need to configure an API key for a supported LLM provider (e.g., DeepSeek).

``` r

# Basic usage
# 'edo' is your enrichment result object
res <- interpret(edo)
print(res)
```

### Tasks and Inputs

`interpret()` is not just for explaining enrichment results. It breaks down LLM capabilities into three distinct tasks:

-   **`task = "interpretation"`**: (Default) Converts enrichment results into a mechanistic narrative suitable for publication (What -\> So What).
-   **`task = "annotation"`**: Performs cell type annotation for single-cell clusters using both marker genes and enrichment terms as evidence.
-   **`task = "phenotyping"`**: Assigns a "state/phenotype label" to a group (e.g., "Pro-inflammatory" or "Senescent-like").

To strengthen the evidence, `interpret()` supports "evidence synthesis" from multiple sources:

-   **Single Object**: `enrichResult`, `gseaResult`, or `compareClusterResult`.
-   **Multiple Objects**: A `list()` of results (e.g., `list(kegg_res, go_res)` or `list(cellmarker_res, go_res)`).
-   **Batch Processing**: If the input is a `compareCluster` result, it automatically splits by cluster and generates a report for each.

The key features of `interpret()` include:

-   **Prompt Skeleton**: A fixed structure to guide the LLM.
-   **Structured Output**: Enforced structure for parsing, comparison, and batch processing.
-   **Reasoning First**: Encourages "deduction before writing" to avoid merely listing pathway names.

### Cell Type Annotation

For example, we can use `Seurat` to identify marker genes for each cluster in a single-cell RNA-seq dataset. Then we can use `compareCluster` to perform enrichment analysis for each cluster. Finally, we can use `interpret` to annotate cell types based on the enrichment results and marker genes.

libray(dplyr) topN_marker \<- function(markers, n) { markers %\>% group_by(cluster) %\>% dplyr::filter(avg_log2FC \> 1) %\>% slice_head(n = n) %\>% ungroup() } top20 \<- topN_marker(pbmc.markers, 20) \# downloaded from: http://www.bio-bigdata.center/CellMarker_download_files/file/Cell_marker_Human.xlsx cm \<- rio::import("Cell_marker_Human.xlsx") x \<- compareCluster(gene\~cluster, data=top10, fun=enricher, TERM2GENE=cm\[,c("cell_name", "marker")\]) y \<- interpret(x, task="annotation")

    The output `y` is a list of interpretation results, one for each cluster. We can extract the inferred cell types.

> sapply(y, (x) x\$cell_type) 0 "Naive T Cell" 1 "Classical Monocyte" 2 "CD4+ T cell" 3 "Follicular B cell" 4 "CD8+ Cytotoxic T Cell" 5 "CD16+ monocyte (Non-classical monocyte)" 6 "Natural Killer (NK) cell" 7 "Plasmacytoid Dendritic Cell (pDC)" 8 "Megakaryocyte"

    This result is highly consistent with the manual annotation from the [Seurat pbmc3k tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#assigning-cell-type-identity-to-clusters):

    | Cluster ID | Markers | Cell Type |
    | :--- | :--- | :--- |
    | 0 | IL7R, CCR7 | Naive CD4+ T |
    | 1 | CD14, LYZ | CD14+ Mono |
    | 2 | IL7R, S100A4 | Memory CD4+ |
    | 3 | MS4A1 | B |
    | 4 | CD8A | CD8+ T |
    | 5 | FCGR3A, MS4A7 | FCGR3A+ Mono |
    | 6 | GNLY, NKG7 | NK |
    | 7 | FCER1A, CST3 | DC |
    | 8 | PPBP | Platelet |

    The full report provides detailed reasoning, confidence levels, and supporting evidence (markers/pathways) for each cluster assignment, offering transparency and explainability that simple label transfer methods lack.

    print(y)

# Enrichment Interpretation / Annotation Report

## Cell Type Annotation

### Cluster: 0

**Cell Type:** Naive T Cell **Confidence:** High

**Reasoning:** The cluster is definitively identified as a naive T cell based on the co-expression of canonical pan-T cell markers (CD3D, CD3E) and the master regulator of naive T cell identity, TCF7 (cited in top terms: Naive CD8+ T cell, Naive CD4+ T cell, etc.). The high expression of CCR7, a critical homing receptor for naive and central memory T cells, and LEF1, another Wnt-pathway TF co-operating with TCF7, further solidifies this identity. The enrichment list is dominated by naive and central memory T cell subtypes, with effector/cytotoxic terms ranking lower and lacking their specific markers (e.g., GZMB, PRF1). The presence of both CD4+ and CD8+ associated terms suggests a mixed population or a shared naive state before lineage commitment, but the core identity is naive T cell.

**Supporting Markers/Pathways:** - CD3D - CD3E - TCF7 - CCR7 - LEF1 - NOSIP - MAL

------------------------------------------------------------------------

## Cell Type Annotation

### Cluster: 1

**Cell Type:** Classical Monocyte **Confidence:** High

**Reasoning:** The enrichment list contains many related myeloid cell types, but the specific marker gene profile is definitive. The cluster expresses the core classical monocyte signature: high expression of CD14, S100A8, S100A9, FCN1, and LYZ (Top Specific/Marker Genes). While 'Myeloid cell' and 'Macrophage' are top-ranked by p-value, they are broad categories. The specific 'Classical monocyte' term (GeneRatio: 5/20, p.adjust: 1.637326e-07) is strongly supported by its gene list (S100A9/FCN1/CD14/S100A8/LYZ), which perfectly matches the top markers. The cluster lacks definitive markers to distinguish it as a Dendritic Cell (e.g., no FLT3, CD1C, CLEC9A), Macrophage (e.g., low/absent MRC1/CD163), or Neutrophil (e.g., absent MPO, ELANE). The presence of FCN1 and CD14 together is a hallmark of classical monocytes, and the absence of FCGR3A (CD16) argues against non-classical monocytes.

**Supporting Markers/Pathways:** - CD14 - S100A8 - S100A9 - FCN1 - LYZ - CST3 - TYROBP - MS4A6A

------------------------------------------------------------------------

## Cell Type Annotation

### Cluster: 2

**Cell Type:** CD4+ T cell **Confidence:** High

**Reasoning:** The cluster is definitively a T cell, as the top enriched term is 'T cell' (p.adjust: 2.86e-17) and the marker list includes the core T-cell receptor complex genes CD3D, CD3E, CD3G, and CD247 (LAT). Among T-cell subtypes, the evidence strongly favors a CD4+ lineage over CD8+. The second most significant term is 'CD4+ T cell' (p.adjust: 2.40e-14), and its gene list (IL32, CD3E, IL7R, CD27, CD3D, TNFRSF4, MAL, CD2, LTB, CD40LG, CD3G) is almost entirely contained within the top 'T cell' markers. Key CD4+ T-cell markers IL7R and CD27 are among the top specific genes. While 'CD8+ T cell' is also enriched, its signature genes (like AQP3) are present but lower in the marker list, and definitive cytotoxic CD8+ markers (e.g., GZMB, PRF1) are absent. The presence of CD40LG and TNFRSF4 (OX40), which are associated with CD4+ T helper and regulatory functions, further supports this assignment. The cluster lacks exclusive markers for NK cells (e.g., NCAM1, KLR genes) or Tregs (FOXP3), though it shows some regulatory association.

**Supporting Markers/Pathways:** - CD3D - CD3E - CD3G - IL7R - CD27 - CD40LG - IL32 - LTB - TNFRSF4

------------------------------------------------------------------------

## Cell Type Annotation

### Cluster: 3

**Cell Type:** Follicular B cell **Confidence:** High

**Reasoning:** The top enriched term is 'Follicular B cell' (p.adjust: 2.95e-16), and its gene list contains definitive B cell lineage markers (CD79A, CD79B, MS4A1, BANK1, FCER2, TCL1A) that are also present in the cluster's top specific genes. While other top terms like 'Secretory cell' and 'Classical monocyte' are enriched, they are driven almost exclusively by MHC Class II genes (HLA-DRA, HLA-DRB1, etc.), which are not specific to those cell types but are also expressed by antigen-presenting B cells. The presence of core B cell receptor components (CD79A/B) and mature B cell markers (MS4A1, FCER2, TCL1A) that are absent from monocyte/dendritic cell definitions, combined with the lack of specific monocyte (e.g., CD14, FCGR3A) or secretory cell markers, confirms the identity as a Follicular B cell.

**Supporting Markers/Pathways:** - CD79A - MS4A1 - CD79B - TCL1A - FCER2 - BANK1 - CD37 - HLA-DRA - HLA-DRB1 - CD74

------------------------------------------------------------------------

## Cell Type Annotation

### Cluster: 4

**Cell Type:** CD8+ Cytotoxic T Cell **Confidence:** High

**Reasoning:** The top enriched terms are a mixture of 'Natural killer cell' and various T cell subtypes, indicating shared cytotoxic function. However, the specific marker gene list is definitive. It includes the core T cell receptor complex genes CD3D, CD8A, and CD8B (present in 'CD8+ T cell' and 'Cytotoxic T cell' enrichments), which are lineage-defining for CD8+ T cells and absent from NK cells. While NKG7, PRF1, GZMA, GZMK, and GZMH are shared cytotoxic molecules, the co-expression of CD3D with CD8A/CD8B specifically identifies a cytotoxic T cell lineage. The absence of definitive NK-specific markers (e.g., NCAM1/CD56, KLRD1/CD94, FCGR3A/CD16) from the top marker list, and the presence of the T cell-specific signaling adaptor HCST (DAP10), supports a T cell identity. The 'Cytotoxic CD8+ T cell' enrichment term (GeneRatio: 8/20) provides the most precise functional and lineage match.

**Supporting Markers/Pathways:** - CD3D - CD8A - CD8B - GZMK - NKG7 - CCL5 - PRF1 - GZMA - GZMH - CST7 - LAG3 - KLRG1

------------------------------------------------------------------------

## Cell Type Annotation

### Cluster: 5

**Cell Type:** CD16+ monocyte (Non-classical monocyte) **Confidence:** High

**Reasoning:** The top enriched term is 'CD1C-CD141- dendritic cell' (p.adjust: 3.68e-23), but this is likely a misannotation due to shared myeloid markers. The gene list for this term (LST1, CKB, HCK, CSF1R, IFITM3, MS4A7, SERPINA1, LILRB1, CDKN1C, PILRA, FCGR3A, HMOX1, RHOC, LRRC25, SIGLEC10, MS4A4A) is a composite of pan-myeloid and monocyte-specific genes, and lacks definitive dendritic cell markers (e.g., CD1C, CLEC9A, BATF3). The cluster's specific marker list is dominated by canonical markers for non-classical monocytes: FCGR3A (CD16) is the defining marker, supported by MS4A7, CSF1R, and HES4 (enriched in 'CD16+ monocyte' and 'Non-classical monocyte' terms). The absence of T/NK/B cell markers and the presence of macrophage/pan-myeloid genes (e.g., CST3, CTSL) confirm a myeloid lineage, while the high expression of FCGR3A, MS4A7, and HES4 specifically distinguishes the non-classical monocyte subset from classical monocytes, macrophages, and dendritic cells.

**Supporting Markers/Pathways:** - FCGR3A (CD16) - MS4A7 - HES4 - CSF1R - LST1 - IFITM3 - SIGLEC10 - PILRA

------------------------------------------------------------------------

## Cell Type Annotation

### Cluster: 6

**Cell Type:** Natural Killer (NK) cell **Confidence:** High

**Reasoning:** The top enriched terms include both 'Natural killer cell' (17/20 genes, p=2.79e-27) and 'Cytotoxic CD4+ T cell' (14/20 genes, p=5.24e-28). While the latter has a slightly better p-value, the discriminatory marker analysis strongly favors NK cells. The cluster's top specific genes include canonical NK markers FGFBP2, SPON2, XCL2, XCL1, SH2D1B, and FCGR3A (CD16a), which are not specific to T cells. Critically, the cluster lacks definitive T-cell lineage markers (e.g., CD3D, CD3E, CD4, CD8A). The shared cytotoxic genes (PRF1, GNLY, GZMB, GZMA, NKG7) are expressed by both NK cells and cytotoxic T cells, but the presence of NK-specific markers and absence of T-cell receptor complex genes confirms NK cell identity.

**Supporting Markers/Pathways:** - FGFBP2 - SPON2 - XCL2 - XCL1 - SH2D1B - FCGR3A (CD16a) - KLRD1 (CD94) - GNLY - PRF1 - GZMB - GZMA - NKG7

------------------------------------------------------------------------

## Cell Type Annotation

### Cluster: 7

**Cell Type:** Plasmacytoid Dendritic Cell (pDC) **Confidence:** High

**Reasoning:** The enrichment list contains two distinct dendritic cell lineages: conventional/myeloid DCs (cDC2, CD1C+ DCs) and plasmacytoid DCs (pDC). While the top-ranked term is 'CD1C+\_A dendritic cell' (a cDC2 subtype), the cluster's specific marker genes are definitive for pDC identity. The pDC-specific markers LILRA4, CLEC4C (BDCA-2), SERPINF1, and P2RY6 are present in the top marker list and are the defining genes for the 'Plasmacytoid dendritic cell(pDC)' and 'Plasmacytoid dendritic cell' enrichment terms. Critically, the cluster lacks the core, non-overlapping markers for cDC2s: while it expresses CD1C and FCER1A (which can be expressed at low levels in some pDCs), it does NOT express the definitive cDC2 markers CLEC10A (in the cDC2 enrichment term but not in the pDC-specific marker set from the top genes) and CD1C at high specificity relative to pDC markers. The rule of exclusion applies: the top term is a cDC2 type, but the specific marker gene list is dominated by pDC markers and lacks exclusive cDC2 commitment.

**Supporting Markers/Pathways:** - LILRA4 (ILT7) - CLEC4C (BDCA-2) - SERPINF1 - P2RY6 - CLIC2 - SCT - LRRC26

------------------------------------------------------------------------

## Cell Type Annotation

### Cluster: 8

**Cell Type:** Megakaryocyte **Confidence:** High

**Reasoning:** The assignment is based on the definitive convergence of enrichment terms and marker genes. The top enriched term is 'Megakaryocyte' (p.adjust: 3.58e-10) with genes SPARC, GNG11, PF4, GP9, ITGA2B, GP1BA. The related terms 'Platelet' and 'Progenitor cell' are lower-ranked and share subsets of these genes (e.g., PF4, GP9, ITGA2B), which is expected as platelets are anucleate fragments of megakaryocytes. The top specific marker gene list is dominated by canonical megakaryocyte/platelet markers (GP9, ITGA2B, GP1BA, PF4, ITGB3, SPARC, GNG11) and lacks definitive markers for other hematopoietic lineages that could challenge this identity (e.g., no CD3, CD19, CD14, ELANE). The 'Progenitor cell' term is likely reflective of the shared SPARC and PF4 expression in some progenitor states, but the presence of terminal differentiation markers like GP1BA and ITGA2B confirms a mature megakaryocyte identity.

**Supporting Markers/Pathways:** - GP9 - ITGA2B - GP1BA - PF4 - SPARC - GNG11 - ITGB3 - CLDN5 - CMTM5 - SDPR

------------------------------------------------------------------------

## The Multi-Agent System

Instead of relying on a single prompt, `clusterProfiler` employs a **Multi-Agent System (MAS)** to ensure accuracy and depth. This system consists of three specialized agents that work in a pipeline:

1.  **Agent Cleaner**: Acts as a curator. It filters out "housekeeping" pathways (e.g., Ribosome, Spliceosome) that may be statistically significant but irrelevant to the specific biological context (e.g., tumor immunology), reducing noise.
2.  **Agent Detective**: Acts as a systems biologist. It analyzes the gene list, looks for **Hub Genes** in Protein-Protein Interaction (PPI) networks, and combines this with Fold Change data to identify **Key Drivers** and infer regulatory mechanisms.
3.  **Agent Storyteller**: Acts as a scientific writer. It synthesizes the findings from the Cleaner and Detective into a logical narrative, distinguishing between observations ("What"), mechanisms ("How"), and implications ("So What").

You can activate this deep mode using the `interpret_agent()` function.

res \<- interpret_agent(edo, context = context)

    ## Knowledge-Guided Interpretation

    Enrichment analysis often treats genes as a "bag of words," ignoring their interactions and expression changes. The **Knowledge-Guided Interpretation** mode injects external knowledge to empower the AI's reasoning.

    ### 1. PPI Networks and Hub Genes
    By setting `add_ppi = TRUE`, the system fetches protein-protein interaction data (from STRING). The AI can then identify functional modules (e.g., a TCR signaling complex) rather than just isolated genes.

    ### 2. Expression Trends
    By providing `gene_fold_change`, the AI can infer the direction of pathway activity. For example, if the "Apoptosis" pathway is enriched but pro-apoptotic genes are downregulated while anti-apoptotic genes are upregulated, the AI will correctly interpret this as "Apoptosis Resistance."

    ```r
    # Prepare a named vector of fold changes
    gene_list <- c("CD8A" = 2.5, "PDCD1" = 1.8, "GZMB" = 3.2)

    res <- interpret(edo,
                     task = "cell_type",
                     add_ppi = TRUE,                # Enable PPI network analysis
                     gene_fold_change = gene_list   # Inject expression data
    )

## Reference-Guided Interpretation

For cell type annotation, LLMs can sometimes hallucinate. To prevent this, `clusterProfiler` supports **Reference-Guided Interpretation**.

### Prior Knowledge Injection

You can provide "prior knowledge" (e.g., results from SingleR, scGPT, or manual rough annotation) to the AI. The AI acts as a validator and refiner: \* **Validation**: Checks if pathway evidence supports the prior label. \* **Refinement**: Refines a broad label (e.g., "T cell") into a specific state (e.g., "Proliferating CD8+ T cell") based on pathway activity. \* **Correction**: Flags potential misannotations if the evidence contradicts the prior.

res \<- interpret(edo, prior = my_priors, task = "cell_type")

    ### Hierarchical Interpretation
    For complex datasets, `interpret_hierarchical()` mimics the human thought process of annotating major lineages first (e.g., Myeloid) and then subtypes (e.g., M1 Macrophage). It enforces lineage constraints to prevent impossible annotations (e.g., a T cell subtype appearing within a Myeloid cluster).

    This approach is also highly applicable to **Single-cell Trajectory Inference**. Developmental processes inherently follow a hierarchical structure (e.g., Stem Cell -> Progenitor -> Terminally Differentiated Cell). By utilizing this hierarchical relationship, `interpret_hierarchical()` can provide context-aware interpretations that respect the biological differentiation path, ensuring that downstream states are interpreted within the context of their upstream progenitors.

    ```r
    # Mapping between minor and major clusters
    cluster_mapping <- c(
      "SubCluster1_1" = "MajorCluster1",
      "SubCluster1_2" = "MajorCluster1",
      "SubCluster2_1" = "MajorCluster2"
    )

    # Hierarchical interpretation
    res_hier <- interpret_hierarchical(
        x_minor = enrich_minor,
        x_major = enrich_major,
        mapping = cluster_mapping
    )

## Gene-Based Fallback Mode

In real-world research, enrichment analysis sometimes fails to return significant pathways due to small gene sets or background noise.

`clusterProfiler` introduces a **Gene-Based Fallback Mode**. When no enriched pathways are found, the agent does not simply error out. Instead, it: 1. Directly analyzes the function of the input genes. 2. Retrieves PPI networks for these genes. 3. Infers biological function based on gene connectivity and function, providing a "Medium Confidence" report instead of an empty result.

``` r
# When enrichment fails, this still works
res <- interpret_agent(tough_genes, add_ppi = TRUE)
```

## Visualization: From Story to Figure

Text reports are great, but graphical representations are often preferred for presentations. The `plot()` method for interpretation results uses `ggtangle` (a grammar of graphics for networks) to visualize the AI-inferred regulatory network.

The resulting plot highlights: \* **Key Drivers**: Central nodes identified by the AI. \* **Activation (Green)** / **Inhibition (Red)**: Regulatory relationships inferred from data. \* **Interactions (Grey)**: Physical associations.

``` r
# Visualize the interpretation result
plot(res)
```

This feature creates a closed loop: from **Enrichment** (Statistics) to **Interpretation** (Insight) and finally to **Visualization** (Communication).

---

# Useful utilities

    ## Translating Biological ID

    ### `bitr`: Biological Id TranslatoR

    The  package provides the `bitr()` and `bitr_kegg()` functions for converting ID types. Both `bitr()` and `bitr_kegg()` support many species including model and many non-model organisms.

    x <- c("GPX3",  "GLRX",   "LBP",   "CRYAB", "DEFB1", "HCLS1",   "SOD2",   "HSPA2",
           "ORM1",  "IGFBP1", "PTHLH", "GPC3",  "IGFBP3","TOB1",    "MITF",   "NDRG1",
           "NR1H4", "FGFR3",  "PVR",   "IL6",   "PTPRM", "ERBB2",   "NID2",   "LAMB1",
           "COMP",  "PLS3",   "MCAM",  "SPP1",  "LAMC1", "COL4A2",  "COL4A1", "MYOC",
           "ANXA4", "TFPI2",  "CST6",  "SLPI",  "TIMP2", "CPM",     "GGT1",   "NNMT",
           "MAL",   "EEF1A2", "HGD",   "TCN2",  "CDA",   "PCCA",    "CRYM",   "PDXK",
           "STC1",  "WARS",  "HMOX1", "FXYD2", "RBP4",   "SLC6A12", "KDELR3", "ITM2B")
    eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    head(eg)

Users should provide an annotation package, and both `fromType` and `toType` can accept any supported types.

Users can use `keytypes` to list all supported types.

`{r} library(org.Hs.eg.db) keytypes(org.Hs.eg.db)`

We can translate from one type to other types. `{r message=FALSE, warning=TRUE} library(clusterProfiler) ids <- bitr(x, fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL"), OrgDb="org.Hs.eg.db") head(ids)`

For GO analysis, users don't need to convert IDs, as all ID types provided by `OrgDb` can be used in `groupGO`, `enrichGO`, and `gseGO` by specifying the `keyType` parameter.

Users can use the `bitr` function to convert IDs using any ID types available in the `OrgDb` object. For example, users may want to know the genes in the background that belong to a specific GO term. Such information can be easily accessed using `bitr`.

`{r} m <- bitr("GO:0006805", fromType="GO", toType = "SYMBOL", OrgDb=org.Hs.eg.db) dim(m) head(m)`

Note: If you want to extract genes in your input gene list that are belong to a specific term/pathway, you can use the [geneInCategory](#how-to-extract-genes-of-a-specific-termpathway) function.

### `bitr_kegg`: converting biological IDs using KEGG API

The `bitr_kegg()` function supports ID conversion via the KEGG API, which is useful for species supported by KEGG but lacking an `OrgDb` object. For detailed usage and examples, please refer to the [KEGG ID Conversion](#kegg-id-conversion) section in the KEGG analysis chapter.

## `setReadable`: translating gene IDs to human readable symbols

Some of the functions, especially those internally supported for [DO](#dose-enrichment), [GO](#clusterprofiler-go), and [Reactome Pathway](#reactomepa), support a parameter, `readable`. If `readable = TRUE`, all the gene IDs will be translated to gene symbols. The `readable` parameter is not available for enrichment analysis of KEGG or using user's own annotation. KEGG analysis using `enrichKEGG` and `gseKEGG`, which internally query annotation information from KEEGG database and thus support all species if it is available in the KEGG database. However, KEGG database doesn't provide gene ID to symbol mapping information. For analysis using user's own annotation data, we even don't know what species is in analyzed. Translating gene IDs to gene symbols is partly supported using the `setReadable` function if and only if there is an `OrgDb` available. The `setReadable()` function also works with `compareCluster()` output.

y \<- setReadable(x, OrgDb = org.Hs.eg.db, keyType="ENTREZID") \## The geneID column is translated to symbol head(y, 3)

    For those functions that internally support `readable` parameter, user can also use `setReadable` for translating gene IDs.

    ## Parsing GMT files

    The GMT (Gene Matrix Transposed) file format is a tab delimited file format that is widely used to describe gene sets. Each row in the GMT format represents one gene set and each gene set is described by a name (or ID), a description and the genes in the gene set as illustrated in @fig-gmtFormat.

The  package implemented a function, `read.gmt()`, to parse GMT file into a `data.frame`. The WikiPathway GMT file encodes information of `version`, `wpid` and `species` into the `Description` column. The  provides the `read.gmt.wp()` function to parse WikiPathway GMT file and supports parsing information that encoded in the `Description` column.

\`\``{r} # use`wget -c`in`download.file\` wget::wget_set() url \<- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=COVID-19_Related_Gene_Sets" download.file(url, destfile = "COVID19_GeneSets.gmt")

covid19_gs \<- read.gmt("COVID19_GeneSets.gmt") head(covid19_gs)

    Gene set resources:

    + <https://maayanlab.cloud/Enrichr/#stats>
    + <http://ge-lab.org/#/data>

    There are many gene sets available online. After parsing by the `read.gmt()` function, the data can be directly used to perform enrichment analysis using `enricher()` and `GSEA()` functions (see also [Chapter 12](#universal-api)).

    ## Data frame interface for accessing enriched results

    Enrichment result objects in `clusterProfiler` support standard data frame operations including `[`, `[[`, and `$` operators, allowing users to subset and access results using familiar data frame syntax.

    ### Subsetting with `[` operator

    The `[` operator can be used to subset enrichment results based on row indices or logical conditions:

    # Subset first 5 rows
    x_first5 <- x[1:5, ]

    # Subset based on p-value threshold
    x_sig <- x[x$pvalue < 0.05, ]

    # Subset specific terms by ID
    x_specific <- x[x$ID %in% c("hsa04110", "hsa04218"), ]

By default, the `[` operator returns a `data.frame`. To preserve the enrichment result object structure for further analysis and visualization, use the `asis = TRUE` parameter:

`{r} # Preserve as enrichResult object x_filtered <- x[x$pvalue < 0.05, asis = TRUE] class(x_filtered)`

### Accessing columns with `[[` and `$` operators

The `$` operator provides convenient access to specific columns:

# Access p-value column using

pvalues \<- x\$pvalue head(pvalues)

# Access geneID column

gene_ids \<- x\$geneID head(gene_ids)

    ### Accessing genes of a specific term using `[[` operator

    When applied with a term/pathway ID, the `[[` operator returns the vector of genes that belong to that specific category in the enrichment result.

    # Genes annotated in a single term (e.g., hsa04218)
    genes_in_term <- x[["hsa04218"]]
    head(genes_in_term)

    # Retrieve genes for multiple terms
    ids <- c("hsa04218", "hsa04814")
    genes_by_term <- geneInCategory(x)[ids]
    str(genes_by_term)

### Practical filtering examples

Users can filter enrichment results based on various criteria:

# Filter by gene count (minimum 5 genes per term)

x_min5 \<- x\[x\$Count \>= 5, asis = TRUE\]

# Filter specific biological processes

x_bp \<- x\[grepl("process", x\$Description, ignore.case = TRUE), asis = TRUE\]

# Filter terms containing specific genes

x_containing_gene \<- x\[grepl("6280", x\$geneID), asis = TRUE\] \`\`\`

This data frame interface provides flexible filtering capabilities while maintaining compatibility with `enrichplot` visualization functions when using `asis = TRUE`.

## Alternative filtering approaches with dplyr

While base R operators provide familiar syntax for filtering enrichment results, users may prefer the more readable and chainable approach offered by `dplyr` verbs. For comprehensive examples of using `filter`, `arrange`, `select`, `mutate`, and other dplyr functions with enrichment result objects, please refer to the [dplyr verbs for manipulating enrichment result](#clusterProfiler-dplyr) chapter.

The `dplyr` approach is particularly useful for complex filtering operations and pipeline-style data manipulation.

## References

---

# Frequently asked questions

    ## How to prepare your own geneList

    GSEA analysis requires a ranked gene list, which contains three features:

    + numeric vector: fold change or other type of numerical variable
    + named vector: every number has a name, the corresponding gene ID
    + sorted vector: number should be sorted in decreasing order

    If you import your data from a `csv` file, the file should contains two columns, one for gene ID (no duplicated ID allowed) and another one for fold change. You can prepare your own `geneList` via the following command:

    ```r
    d = read.csv(your_csv_file)
    ## assume 1st column is ID
    ## 2nd column is FC

    ## feature 1: numeric vector
    geneList = d[,2]

    ## feature 2: named vector
    names(geneList) = as.character(d[,1])

    ## feature 3: decreasing order
    geneList = sort(geneList, decreasing = TRUE)

## No gene can be mapped

-   <https://www.biostars.org/p/431270/>
-   <https://github.com/YuLab-SMU/clusterProfiler/issues/280>

## Showing specific pathways

By default, all the visualization methods provided by  display most significant pathways. If users are interested to show some specific pathways (e.g. excluding some unimportant pathways among the top categories), users can pass a vector of selected pathways to the `showCategory` parameter in `dotplot()`, `barplot()`, `treeplot()`, `cnetplot()` and `emapplot()` etc.

See @fig-selectedPaths for a side-by-side comparison of default top categories versus user-selected pathways.

x \<- enrichKEGG(de)

## show top 10 most significant pathways and want to exclude the second one

## dotplot(x, showCategory = x\$Description\[1:10\]\[-2\])

set.seed(2020-10-27) selected_pathways \<- sample(x\$Description, 5) selected_pathways

p1 \<- dotplot(x, showCategory = 10, font.size=14) p2 \<- dotplot(x, showCategory = selected_pathways, font.size=14)

aplot::plot_list(p1, p2, tag_levels = "A")

    Note: Another solution is using the `filter` verb to extract a subset of the result as described in [Chapter 16](#clusterProfiler-dplyr).

    ## How to extract genes of a specific term/pathway

    id <- x$ID[1:3]
    id
    x[[id[1]]]
    geneInCategory(x)[id]

## Wrap long axis labels

Most of the functions in  can automatically split long labels across multiple lines. Users can passed a line width to the `label_format` parameter (default is 30). It also supports user defined function to format label strings.

See @fig-formatLabel for examples of wrapping long axis labels using numeric width and a custom labeller.

p1 \<- dotplot(y, label_format = 20) p2 \<- dotplot(y, label_format = function(x) stringr::str_wrap(x, width=20)) cowplot::plot_grid(p1, p2, ncol=2, labels=c("A", "B"))

    The `label_format` option works with `barplot()`, `dotplot()`, `heatplot()`, `treeplot` and `ridgeplot()`.

    ## Why many genes are not annotated in KEGG?

    Users often find that a significant portion of their gene list is dropped in the KEGG enrichment analysis. For example, a user reported that only 281 genes were retained from a list of over 700 genes.

    This is not a software issue but rather a reflection of the database coverage. KEGG pathways mainly focus on metabolic and signaling pathways, and many genes are not annotated in any KEGG pathway.

    We can verify the annotation coverage using `bitr_kegg()`:

    # gene is a vector of gene IDs (e.g., Entrez IDs)
    kk <- bitr_kegg(geneID = gene, fromType='ncbi-geneid', toType='Path', organism='hsa')

The warning message will indicate the percentage of genes that failed to map.

To manually verify a specific gene, you can visit the KEGG website using a URL constructed with the organism code and gene ID, for example: `http://www.genome.jp/dbget-bin/www_bget?hsa:100506775`. If the gene page does not list any pathway, it means the gene has no KEGG pathway annotation.

---
