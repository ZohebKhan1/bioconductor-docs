# bioconductor-docs

LLM-parsable markdown-formatted Bioconductor documentation
and associated papers for bulk RNA-seq packages. Designed
to be indexed by [btca](https://btca.dev) and consumed as
grounded context by AI coding agents.

## Clone

```sh
git clone https://github.com/ZohebKhan1/bioconductor-docs.git
```

## Setup with btca

[btca](https://btca.dev) indexes repository content and
exposes it to LLM agents via MCP. Install, add this repo,
then query package documentation directly from your editor.

```sh
# install
bun add -g btca

# initialize in your project
btca init

# add from GitHub
btca add -n bioc-docs \
  https://github.com/ZohebKhan1/bioconductor-docs

# or add a local clone
btca add -n bioc-docs /path/to/bioconductor-docs

# query
btca ask -r bioc-docs \
  -q "DESeq2 shrinkage estimator parameters"
```

Full btca documentation: https://docs.btca.dev

## AGENTS.md integration

Add the following to your project's `AGENTS.md` or
`CLAUDE.md` to instruct coding agents to consult package
documentation before writing analysis code:

```
When using any of the following R/Bioconductor packages,
query bioconductor-docs via btca to verify function
signatures, parameter defaults, and statistical
assumptions before generating code: DESeq2, edgeR, limma,
apeglm, tximport, tximeta, clusterProfiler, fgsea, GSVA,
singscore, maSigPro, ImpulseDE2, RAPToR, ComplexHeatmap,
PCAtools.
```

## Packages

### 01 Data import and annotation

| Package       | Files              |
| ------------- | ------------------ |
| tximport      | vignette           |
| tximeta       | vignette           |
| biomaRt       | vignette, tutorial |
| AnnotationDbi | vignette           |
| GO.db         | vignette           |

### 02 Differential expression

| Package | Files                            |
| ------- | -------------------------------- |
| DESeq2  | vignette, paper, workflow        |
| edgeR   | vignette, paper, tutorial        |
| limma   | vignette, paper, tutorial, intro |
| apeglm  | vignette, paper                  |

### 03 Time series

| Package    | Files                     |
| ---------- | ------------------------- |
| maSigPro   | vignette, paper           |
| ImpulseDE2 | vignette, paper           |
| RAPToR     | vignette, paper, tutorial |

### 04 Pathway and gene set analysis

| Package         | Files                     |
| --------------- | ------------------------- |
| clusterProfiler | vignette, book            |
| enrichplot      | vignette                  |
| fgsea           | vignette, paper, tutorial |
| GSVA            | vignette, paper, tutorial |
| singscore       | vignette, paper           |
| msigdbr         | vignette                  |
| ReactomePA      | vignette                  |

### 05 Visualization

| Package        | Files                    |
| -------------- | ------------------------ |
| ComplexHeatmap | vignette, book, tutorial |
| PCAtools       | vignette                 |

## File types

- `vignette.md` -- official Bioconductor/CRAN vignette
- `paper.md` -- primary publication
- `tutorial.md` -- supplementary tutorials
- `workflow.md` -- end-to-end analysis workflows
- `book.md` -- online book content (bookdown/quarto)
- `_metadata.yml` -- source URLs and conversion metadata

## Conversion methods

- HTML vignettes: pandoc
- PDF vignettes and papers: Docling (CPU)
- Rmd/qmd book sources: pandoc
- RAPToR vignettes: pre-rendered HTML from repo

## Contact

**Author:** Zoheb Khan

**Affiliation:** Bioinformatician @ Moskowitz Lab at the University of Chicago Department of Pathology, Pediatrics, and Human Genetics

**Email:** zohebkhan600@gmail.com

**Website:** https://zohebkhan1.github.io/
