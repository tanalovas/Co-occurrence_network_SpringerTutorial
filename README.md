# Microbial Co-occurrence Network Reconstruction: a SpringerTutorial

R tutorial for reconstructing microbial co-occurrence networks from 16S rRNA amplicon data using SPIEC-EASI. It accompanies the book chapter Microbial Co-Occurrence Network Reconstruction (Methods in Molecular Biology, Springer Nature).

The workflow is demonstrated on two functionally distinct tanks (Liquid waste and Denitrification) of a wastewater treatment plant, restricted to a single season so that the per-tank sample size remains adequate for StARS stability selection.

## Contents

| File | Description |
|---|---|
| `network.R` | Complete analysis script |
| `network_tutorial_phyloseq_object.RData` | Demonstration dataset: phyloseq object `ps_tutorial` (ASV table, taxonomy, sample metadata with `Tank` and `Season`) |

## Requirements

R (>= 4.2) and the following packages:

```r
# CRAN
install.packages(c("tidyverse", "igraph", "ggraph", "tidygraph", "patchwork",
                   "scales", "ggrepel", "pulsar", "huge", "tictoc"))

# Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("phyloseq")

# GitHub
remotes::install_github("zdk123/SpiecEasi")
```

The script runs end to end: subsetting, genus-level agglomeration, prevalence filtering, SPIEC-EASI fitting per tank, degree extraction, Spearman correlation and the three figures.

## Parameters

All user-adjustable settings are grouped in Section 2 of the script:

| Parameter | Default | Meaning |
|---|---|---|
| `SEASON` | `"Winter"` | Season to analyse |
| `tanks` | Liquid waste, Denitrification | Tanks to compare |
| `PREVALENCE_THRESHOLD` | `0.10` | Keep genera present in at least 10% of samples |
| `MIN_SAMPLES` | `10` | Minimum samples per tank; smaller tanks are skipped rather than fitted |
| `SE_METHOD` | `"mb"` | Meinshausen-Buhlmann neighbourhood selection |
| `SE_NLAMBDA` | `20` | Number of lambda values along the regularisation path |
| `SE_REP_NUM` | `50` | StARS subsampling replicates |
| `FOCAL_ONLY` | `TRUE` | Restrict Figures 1 and 2 to the dominant genera |
| `AB_CUMULATIVE_THRESHOLD` | `0.6` | Cumulative abundance fraction defining "dominant" |

Increasing `SE_NLAMBDA` and `SE_REP_NUM` gives a more stable network selection at the cost of longer runtime.

## Output

Three figures:

1. Node centrality: relative degree of each dominant genus, one panel per tank.
2. Co-occurrence network topology, one panel per tank; node size is proportional to relative degree and edges are SPIEC-EASI inferred associations.
3. Relative network degree against mean relative abundance for all genera retained after filtering, with the per-tank Spearman correlation (BH-adjusted).

Two `write.csv()` calls for the degree/abundance table and the correlation results are present but commented out.

## Citation

If you use this workflow, please cite the chapter and SPIEC-EASI:

> Citation of the chapter to be added upon publication.

> Kurtz ZD, Muller CL, Miraldi ER, Littman DR, Blaser MJ, Bonneau RA (2015) Sparse and compositionally robust inference of microbial ecological networks. PLoS Comput Biol 11(5):e1004226.


