# Analysis Code for *"Uncertainty Modeling Outperforms Machine Learning for Microbiome Data Analysis"*

**Manuscript:** [PubMed: 41000811](https://pubmed.ncbi.nlm.nih.gov/41000811/)  
**Data Repository:** [Silverman-Lab/mutt](https://github.com/Silverman-Lab/mutt)

---

## Overview

This repository contains the R analysis code, helper functions, and R Markdown workflows used in the manuscript *"Uncertainty Modeling Outperforms Machine Learning for Microbiome Data Analysis."*  
It supports the full analytical pipeline — from data preprocessing and model comparison to statistical uncertainty estimation — used to evaluate microbiome datasets under various compositional and technical noise assumptions.

---

## Repository Structure

| Path | Description |
|------|--------------|
| **ALDEx3/** | Scripts and resources for ALDEx3 differential abundance analyses. |
| **functions.R** | Core analysis functions shared across scripts. |
| **helperfunctions.R** | Utility and plotting functions for reproducibility and report generation. |
| **microbial_load_predictor/** | Main microbial load prediction workflow (primary implementation). |
| **microbial_load_predictor1/** | Alternate version or experimental branch of the microbial load predictor pipeline. |
| **mlpanalysis_pub.Rmd** | Primary Analysis Script for Model Evaluation -- Machine learning and uncertainty model performance analyses (publication-ready version). |
| **studycharacteristics.R** | Script summarizing dataset characteristics and metadata harmonization. |
| **wirbel/** | Supporting files for Wirbel et al. datasets and reproducibility checks. |

---

Refer to the top of each `.Rmd` file for full session information and reproducibility environment.

---

## mutt

Clone this repository:
 ```bash
 git clone https://github.com/Silverman-Lab/mutt.git
```

## Citation
If you use this code or data (and more specifically the mutt repository), please cite:

> Konnaris, MA *et al.*  
> **Uncertainty Modeling Outperforms Machine Learning for Microbiome Data Analysis.**  
> [PubMed: 41000811](https://pubmed.ncbi.nlm.nih.gov/41000811/)

