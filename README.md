# Breast Cancer Subtype Classification Using Relative Gene Expression Order Representations

This repository contains the code files associated with the manuscript 'A Robust Machine Learning Approach for Breast
Cancer Subtype Classification Using Relative Gene
Expression Order Representations', which is currently under review.

## Problem Statement
Breast cancer subtype classification relies on RNA-seq expression profiles to aid personalized treatment. However, existing approaches lack robustness due to technical noise across different RNA-seq datasets. 

## Methods
The manuscript proposes a machine learning approach that uses relative gene expression order representations to capture biological variation while minimizing cross-sample dependence for improving the robustness of subtype classification across various technical settings. Two types of representations are derived based on relative gene expression order within each sample:
 - Rank-normalized representation: replaces gene expression values with their within-sample ranks
 - Word2vec embeddings: converts sorted gene identifiers into textual sequences, then encodes them into dense vectors using the word2vec model.

Both representations are used to train machine learning models and compared with models trained on cross-sample normalized raw counts which were normalized as follows: PyDESeq2's median of ratios normalization --> Log2 transformation --> Z-score standardization.

### Datasets
Two publicly available bulk RNA-seq datasets are used:
- SCAN-B (downloaded from the NCBI Gene Expression Omnibus (GEO) website under accession ID GSE202203)
  - SCAN-B was further divided into SCAN-B HiSeq (80% training and 20% test sets) and SCAN-B NextSeq as semi-external validation set
- TCGA-BRCA (downloaded from GDC Data Portal using TCGA biolinks R package)

The raw datasets are not provided in this repo. The [`datasets`](https://github.com/shamita98/gene-exp-ml-breast-cancer-subtyping/tree/main/datasets) folder contain the code files for:
- retrieval of TGCA-BRCA gene expression dataset
- cleaning and organization of both SCAN-B and TCGA-BRCA datasets
- train/test split of SCAN-B HiSeq dataset
- survival data extraction for survival analysis
- PAM50 gene list with gene symbols and ENSEMBL gene ids

### Exploratory Data Analysis
The [`exploratory_data_analysis`](https://github.com/shamita98/gene-exp-ml-breast-cancer-subtyping/tree/main/exploratory_data_analysis) folder demonstrates visualization and statistical tests conducted to assess gene expression distribution shifts across datasets.

### Model Development and Evaluation
The [`counts_classification`](https://github.com/shamita98/gene-exp-ml-breast-cancer-subtyping/tree/main/counts_classification), [`ranks_classification`](https://github.com/shamita98/gene-exp-ml-breast-cancer-subtyping/tree/main/ranks_classification) and [`word2vec_classification`](https://github.com/shamita98/gene-exp-ml-breast-cancer-subtyping/tree/main/word2vec_classification) folders contain the source code files for the development, cross-validation and test set evaluation of counts-, rank- and word2vec embedding-based machine learning models.

### Downstream analyses
The best machine learning models from each classification category underwent further downstream analyses which are demonstrated in the following folders:
- [`survival_analysis`](https://github.com/shamita98/gene-exp-ml-breast-cancer-subtyping/tree/main/survival_analysis): compares survival distribution of actual PAM50 subtype groups with those based on model predictions
- [`pca_gene_exp_representations`](https://github.com/shamita98/gene-exp-ml-breast-cancer-subtyping/tree/main/pca_gene_exp_representations): compares sample segregation in each of the datasets based on gene expression representation type
- [`permutation_feature_importance`](https://github.com/shamita98/gene-exp-ml-breast-cancer-subtyping/tree/main/permutation_feature_importance): identifies top 10 most important genes in all datasets for rank-based model
- [`word2vec_gene_embedding_clustering_functional_analysis`](https://github.com/shamita98/gene-exp-ml-breast-cancer-subtyping/tree/main/word2vec_gene_embedding_clustering_functional_analysis): analyzes gene clusters derived from gene-level embeddings of PAM50 genes; studies rank distribution within each gene cluster for co-occurrence analysis; retrieves non-redundant GO biological process terms enriched in each cluster

## Requirements
- Python dependencies are listed in [`environment.yml`](environment.yml). To create the conda environment, run:
```bash
conda env create -f environment.yml
conda activate proj_env
```
- Jupyter Notebook
- RStudio

