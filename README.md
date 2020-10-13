# SARS-CoV-2 Signature

This repository includes code for processing data and the analyses done for the paper "SARS-CoV-2 Early Infection Signature Identified Potential Key Infection Mechanisms and Drug Targets." `covid19_total.Rmd` script has code for the results described in the paper.


# What is this repository for?

* 'DESeq2' R package was used to normalization the reads count and find out the differetially expressed genes (DEGs). These DEGs were used in Ingenuity Pathway Analysis. 'ASSIGN' R package was used to generate a 25-gene signature. The R scripts we provide here can be used to varify the signatures in RNA-seq and single cell data. Also, it can be used to identity potential drug in ConnectivityMap database. 
* We have provided the code and various intermediate data files that we produced in performing the analyses we describe in the manuscript.


# How to normalize raw RNA-Seq data and generate the 25-gene signature

This pipeline is designed to be used in R environment.

1. install the R statistical package. We used version 3.6.1.



All the datasets are deposited in [covid19_inputdata](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx)

1. Generate Table 1 and Figure 1A, B of 25 gene expression signature: input file is [24celllines_4patients_norm.txt](https://github.com/yueli711/COVID-19/blob/master/Table1_Figure1AB_cell24_4/24celllines_4patients_norm.txt), it will generate [cell24_4.csv](https://github.com/yueli711/COVID-19/blob/master/Table1_Figure1AB_cell24_4/cell24_4.csv) during the process.  
2. Generate Figure 1C left: input files are [signature_gene_list_prior_25yueli.csv](https://github.com/yueli711/COVID-19/blob/master/Figure1C_left_5_6_7_16_15positive/signature_gene_list_prior_25yueli.csv) and [56716_15positive.txt](https://github.com/yueli711/COVID-19/blob/master/Figure1C_left_5_6_7_16_15positive/56716_15positive.txt), it will generate [cell5_6_7_16_15.csv](https://github.com/yueli711/COVID-19/blob/master/Figure1C_left_5_6_7_16_15positive/cell5_6_7_16_15.csv) during the process.
3. Generate Figure 1C right: input files are [signature_gene_list_prior_25yueli.csv](https://github.com/yueli711/COVID-19/blob/master/Figure1C_right_5_6_7_16_2negative/signature_gene_list_prior_25yueli.csv) and [56716_2negative.txt](https://github.com/yueli711/COVID-19/blob/master/Figure1C_right_5_6_7_16_2negative/56716_2negative.txt), it will generate [cell5_6_7_16_2.csv](https://github.com/yueli711/COVID-19/blob/master/Figure1C_right_5_6_7_16_2negative/cell5_6_7_16_2.csv) during the process.
4. Generate 1D BALF: input files are [signature_gene_list_prior_25yueli.csv](https://github.com/yueli711/COVID-19/blob/master/Figure1D_BALF/signature_gene_list_prior_25yueli.csv) and [cell_5_6_7_16_BALF_norm.txt](https://github.com/yueli711/COVID-19/blob/master/Figure1D_BALF/cell_5_6_7_16_BALF_norm.txt), it will generate [cell_5_6_7_16_BALF.csv](https://github.com/yueli711/COVID-19/blob/master/Figure1D_BALF/cell_5_6_7_16_BALF.csv) and [combat_cell56716_BALF.csv](https://github.com/yueli711/COVID-19/blob/master/Figure1D_BALF/combat_cell56716_BALF.csv) during the process.
5. Generate 1D PBMC: input files are[signature_gene_list_prior_25yueli.csv](https://github.com/yueli711/COVID-19/blob/master/Figure1D_PBMC/signature_gene_list_prior_25yueli.csv) and [cell_5_6_7_16_PBMC_norm.txt](https://github.com/yueli711/COVID-19/blob/master/Figure1D_PBMC/cell_5_6_7_16_PBMC_norm.txt), it will generate [cell_5_6_7_16_PBMC.csv](https://github.com/yueli711/COVID-19/blob/master/Figure1D_PBMC/cell_5_6_7_16_PBMC.csv).
6. integrate and figures of single cell: [hms_individual_integrated_OK.rds](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx).To download all the required files click [single_cell](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx).


# Required R packages

1. Seurat
2. ggplots
3. ggplot2
4. cowplot
5. scater
6. scran
7. BioParallel
8. BiocNeighbors
9. data.table
10. ASSIGN
11. sva
12. stringr
13. DESeq2
14. DESeq
15. pamr
16. readxl
17. ggpubr


# Explain all the R files:

`covid19_total.R` A script that combined all the scripts.

`cell24_4patients_25yueli.R` Use 24 celllines as training set, 4 patients as test set, generate 25 gene expression signature.

`5_6_7_16_15positive.R` Use "signature_gene_list_prior_25yueli.csv" to test series15 as positive control.

`5_6_7_16_2negative.R` Use "signature_gene_list_prior_25yueli.csv" to test series2 as negative control.

`cell_5_6_7_16_BALF_assign.R` Use "signature_gene_list_prior_25yueli.csv" to test BALF data.

`cell_5_6_7_16_PBMC_assign.R`  Use "signature_gene_list_prior_25yueli.csv" to test data.

`single_cell_integrate.R` Input single cell files, filter, CreateSeuratObject, normalize, then integrat 12 samples.

`single_cell.R` Use hms_individual_integrated_OK.rds to process single cell data, draw the figures.

`DEseq2_Norm_remove_batch_volcano_finalized.R` Excluded all the fles of not infected with SARS-CoV-2, not from homo sapiens or another file Series16_A549-ACE2_SARS-CoV-2_Rux_2. Normalize,remove batch effect of 36 and 24 celllines, DESeq, differetially expressed genes, volcano plot.

`complex_heatmap.R` Draw complex heatmap to test remove batch effect.
