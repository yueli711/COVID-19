# SARS-CoV-2 Signature

This repository includes code for processing data and the analyses done for the paper "SARS-CoV-2 Early Infection Signature Identified Potential Key Infection Mechanisms and Drug Targets." [covid19_total.Rmd](https://github.com/yueli8/COVID-19/blob/master/covid19_total.Rmd) script has code for the results described in the paper.


## What is this repository for?

* 'DESeq2' R package was used to normalization the reads count and find out the differetially expressed genes (DEGs). These DEGs were used in Ingenuity Pathway Analysis. 'ASSIGN' R package was used to generate a 25-gene signature. The R scripts we provide here can be used to varify the signatures in RNA-seq and single cell data. Also, it can be used to identity potential drug in ConnectivityMap database. 
* We have provided the code and various intermediate data files that we produced in performing the analyses we describe in the manuscript.


### How to normalize and find out the differentially expressed genes of the RNA-Seq data

This pipeline is designed to be used in R environment.

1. install the [R statistical package](https://www.r-project.org/). We used version 3.6.1.

2. Install the following R packages, which can be obtained using either the install.packages function in R or via the [Bioconductor framework](http://www.bioconductor.org/):

* Seurat
* ggplots
* ggplot2
* cowplot
* scater
* scran
* BioParallel
* BiocNeighbors
* data.table
* ASSIGN
* sva
* stringr
* DESeq2
* DESeq
* pamr
* readxl
* ggpub

3. Clone this git repository to your computer.

4. Download all the data in [covid19_inputdata](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx). Set the directory of ```~/covid19``` and store all the data in that directory.

5. Run the R script at [Normalize_differentially expressed genes.R](https://github.com/yueli8/COVID-19/blob/master/Normalize_differentially%20expressed%20genes/Normalize_differentially%20expressed%20genes.R). The input file is [gse147507_counts](https://github.com/yueli8/COVID-19/blob/master/input_files/gse147507_counts), the nori8/COVID-19/blob/master/input_files/combat_cell56716_BALF.csv) during the process.


6. Verify the 25 gene expresson signature in PBMC: the input file is [cell_5_6_7_16_PBMC_norm.txt](https://github.com/yueli8/COVID-19/blob/master/input_files/cell_5_6_7_16_PBMC_norm.txt). The code is [cell_5_6_7_16_PBMC_assign.R](https://github.com/yueli8/COVID-19/blob/master/PBMC/cell_5_6_7_16_PBMC_assign.R). It will generate [cell_5_6_7_16_PBMC.csv](https://github.com/yueli8/COVID-19/blob/master/input_files/cell_5_6_7_16_PBMC.csv).

7. All the results can been shown in Barplots. There are two ways to create barplot: one is without excel file,and the other is with excel file. If without excel file, we have to create an excel file, so the code is [bar_code_yueli.R](https://github.com/yueli8/COVID-19/blob/master/bar_code_yueli/bar_code_yueli.R). The input file should be the result of pathway_activity_testset in samples Series2, 5 , BALF and PBMC. If we can creat an excel file by our own: [DATA.xlsx](https://github.com/yueli8/COVID-19/blob/master/input_files/DATA.xlsx), then we can directly use [Barplot_verify.R](https://github.com/yueli8/COVID-19/blob/master/Barplot_verify/Barplot_verify.R). 

### Expression patterns of the 25 gene signature in single cell RNA-seq data

1. The single sell RNA-seq files [hc_51_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [hc_52_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [hc_100_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [mild_141_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [mild_142_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [mild_144_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [severe_143_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [severe_145_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [severe_146_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [severe_148_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [severe_149_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [severe_152_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx) have to be stored in the directory ```~/covid19```.

2. Runnung the code [single_cell_integrate.R](https://github.com/yueli8/COVID-19/blob/master/Single_cell_integrate/single_cell_integrate.R). After filter, and integration, it will generate the [hms_individual_integrated_OK.rds](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx) file during the process.

3. Running the code [single_cell.R](https://github.com/yueli8/COVID-19/blob/master/Single_cell_cluster_figure/single_cell.R). After name each clustered, plotted, we can have the Figure 2 and Extended Data Fig. 4-11.

### Pharmacologic signature connections identified in the ConnectivityMap (CMAP) database.

The running code is [connectivity_map.R](https://github.com/yueli8/COVID-19/blob/master/Connectivity_map/connectivity_map.R). Input file is [SARS-cov2](https://github.com/yueli8/COVID-19/blob/master/Connectivity_map/SARS-cov2.csv), It will come out connectivity map figures after running the code.

10. 
## Explain all the R files:

1. [covid19_total.R](https://github.com/yueli8/COVID-19/blob/master/covid19_total.Rmd) A script that combined all the scripts.

2. [cell24_4patients_25yueli.R](https://github.com/yueli8/COVID-19/blob/master/25_gene_expression_signature/cell24_4patients_25yueli.R) Use 24 celllines [24celllines_4patients_norm.txt](https://github.com/yueli8/COVID-19/blob/master/25_gene_expression_signature/24celllines_4patients_norm.txt)as training set, 4 patients as test set, generate 25 gene expression signature.

3. [5_6_7_16_15positive.R](https://github.com/yueli8/COVID-19/blob/master/Series15/5_6_7_16_15positive.R) Use [signature_gene_list_prior_25yueli.csv](https://github.com/yueli8/COVID-19/blob/master/input_files/signature_gene_list_prior_25yueli.csv) to test series15 as positive control.

4. [5_6_7_16_2negative.R](https://github.com/yueli8/COVID-19/blob/master/Series15/5_6_7_16_15positive.R) Use [signature_gene_list_prior_25yueli.csv](https://github.com/yueli8/COVID-19/blob/master/input_files/signature_gene_list_prior_25yueli.csv) to test series2 as negative control.

5. [cell_5_6_7_16_BALF_assign.R](https://github.com/yueli8/COVID-19/blob/master/BALF/cell_5_6_7_16_BALF_assign.R) Use [signature_gene_list_prior_25yueli.csv](https://github.com/yueli8/COVID-19/blob/master/input_files/signature_gene_list_prior_25yueli.csv) to test BALF data.

6. [cell_5_6_7_16_PBMC_assign.R](https://github.com/yueli8/COVID-19/blob/master/PBMC/cell_5_6_7_16_PBMC_assign.R)  Use [signature_gene_list_prior_25yueli.csv](https://github.com/yueli8/COVID-19/blob/master/input_files/signature_gene_list_prior_25yueli.csv) to test data.

7. [single_cell_integrate.R](https://github.com/yueli8/COVID-19/blob/master/Single_cell_integrate/single_cell_integrate.R) Input single cell files, filter, CreateSeuratObject, normalize, then integrat 12 samples.

8. [single_cell.R](https://github.com/yueli8/COVID-19/blob/master/Single_cell_integrate/single_cell_integrate.R) Use [hms_individual_integrated_OK.rds](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx) to process single cell data, draw the figures.

9. [Normalize_differentially expressed genes.R](https://github.com/yueli8/COVID-19/blob/master/script/complex_heatmap.R) Excluded all the fles of not infected with SARS-CoV-2, not from homo sapiens or another file Series16_A549-ACE2_SARS-CoV-2_Rux_2. Normalize,remove batch effect of 36 and 24 celllines, DESeq, differetially expressed genes, volcano plot.

10. [complex_heatmap.R](https://github.com/yueli8/COVID-19/blob/master/script/complex_heatmap.R) Draw complex heatmap to test remove batch effect.


## Contact information

* Moom R. Roosam. [roosan@chapman.edu](mailto:roosan@chapman.edu)
* Yue Li. [yli1@chapman.edu](mailto:yli1@chapman.edu)
