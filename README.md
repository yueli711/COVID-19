# SARS-CoV-2 Signature

This repository includes code for processing data and the analyses done for the paper "SARS-CoV-2 Early Infection Signature Identified Potential Key Infection Mechanisms and Drug Targets." [covid19_total.Rmd](https://github.com/yueli8/COVID-19/blob/master/covid19_total.Rmd) script has code for the results described in the paper.


## What is this repository for?

* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) R package was used to normalization the read count and find out the differetially expressed genes (DEGs). These DEGs were used in Ingenuity Pathway Analysis. [ASSIGN](https://bioconductor.org/packages/release/bioc/html/ASSIGN.html) R package was used to generate a 25-gene signature. The R scripts we provide here can be used to verify the signatures in bulk and single cell RNA sequencing data. Also, it can be used to identity potential drug in ConnectivityMap database. 
* We have provided the code and various intermediate data files that we produced in performing the analyses described in the manuscript.


### Install the software, download the data and set up the directory

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
* ggpubr

3. Clone this github repository to your computer.

4. Download all the data in [covid19_inputdata](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx). Set the directory of ```~/covid19``` and store all the data in that directory.

### Normalize the read count, remove batch effect and find out the differentially expresssed genes

1. Run the R script of [Normalize_differentially expressed genes.R](https://github.com/yueli8/COVID-19/blob/master/script/Normalize_differentially_expressed_genes.R). The input file is [gse147507_counts](https://github.com/yueli8/COVID-19/blob/master/input_files/gse147507_counts).

2. After remove the batch effect, it will come out two files: [ComBat_data(all)(counts).txt_24](https://github.com/yueli8/COVID-19/blob/master/Normalize_differentially_expressed_genes/ComBat_data(all)(counts).txt_24) and [ComBat_data(all)(counts).txt_36](https://github.com/yueli8/COVID-19/blob/master/Normalize_differentially_expressed_genes/ComBat_data(all)(counts).txt_36). The first one is for series 5, 6, 7, 16, and the second is series 2, 5, 6, 7, 15, 16. Series 5, 6, 7, 16 were used to generate gene expression signatures, series 2 is used as negative control, and series 15 is positive control. 

3. The normalized file is [gse147507_norm.csv](https://github.com/yueli8/COVID-19/blob/master/Normalize_differentially_expressed_genes/gse147507_norm.csv), and the differentially expressed genes file is [gse147507_deg.csv](https://github.com/yueli8/COVID-19/blob/master/Normalize_differentially_expressed_genes/gse147507_deg.csv).



### Generate a 25 gene expression signature

1. Store the file [24celllines_4patients_norm.txt](https://github.com/yueli8/COVID-19/blob/master/25_gene_expression_signature/24celllines_4patients_norm.txt) in the directory ```~/covid19```, run the code [cell24_4patients_25yueli.R](https://github.com/yueli8/COVID-19/blob/master/script/cell24_4patients_25yueli.R). It will generate [cell24_4.csv](https://github.com/yueli8/COVID-19/blob/master/25_gene_expression_signature/cell24_4.csv) during the process.

2. The 25 gene expression signature and Figure 1a, Table 1, Extended Data Fig. 2 will generate after running the code.

### Verify the 25 gene expression signature in the test data (cell lines and clinical patients)

1. Four test sets were used. (1) Series 2 is A549 cell line infeted with SARS-CoV-2 compared with mock; (2) Series 15 is postmortem COVID-19 patients compared with healthy lung biospy; (3) Bronchoalveolar lavage fluid (BALF) cells of COVID-19 patients compared with healthy control; (4) Peripheral blood mononuclear cells (PBMC) of COVID-19 patients compared with healthy controls.

2. We have to always put the [signature_gene_list_prior_25yueli.csv](https://github.com/yueli8/COVID-19/blob/master/input_files/signature_gene_list_prior_25yueli.csv) file in the directory of ```~/covid19```.

3. Series 15 used as positive control: the input files is [56716_15positive.txt](https://github.com/yueli8/COVID-19/blob/master/Series15/56716_15positive.txt). The code is [5_6_7_16_15positive.R](https://github.com/yueli8/COVID-19/blob/master/Series15/5_6_7_16_15positive.R). It will generate [cell5_6_7_16_15.csv](https://github.com/yueli8/COVID-19/blob/master/Series15/cell5_6_7_16_15.csv) during the process.

4. Series 2 used as negative control: the input fiel is [56716_2negative.txt](https://github.com/yueli8/COVID-19/blob/master/Series2/56716_2negative.txt). The code is [5_6_7_16_2negative.R](https://github.com/yueli8/COVID-19/blob/master/Series2/5_6_7_16_2negative.R). It will generate [cell5_6_7_16_2.csv](https://github.com/yueli8/COVID-19/blob/master/Series2/cell5_6_7_16_2.csv) during the process.

5. Verify the 25 gene expression signature in BALF: the input file is [cell_5_6_7_16_BALF_norm.txt](https://github.com/yueli8/COVID-19/blob/master/BALF/cell_5_6_7_16_BALF_norm.txt). The code is [cell_5_6_7_16_BALF_assign.R](https://github.com/yueli8/COVID-19/blob/master/BALF/cell_5_6_7_16_BALF_assign.R). It will generate [cell_5_6_7_16_BALF.csv](https://github.com/yueli8/COVID-19/blob/master/BALF/cell_5_6_7_16_BALF.csv) and [combat_cell56716_BALF.csv](https://github.com/yueli8/COVID-19/blob/master/BALF/combat_cell56716_BALF.csv) during the process.

6. Verify the 25 gene expresson signature in PBMC: the input file is [cell_5_6_7_16_PBMC_norm.txt](https://github.com/yueli8/COVID-19/blob/master/PBMC/cell_5_6_7_16_PBMC_norm.txt). The code is [cell_5_6_7_16_PBMC_assign.R](https://github.com/yueli8/COVID-19/blob/master/PBMC/cell_5_6_7_16_PBMC_assign.R). It will generate [cell_5_6_7_16_PBMC.csv](https://github.com/yueli8/COVID-19/blob/master/PBMC/cell_5_6_7_16_PBMC.csv).

7. All the results can been shown in Barplots. There are two ways to create barplot: one is without excel file,and the other is with excel file. If without excel file, we have to create an excel file, so the code is [bar_code_yueli.R](https://github.com/yueli8/COVID-19/blob/master/script/bar_code_yueli.R). The input file should be the result of pathway_activity_testset in samples Series2, 5 , BALF and PBMC. If we can creat an excel file by our own: [DATA.xlsx](https://github.com/yueli8/COVID-19/blob/master/input_files/DATA.xlsx), then we can directly use [Barplot_verify.R](https://github.com/yueli8/COVID-19/blob/master/Barplot_verify/Barplot_verify.R).
 
 
### Expression patterns of the 25 gene signature in single cell RNA-seq data

1. The single sell RNA-seq files [hc_51_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [hc_52_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [hc_100_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [mild_141_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [mild_142_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [mild_144_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [severe_143_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [severe_145_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [severe_146_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [severe_148_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [severe_149_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx), [severe_152_02.csv](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx) have to be stored in the directory ```~/covid19```.

2. Runnung the code [single_cell_integrate.R](https://github.com/yueli8/COVID-19/blob/master/Single_cell_integrate/single_cell_integrate.R). After filter, and integration, it will generate the [hms_individual_integrated_OK.rds](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx) file during the process.

3. Running the code [single_cell.R](https://github.com/yueli8/COVID-19/blob/master/Single_cell_cluster_figure/single_cell.R). After name each clustered, plotted, we can have the Figure 2 and Extended Data Fig. 4-11.

### Pharmacologic signature connections identified in the ConnectivityMap (CMAP) database.

The running code is [connectivity_map.R](https://github.com/yueli8/COVID-19/blob/master/Connectivity_map/connectivity_map.R). Input file is [SARS-cov2](https://github.com/yueli8/COVID-19/blob/master/Connectivity_map/SARS-cov2.csv), It will come out connectivity map figures after running the code.


## Explain all the R files:

1. [covid19_total.R](https://github.com/yueli8/COVID-19/blob/master/covid19_total.Rmd) A script that combined all the scripts.

2. [cell24_4patients_25yueli.R](https://github.com/yueli8/COVID-19/blob/master/25_gene_expression_signature/cell24_4patients_25yueli.R) Use 24 celllines [24celllines_4patients_norm.txt](https://github.com/yueli8/COVID-19/blob/master/25_gene_expression_signature/24celllines_4patients_norm.txt) as training set, 4 patients as test set, generate 25 gene expression signature.

3. [5_6_7_16_15positive.R](https://github.com/yueli8/COVID-19/blob/master/Series15/5_6_7_16_15positive.R) Use [signature_gene_list_prior_25yueli.csv](https://github.com/yueli8/COVID-19/blob/master/input_files/signature_gene_list_prior_25yueli.csv) to test series15 as positive control.

4. [5_6_7_16_2negative.R](https://github.com/yueli8/COVID-19/blob/master/Series2/5_6_7_16_2negative.R) Use [signature_gene_list_prior_25yueli.csv](https://github.com/yueli8/COVID-19/blob/master/input_files/signature_gene_list_prior_25yueli.csv) to test series2 as negative control.

5. [cell_5_6_7_16_BALF_assign.R](https://github.com/yueli8/COVID-19/blob/master/BALF/cell_5_6_7_16_BALF_assign.R) Use [signature_gene_list_prior_25yueli.csv](https://github.com/yueli8/COVID-19/blob/master/input_files/signature_gene_list_prior_25yueli.csv) to test BALF data.

6. [cell_5_6_7_16_PBMC_assign.R](https://github.com/yueli8/COVID-19/blob/master/PBMC/cell_5_6_7_16_PBMC_assign.R)  Use [signature_gene_list_prior_25yueli.csv](https://github.com/yueli8/COVID-19/blob/master/input_files/signature_gene_list_prior_25yueli.csv) to test data.

7. [single_cell_integrate.R](https://github.com/yueli8/COVID-19/blob/master/Single_cell_integrate/single_cell_integrate.R) Input single cell files, filter, CreateSeuratObject, normalize, then integrat 12 samples.

8. [single_cell.R](https://github.com/yueli8/COVID-19/blob/master/Single_cell_cluster_figure/single_cell.R) Use [hms_individual_integrated_OK.rds](https://drive.google.com/drive/folders/1mIFiEcPm3o5FEkeBD4v3MaGjM0xazGbx) to process single cell data, draw the figures.

9. [Normalize_differentially_expressed_genes.R](https://github.com/yueli8/COVID-19/blob/master/Normalize_differentially_expressed_genes/Normalize_differentially_expressed_genes.R) Excluded all the fles of not infected with SARS-CoV-2, not from homo sapiens or another file Series16_A549-ACE2_SARS-CoV-2_Rux_2. Normalize,remove batch effect of 36 and 24 celllines, DESeq, differetially expressed genes, volcano plot.

10. [complex_heatmap.R](https://github.com/yueli8/COVID-19/blob/master/script/complex_heatmap.R) Draw complex heatmap to test remove batch effect.

11. [Barplot_verify.R](https://github.com/yueli8/COVID-19/blob/master/Barplot_verify/Barplot_verify.R) Used to draw bar plot of the Series 2, 15, BALF and PBMC.

12. [connectivity_map.R](https://github.com/yueli8/COVID-19/blob/master/Connectivity_map/connectivity_map.R) Used the [SARS-cov2.csv](https://github.com/yueli8/COVID-19/blob/master/Connectivity_map/SARS-cov2.csv) to draw the connectivity map. 

13. [bar_code_yueli.R](https://github.com/yueli8/COVID-19/blob/master/script/bar_code_yueli.R) Create an excel file of all the cell lines and clinical patients, then draw the bar code. 


## Contact information

* Moom R. Roosam. [roosan@chapman.edu](mailto:roosan@chapman.edu)
* Yue Li. [yli1@chapman.edu](mailto:yli1@chapman.edu)
