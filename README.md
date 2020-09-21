# COVID-19

This repository presents all analyses done for the paper "Gene Expression Signature of SARS-CoV-2 Drives Development of COVID-19". Each folder contains the codes for the each figures shown in the paper.


# Datasets

You will need the following datasets for running all the codes:

1. Generate Table 1 and Figure 1A, B of 25 gene expression signature: input file is [24celllines_4patients_norm.txt](https://github.com/yueli711/COVID-19/blob/master/Table1_Figure1AB_cell24_4/24celllines_4patients_norm.txt), it will generate [cell24_4.csv](https://github.com/yueli711/COVID-19/blob/master/Table1_Figure1AB_cell24_4/cell24_4.csv) during the process.  
2. Generate Figure 1C left: input files are [signature_gene_list_prior_25yueli.csv](https://github.com/yueli711/COVID-19/blob/master/Figure1C_left_5_6_7_16_15positive/signature_gene_list_prior_25yueli.csv), [56716_15positive.txt](https://github.com/yueli711/COVID-19/blob/master/Figure1C_left_5_6_7_16_15positive/56716_15positive.txt), it will generate [cell5_6_7_16_15.csv](https://github.com/yueli711/COVID-19/blob/master/Figure1C_left_5_6_7_16_15positive/cell5_6_7_16_15.csv) during the process.
3. Generate Figure 1C right: cell5_6_7_16_2.csv, signature_gene_list_prior_25yueli.csv, it will generate 56716_2negative.txt during the process.
4. Generate 1D BALF: cell_5_6_7_16_BALF_norm.txt, signature_gene_list_prior_25yueli.csv, it will generate cell_5_6_7_16_balf.csv and combat_cell56716_BALF.csv during the process.
5. Generate 

cell24_4patients_25yueli.R

Use 24 celllines as training set, 4 patients as test set, generate 25 gene expression signature.

5_6_7_16_15positive.R

Use "signature_gene_list_prior_25yueli.csv" to test series15 as positive control.

5_6_7_16_2negative.R

Use "signature_gene_list_prior_25yueli.csv" to test series2 as negative control.

cell_5_6_7_16_BALF_assign.R

Use "signature_gene_list_prior_25yueli.csv" to test BALF data.

cell_5_6_7_16_PBMC_assign.R

Use "signature_gene_list_prior_25yueli.csv" to test  data.

single_cell_integrate.R

input single cell files, filter, CreateSeuratObject, normalize, then integrat 12 samples 

single_cell.R

Use hms_individual_integrated_OK.rds to process single cell data, draw the figures. 

DEseq2_Norm_remove_batch_volcano_finalized.R 

Downloaded the GSE147507_RawReadCounts_Human.tsv from gse147505. 
Delete RSV, IAV and Series16_A549-ACE2_SARS-CoV-2_Rux_2. Normalize,remove batch effect of 36 and 24 celllines, DESeq, differetially expressed genes, volcano plot.

complex_heatmap.R

Draw complex heatmap to test remove batch effect.

All the single cell data were uploaded in the https://www.dropbox.com/home/Research.Data/Yue_Li/singlecell_data.
