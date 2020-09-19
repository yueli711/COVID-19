library(stringr)
library(sva)
library(ASSIGN)
library(data.table)
setwd("/home/li/covid19/result01/5_6_7_16_PBMC")
cell_5_6_7_16_PBMC<-read.table("cell_5_6_7_16_PBMC_norm.txt",head=TRUE)
write.csv(cell_5_6_7_16_PBMC,"cell_5_6_7_16_PBMC.csv")
cells_4_PBMC<-read.csv("cell_5_6_7_16_PBMC.csv",header=TRUE)
cells_4_PBMC<-(log2(data.frame(cells_4_PBMC,check.names = F, row.names=1)+1))
colnames(cells_4_PBMC)
#rownames(cells_4_PBMC)
plot(hclust(dist(t(cells_4_PBMC)),method="complete"))##samples are rather clustering by type than the COV2 infection status
#Series 15 is furthest from the others and that makes sense since these are from patients rather than cell lines. We should exclude Series 15 from the signature generation dataset.

##filter all zeroes
cells_4_PBMC_filt<-(cells_4_PBMC[apply(cells_4_PBMC==0,1,mean)<0.1,])
#precombat PCA
pca<-prcomp(t(cells_4_PBMC_filt))
plot(pca)
{plot(pca$x[,1],pca$x[,2])
  points(pca$x[1:15,1],pca$x[1:15,2],col=2,pch=2)
  points(pca$x[16:30,1],pca$x[16:30,2],col=3,pch=2)
}
which(pca$x[,1]< -100)
which(pca$x[,1]< -20)
which(pca$x[,1]< 0)
###we need to batch adjust based on the cell types
##Series5_A549-1 Series6_A549.ACE2-2 Series7_Calu3-3 Series16_A549.ACE2-4, patients-5
bat <-as.data.frame(cbind(
    colnames(cells_4_PBMC_filt),
c(rep(1, 3),rep(2, 3),rep(3, 3),rep(4, 3),rep(5, 3),rep(1, 3),rep(2, 3),rep(3, 3),rep(4, 3),rep(6, 3)),
 c(rep(1, 15), rep(2, 15))))
mod<-model.matrix(~as.factor(bat[,3]), data=bat)
combat_mock_cov<-ComBat(dat = as.matrix(cells_4_PBMC_filt,),batch = (bat[,2]),mod=mod,par.prior = T)
##PCA post combat
pca<-prcomp(t(combat_mock_cov))
plot(pca)
{plot(pca$x[,1],pca$x[,2])
  points(pca$x[1:15,1],pca$x[1:15,2], main="Top 2 PCs",col=2)
  points(pca$x[16:30,1],pca$x[16:30,2], main="Top 2 PCs",col=3)}
which(pca$x[,1]< -50)
which(pca$x[,2]< -50)

plot(hclust(dist(t(combat_mock_cov)),method="complete"))

########running assign with the best 25 genes found in the cell line data########
#combat_mock_cov<-read.table("combat_cell56716_PBMC.txt",header=TRUE)
c_mock<-as.matrix(combat_mock_cov[,c(1:12)])
c_cov<-as.matrix(combat_mock_cov[,c(16:27)])
test<-as.matrix(combat_mock_cov[,c(13:15,28:30)])
trainingLabela <- list(control=list(mock=1:12),cov=13:24)
genelist_25<-read.csv("signature_gene_list_prior_25yueli.csv")
basedir<-getwd()
sub_dir <- paste(basedir,paste("cov", 25, sep=""),sep='/')
dir.create(sub_dir)
set.seed(1220)
assign.wrapper(
  trainingData = cbind(c_mock, c_cov),
  testData = test,
  trainingLabel = trainingLabela,
  geneList = list(genelist_25$X),
  n_sigGene = 25,
  adaptive_B = T,
  adaptive_S = F,
  outputDir = sub_dir,
  p_beta = 0.01,
  theta0 = 0.05,
  theta1 = 0.9,
  iter = 2000,
  burn_in = 1000)
