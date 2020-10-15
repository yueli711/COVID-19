library(stringr)
library(sva)
library(ASSIGN)
library(data.table)
setwd("~/covid19/result01/cell24_4")
cell24_4<-read.table("24celllines_4patients_norm.txt")
write.csv(cell24_4,"cell24_4.csv")

cell24_4<-read.csv("cell24_4.csv",header=TRUE)
cell24_4<-(log2(data.frame(cell24_4,check.names = F, row.names=1)+1))
colnames(cell24_4)
rownames(cell24_4)
plot(hclust(dist(t(cell24_4)),method="complete"))
cell24_4_filt<-(cell24_4[apply(cell24_4==0,1,mean)<0.1,])
pca<-prcomp(t(cell24_4_filt))
plot(pca)
{plot(pca$x[,1],pca$x[,2])
  points(pca$x[1:14,1],pca$x[1:14,2],col=2,pch=2)
  points(pca$x[15:28,1],pca$x[15:28,2],col=3,pch=2)
}
which(pca$x[,1]< -100)
which(pca$x[,1]< -20)
which(pca$x[,1]< 0)

bat <-
  as.data.frame(cbind(
    colnames(cell24_4_filt),
    c(
      rep(1, 3),
      rep(2, 3),
      rep(3, 3),
      rep(4, 3),
      rep(5, 2),
      rep(1, 3),
      rep(2, 3),
      rep(3, 3),
      rep(4, 3),
      rep(5, 2)
    ),
    c(rep(1, 14), rep(2, 14))
  ))

mod<-model.matrix(~as.factor(bat[,3]), data=bat)

combat_cell24_4<-ComBat(dat = as.matrix(cell24_4_filt,),batch = (bat[,2]),mod=mod,par.prior = T)
##PCA post combat
pca<-prcomp(t(combat_cell24_4))
plot(pca)
{plot(pca$x[,1],pca$x[,2])
  points(pca$x[1:14,1],pca$x[1:14,2], main="Top 2 PCs",col=2)
  points(pca$x[15:28,1],pca$x[15:28,2], main="Top 2 PCs",col=3)
}

which(pca$x[,1]< -50)
which(pca$x[,2]< -50)

c_mock<-as.matrix(combat_cell24_4[,c(1:12)])
c_cov<-as.matrix(combat_cell24_4[,c(15:26)])
test<-as.matrix(combat_cell24_4[,c(13:14,27:28)])
trainingLabela <- list(control=list(mock=1:12),cov=13:24)
basedir<-getwd()
sub_dir <- paste(basedir,paste("cov", 25, sep=""),sep='/')
dir.create(sub_dir)
set.seed(1220)
assign.wrapper(
  trainingData = cbind(c_mock, c_cov),
  testData = test,
  trainingLabel = trainingLabela,
  geneList = NULL,
  n_sigGene = 25,
  adaptive_B = T,
  adaptive_S = F,
  outputDir = sub_dir,
  p_beta = 0.01,
  theta0 = 0.05,
  theta1 = 0.9,
  iter = 2000,
  burn_in = 1000)
