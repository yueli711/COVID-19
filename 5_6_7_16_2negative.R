library(stringr)
library(sva)
library(ASSIGN)
library(data.table)
setwd("~/covid19/result01/5_6_7_16_2negative")
cell5_6_7_16_2<-read.table("56716_2negative.txt",header=TRUE)
write.csv(cell5_6_7_16_2,"cell5_6_7_16_2.csv")
mock_cov<-read.csv("cell5_6_7_16_2.csv", header = TRUE)
mock_cov<-(log2(data.frame(mock_cov,check.names = F, row.names=1)+1))
mock_cov_filt<-(mock_cov[apply(mock_cov==0,1,mean)<0.1,])
bat <-as.data.frame(cbind(
  colnames(cell5_6_7_16_2),
  c(rep(1, 3),rep(2, 3),rep(3, 3),rep(4, 3),rep(5, 3),rep(1, 3),rep(2, 3),rep(3, 3),rep(4, 3),rep(5, 3)),
  c(rep(1, 15), rep(2, 15))))
mod<-model.matrix(~as.factor(bat[,3]), data=bat)
combat_mock_cov<-ComBat(dat = as.matrix(mock_cov_filt,),batch = (bat[,2]),mod=mod,par.prior = T)

c_mock<-as.matrix(combat_mock_cov[,c(1:12)])
c_cov<-as.matrix(combat_mock_cov[,c(16:27)])
test<-as.matrix(combat_mock_cov[,c(13:15,28:30)])

trainingLabela <- list(control=list(mock=1:12),cov=13:24)
genelist_25<-read.csv("signature_gene_list_prior_25yueli.csv")

basedir<-getwd()
sub_dir <- paste(basedir,paste("cov_negative", 25, sep=""),sep='/')
dir.create(sub_dir)
set.seed(1220)
assign.wrapper(
  trainingData = cbind(c_mock, c_cov),
  testData = test,
  trainingLabel = trainingLabela,
  geneList = list(genelist_25$X),
  adaptive_B = T,
  adaptive_S = F,
  outputDir = sub_dir,
  p_beta = 0.01,
  theta0 = 0.05,
  theta1 = 0.9,
  iter = 2000,
  burn_in = 1000)

