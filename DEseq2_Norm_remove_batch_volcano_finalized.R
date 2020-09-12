library(DESeq2)
library(data.table)
library(DESeq)
library(limma)
library(DESeq)
setwd("/home/li/covid19/result01/DEseq2_Norm_remove_batch_volcano_finalized")
cts<-read.table("gse147507_counts",head=TRUE)
write.csv(cts,"gse147507_counts.csv")
#setup DESeqDataSetFromMatrix
countData<-cts
condition <- factor(c(rep("mock",20),rep("CoV",20)), levels = c("mock","CoV"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("mock","CoV"))
dds$condition <- relevel(dds$condition, ref = "mock")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"gse147507_norm.csv")

#remove batch effect_36celllines
library(sva) #blog.sina.com.cn/s/blog_70a5f5210102wibx.html
library(pamr)
library(limma)
cdata <- read.csv("gse147507_norm_36celllines.csv", header = T,  row.names = 1)
cdata_0<-cdata[apply(cdata==0,1,sum)!=ncol(cdata),]
dim(cdata_0)
cdata_0 <- as.matrix(cdata_0)
csif <- read.table("pheno_gse147507_1_2_5_6_7_16.txt", header = T, row.names = 1)
modcombat = model.matrix(~1, data = csif)
batch = csif$batch
combat_edata = ComBat(dat=cdata_0, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=F)
write.table(combat_edata, "ComBat_data(all)(counts).txt_36", sep = "\t", quote = F)

#remove batch effect_24celllines
cdata <- read.table("gse147507_from40_5_6_7_16_norm.txt", header = T,  row.names = 1)
cdata_0<-cdata[apply(cdata==0,1,sum)!=ncol(cdata),]
dim(cdata_0)
cdata_0 <- as.matrix(cdata_0)
csif <- read.table("pheno_gse147507_5_6_7_16.txt", header = T, row.names = 1)
modcombat = model.matrix(~1, data = csif)
batch = csif$batch
combat_edata = ComBat(dat=cdata_0, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=F)
write.table(combat_edata, "ComBat_data(all)(counts).txt_24", sep = "\t", quote = F)

#volcano plot
cts<-read.table("gse147507_from40_5_6_7_16_norm.txt",head=TRUE)
countData<-round(cts)
condition <- factor(c(rep("mock",12),rep("CoV",12)), levels = c("mock","CoV"))
#coldata<-data.frame(row.names=colnames(countData),condition)
condition
#coldata
featureData <- data.frame(gene=rownames(cts))
name01<-read.table("name.txt",header=FALSE)
name02<-factor(c(t(name01)))
coldata<-data.frame(name=name02, condition)
coldata
dds<-DESeqDataSetFromMatrix(countData=countData,colData=coldata, design=~condition)

dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res
res <- results(dds2, name="condition_CoV_vs_mock")
res <- results(dds2, contrast=c("condition","CoV","mock"))
resultsNames(dds2)
write.csv(res,"gse147507_deg.csv")

DEG1=as.data.frame(res)
DEG1=na.omit(DEG1)
logFC_cutoff=with(DEG1,mean(abs(log2FoldChange))+2*sd(abs(log2FoldChange)))
logFC_cutoff
2^logFC_cutoff   
# 返回结果为 3.67
DEG1$change=as.factor(ifelse(DEG1$pvalue<0.05 &  abs(DEG1$log2FoldChange)>logFC_cutoff,
                             ifelse(DEG1$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))
this_tile=paste('Cutoff for logFC is ',round(logFC_cutoff,3),
                '\nThe number of up gene is ',nrow(DEG1[DEG1$change=='UP',]),
                '\nThe number of down gene is ',nrow(DEG1[DEG1$change=='DOWN',]))
library(ggplot2)
g=ggplot(data=DEG1,
         aes(x=log2FoldChange,y=-log10(pvalue),   #这里将pvalue取负对数
             color=change)) +
  geom_point(alpha=0.4,size=1.75) +     #绘制点图
  theme_set(theme_set(theme_bw(base_size=20))) +
  xlab("log2 fold change")+ylab("-log10 pvalue") +    #轴标签
  ggtitle(this_tile)+theme(plot.title=element_text(size=15,hjust=0.5)) +
  scale_color_manual(values=c('blue','black','red'))   #设定颜色
g

