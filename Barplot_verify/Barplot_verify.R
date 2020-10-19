library(stringr)
library(sva)
library(ASSIGN)
library(data.table)
library(readxl)
library(ggpubr)
suppressMessages(require(Seurat))
suppressMessages(require(ggplot2))
suppressMessages(require(cowplot))
suppressMessages(require(scater))
suppressMessages(require(scran))
suppressMessages(require(BiocParallel))


###Barplots
setwd("~/bar_code")
test_validation<- as.data.frame(read_excel("DATA.xlsx", sheet= "Combined"))
test_validation$Sample<- factor(test_validation$Sample, levels= c( "Lung Biopsy","A549","BALF/PBMC"))
test_validation$Point<- factor(test_validation$Point, levels= c("Lung Biopsy","A549", "BALF", "PBMC"))
t.test(test_validation$Activity[c(11:13,18:20)],test_validation$Activity[c(14:16,21:23)])
DATA_Graph2<- test_validation %>% ggplot(aes(x= Description, y= Activity,fill=Description)) + 
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9","brown2"))+
  stat_summary(geom= 'bar', fun= 'mean') + geom_point(position= "jitter", aes(shape= Point), size= 3) +
  facet_wrap(~Sample, scales= "free_x") + theme_classic()+ ylim(NA,1)+ labs(y= "Predicted \nSARS-CoV-2 infection activity",x= element_blank())+
  theme(axis.title= element_text(face= "bold", size= 12), axis.text= element_text(face= "bold", size= 12), strip.text.x= element_text(size= 12, face= "bold.italic") ,plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ann_text<-data.frame(Description="Healthy",Activity=0.995,lab="p<0.0001",Sample = factor("BALF/PBMC",levels = c("Lung Biopsy","A549","BALF/PBMC")))
DATA_Graph2+geom_text(data=ann_text, label=ann_text$lab, fontface="bold")
