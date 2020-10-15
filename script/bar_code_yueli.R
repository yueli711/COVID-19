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
suppressMessages(require(BiocNeighbors))

setwd("~/cov25")
Series15<-read.csv("pathway_activity_testset.csv",header=TRUE)
Series15_cov<-Series15[,2]

setwd("~/cov_negative25")
Series2<-read.csv("pathway_activity_testset.csv",header=TRUE)
Series2_V1<-Series2[,2]

setwd("~/cov_BALF25")
BALF<-read.csv("pathway_activity_testset.csv",header=TRUE)
BALF_V1<-BALF[,2]

setwd("~/cov_PBMC25")
PBMC<-read.csv("pathway_activity_testset.csv",header=TRUE)
PBMC_V1<-PBMC[,2]

Name<-t(c("Series15_HealthyLungBiopsy_2","Series15_HealthyLungBiopsy_1","Series15_COVID19Lung_2","Series15_COVID19Lung_1",
            "Series2_A549_Mock_1","Series2_A549_Mock_2","Series2_A549_Mock_3","Series2_A549_SARS.CoV.2_1", 
            "Series2_A549_SARS.CoV.2_2","Series2_A549_SARS.CoV.2_3",
            "SRR10571724","SRR10571730","SRR10571732","CRR119894","CRR119895","CRR119896",
            "CRR119897","CRR119890","CRR125445","CRR125446","CRR119891","CRR119892","CRR119893"))

Activity<-c(Series15_cov, Series2_V1, BALF_V1,PBMC_V1)


Description<-t(c("Control","Control", "SARS-CoV2 Infected","SARS-CoV2 Infected","Mock","Mock","Mock","SARS-CoV2 Infected","SARS-CoV2 Infected",
            "SARS-CoV2 Infected", "Healthy", "Healthy", "Healthy", "SARS-CoV2 Infected","SARS-CoV2 Infected","SARS-CoV2 Infected","SARS-CoV2 Infected",
            "Healthy", "Healthy", "Healthy", "SARS-CoV2 Infected","SARS-CoV2 Infected","SARS-CoV2 Infected"))
            

Sample<-t(c("Lung Biopsy", "Lung Biopsy","Lung Biopsy","Lung Biopsy","A549","A549", "A549", "A549", "A549", "A549","BALF/PBMC", "BALF/PBMC", "BALF/PBMC", "BALF/PBMC", "BALF/PBMC", "BALF/PBMC","BALF/PBMC",
                         "BALF/PBMC","BALF/PBMC","BALF/PBMC","BALF/PBMC","BALF/PBMC","BALF/PBMC"))

Point<-t(c("Lung_Biopsy/A549","Lung_Biopsy/A549","Lung_Biopsy/A549","Lung_Biopsy/A549","Lung_Biopsy/A549","Lung_Biopsy/A549",
                   "Lung_Biopsy/A549", "Lung_Biopsy/A549", "Lung_Biopsy/A549","Lung_Biopsy/A549","BALF","BALF","BALF","BALF","BALF","BALF","BALF",
                   "PBMC", "PBMC", "PBMC", "PBMC", "PBMC", "PBMC"))

DATA01<-t(rbind(Name,Activity,Description,Sample,Point))

cnames=c("Name","Activity","Description","Sample","Point")

colnames(DATA01)=cnames

write.csv(DATA01,"DATA01.csv")

DATA01<-read.csv("DATA01.csv", header = TRUE)


###Barplots
setwd("/home/li/covid19/result01/bar_code")
DATA<- data.frame(DATA01, check.names=F)
#DATA$Sample<- factor(DATA$Sample, levels= c( "Lung Biopsy","A549","BALF/PMBC"))
#DATA$Point<- factor(DATA$Point, levels= c("Lung Biopsy/A549", "BALF", "PMBC"))

DATA_Graph1<- DATA %>% ggplot(aes(x= Description, y= Activity)) + stat_summary(geom= 'bar', fun= 'mean', fill= 'grey60') + geom_point(position= "jitter", aes(shape= Point), size= 4) +
  facet_wrap(~Sample, scales= "free_x") + theme_classic()+ ylim(NA,1)+ labs(y= "SARS-CoV2 infection activity", x= element_blank())+
  theme_minimal() + theme(axis.title= element_text(face= "bold", size= 16), axis.text= element_text(face= "bold", size= 12), strip.text.x= element_text(size= 16, face= "bold.italic") ,plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  

DATA_Graph2<- DATA %>% ggplot(aes(x= Description, y= Activity)) + stat_summary(geom= 'bar', fun= 'mean', fill= 'grey60') + geom_point(position= "jitter", aes(shape= Point), size= 4) +
  facet_wrap(~Sample, scales= "free_x") + theme_classic()+ ylim(NA,1)+ labs(y= "SARS-CoV2 infection activity",x= element_blank())+
  theme(axis.title= element_text(face= "bold", size= 16), axis.text= element_text(face= "bold", size= 12), strip.text.x= element_text(size= 16, face= "bold.italic") ,plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

DATA_Graph1
                       
DATA_Graph2

###Barplots
setwd("/home/li/covid19/result01/bar_code")
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


##Figure 3
setwd("/home/li/covid19/result01/bar_code")
cov2<-read.csv("SARS-cov2.csv", header=TRUE)
strong_cov2<-droplevels.data.frame(subset(cov2,cov2$Score< -90|cov2$Score > 90))
dim(strong_cov2)
drugs<-c("ivermectin","ibuprofen","irbesartan","olmesartan","losartan","chloroquine","dexamethasone","fluticasone","perindopril","lopinavir","ribavirin","ritonavir","ramipril","tamoxifen","atorvastatin","ketoconazole")
strong_cov2_interesting<-rbind(strong_cov2,cov2[(cov2$Name%in%drugs),])
strong_cov2_interesting<-strong_cov2_interesting[duplicated(strong_cov2_interesting)==F,]
#write.csv(strong_cov2_interesting,"~/Desktop/COVID19/Manuscript/Source Data/interesting_CMAP_cs.csv")
dim(droplevels.data.frame(subset(cov2,cov2$Score< -90&cov2$Type=="cp")))##45 potential targets
#cp
strong_cov2_class<-droplevels(strong_cov2[strong_cov2$Score< -90 &strong_cov2$Description!="-"&strong_cov2$Type=="cp",])
classes<-names(table(strong_cov2_class[,6])[table(strong_cov2_class[,6])>2])
strong_cov2_class_filt<-droplevels(strong_cov2_class[strong_cov2_class$Description%in%classes,])
all_classes<-names(table(strong_cov2[,6])[table(strong_cov2[,6])>3])
all_classes<-all_classes[!all_classes%in%c("-","CD molecules" ,"Mitochondrial respiratory chain complex / Complex I","RNA binding motif (RRM) containing","Zinc fingers, C2H2-type" )]
all_strong_cov2_class_filt<-droplevels(strong_cov2[strong_cov2$Description%in%all_classes,])
# ggplot(all_strong_cov2_class_filt, aes(x=Type, y=Score, fill=Description)) + 
#   geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   theme_classic()
cov2$Description[cov2$Description=="Angiotensin antagonist"]<-"Angiotensin receptor antagonist"
cov2$Description[cov2$Description=="HIV protease inhibitor"]<-"Antiviral"
fig3a<-ggplot(strong_cov2_class_filt, aes(x=Description, y=Score, fill=Description)) + 
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(y= "CS",x= element_blank())+ylim(-100,-90)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
fig3b<-ggplot(cov2[cov2$Name%in%drugs,],aes(x =Name, y=Score, fill=Description))+
  geom_bar(stat="identity", position=position_dodge())+
  labs(y= "CS",x= element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggarrange(fig3a,fig3b,nrow = 2,heights = c(1.5,1),labels = c("a","b"))



