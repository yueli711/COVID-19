suppressMessages(require(Seurat))
suppressMessages(require(ggplot2))
suppressMessages(require(cowplot))
suppressMessages(require(scater))
suppressMessages(require(scran))
suppressMessages(require(BiocParallel))
suppressMessages(require(BiocNeighbors))
library(data.table)
setwd("/home/li/Transcriptomic_patients/gse145926/test")

#input, CreateSeuratObject, filter, save
a<-"hc_51_02.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
#构建Seurat对象，这里会有个初筛，保证所有基因在至少3个细胞中表达（0.1%细胞数），保证每个细胞至少能检测到200个基因。
pbmc<- CreateSeuratObject(counts = a1, project = "h1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
#用subset函数，质控：筛选检测到基因数目超过2500或低于200的细胞，单个细胞中线粒体基因数目占比超过>5%
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#数据标准化：默认使用数据标准化方法是LogNormalize, 每个细胞总的表达量都标准化到10000，然后log取对数
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#变化基因鉴定：鉴定在细胞间表达高度变化的基因，后续研究需要集中于这部分基因，首先计算每一个基因的均值和方差，并且直接模拟其关系。默认返回2000个基因。
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#数据缩放 线性转换缩放数据，ScaleData()函数可以实现此功能。最终每个基因均值为0，方差为1。结果存放于pbmc[["RNA"]]@scale.data
all.genes <- rownames(pbmc)
#设置参数features是因为ScaleData默认处理前面鉴定的差异基因。这一步怎么做都不会影响到后续pca和聚类，但是会影响做热图
pbmc <- ScaleData(pbmc, features = all.genes)
#移除影响方差的因素
#pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
#对缩放后的数据进行PCA分析，默认使用前面鉴定表达变化大的基因。使用features参数可以重新定义数据集。
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "hc51.rds")

a<-"hc_52_02.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
pbmc<- CreateSeuratObject(counts = a1, project = "h2", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "hc52.rds")

a<-"hc_100_02.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
pbmc<- CreateSeuratObject(counts = a1, project = "h3", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "hc100.rds")

a<-"mild_141_02.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
pbmc<- CreateSeuratObject(counts = a1, project = "m1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "mild141.rds")


a<-"mild_142_02.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
pbmc<- CreateSeuratObject(counts = a1, project = "m2", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "mild142.rds")

a<-"mild_144_02.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
pbmc<- CreateSeuratObject(counts = a1, project = "m3", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "mild144.rds")

a<-"severe_143_02.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
pbmc<- CreateSeuratObject(counts = a1, project = "s1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "severe143.rds")

a<-"severe_145_02.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
pbmc<- CreateSeuratObject(counts = a1, project = "s2", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "severe145.rds")

a<-"severe_146_02.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
pbmc<- CreateSeuratObject(counts = a1, project = "s3", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "severe146.rds")

a<-"severe_148_02.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
pbmc<- CreateSeuratObject(counts = a1, project = "s4", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "severe148.rds")

a<-"severe_149_02.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
pbmc<- CreateSeuratObject(counts = a1, project = "s5", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "severe149.rds")

a<-"severe_152_02.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
pbmc<- CreateSeuratObject(counts = a1, project = "s6", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "severe152.rds")
#input readRDS
hc51<-readRDS(file="hc51.rds")
hc52<-readRDS(file="hc52.rds")
hc100<-readRDS(file="hc100.rds")
mild141<-readRDS(file="mild141.rds")
mild142<-readRDS(file="mild142.rds")
mild144<-readRDS(file="mild144.rds")
severe143<-readRDS(file="severe143.rds")
severe145<-readRDS(file="severe145.rds")
severe146<-readRDS(file="severe146.rds")
severe148<-readRDS(file="severe148.rds")
severe149<-readRDS(file="severe149.rds")
severe152<-readRDS(file="severe152.rds")
#setup tech and celltype
hc51<-RenameCells(hc51,add.cell.id="hc51",for.merge=T)
hc51@meta.data$tech<-"healthy_control"
hc51@meta.data$celltype<-"healthy_control_51"

hc52<-RenameCells(hc52,add.cell.id="hc52",for.merge=T)
hc52@meta.data$tech<-"healthy_control"
hc52@meta.data$celltype<-"healthy_control_52"

hc100<-RenameCells(hc100,add.cell.id="hc100",for.merge=T)
hc100@meta.data$tech<-"healthy_control"
hc100@meta.data$celltype<-"healthy_control_100"

mild141<-RenameCells(mild141,add.cell.id="mild141",for.merge=T)
mild141@meta.data$tech<-"mild"
mild141@meta.data$celltype<-"mild_141"

mild142<-RenameCells(mild142,add.cell.id="mild142",for.merge=T)
mild142@meta.data$tech<-"mild"
mild142@meta.data$celltype<-"mild_142"

mild144<-RenameCells(mild144,add.cell.id="mild144",for.merge=T)
mild144@meta.data$tech<-"mild"
mild144@meta.data$celltype<-"mild_144"

severe143<-RenameCells(severe143,add.cell.id="severe143",for.merge=T)
severe143@meta.data$tech<-"severe"
severe143@meta.data$celltype<-"severe_143"

severe145<-RenameCells(severe145,add.cell.id="severe145",for.merge=T)
severe145@meta.data$tech<-"severe"
severe145@meta.data$celltype<-"severe_145"

severe146<-RenameCells(severe146,add.cell.id="severe146",for.merge=T)
severe146@meta.data$tech<-"severe"
severe146@meta.data$celltype<-"severe_146"

severe148<-RenameCells(severe148,add.cell.id="severe148",for.merge=T)
severe148@meta.data$tech<-"severe"
severe148@meta.data$celltype<-"severe_148"

severe149<-RenameCells(severe149,add.cell.id="severe149",for.merge=T)
severe149@meta.data$tech<-"severe"
severe149@meta.data$celltype<-"severe_149"

severe152<-RenameCells(severe152,add.cell.id="severe152",for.merge=T)
severe152@meta.data$tech<-"severe"
severe152@meta.data$celltype<-"severe_152" 
#merge
h1_2<-merge(hc51,hc52)
h123<-merge(h1_2,hc100)
m1_2<-merge(mild141,mild142)
m123<-merge(m1_2,mild144)
s1_2<-merge(severe143,severe145)
s123<-merge(s1_2,severe146)
s4_5<-merge(severe148,severe149)
s456<-merge(s4_5,severe152)
s123456<-merge(s123,s456)
hm<-merge(h123,m123)
hms<-merge(hm,s123456)
saveRDS(hms, file="hms_before_integrate.rds")
#before integrate
hms[["percent.mt"]] <- PercentageFeatureSet(hms, pattern = "^Mt-")
VlnPlot(hms, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = hms, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + 
  NoLegend()
plot_grid(p1,p2)
#integrate
pancreas.list <- SplitObject(pancreas, split.by = "celltype")
for (i in 1: length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)
}
reference.list <- pancreas.list[c("healthy_control_51","healthy_control_52","healthy_control_100","mild_141",
                                  "mild_142","mild_144", "severe_143", "severe_145", "severe_146", "severe_148",
                                  "severe_149", "severe_152")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "hms_after_integrated.rds")
