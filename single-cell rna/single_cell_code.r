library(dplyr)
library(Seurat)
library(patchwork) 
library(ggplot2)
library(harmony)



setwd("/share/home/gjzzzy/lixin/single_cell")

#产生单细胞数据集
sc.data <- Read10X(data.dir = "/share/home/gjzzzy/lixin/Matrix")
sc <- CreateSeuratObject(counts = sc.data, project = "GZ", min.cells = 3, min.features = 200)

#计算单细胞线粒体基因表达比例
mt_gene<-read.table(file="./MT_genes.txt",header=F,sep="\t",quote="",row.names=2) #extract from gtf file
mt_select_gene<-intersect(rownames(mt_gene),rownames(sc))
sc[["percent.mt"]] <-PercentageFeatureSet(sc, features = mt_select_gene)

#质控与可视化
pdf(file="01_featureViolin_1.pdf",width=14,height=6)
VlnPlot(object = sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster=FALSE)
dev.off()
sc<- subset(sc,subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 20 & nCount_RNA < 10000)#进行子集筛选，根据01图进行调整)
pdf(file="01_featureViolin_2.pdf",width=14,height=6)
VlnPlot(object = sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster=FALSE)
dev.off()

sc<- NormalizeData(sc) #归一化

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000) #查找高变基因

top10 <- head(VariableFeatures(sc), 10)

pdf(file="02_variable_gene.pdf",width=12,height=6)
plot1 <- VariableFeaturePlot(sc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()


all.genes <- rownames(sc)
sc <- ScaleData(sc, features = all.genes) #标准化
sc <- RunPCA(sc, features = VariableFeatures(object = sc))

pdf(file="03.DimHeatmap.pdf",width=8,height=6)
DimHeatmap(sc, dims = 1, cells = 500, balanced = TRUE)
dev.off()


pdf(file="04.ElbowPlot.pdf",width=10,height=6)
ElbowPlot(sc, ndims = 50)
dev.off()

#sc$batch <- "1"
#sc <- RunHarmony(object = sc, group.by.vars=c("batch"), plot_convergence = TRUE, max.iter = 100) #去批次


sc <- FindNeighbors(sc, dims = 1:30)
sc <- FindClusters(sc, resolution = 0.5)
sc <- RunUMAP(sc, dims = 1:30) 
DimPlot(sc, reduction = "umap", pt.size=1, label = TRUE , label.size = 3,raster=FALSE)
name <- paste("05_umap_zong_0.5.pdf",sep="")
ggsave(name,width=10,height=8)
sc <- RunTSNE(sc, dims = 1:30)
TSNEPlot(sc, reduction = "tsne", pt.size=1, label = TRUE , label.size = 3,raster=FALSE)
name <-paste("06_tsne_zong_0.5.pdf",sep="")
ggsave(name,width=10,height=8)

name <-paste("07_zong.rds",sep="")
saveRDS(sc, file = name)

sc.markers <- FindAllMarkers(sc, min.pct = 0.25,logfc.threshold = 0.25)
name <-paste("08_marker_0.5.csv",sep="")
write.csv(sc.markers, file= name)

pdf(file="09.FeaturePlot.pdf",width=15,height=10)
FeaturePlot(sc, features = c("ALB", "TTR.1", "PECAM1","VWF", "DCN","HBBA", "HBA1","MARCO", "TOP2A","CENPF","CENPE"))
dev.off()


re.cluster.ids <- c(
"Hepatocyte",
"Hepatocyte",
"Hepatocyte",
"Hepatocyte",
"Endothelial cell",
"Hepatocyte",
"Hepatocyte",
"Proliferating cell",
"Erythrocyte",
"Hepatic stellate cell",
"Cholangiocyte",
"Kupffer cell",
"Platelet")

names(re.cluster.ids)<- levels(sc)
sc <-RenameIdents(sc,re.cluster.ids)
sc$cell_type <- Idents(sc)

levels(sc$cell_type)


pdf(file="10.umap_labels.pdf",width=10,height=8)
DimPlot(sc, reduction = "umap", group.by = "cell_type", label = TRUE)
dev.off()

name <-paste("11_zong_labels.rds",sep="")
saveRDS(sc, file = name)