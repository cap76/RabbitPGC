library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)
library("pheatmap")
library(scran)
set.seed(1)

#Tediously add all the folders can shortcut this later when we have finalised everything.
saveext = "~/Desktop/Thorsten/FINAL/Rabbit_Analysis_with_new_datav3/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))
dir.create(paste(saveext,"/Docs/",sep=""))
dir.create(paste(saveext,"/Monocle/",sep=""))


oldrabbit <- readRDS(paste("~/Desktop/Thorsten/FINAL/Rabbit_Analysis_Fin5/OldAnnotations.rds",sep=""))


#Rabbit 2
RabbitBS2<-read.table("/Users/christopherpenfold/Desktop/Ara/Rabbit_meta_2.csv",sep=",", header = T, row.names=1)
#R2 <- read.table("/Users/christopherpenfold/Desktop/Ara/featurecountsAll_base_rabbit_v2.csv",sep=",",header = T, row.names=1)
R2 <- read.table("/Users/christopherpenfold/Desktop/Ara/featurecountsAll_base_test.txt",sep="\t",header = T, row.names=1)
RabitSCE <- SingleCellExperiment( list(counts = as.matrix(R2) ) )
clusters <- quickCluster(RabitSCE, min.size=30)
Rsce <- computeSumFactors(RabitSCE, clusters=clusters)
Rsce <- assay(Rsce, "counts")
#sce <- logNormCounts(Rsce)
Genes<-read.table("/Users/christopherpenfold/Desktop/Misc/RabbitGenes2.csv",sep=",",header = T)
counts1 <- assay(RabitSCE, "counts")
rownames(counts1) <- Genes$Gene
rabbit_data2 <- CreateSeuratObject(counts = counts1, assay = "RNA",min.cells = 0, min.features = 0)
rabbit_data2$species <- "1) Rabbit"
rabbit_data2$species2 <- "1) Rabbit_2"
Idents(rabbit_data2) <- RabbitBS2$Type2
rabbit_data2 <- NormalizeData(rabbit_data2, verbose = FALSE)
rabbit_data2 <- FindVariableFeatures(rabbit_data2, selection.method = "vst", nfeatures = 20000)


#ENSOCUG00000012651 <- "SOX17"
#ENSOCUG00000010996 <- "TFAP2C"
#ENSOCUG00000013783 <- "NANOG"
#ENSOCUG00000014595 <- "TBX3"
#ENSOCUG00000012254 <- "SNAI2"
#ENSOCUG00000004992 <- "TDGF1"
#rownames(RabitSCE)[which(rownames(RabitSCE)=="ENSOCUG00000012254")] <- "SNAI2"
#rownames(RabitSCE)[which(rownames(RabitSCE)=="ENSOCUG00000004992")] <- "SNAI2"

#rabbit_data2 <- CreateSeuratObject(counts = Rsce, assay = "RNA",min.cells = 0, min.features = 0)
#rabbit_data2$species <- "1) Rabbit"
#rabbit_data2$species2 <- "1) Rabbit_2"
#rabbit_data2 <- subset(rabbit_data2, subset = nFeature_RNA > 0)
#rabbit_data2 <- NormalizeData(rabbit_data2, verbose = FALSE)
#rabbit_data2 <- FindVariableFeatures(rabbit_data2, selection.method = "vst", nfeatures = 20000)
#Idents(rabbit_data2) <- RabbitBS2$Type#

#Human
hBS<-read.table("/Users/christopherpenfold/Desktop/Aracely/humanIDs.csv",sep=",",header = T, row.names=1)
raw_countsA1<-read.table("/Users/christopherpenfold/Desktop/Aracely/featurecountsLi.csv",sep=",",header = T, row.names=1)
raw_countsA1 <- raw_countsA1[,2:dim(raw_countsA1)[2]]
human_dataA1 <- CreateSeuratObject(counts = na.omit(raw_countsA1), assay = "RNA",min.cells = 0, min.features = 0)
human_dataA1$species <- "4) Human (in vitro)"
human_dataA1$species2 <- "4) Human (in vitro)"
Idents(human_dataA1) <- hBS$Type5
human_dataA1 <- subset(human_dataA1, subset = nFeature_RNA > 0)
human_dataA1 <- NormalizeData(human_dataA1, verbose = FALSE)
human_dataA1 <- FindVariableFeatures(human_dataA1, selection.method = "vst", nfeatures = 20000)
human_dataA1$Cells <- Idents(human_dataA1)

saveext1 = "~/Desktop/Thorsten/FINAL/AllPlatesCCA_redooriginal/"
saveext2 = "~/Desktop/Thorsten/FINAL/AllPlatesCCA_3_dataset/"

Am <- readRDS(file = paste(saveext1,"/AmLC.rds",sep="") )
Am <- intersect(Am,colnames(human_dataA1))
ICM <- readRDS(paste(saveext2,"/pICM",".rds",sep="") )
Tb <- readRDS(paste(saveext2,"/pTb",".rds",sep="") )
Hyp <- readRDS(paste(saveext2,"/pHyp",".rds",sep="") )
TbCS5 <- readRDS(paste(saveext2,"/pTbCS5",".rds",sep="") )

Idents(human_dataA1,cells=TbCS5) <- "Tb_CS5"
Idents(human_dataA1,cells=ICM) <- "ICM_CS3"
Idents(human_dataA1,cells=Hyp) <- "Hyp_CS3"
Idents(human_dataA1,cells=Tb) <- "Tb_CS4"
Idents(human_dataA1,cells=Am) <- "Am"
Am_CS5 <- colnames(human_dataA1)[which( ( Idents(human_dataA1)=="Am") &  human_dataA1$Cells=="EmDisc_CS5" )]
Am_CS6 <- colnames(human_dataA1)[which( ( Idents(human_dataA1)=="Am") &  human_dataA1$Cells=="EmDisc_CS6" )]
Am_CS52 <- colnames(human_dataA1)[which( ( Idents(human_dataA1)=="Am") &  human_dataA1$Cells=="EmDiscPS_CS5" )]
Am_CS62 <- colnames(human_dataA1)[which( ( Idents(human_dataA1)=="Am") &  human_dataA1$Cells=="EmDiscPS_CS6" )]
Idents(human_dataA1,cells=Am_CS5) <- "Am_CS5"
Idents(human_dataA1,cells=Am_CS6) <- "Am_CS6"
Idents(human_dataA1,cells=Am_CS62) <- "Am_CS6"


#Cyno
cyBS<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/cynomolgous/cyKey.csv",sep=",", header = T, row.names=1)
raw_counts1<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/cynomolgous/featurecountsNakamura3.csv",sep=",",header = T, row.names=1)
cynomolgous_data <- CreateSeuratObject(counts = raw_counts1[,which(cyBS$Emb5>0)], assay = "RNA",min.cells = 0, min.features = 0)
Idents(cynomolgous_data) <- cyBS$LABEL_4[which(cyBS$Emb5>0)]
cynomolgous_data$species <- "2) Cynomolgous"
cynomolgous_data$species2 <- "2) Cynomolgous"
cynomolgous_data <- subset(cynomolgous_data, subset = nFeature_RNA > 0)
cynomolgous_data <- NormalizeData(cynomolgous_data, verbose = FALSE)
cynomolgous_data <- FindVariableFeatures(cynomolgous_data, selection.method = "vst", nfeatures = 20000)


#Load marmoset data key
RabbitBS<-read.table("/Users/christopherpenfold/Desktop/Ara/Rabbit_meta.tsv",sep="\t", header = T, row.names=1)

D1 <-readRDS('/Users/christopherpenfold/Desktop/Ara/Rabbit_SCE.RDS')
Genes<-read.table("/Users/christopherpenfold/Desktop/Misc/RabbitGenes1.csv",sep=",",header = T)
counts1 <- assay(D1, "counts")
rownames(counts1) <- Genes$Gene
rabbit_data <- CreateSeuratObject(counts = counts1, assay = "RNA",min.cells = 0, min.features = 0)
rabbit_data$species <- "1) Rabbit"
rabbit_data$species2 <- "1) Rabbit"
Idents(rabbit_data) <- as.character(RabbitBS$TimePoint)
rabbit_data$ID1 <- as.character(RabbitBS$TimePoint)
rabbit_data$ID2 <- as.character(RabbitBS$TimePoint)
rabbit_data$ID3 <- as.character(RabbitBS$TimePoint)
rabbit_data <- subset(rabbit_data, idents = c("E6","E7"))
rabbit_data <- NormalizeData(rabbit_data, verbose = FALSE)
rabbit_data <- FindVariableFeatures(rabbit_data, selection.method = "vst", nfeatures = 20000)


TF<-read.table("/Users/christopherpenfold/Desktop/Thorsten/Leaving_package/Dimensionality\ reduction\ techniques\ smart-seq2\ -\ Boroviaklab\ data/Human_TF_MasterList_v1_02.csv",sep=",",header = F, row.names=2)


Markers <- c("SALL4",
             "ESRRB",
             "KLF5",
             "KLF4",
             "DPPA5",
             "TBX3",
             "POU5F1",
             "NANOG",
             "SOX2",
             "PRDM14",
             "T",
             "MIXL1",
             "EVX1",
             "TFAP2C",
             "TFAP2A",
             "GATA3",
             "HAND1",
             "CD36",
             "APOA1",
             "APOB",
             "FOXA2",
             "SOX17",
             "PRDM1",
             "EOMES",
             "LHX1",
             "PITX2",
             "NANOS3",
             "PDPN",
             "ALPL")

nrabbit_data3 <- rabbit_data2 #merge(rabbit_data, y = c(rabbit_data2), project = "merged")
nrabbit_data3 <- FindVariableFeatures(nrabbit_data3, selection.method = "vst", nfeatures = 3000)
nrabbit_data3 <- ScaleData(nrabbit_data3, verbose = FALSE)
nrabbit_data3 <- RunPCA(nrabbit_data3, npcs = 30, verbose = FALSE)
nrabbit_data3 <- RunUMAP(nrabbit_data3, reduction = "pca", dims = 1:20)
nrabbit_data3 <- RunTSNE(nrabbit_data3, reduction = "pca", dims = 1:20)
nrabbit_data3 <- FindNeighbors(nrabbit_data3, reduction = "pca", dims = 1:20)

DimPlot(nrabbit_data3, reduction = "umap",  label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_nrab3",".pdf",sep=""),width = 8, height = 8)


FeaturePlot(nrabbit_data3, features = "SOX17", split.by = "species")
ggsave(filename=paste(saveext,"/RabbitSOX17.pdf",sep=""),width = 8, height = 8,limitsize = FALSE)
 

DimPlot(nrabbit_data3, reduction = "tsne",  label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_nrab3",".pdf",sep=""),width = 8, height = 8)


#PGC with #2F3C56 lightsalmon1, PGCLC with orangered #ff4500, PSC with lightskyblue4 #607b8b, and PSC FGF low with lightskyblue2 ## a4d3ee, if possible.


#cType <- c("cPSC","rPSC_low","rPSC","rPGCLC","rPGC","cPGC","cExMes","cGastEnd","cEmDisc3","cGast","cEmDisc2","cHyp","cICM","cPTb","cEmDisc1","cEpi","cTb","cVE","cSYS","rVE","rEmDisc3","rEmDisc1","rHyp","rPTb","rTb","rEmDisc2","rGast","rICM","hPTb","hExMes","hGast","hEmDisc1","hEmDisc3","hGastEndo","hSYS","hEmDisc2","hVE","hEpi","hTb","hICM","hHyp")  
#BaseCol<-c("#607b8b","#a4d3ee","#607b8b","#ff4500","#2F3C56","#2F3C56","#01665e","#cb181d","#abd9e9","#cb181d","#74add1","#fd8d3c","#2d004b","#1b7837","#4575b4","#313695","#7fbc41","#feb24c","#feb24c","#fe9929","#fe9929","#2166ac","#ec7014","#00441b","#4d9221","#4393c3","#a50f15","#053061","#5aae61","#35978f","#ef3b2c","#045a8d","#3690c0","#ef3b2c","#fdd0a2","#0570b0","#fdd0a2","#023858","#b8e186","#542788","#fdae6b")



cType <- c("rEpi","rEmDisc1","rEmDisc2","rEmDisc3","rEmDisc4","rGast1","rGast2","rTb","rPTb","rHyp","rVE","rPGC","rPGCLC","rPSC","rPSC_low",
"cEpi","cEmDisc1","cEmDisc2","cEmDisc3","cGast1","cGast2","cTb","cPTb","cHyp","cVE","cExMes","cPGC","cICM","cPSC",
"hEpi","hEmDisc1","hEmDisc2","hEmDisc3","hEmDisc4","hEmDiscEnd","hGast1","hGast2","hTb","hPTb","hHyp","hVE","hExMes","hPGC","hICM","hAm")
BaseCol<-c("#053061","#2166ac","#4393c3","#92c5de","#a5cfe4","#db242e","#a50f15","#00441b","#4d9221","#ec7014","#fe9929","#e43a9c","#cf324a","#607b8b","#a4d3ee",
"#313695","#4575b4","#74add1","#abd9e9","#d72830","#cb181d","#1b7837","#7fbc41","#fd8d3c","#feb24c","#01665e","#c51b7d","#2d004b","#607b8b",
"#023858","#045a8d","#0570b0","#3690c0","#61abd3","#69afcf","#cd4832","#ef3b2c","#5aae61","#b8e186","#fdae6b","#fdd0a2","#35978f","#de77ae","#542788","#5a535e")

#cType <- c("cPSC","rPSC_low","rPSC","rPGCLC","rPGC","cPGC","cPGC","cExMes","cGastEnd","cEmDisc3","cGast","cEmDisc2","cHyp","cICM","cPTb","cEmDisc1","cEpi","cTb","cVE","cSYS","rVE","rEmDisc3","rEmDisc1","rHyp","rPTb","rTb","rEmDisc2","rGast","rICM","hPTb","hExMes","hGast","hEmDisc1","hEmDisc3","hGastEnd","hSYS","hEmDisc2","hVE","hEpi","hTb","hICM","hHyp")  
#BaseCol<-c("#607b8b","#a4d3ee","#607b8b","#ff4500","#2F3C56","#2F3C56","#c51b7d","#01665e","#cb181d","#abd9e9","#cb181d","#74add1","#fd8d3c","#2d004b","#1b7837","#4575b4","#313695","#7fbc41","#feb24c","#feb24c","#fe9929","#fe9929","#2166ac","#ec7014","#00441b","#4d9221","#4393c3","#a50f15","#053061","#5aae61","#35978f","#ef3b2c","#045a8d","#3690c0","#ef3b2c","#fdd0a2","#0570b0","#fdd0a2","#023858","#b8e186","#542788","#fdae6b")


#DimPlot(nrabbit_data3, reduction = "umap",  label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_nrab3",".pdf",sep=""),width = 8, height = 8)

#DimPlot(nrabbit_data3, reduction = "pca",  label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/PCA_split_nrab3",".pdf",sep=""),width = 8, height = 8)

#nrabbit_data3 <- rabbit_data2 #merge(rabbit_data, y = c(rabbit_data2), project = "merged")
#nrabbit_data3 <- FindVariableFeatures(nrabbit_data3, selection.method = "vst", nfeatures = 20000)
#nrabbit_data3 <- ScaleData(nrabbit_data3, verbose = FALSE)
#nrabbit_data3 <- RunPCA(nrabbit_data3, npcs = 30, verbose = FALSE)
#nrabbit_data3 <- RunUMAP(nrabbit_data3, reduction = "pca", dims = 1:20)
#nrabbit_data3 <- RunTSNE(nrabbit_data3, reduction = "pca", dims = 1:20)
#nrabbit_data3 <- FindNeighbors(nrabbit_data3, reduction = "pca", dims = 1:20)

#DimPlot(nrabbit_data3, reduction = "tsne",  label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/TSNE_nrab4",".pdf",sep=""),width = 8, height = 8)

#DimPlot(nrabbit_data3, reduction = "umap",  label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_nrab4",".pdf",sep=""),width = 8, height = 8)

#DimPlot(nrabbit_data3, reduction = "pca",  label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/PCA_split_nrab4",".pdf",sep=""),width = 8, height = 8)

#nrabbit_data3 <- FindClusters(nrabbit_data3, resolution = 1.5)

#DimPlot(nrabbit_data3, reduction = "tsne",  label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/TSNE_rab3","_cl.pdf",sep=""),width = 8, height = 8)

#DimPlot(nrabbit_data3, reduction = "umap",  label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_rab3","_cl.pdf",sep=""),width = 8, height = 8)

#DimPlot(nrabbit_data3, reduction = "pca",  label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/PCA_split_rab3","_cl.pdf",sep=""),width = 8, height = 8)


#library("pheatmap")
#avexp  <- AverageExpression(object = nrabbit_data3, slot = "data")

#a <- avexp$RNA
#a <- a[Markers,c("1","0","2","3")]



#intgenes <- rownames(avexp$RNA)
#b <- avexp$RNA[intgenes,c("1","0","2","3")]

#C1 <- cor(b, method = "pearson")

#redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
#mat_breaks <- seq(0, 2, length.out = 20)
#pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE,  filename = "~/Desktop/HM_ara_logged.pdf",,width=10,height=16)

#mat_breaks <- seq(-1, 1, length.out = 20)
#pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, scale="row", filename = "~/Desktop/HM_ara_logged_scaled.pdf",width=10,height=16)

#mat_breaks <- seq(0, 2, length.out = 20)
#pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks,gaps_row=c(6,11,16,21,29), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/HM_ara_logged",".pdf",sep=""),width=8,height=16)

#mat_breaks <- seq(-1, 1, length.out = 20)
#pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks,gaps_row=c(6,11,16,21,29), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, scale="row", filename = paste(saveext,"/DimRed/HM_ara_logged_scaled",".pdf",sep=""),width=8,height=16)


#mat_breaks <- seq(0.8, 0.92, length.out = 20)
#pheatmap(C1[c("1","0","2","3"),c("1","0","2","3")],color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/Sample_corr",".pdf",sep=""),width=4,height=4)



#FeaturePlot(nrabbit_data3, features = "SOX17", cols =  c("lightgrey", "darkblue"), pt.size = 5)
#ggsave(filename=paste(saveext,"/Markers/Marker_full_SOX17.pdf",sep=""),width = 10, height = 10,limitsize = FALSE)


#FeaturePlot(nrabbit_data3, features = "NANOS3", cols =  c("lightgrey", "darkblue"), pt.size = 5)
#ggsave(filename=paste(saveext,"/Markers/Marker_full_NANOS3.pdf",sep=""),width = 10, height = 10,limitsize = FALSE)


#FeaturePlot(nrabbit_data3, features = "SOX2", cols =  c("lightgrey", "darkblue"), pt.size = 5)
#ggsave(filename=paste(saveext,"/Markers/Marker_full_SOX2.pdf",sep=""),width = 10, height = 10,limitsize = FALSE)



#saveRDS(rabbit_data, file = paste(saveext,"rabbit_data.rds",sep=""))
#Already regressed out cell cycle? Can't help to redo to be sure.
#humaninvit2_data <- CellCycleScoring(humaninvit2_data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#RidgePlot(humaninvit2_data, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
#ggsave(filename=paste(saveext,"Human_invitro2_ridge",".pdf",sep=""),width = 8, height = 8)
#humainvit2_data <- ScaleData(humaninvit2_data, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(humaninvit2_data))
#humainvit2_data <- ScaleData(humaninvit2_data) #, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(humaninvit2_data))
#saveRDS(humaninvit2_data, file = paste(saveext,"humaninvit2_data.rds",sep=""))


#rabbit_data3 <- merge(rabbit_data, y = c(rabbit_data2), project = "merged")
#rabbit_data3 <- FindVariableFeatures(rabbit_data3, selection.method = "vst", nfeatures = 3000)
#rabbit_data3 <- ScaleData(rabbit_data3, verbose = FALSE)
#rabbit_data3 <- RunPCA(rabbit_data3, npcs = 30, verbose = FALSE)
#rabbit_data3 <- RunUMAP(rabbit_data3, reduction = "pca", dims = 1:20)
#rabbit_data3 <- RunTSNE(rabbit_data3, reduction = "pca", dims = 1:20)
#rabbit_data3 <- FindNeighbors(rabbit_data3, reduction = "pca", dims = 1:20)

#DimPlot(rabbit_data3, reduction = "tsne",  label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/TSNE_rab3",".pdf",sep=""),width = 8, height = 8)

#DimPlot(rabbit_data3, reduction = "umap",  label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_rab3",".pdf",sep=""),width = 8, height = 8)

#DimPlot(rabbit_data3, reduction = "pca",  label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/PCA_split_rab3",".pdf",sep=""),width = 8, height = 8)


##################################################################
#Now do preliminary alignment
#cyno <- subset(cynomolgous_data,idents="ESC_primed")
#cyno <- subset(cynomolgous_data,idents="ESC_primed")
#rabbit <- subset(rabbit_data2,idents=c("Soma","PGC") )
#rabbit_data3 <- merge(rabbit_data, y = c(rabbit), project = "merged")


#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data2, rabbit_data, human_dataA1 ), dims = 1:20, anchor.features = 4000, k.filter = 20)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
DimPlot(mammal.combined, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_20",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_20",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_20",".pdf",sep=""),width = 16, height = 4)

mammal.combined2 <- subset(mammal.combined,idents = c("ESC_primed","PGCLC","PSC","PSC_low"),invert = TRUE)

DimPlot(mammal.combined2, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_20","_subset.pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined2, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_20","_subset.pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined2, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_20","_subset.pdf",sep=""),width = 16, height = 4)



#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data2, rabbit_data, human_dataA1 ), dims = 1:20, anchor.features = 4000, k.filter = 30)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
DimPlot(mammal.combined, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_30",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_30",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_30",".pdf",sep=""),width = 16, height = 4)

mammal.combined2 <- subset(mammal.combined,idents = c("ESC_primed","PGCLC","PSC","PSC_low"),invert = TRUE)

DimPlot(mammal.combined2, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_30","_subset.pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined2, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_30","_subset.pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined2, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_30","_subset.pdf",sep=""),width = 16, height = 4)




#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data2, rabbit_data, human_dataA1 ), dims = 1:20, anchor.features = 4000, k.filter = 40)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
DimPlot(mammal.combined, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_40",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_40",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_40",".pdf",sep=""),width = 16, height = 4)

mammal.combined2 <- subset(mammal.combined,idents = c("ESC_primed","PGCLC","PSC","PSC_low"),invert = TRUE)

DimPlot(mammal.combined2, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_40","_subset.pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined2, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_40","_subset.pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined2, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_40","_subset.pdf",sep=""),width = 16, height = 4)





#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data2, rabbit_data, human_dataA1 ), dims = 1:20, anchor.features = 4000, k.filter = 50)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
DimPlot(mammal.combined, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_50",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_50",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_50",".pdf",sep=""),width = 16, height = 4)

mammal.combined2 <- subset(mammal.combined,idents = c("ESC_primed","PGCLC","PSC","PSC_low"),invert = TRUE)

DimPlot(mammal.combined2, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_50","_subset.pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined2, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_50","_subset.pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined2, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_50","_subset.pdf",sep=""),width = 16, height = 4)






#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data2, rabbit_data, human_dataA1 ), dims = 1:20, anchor.features = 4000, k.filter = 25)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
DimPlot(mammal.combined, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_25",".pdf",sep=""),width = 15, height = 4)
DimPlot(mammal.combined, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_25",".pdf",sep=""),width = 15, height = 4)
DimPlot(mammal.combined, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_25",".pdf",sep=""),width = 15, height = 4)

mammal.combined2 <- subset(mammal.combined,idents = c("ESC_primed","PGCLC","PSC","PSC_low"),invert = TRUE)

DimPlot(mammal.combined2, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_25","_subset.pdf",sep=""),width = 15, height = 4)
DimPlot(mammal.combined2, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_25","_subset.pdf",sep=""),width = 15, height = 4)
DimPlot(mammal.combined2, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_25","_subset.pdf",sep=""),width = 15, height = 4)

#Do various plots
DimPlot(mammal.combined, reduction = "pca", group.by = "species", label = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA",".pdf",sep=""),width = 8, height = 4)
ElbowPlot(mammal.combined)
ggsave(filename=paste(saveext,"/DimRed/PCA_var",".pdf",sep=""),width = 8, height = 8)
VizDimLoadings(mammal.combined, dims = 1:6, nfeatures = 30, col = "blue",
               reduction = "pca", projected = FALSE, balanced = FALSE,
               ncol = NULL, combine = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_loadings",".pdf",sep=""),width = 8, height = 20)

panel1 <- c("SALL4",
             "ESRRB",
             "KLF5",
             "KLF4",
             "DPPA5",
             "TBX3",
             "POU5F1",
             "SOX2",
             "PRDM14",
             "NANOG",
             "SOX2",
             "SOX3",
             "SALL2",
             "OTX2")
panel2 <- c("T",
            "MIXL1",
            "EVX1",
            "CDH1",
            "FST",
            "CDH2",
            "GSC",
            "HAND1",
            "TBX6")
panel3 <- c(
"TFAP2C",
"TFAP2A",
"GATA3",
"HAND1",
"CD36")

panel4 <- c("APOA1",
            "APOB",
            "FOXA2",
            "SOX17",
            "PRDM1",
            "EOMES",
            "LHX1",
            "PITX2")

panel5 <- c("SOX17",
"PRDM1",
"TFAP2C",
"PRDM14",
"TFCP2L1",
"CD38",
"PDPN",
"NANOS3",
"UTF1",
"ARID5B",
"FGFR3",
"DND1",
"DPPA3",
"SOX15",
"TFAP2A",
"LAMA4",
"KIT",
"DPPA5",
"KLF4")

panel6 <- c("DAZL",
             "DDX4",
             "MAEL",
             "SYCP3",
             "RNF17",
             "TDRD5",
             "TDRD9",
             "TDRD10",
             "PIWIL1",
             "PIWIL2",
             "PIWIL3",
             "PRAME",
             "DND1",
             "DMRT1",
             "DNMT3B",
             "MAGEB2")

panel7 <- c("ZNF597",
            "ZNF99",
            "ZNF208",
            "ZNF354C",
            "ZNF534",
            "ZNF560",
            "ZNF662",
            "ZNF726",
            "ZNF728",
            "ZNF729")

panel8 <- c("TRIM28",
            "KAP1",
            "UHRF1",
            "DNMT3A",
            "DNMT3B",
            "DNMT3L",
            "TET1",
            "TET2",
            "TET3",
            "TDG")
panel9 <- c("SETDB1",
            "PRMT5",
            "KAP1",
            "KDM3A",
            "EHMT1",
            "EZH2",
            "KMDM4B",
            "EHMT2")


panel10 <- c("NANOG",
"FAM162",
"S100A10",
"LAMA4",
"DPPA5",
"IFI16",
"NANOS3",
"PCSK1N",
"RASD1",
"SDCBP",
"SAT1",
"FAM122C",
"PDPN",
"ASRGL1",
"HERC5",
"MKRN1",
"CDK2AP1",
"SOX15",
"IFITM1",
"RMND1")

panel11 <- c("HMGN2",
             "MDK",
             "SRP9",
             "GNAS",
             "PPIA",
             "BEX3",
             "TSTD1",
             "PTMA",
             "HMGB1",
             "HSBP1",
             "DEK",
             "FTL",
             "GSTO1",
             "PCLAF",
             "LDHB",
             "RANBP1",
             "UBE2V2",
             "SINHCAF",
             "PRDX2",
             "TMSB15A")

#testrabbit <- readRDS("/Users/christopherpenfold/Desktop/Thorsten/FINAL/Rabbit_Analysis_Fin5/mammal.combined.rds")
#testrabbit <- FindClusters(testrabbit, resolution = 1.2)


mammal.combined$ID1 <- Idents(mammal.combined)
Idents(mammal.combined,WhichCells(oldrabbit,idents="PGC") ) <- "PGC"
Idents(mammal.combined,WhichCells(oldrabbit,idents="EmDisc1") ) <- "EmDisc1"
Idents(mammal.combined,WhichCells(oldrabbit,idents="EmDisc2") ) <- "EmDisc2"
Idents(mammal.combined,WhichCells(oldrabbit,idents="EmDisc3") ) <- "EmDisc3"
Idents(mammal.combined,WhichCells(oldrabbit,idents="Epi") ) <- "Epi"
Idents(mammal.combined,WhichCells(oldrabbit,idents="Gast1") ) <- "Gast1"
Idents(mammal.combined,WhichCells(oldrabbit,idents="Gast2") ) <- "Gast2"
Idents(mammal.combined,WhichCells(oldrabbit,idents="ExMes") ) <- "ExMes"
Idents(mammal.combined,WhichCells(oldrabbit,idents="Hyp") ) <- "Hyp"
Idents(mammal.combined,WhichCells(oldrabbit,idents="ICM") ) <- "ICM"
Idents(mammal.combined,WhichCells(oldrabbit,idents="PTb") ) <- "PTb"
Idents(mammal.combined,WhichCells(oldrabbit,idents="SYS") ) <- "SYS"
Idents(mammal.combined,WhichCells(oldrabbit,idents="Tb") ) <- "Tb"
Idents(mammal.combined,WhichCells(oldrabbit,idents="VE") ) <- "VE"
Idents(mammal.combined,WhichCells(oldrabbit,idents="EmDiscEnd") ) <- "EmDiscEnd"
DimPlot(mammal.combined, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TestRab",".pdf",sep=""),width = 15, height = 4)
DimPlot(mammal.combined, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TestRabUMAP",".pdf",sep=""),width = 15, height = 4)


mammal.combined$OldIDs <- Idents(mammal.combined)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)
DimPlot(mammal.combined, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_Cl5",".pdf",sep=""),width = 15, height = 4)
DimPlot(mammal.combined, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_Cl5",".pdf",sep=""),width = 15, height = 4)
DimPlot(mammal.combined, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_Cl5",".pdf",sep=""),width = 15, height = 4)


mammal.combined <- FindClusters(mammal.combined, resolution = 1.2)
mammal.combined$Cl12 <- Idents(mammal.combined)
DimPlot(mammal.combined, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_Cl12",".pdf",sep=""),width = 15, height = 4)
DimPlot(mammal.combined, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_Cl12",".pdf",sep=""),width = 15, height = 4)
DimPlot(mammal.combined, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_Cl12",".pdf",sep=""),width = 15, height = 4)

Idents(mammal.combined) <- mammal.combined$ID1
mammal.combined2 <- subset(mammal.combined,idents = c("ESC_primed","PGCLC","PSC","PSC_low"),invert = TRUE)
Idents(mammal.combined2) <- mammal.combined2$Cl12

#Highlight by stage and by new e.g., 1,2,3

FinalIDs <- as.character(Idents(mammal.combined))

#Rabbit
FinalIDs[which(mammal.combined$Cl12=="5" & (mammal.combined$species2=="1) Rabbit" | mammal.combined$species2=="1) Rabbit_2") )]<-"rGast2"
FinalIDs[which(mammal.combined$Cl12=="13" & (mammal.combined$species2=="1) Rabbit" | mammal.combined$species2=="1) Rabbit_2") )]<-"rGast1"
FinalIDs[which(mammal.combined$Cl12=="9" & (mammal.combined$species2=="1) Rabbit" | mammal.combined$species2=="1) Rabbit_2") )]<-"rPGC"
FinalIDs[which(mammal.combined$Cl12=="6" & (mammal.combined$species2=="1) Rabbit" | mammal.combined$species2=="1) Rabbit_2") )]<-"rHyp"
FinalIDs[which(mammal.combined$Cl12=="7" & (mammal.combined$species2=="1) Rabbit" | mammal.combined$species2=="1) Rabbit_2") )]<-"rVE"
FinalIDs[which(mammal.combined$Cl12=="0" & (mammal.combined$species2=="1) Rabbit" | mammal.combined$species2=="1) Rabbit_2") )]<-"rEmDisc1"
FinalIDs[which(mammal.combined$Cl12=="1" & (mammal.combined$species2=="1) Rabbit" | mammal.combined$species2=="1) Rabbit_2") )]<-"rEmDisc2"
FinalIDs[which(mammal.combined$Cl12=="11" & (mammal.combined$species2=="1) Rabbit" | mammal.combined$species2=="1) Rabbit_2") )]<-"rEmDisc3"
FinalIDs[which(mammal.combined$Cl12=="8" & (mammal.combined$species2=="1) Rabbit" | mammal.combined$species2=="1) Rabbit_2") )]<-"rEmDisc4"
FinalIDs[which(mammal.combined$Cl12=="10" & (mammal.combined$species2=="1) Rabbit" | mammal.combined$species2=="1) Rabbit_2") )]<-"rEmDisc4"
FinalIDs[which(mammal.combined$Cl12=="3" & (mammal.combined$species2=="1) Rabbit" | mammal.combined$species2=="1) Rabbit_2") )]<-"rTb"
FinalIDs[which(mammal.combined$Cl12=="2" & (mammal.combined$species2=="1) Rabbit" | mammal.combined$species2=="1) Rabbit_2") )]<-"rPTb"

FinalIDs[which(mammal.combined$Cl12=="16" & (mammal.combined$species2=="1) Rabbit" | mammal.combined$species2=="1) Rabbit_2") )]<-"rEmDisc1"

FinalIDs[which(mammal.combined$ID1=="PGCLC" & (mammal.combined$species2=="1) Rabbit" | mammal.combined$species2=="1) Rabbit_2") )]<-"rPGCLC"
FinalIDs[which(mammal.combined$ID1=="PSC_low" & (mammal.combined$species2=="1) Rabbit" | mammal.combined$species2=="1) Rabbit_2") )]<-"rPSC_low"
FinalIDs[which(mammal.combined$ID1=="PSC" & (mammal.combined$species2=="1) Rabbit" | mammal.combined$species2=="1) Rabbit_2") )]<-"rPSC"


#Cyno
FinalIDs[which(mammal.combined$Cl12=="5" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cGast2"
FinalIDs[which(mammal.combined$Cl12=="13" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cGast1"
FinalIDs[which(mammal.combined$Cl12=="9" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cPGC"
FinalIDs[which(mammal.combined$Cl12=="6" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cHyp"
FinalIDs[which(mammal.combined$Cl12=="7" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cVE"
FinalIDs[which(mammal.combined$Cl12=="0" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cEmDisc1"
FinalIDs[which(mammal.combined$Cl12=="1" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cEmDisc2"
FinalIDs[which(mammal.combined$Cl12=="11" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cEmDisc3"
FinalIDs[which(mammal.combined$Cl12=="8" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cEmDisc4"
FinalIDs[which(mammal.combined$Cl12=="10" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cEmDisc4"
FinalIDs[which(mammal.combined$Cl12=="3" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cTb"
FinalIDs[which(mammal.combined$Cl12=="15" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cExMes"
FinalIDs[which(mammal.combined$Cl12=="4" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cPTb"
FinalIDs[which(mammal.combined$Cl12=="2" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cPTb"
FinalIDs[which(mammal.combined$Cl12=="14" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cPTb"


FinalIDs[which(mammal.combined$ID1=="Hyp_CS3" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cHyp"
FinalIDs[which(mammal.combined$ID1=="ICM_CS3" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cICM"
FinalIDs[which(mammal.combined$ID1=="Epi_CS3" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cEpi"
FinalIDs[which(mammal.combined$ID1=="Tb_CS3" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cTb"
FinalIDs[which(mammal.combined$ID1=="Tb_CS3" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cTb"
FinalIDs[which(mammal.combined$ID1=="ESC_primed" & (mammal.combined$species2=="2) Cynomolgous") )]<-"cPSC"



#Hyp_CS3 , ESC_primed, Hyp_CS3, ICM_CS3, Epi_CS3,  Tb_CS3

#Human
#
FinalIDs[which(mammal.combined$Cl12=="8" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hEmDisc4"

FinalIDs[which(mammal.combined$Cl12=="5" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hGast2"
FinalIDs[which(mammal.combined$Cl12=="6" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hHyp"
FinalIDs[which(mammal.combined$Cl12=="7" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hVE"
FinalIDs[which(mammal.combined$Cl12=="0" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hEmDisc1"
FinalIDs[which(mammal.combined$Cl12=="1" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hEmDisc2"
FinalIDs[which(mammal.combined$Cl12=="3" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hTb"
FinalIDs[which(mammal.combined$Cl12=="15" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hExMes"
FinalIDs[which(mammal.combined$Cl12=="4" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hPTb"
FinalIDs[which(mammal.combined$Cl12=="2" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hPTb"
FinalIDs[which(mammal.combined$Cl12=="14" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hPTb"
FinalIDs[which(mammal.combined$Cl12=="12" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hPTb"

FinalIDs[which(mammal.combined$ID1=="Tb_CS4" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hTb"
FinalIDs[which(mammal.combined$ID1=="VE_CS4" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hVE"
FinalIDs[which(mammal.combined$ID1=="Epi_CS3" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hEpi"
FinalIDs[which(mammal.combined$ID1=="Epi_CS4" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hEpi"
FinalIDs[which(mammal.combined$ID1=="ICM_CS3" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hICM"
FinalIDs[which(mammal.combined$ID1=="Hyp_CS3" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hHyp"
FinalIDs[which(mammal.combined$ID1=="Am_CS5" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hAm"
FinalIDs[which(mammal.combined$ID1=="Am_CS6" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hAm"
FinalIDs[which(mammal.combined$ID1=="Am_CS7" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hAm"

FinalIDs[which(mammal.combined$ID1=="Hyp_CS3" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hHyp"
FinalIDs[which(mammal.combined$ID1=="VE_CS3" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hVE"
FinalIDs[which(FinalIDs=="hVE" & mammal.combined$ID1=="EmDiscPS_CS6" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hEmDiscEnd"
FinalIDs[which(FinalIDs=="hExMes" & mammal.combined$ID1=="EmDiscPS_CS6" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hGast2"
FinalIDs[which(FinalIDs=="hExMes" & mammal.combined$ID1=="EmDisc_CS5" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hGast2"
FinalIDs[which(FinalIDs=="hExMes" & mammal.combined$ID1=="EmDisc_CS6" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hGast2"
FinalIDs[which(FinalIDs=="hExMes" & mammal.combined$ID1=="Tb_CS6" & (mammal.combined$species2=="4) Human (in vitro)") )]<-"hPTb"

cluster_letters <- LETTERS[mammal.combined$orig.ident]
cluster_letters[1:length(cluster_letters)] <- 0
cluster_letters[which(mammal.combined$ID2=="E6")] <- 0
cluster_letters[which(mammal.combined$ID2=="E7")] <- 1
cluster_letters[which(mammal.combined$ID1=="Soma")] <- 2
cluster_letters[which(mammal.combined$ID1=="PGC")] <- 2
cluster_letters[which(mammal.combined$ID1=="PGCLC")] <- 3
cluster_letters[which(mammal.combined$ID1=="PSC")] <- 3
cluster_letters[which(mammal.combined$ID1=="PSC_low")] <- 3
cluster_letters[which(mammal.combined$ID1=="ESC_primed")] <- 3
cluster_letters <- as.factor(cluster_letters)
names(cluster_letters)=rownames(mammal.combined@meta.data)
mammal.combined <- AddMetaData(object = mammal.combined,metadata = cluster_letters,col.name = 'cell.orig')


Idents(mammal.combined) <- FinalIDs
mammal.combined$FinalIDs <- FinalIDs

colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]


DimPlot(mammal.combined, cols = coluse, shape.by = 'cell.orig', reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_newLab",".pdf",sep=""),width = 15, height = 4)

DimPlot(mammal.combined, cols = coluse, shape.by = 'cell.orig', reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_newLab",".pdf",sep=""),width = 15, height = 4)

DimPlot(mammal.combined, cols = coluse, shape.by = 'cell.orig', reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_newLab",".pdf",sep=""),width = 15, height = 4)

Idents(mammal.combined) <- mammal.combined$ID1
mammal.combined2 <- subset(mammal.combined,idents = c("ESC_primed","PGCLC","PSC","PSC_low"),invert = TRUE)
Idents(mammal.combined2) <- mammal.combined2$FinalIDs


colind <- integer( length( levels(Idents(mammal.combined2)) )  )
for (i in 1:length( levels(Idents(mammal.combined2)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined2))[i])
}
coluse <- BaseCol[colind]


DimPlot(mammal.combined2, cols = coluse, shape.by = 'cell.orig', reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_subnewLab",".pdf",sep=""),width = 15, height = 4)

DimPlot(mammal.combined2, cols = coluse, shape.by = 'cell.orig', reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_subnewLab",".pdf",sep=""),width = 15, height = 4)

DimPlot(mammal.combined2, cols = coluse, shape.by = 'cell.orig', reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_suubnewLab",".pdf",sep=""),width = 15, height = 4)



Idents(mammal.combined) <- mammal.combined$ID1
mammal.combined4 <- subset(mammal.combined,idents = c("ESC_primed","PGCLC","PSC","PSC_low","Soma","PGC"),invert = TRUE)
Idents(mammal.combined4) <- mammal.combined4$FinalIDs


colind <- integer( length( levels(Idents(mammal.combined4)) )  )
for (i in 1:length( levels(Idents(mammal.combined4)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined4))[i])
}
coluse <- BaseCol[colind]


DimPlot(mammal.combined4, cols = coluse, shape.by = 'cell.orig', reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_subnewLab_orig",".pdf",sep=""),width = 15, height = 4)

DimPlot(mammal.combined4, cols = coluse, shape.by = 'cell.orig', reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_subnewLab_orig",".pdf",sep=""),width = 15, height = 4)

DimPlot(mammal.combined4, cols = coluse, shape.by = 'cell.orig', reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_suubnewLab_orig",".pdf",sep=""),width = 15, height = 4)

#mammal.combined4 == original


Idents(mammal.combined) <- mammal.combined$FinalIDs
mammal.combined3 <- subset(mammal.combined,idents = c("cEpi","cEmDisc1","cEmDisc3","cGast1","cGast2","cPGC","cPSC",
  "hEmDisc1","hEmDisc2","hEmDisc4","hGast2", "hEpi",
  "rEmDisc1","rEmDisc2","rEmDisc3","rEmDisc4","rGast1","rGast2",
  "rPGC","rPGCLC","rPSC","rPSC_low"),invert = FALSE)


mammal.combined3 <- RunPCA(mammal.combined3, npcs = 30, verbose = FALSE)
mammal.combined3 <- RunUMAP(mammal.combined3, reduction = "pca", dims = 1:20)
mammal.combined3 <- FindNeighbors(mammal.combined3, reduction = "pca", dims = 1:20)


colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- BaseCol[colind]

DimPlot(mammal.combined3, cols = coluse,shape.by = 'cell.orig', reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_EmbsubnewLab",".pdf",sep=""),width = 15, height = 4)

DimPlot(mammal.combined3, cols = coluse, shape.by = 'cell.orig', reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_EmbsubnewLab",".pdf",sep=""),width = 15, height = 4)

DimPlot(mammal.combined3, cols = coluse, shape.by = 'cell.orig', reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_EmbsuubnewLab",".pdf",sep=""),width = 15, height = 4)


mammal.combined3 <- subset(mammal.combined,idents = c("cEmDisc1","cEmDisc3","cGast1","cGast2","cPGC","cPSC",
                                                      "hEmDisc1","hEmDisc2","hEmDisc4","hGast2",
                                                      "rEmDisc1","rEmDisc2","rEmDisc3","rGast1","rGast2",
                                                      "rPGC","rPGCLC","rPSC","rPSC_low"),invert = FALSE)



colind <- integer( length( levels(Idents(mammal.combined3)) )  )
for (i in 1:length( levels(Idents(mammal.combined3)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined3))[i])
}
coluse <- BaseCol[colind]

mammal.combined3 <- RunPCA(mammal.combined3, npcs = 30, verbose = FALSE)
mammal.combined3 <- RunUMAP(mammal.combined3, reduction = "pca", dims = 1:20)
mammal.combined3 <- FindNeighbors(mammal.combined3, reduction = "pca", dims = 1:20)


DimPlot(mammal.combined3, cols = coluse, shape.by = 'cell.orig', reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_EmbsubnewLab2",".pdf",sep=""),width = 15, height = 4)

DimPlot(mammal.combined3, cols = coluse, shape.by = 'cell.orig', reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_EmbsubnewLab2",".pdf",sep=""),width = 15, height = 4)

DimPlot(mammal.combined3, cols = coluse, shape.by = 'cell.orig', reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_EmbsuubnewLab2",".pdf",sep=""),width = 15, height = 4)

###################################################
#Now we do UMAP feature plots
Data <- GetAssayData(mammal.combined,assay = "RNA")
Data2 <- GetAssayData(mammal.combined3,assay = "RNA")
Data3 <- GetAssayData(mammal.combined4,assay = "RNA")


DefaultAssay(mammal.combined) <- "RNA"
FeaturePlot(mammal.combined, reduction = "umap", features = "BMP2", split.by = "species", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_BMP2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "BMP4", split.by = "species", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_BMP4.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "BMP8", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_BMP8.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "WNT3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_WNT3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "WNT5A", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_WNT5A.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "WNT5B", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_WNT5B.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "FGF2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_FGF2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "FGF4", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_FGF4.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "FGF5", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_FGF5.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "FGF10", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_FGF10.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "NODAL", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_NODAL.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "TGFB1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TGFB1.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "TGFB2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TGFB2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "TGFB3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TGFB3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)


FeaturePlot(mammal.combined, reduction = "umap", features = "POU5F1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_POU5F1.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "NANOG", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_NANOG.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "SOX2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_SOX2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "PRDM14", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_PRDM14.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "T", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_T.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "MIXL1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_MIXL1.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "TFAP2C", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TFAP2C.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "NANOS3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_NANOS3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "PRDM1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_PRDM1.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "SOX17", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_SOX17.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "GATA3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_GATA3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "TFAP2A", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TFAP2A.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "APOA1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_APOA1.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "LHX1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_LHX1.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)


FeaturePlot(mammal.combined, reduction = "pca", features = "BMP2", split.by = "species", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_BMP2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "BMP4", split.by = "species", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_BMP4.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "BMP8", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_BMP8.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "WNT3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_WNT3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "WNT5A", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_WNT5A.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "WNT5B", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_WNT5B.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "FGF2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_FGF2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "FGF4", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_FGF4.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "FGF5", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_FGF5.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "FGF10", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_FGF10.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "NODAL", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_NODAL.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "TGFB1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_TGFB1.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "TGFB2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_TGFB2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "TGFB3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_TGFB3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)


FeaturePlot(mammal.combined, reduction = "pca", features = "POU5F1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_POU5F1.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "NANOG", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_NANOG.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "SOX2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_SOX2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "PRDM14", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_PRDM14.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "T", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_T.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "MIXL1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_MIXL1.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "TFAP2C", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_TFAP2C.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "NANOS3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_NANOS3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "PRDM1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_PRDM1.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "SOX17", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_SOX17.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "GATA3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_GATA3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "TFAP2A", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_TFAP2A.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "APOA1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_APOA1.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "pca", features = "LHX1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_LHX1.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)



DefaultAssay(mammal.combined2) <- "RNA"
FeaturePlot(mammal.combined2, reduction = "umap", features = "BMP2", split.by = "species", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_BMP2_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "BMP4", split.by = "species", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_BMP4_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "BMP8", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_BMP8_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "WNT3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_WNT3_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "WNT5A", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_WNT5A_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "WNT5B", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_WNT5B_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "FGF2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_FGF2_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "FGF4", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_FGF4_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "FGF5", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_FGF5_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "FGF10", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_FGF10_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "NODAL", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_NODAL_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "TGFB1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TGFB1_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "TGFB2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TGFB2_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "TGFB3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TGFB3_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)


FeaturePlot(mammal.combined2, reduction = "umap", features = "POU5F1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_POU5F1_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "NANOG", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_NANOG_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "SOX2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_SOX2_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "PRDM14", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_PRDM14_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "T", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_T_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "MIXL1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_MIXL1_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "TFAP2C", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TFAP2C_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "NANOS3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_NANOS3_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "PRDM1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_PRDM1_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "SOX17", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_SOX17_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "GATA3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_GATA3_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "TFAP2A", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TFAP2A_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "APOA1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_APOA1_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "umap", features = "LHX1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_LHX1_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)



FeaturePlot(mammal.combined2, reduction = "pca", features = "BMP2", split.by = "species", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_BMP2_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "BMP4", split.by = "species", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_BMP4_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "BMP8", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_BMP8_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "WNT3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_WNT3_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "WNT5A", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_WNT5A_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "WNT5B", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_WNT5B_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "FGF2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_FGF2_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "FGF4", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_FGF4_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "FGF5", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_FGF5_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "FGF10", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_FGF10_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "NODAL", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_NODAL_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "TGFB1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_TGFB1_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "TGFB2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_TGFB2_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "TGFB3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_TGFB3_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)


FeaturePlot(mammal.combined2, reduction = "pca", features = "POU5F1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_POU5F1_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "NANOG", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_NANOG_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "SOX2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_SOX2_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "PRDM14", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_PRDM14_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "T", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_T_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "MIXL1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_MIXL1_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "TFAP2C", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_TFAP2C_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "NANOS3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_NANOS3_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "PRDM1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_PRDM1_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "SOX17", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_SOX17_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "GATA3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_GATA3_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "TFAP2A", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_TFAP2A_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "APOA1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_APOA1_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined2, reduction = "pca", features = "LHX1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_LHX1_2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)



DefaultAssay(mammal.combined3) <- "RNA"
FeaturePlot(mammal.combined3, reduction = "umap", features = "BMP2", split.by = "species", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_BMP2_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "BMP4", split.by = "species", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_BMP4_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "BMP8", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_BMP8_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "WNT3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_WNT3_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "WNT5A", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_WNT5A_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "WNT5B", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_WNT5B_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "FGF2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_FGF2_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "FGF4", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_FGF4_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "FGF5", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_FGF5_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "FGF10", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_FGF10_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "NODAL", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_NODAL_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "TGFB1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TGFB1_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "TGFB2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TGFB2_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "TGFB3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TGFB3_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)


FeaturePlot(mammal.combined3, reduction = "umap", features = "POU5F1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_POU5F1_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "NANOG", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_NANOG_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "SOX2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_SOX2_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "PRDM14", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_PRDM14_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "T", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_T_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "MIXL1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_MIXL1_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "TFAP2C", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TFAP2C_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "NANOS3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_NANOS3_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "PRDM1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_PRDM1_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "SOX17", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_SOX17_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "GATA3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_GATA3_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "TFAP2A", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TFAP2A_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "APOA1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_APOA1_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "umap", features = "LHX1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_LHX1_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)



FeaturePlot(mammal.combined3, reduction = "pca", features = "BMP2", split.by = "species", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_BMP2_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "BMP4", split.by = "species", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_BMP4_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "BMP8", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_BMP8_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "WNT3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_WNT3_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "WNT5A", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_WNT5A_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "WNT5B", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_WNT5B_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "FGF2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_FGF2_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "FGF4", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_FGF4_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "FGF5", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_FGF5_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "FGF10", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_FGF10_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "NODAL", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_NODAL_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "TGFB1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_TGFB1_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "TGFB2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_TGFB2_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "TGFB3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_TGFB3_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)



FeaturePlot(mammal.combined3, reduction = "pca", features = "POU5F1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_POU5F1_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "NANOG", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_NANOG_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "SOX2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_SOX2_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "PRDM14", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_PRDM14_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "T", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_T_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "MIXL1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_MIXL1_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "TFAP2C", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_TFAP2C_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "NANOS3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_NANOS3_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "PRDM1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_PRDM1_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "SOX17", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_SOX17_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "GATA3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_GATA3_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "TFAP2A", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_TFAP2A_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "APOA1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_APOA1_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined3, reduction = "pca", features = "LHX1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_pca_LHX1_3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)



#Now we will start looking at markers

finalList <- c("SALL4",
              "ESRRB",
              "KLF4",
              "KLF5",
              "TBX3",
              "POU5F1",
              "NANOG",
              "SOX2",
              "SOX3",
              "PRDM14",
              "TFCP2L1",
              "T",
              "EVX1",
              "GSC",
              "HAND1",
              "TBX6",
              "MIXL1",
              "FST",
              "CDH2",
              "PPIA",
              "HMGB1",
              "TMSB15A",
              "NANOS3",
              "SOX17",
              "SOX15",
              "TFAP2C",
              "PRDM1",
              "KIT",
              "UTF1",
              "ARID5B",
              "LAMA4",
              "MKRN1",
              "DND1",
              "MAEL",
              "TDRD9",
              "DDX4",
              "TDRD5",
              "PIWIL1",
              "PIWIL2",
              "PIWIL4",
              "PRAME",
              "DMRT1",
              "MAGEB2",
              "MAGEB2",
              "TFAP2C",
              "TFAP2A",
              "GATA3",
              "CD36",
              "HAND1",
              "EOMES",
              "APOA1",
              "APOB",
              "FOXA1",
              "FOXA2",
              "LHX1",
              "PITX2",
              "ZNF597",
              "DNMT3A",
              "TET2",
              "TET3",
              "DNMT3B",
              "EHMT2",
              "PRMT5",
              "EHMT1",
              "EZH2",
              "KDM3A")


library("pheatmap")
avexp1  <- AverageExpression(object = mammal.combined, slot = "data")
avexp2  <- AverageExpression(object = mammal.combined2, slot = "data")
avexp3  <- AverageExpression(object = mammal.combined4, slot = "data")

a <- avexp1$RNA
a <- a[finalList ,]

orders <- c("rEmDisc4","rEmDisc2","rEmDisc1","rGast1","rGast2","rPGC","rTb","rHyp","rVE")
a <- a[,orders]

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(-2, 2, length.out =40)
aa <- na.omit( log2(a+1) )
pheatmap(aa,color =  redblue1(40),breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, scale = "row",  filename = paste(saveext, "Rabbit_HM_logged.pdf"),width=4,height=16)

mat_breaks <- seq(0, 100, length.out =40)
pheatmap( (exp(aa)-1)*100,color =  redblue1(40),breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext, "Rabbit_HM_TPM.pdf"),width=4,height=16)

#avexp1  <- AverageExpression(object = mammal.combined, slot = "data")
#avexp2  <- AverageExpression(object = mammal.combined2, slot = "data")

a <- avexp1$RNA
a <- a[finalList ,]

orders <- c("rPGC","rPGCLC","rPSC","rPSC_low")
a <- a[,orders]

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(-2, 2, length.out =40)
aa <- na.omit( log2(a+1) )
pheatmap(aa,color =  redblue1(40),breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, scale = "row",  filename = paste(saveext, "Rabbit_HM_logged_3.pdf"),width=4,height=16)

mat_breaks <- seq(0, 100, length.out =40)
pheatmap( (exp(aa)-1)*100,color =  redblue1(40),breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext, "Rabbit_HM_TPM_3.pdf"),width=4,height=16)



a <- avexp1$RNA
a <- a[finalList ,]

orders <- c("rEmDisc4","rEmDisc2","rEmDisc1","rGast1","rGast2","rPGC","rPSC","rPSC_low","rPGCLC")
a <- a[,orders]

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(-2, 2, length.out =40)
aa <- na.omit( log2(a+1) )
pheatmap(aa,color =  redblue1(40),breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, scale = "row",  filename = paste(saveext, "Rabbit_HM_logged_2.pdf"),width=4,height=16)

mat_breaks <- seq(0, 100, length.out =40)
pheatmap( (exp(aa)-1)*100,color =  redblue1(40),breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext, "Rabbit_HM_TPM_2.pdf"),width=4,height=16)



l1<-grep("rEmDisc4", Idents(mammal.combined))
l2<-grep("rEmDisc2", Idents(mammal.combined))
l4<-grep("rEmDisc1", Idents(mammal.combined))
l5<-grep("rGast1", Idents(mammal.combined))
l6<-grep("rGast2", Idents(mammal.combined))
l7<-grep("rPGC", Idents(mammal.combined))
l10<-grep("rHyp", Idents(mammal.combined))
l11<-grep("rVE", Idents(mammal.combined))
l8<-grep("rTb", Idents(mammal.combined))
#l9<-grep("rPTb", Idents(mammal.combined))


lfull <- c(l1,l2,l4,l5,l6,l7,l10,l11,l8)

DN1 <- GetAssayData(mammal.combined, slot = "data")

mat_breaks <- seq(0, 2, length.out = 20)
pheatmap(DN1[finalList,lfull],color =  redblue1(20),breaks = mat_breaks,gaps_col=c(4,101,145,173,230,261,289,361,386), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/HMrf_ara_logged",".pdf",sep=""),width=20,height=6)

mat_breaks <- seq(-1, 1, length.out = 20)
pheatmap(DN1[finalList,lfull],color =  redblue1(20),breaks = mat_breaks,gaps_col=c(4,101,145,173,230,261,289,361,386), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, scale="row", filename = paste(saveext,"/DimRed/HMrf_ara_logged_scaled",".pdf",sep=""),width=20,height=6)

DN2 <- DN1[finalList,lfull]
annotation_col = data.frame(Lineage = Idents(mammal.combined)[lfull])
rownames(annotation_col) <- colnames(DN2)
mat_breaks <- seq(-1, 1, length.out = 20)
pheatmap(DN1[finalList,lfull],color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,annotation_col = annotation_col,cluster_cols=TRUE, scale="row", filename = paste(saveext,"/DimRed/HMrf_ara_logged_scaled_header",".pdf",sep=""),width=20,height=6)
l2<-grep("rPSC", Idents(mammal.combined))
l3<-grep("rEmDisc4", Idents(mammal.combined))
l4<-grep("rEmDisc2", Idents(mammal.combined))
l5<-grep("rEmDisc1", Idents(mammal.combined))
l6<-grep("rGast1", Idents(mammal.combined))
l7<-grep("rGast2", Idents(mammal.combined))
l81<-grep("rPGC", Idents(mammal.combined))
l9<-grep("rPGCLC", Idents(mammal.combined))
l8<-setdiff(l81,l9)
#l10<-grep("rPGC1", Idents(mammal.combined))
#l12<-grep("rPGC3", Idents(mammal.combined))
##l11<-grep("rPGC2", Idents(mammal.combined))
lfull <- c(l2,l3,l4,l5,l6,l7,l8,l9)
#DN1 <- GetAssayData(mammal.combined, slot = "data")
mat_breaks <- seq(0, 2, length.out = 20)
pheatmap(DN1[finalList,lfull],color =  redblue1(20),breaks = mat_breaks,gaps_col=c(31,35,132,176,204,261,276), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/HMrf_ara_logged_pgclc",".pdf",sep=""),width=20,height=6)
mat_breaks <- seq(-1, 1, length.out = 20)
pheatmap(DN1[finalList,lfull],color =  redblue1(20),breaks = mat_breaks,gaps_col=c(16,20,117,161,189,246,261), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, scale="row", filename = paste(saveext,"/DimRed/HMrf_ara_logged_scaled_pgclc",".pdf",sep=""),width=20,height=6)

DN2 <- DN1[finalList,lfull]
annotation_col = data.frame(Lineage = Idents(mammal.combined)[lfull])
rownames(annotation_col) <- colnames(DN2)

mat_breaks <- seq(-1, 1, length.out = 20)
pheatmap(DN1[finalList,lfull],color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,annotation_col = annotation_col,cluster_cols=TRUE, scale="row", filename = paste(saveext,"/DimRed/HMrf_ara_logged_scaled_pgclc_header",".pdf",sep=""),width=20,height=6)


#Redo heatmaps in Figure 1 ...
a <- avexp3$RNA
a <- a[Markers ,]

orders <- c("rEmDisc4","rEmDisc2","rEmDisc1","rGast1","rGast2","rTb","rPTb","rHyp","rVE")
a <- a[,orders]

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(-2, 2, length.out =40)
aa <- na.omit( log2(a+1) )
pheatmap(aa,color =  redblue1(40),breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, scale = "row",  filename = paste(saveext, "Rabbit_HM_origlogged_2.pdf"),width=4,height=16)

mat_breaks <- seq(0, 100, length.out =40)
pheatmap( (exp(aa)-1)*100,color =  redblue1(40),breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext, "Rabbit_HM_origTPM_2.pdf"),width=4,height=16)


l1<-grep("rEmDisc4", Idents(mammal.combined4))
l2<-grep("rEmDisc2", Idents(mammal.combined4))
l4<-grep("rEmDisc1", Idents(mammal.combined4))
l5<-grep("rGast1", Idents(mammal.combined4))
l6<-grep("rGast2", Idents(mammal.combined4))
l8<-grep("rTb", Idents(mammal.combined4))
l9<-grep("rPTb", Idents(mammal.combined4))
l10<-grep("rHyp", Idents(mammal.combined4))
l11<-grep("rVE", Idents(mammal.combined4))


lfull <- c(l1,l2,l4,l5,l6,l8,l9,l10,l11)
DefaultAssay(mammal.combined4) <- "RNA"
DN1 <- GetAssayData(mammal.combined4, slot = "data")

mat_breaks <- seq(0, 2, length.out = 20)
pheatmap(DN1[Markers,lfull],color =  redblue1(20),breaks = mat_breaks,gaps_col=c(4,101,134,153,199,224,226,254), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/HMrf_ara_logged",".pdf",sep=""),width=20,height=6)

mat_breaks <- seq(-1, 1, length.out = 20)
pheatmap(DN1[Markers,lfull],color =  redblue1(20),breaks = mat_breaks,gaps_col=c(4,101,134,153,199,224,226,254), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, scale="row", filename = paste(saveext,"/DimRed/HMrf_ara_logged_scaled",".pdf",sep=""),width=20,height=6)

DN2 <- DN1[Markers,lfull]
annotation_col = data.frame(Lineage = Idents(mammal.combined)[lfull])
rownames(annotation_col) <- colnames(DN2)
mat_breaks <- seq(-1, 1, length.out = 20)
pheatmap(DN1[Markers,lfull],color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,annotation_col = annotation_col,cluster_cols=TRUE, scale="row", filename = paste(saveext,"/DimRed/HMrf_ara_logged_scaled_header",".pdf",sep=""),width=20,height=6)


#Cross covariance plots
a <- avexp1$RNA
orders1 <- c("rEmDisc4","rEmDisc2","rEmDisc1","rGast1","rGast2","rPGC")
orders2 <- c("rPSC","rPSC_low","rPGCLC")

a1 <- a[,orders1]
a2 <- a[,orders2]

C1 <- cor(log(a1 +1), log(a2+1), method = "pearson")

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.8, 1, length.out = 20)
pheatmap(C1,color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/HM_ara_logged","_cross1.pdf",sep=""),width=10,height=16)


a <- avexp1$integrated

orders1 <- c("rEmDisc4","rEmDisc2","rEmDisc1","rGast1","rGast2","rPGC")
orders2 <- c("cEmDisc3","cEmDisc1","cGast1","cGast2","cPGC")

a1 <- a[,orders1]
a2 <- a[,orders2]


C1 <- cor(log(a1 +1), log(a2+1), method = "pearson")

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.8, 1, length.out = 20)
pheatmap(C1,color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/HM_ara_logged_cross2",".pdf",sep=""),width=10,height=16)


#Old cross covariances
a <- avexp3$integrated
orders1 <- c("rEmDisc4","rEmDisc2","rEmDisc1","rGast1","rGast2","rTb","rPTb","rVE")
orders2 <- c("cICM","cEpi","cEmDisc3","cEmDisc1","cGast1","cGast2","cTb","cPTb","cHyp","cVE")
orders3 <- c("hICM","hEpi","hEmDisc4","hEmDisc2","hEmDisc1","hGast2","hTb","hPTb","hHyp","hVE")
a1 <- a[,orders1]
a2 <- a[,orders2]

C1 <- cor(log(a1 +1), log(a2+1), method = "pearson")

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.8, 1, length.out = 20)
pheatmap(C1,color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/HM_ara_logged_old_cross1",".pdf",sep=""),width=10,height=16)

a3 <- a[,orders3]
C1 <- cor(log(a1 +1), log(a3+1), method = "pearson")
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.8, 1, length.out = 20)
pheatmap(C1,color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/HM_ara_logged_old_cross2",".pdf",sep=""),width=10,height=16)


#Full cross covariances
a <- avexp1$integrated
orders1 <- c("rEmDisc4","rEmDisc2","rEmDisc1","rGast1","rGast2","rPGC","rTb","rPTb","rVE")
orders2 <- c("cICM","cEpi","cEmDisc3","cEmDisc1","cGast1","cGast2","cPGC","cTb","cPTb","cHyp","cVE")
orders3 <- c("hICM","hEpi","hEmDisc4","hEmDisc2","hEmDisc1","hGast2","hTb","hPTb","hHyp","hVE")
a1 <- a[,orders1]
a2 <- a[,orders2]

C1 <- cor(log(a1 +1), log(a2+1), method = "pearson")

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.8, 1, length.out = 20)
pheatmap(C1,color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/HM_ara_logged_old_fullcross1",".pdf",sep=""),width=10,height=16)

a3 <- a[,orders3]
C1 <- cor(log(a1 +1), log(a3+1), method = "pearson")
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.8, 1, length.out = 20)
pheatmap(C1,color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/HM_ara_logged_old_fullcross2",".pdf",sep=""),width=10,height=16)




#
DefaultAssay(mammal.combined)
Cl1 <- FindMarkers(mammal.combined, ident.2 = c("rPGC"), ident.1 = c("rPGCLC"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/PGC_PGCLC_DEseq2.csv",sep=""))

SIGNAL<-read.table("/Users/christopherpenfold/Desktop/LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]
TF<-read.table("/Users/christopherpenfold/Desktop/Toshiaki/UberA/TF.txt",header = F)
TF <- TF$V1

Ae5 <- avexp1$RNA
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.05))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.2)) ,]

p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.05),linetype="dashed") + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"PGCLC_PGC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.05),linetype="dashed") + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"PGCLC_PGC_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()




Cl1 <- FindMarkers(mammal.combined, ident.2 = c("rPSC","rPSC_low"), ident.1 = c("rPGCLC"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
write.csv(as.data.frame(Cl1), file=paste(saveext,"/PGC_PGCLC_DEseq2.csv",sep=""))

SIGNAL<-read.table("/Users/christopherpenfold/Desktop/LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]
TF<-read.table("/Users/christopherpenfold/Desktop/Toshiaki/UberA/TF.txt",header = F)
TF <- TF$V1

Ae5 <- avexp1$RNA
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.05))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.2)) ,]

p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.05),linetype="dashed") + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"PGCLC_ESC_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.05),linetype="dashed") + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"PGCLC_ESC_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()


Cl1 <- FindMarkers(mammal.combined, ident.1 = c("rPSC","rPSC_low"), ident.2 = c("rEmDisc1","rEmDisc2"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
Ae5 <- avexp1$RNA
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.05))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.2)) ,]

p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.05),linetype="dashed") + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"PSC_EmDisc_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.05),linetype="dashed") + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"PSC_EmDisc_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()





Cl1 <- FindMarkers(mammal.combined, ident.1 = c("rPGC"), ident.2 = c("rEmDisc1","rEmDisc2"), verbose = FALSE, test.use = "DESeq2", only.pos = FALSE, logfc.threshold = log(0))
SIGNAL<-read.table("/Users/christopherpenfold/Desktop/LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]
TF<-read.table("/Users/christopherpenfold/Desktop/Toshiaki/UberA/TF.txt",header = F)
TF <- TF$V1

Ae5 <- avexp1$RNA
Ae5$Pval1 <- 1
Ae5$FC1 <- NA
Ae5$Indi <- 0
Ae5[rownames(Cl1),"Pval1"] <- -log(Cl1$p_val_adj)
Ae5[rownames(Cl1),"FC1"] <- Cl1$avg_logFC
pospos1 <- which( abs(Ae5$FC1)>log(1.2) & Ae5$Pval1>-log(0.05))
#Any DE?
Ae5$Indi[pospos1] <- 1
Ae5$Indi <- as.factor(Ae5$Indi)
genes.to.label1 = c(intersect(rownames(Ae5)[pospos1],TF),intersect(rownames(Ae5)[pospos1],TF))
genes.to.label2 = c(intersect(rownames(Ae5)[pospos1],SIGNAL1),intersect(rownames(Ae5)[pospos1],SIGNAL1))
genes.to.label3 = c(intersect(rownames(Ae5)[pospos1],SIGNAL2),intersect(rownames(Ae5)[pospos1],SIGNAL2))
genes.to.label4 = c(intersect(rownames(Ae5)[pospos1],SIGNAL3),intersect(rownames(Ae5)[pospos1],SIGNAL3))
genes.to.label5 <- unique(c(genes.to.label2,genes.to.label3,genes.to.label4))
Ae5 <- Ae5[rownames(Cl1),]
Ae5 <- Ae5[order(Ae5$Pval1),]
Ae5 <- Ae5[which(-Ae5$Pval1<log(0.2)) ,]

p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.05),linetype="dashed") + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label1, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"PGC_EmDisc_TF.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()
p1 <- ggplot(Ae5, aes(FC1,Pval1)) + geom_point(aes(color=Indi)) + theme_classic() +  scale_color_manual(values=c('lightgrey','black')) + geom_hline(yintercept=-log(0.05),linetype="dashed") + geom_vline(xintercept = log(1.2), linetype="dashed") + geom_vline(xintercept = -log(1.2), linetype="dashed")  #+ geom_vline(aes(xintercept = 0)) + geom_hline(aes(yintercept = 0)) + theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
p1 <- LabelPoints(plot = p1, points = genes.to.label5, repel = TRUE, color = 'red')
p1 <- p1 + labs(x = "Delta Expression", y = "-log Pval")
ggsave(filename=paste(saveext,"PGC_EmDisc_lig1.pdf",sep=""),width = 13, height = 13, plot = p1)
dev.off()







mD2 <- read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/Toshi_corrections/Old/Mouse/GSE121650_rna_counts.tsv",sep="\t", header = T, row.names=1)
mName <- read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/Toshi_corrections/Old/Mouse/Genenames.csv",sep=",", header = T, row.names=1)
mMeta <- read.table("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/Toshi_corrections/Old/Mouse/Meta.csv",sep=",", header = T, row.names=1)
rownames(mD) <- mName$ID
rownames(mD)[which(rownames(mD)=="TBXT")] <- "T"
mouse_data2 <- CreateSeuratObject(counts = mD, assay = "RNA",min.cells = 0, min.features = 0)
mouse_data2$species <- "4) Mouse"
mouse_data2$stage <- mMeta$ID
Idents(mouse_data2) <- mouse_data2$stage
mouse_data2 <- subset(mouse_data2,idents = c("E6.5","E7.5"))
mouse_data2 <- NormalizeData(mouse_data2, verbose = FALSE)
mouse_data2 <- FindVariableFeatures(mouse_data2, selection.method = "vst", nfeatures = 3000)
mouse_data2 <- ScaleData(mouse_data2, verbose = FALSE)
mouse_data2 <- RunPCA(mouse_data2, npcs = 30, verbose = FALSE)
mouse_data2 <- RunUMAP(mouse_data2, reduction = "pca", dims = 1:20)
mouse_data2 <- RunTSNE(mouse_data2, reduction = "pca", dims = 1:20)
mouse_data2 <- FindNeighbors(mouse_data2, reduction = "pca", dims = 1:20)
mouse_data2 <- FindClusters(mouse_data2, resolution = 0.5)

DimPlot(mouse_data2, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/Mouse_pca",".pdf",sep=""),width = 8, height = 8)
DimPlot(mouse_data2, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/Mouse_tsne",".pdf",sep=""),width = 8, height = 8)
DimPlot(mouse_data2, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/Mouse_umap",".pdf",sep=""),width = 8, height = 8)
mouse_data2$Cl <- Idents(mouse_data2)
Idents(mouse_data2) <- mouse_data2$stage
DimPlot(mouse_data2, reduction = "pca",  label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/Mouse_pca_stage",".pdf",sep=""),width = 8, height = 8)
DimPlot(mouse_data2, reduction = "tsne",label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/Mouse_tsne_stage",".pdf",sep=""),width = 8, height = 8)
DimPlot(mouse_data2, reduction = "umap",  label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/Mouse_umap_stage",".pdf",sep=""),width = 8, height = 8)
FeaturePlot(mouse_data2, reduction = "umap", features = "PDGFRA",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Mouse_umap_PDGFRA.pdf",sep=""),width = 8, height = 8,limitsize = FALSE)
FeaturePlot(mouse_data2, reduction = "umap", features = "GATA4",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Mouse_umap_GATA4.pdf",sep=""),width = 8, height = 8,limitsize = FALSE)
FeaturePlot(mouse_data2, reduction = "umap", features = "LEFTY2",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Mouse_umap_LEFTY2.pdf",sep=""),width = 8, height = 8,limitsize = FALSE)
FeaturePlot(mouse_data2, reduction = "umap", features = "T",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Mouse_umap_T.pdf",sep=""),width = 8, height = 8,limitsize = FALSE)
FeaturePlot(mouse_data2, reduction = "umap", features = "POU5F1",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Mouse_umap_POU5F1.pdf",sep=""),width = 8, height = 8,limitsize = FALSE)
FeaturePlot(mouse_data2, reduction = "umap", features = "SOX2",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Mouse_umap_SOX2.pdf",sep=""),width = 8, height = 8,limitsize = FALSE)
FeaturePlot(mouse_data2, reduction = "umap", features = "TFAP2C",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Mouse_umap_TFAP2C.pdf",sep=""),width = 8, height = 8,limitsize = FALSE)
FeaturePlot(mouse_data2, reduction = "umap", features = "NANOS3",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Mouse_umap_NANOS3.pdf",sep=""),width = 8, height = 8,limitsize = FALSE)
FeaturePlot(mouse_data2, reduction = "umap", features = "PRDM1", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Mouse_umap_PRDM1.pdf",sep=""),width = 8, height = 8,limitsize = FALSE)
FeaturePlot(mouse_data2, reduction = "umap", features = "NANOG", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Mouse_umap_NANOG.pdf",sep=""),width = 8, height = 8,limitsize = FALSE)
FeaturePlot(mouse_data2, reduction = "umap", features = "PRDM14", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Mouse_umap_PRDM14.pdf",sep=""),width = 8, height = 8,limitsize = FALSE)

Idents(mouse_data2) <- mouse_data2$Cl
MD <- merge(mouse_data2, y = c(mammal.combined), project = "merged")

avexp1  <- AverageExpression(object = MD, slot = "data")


a <- avexp1$RNA

orders1 <- c("2","0","3","1")
orders2 <- c("rEmDisc4","rEmDisc2","rEmDisc1","rGast1","rGast2","rPGC","rHyp","rVE","rTb","rPGCLC","rPSC")

a1 <- a[intersect(rownames(mouse_data2),rownames(rabbit_data2)),orders1]
a2 <- a[intersect(rownames(mouse_data2),rownames(rabbit_data2)),orders2]

C1 <- cor(log(a1 +1), log(a2+1), method = "pearson")

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
#mat_breaks <- seq(0.8, 1, length.out = 20)
pheatmap(C1,color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Mouse_cross2",".pdf",sep=""),width=10,height=5)



orders1 <- c("2","0","3","1")
orders2 <- c("rEmDisc1","rGast1","rGast2")

a1 <- a[intersect(rownames(mouse_data2),rownames(rabbit_data2)),orders1]
a2 <- a[intersect(rownames(mouse_data2),rownames(rabbit_data2)),orders2]

C1 <- cor(log(a1 +1), log(a2+1), method = "pearson")

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
#mat_breaks <- seq(0.8, 1, length.out = 20)
pheatmap(C1,color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/Mouse_cross3",".pdf",sep=""),width=10,height=5)



#try 2
#rabbit_data3 <- merge(rabbit_data, y = c(rabbit_data2), project = "merged")
#Rename to keep variable convention for the join species modelling
#mammal.anchors <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data3, human_dataA1 ), dims = 1:20, anchor.features = 4000, k.filter = 20)
#mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
#DefaultAssay(mammal.combined) <- "integrated"
#mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
#mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE)
#mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
#DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/PCA_split","_v2_1.pdf",sep=""),width = 16, height = 4)
#DimPlot(mammal.combined, reduction = "tsne", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/TSNE_split","_v2_1.pdf",sep=""),width = 16, height = 4)
#DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_split","_v2_1.pdf",sep=""),width = 16, height = 4)

#Rename to keep variable convention for the join species modelling
#mammal.anchors <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data3, human_dataA1 ), dims = 1:20, anchor.features = 4000, k.filter = 40)
#mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
#DefaultAssay(mammal.combined) <- "integrated"
#mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
#mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE)
#mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
#DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/PCA_split","_v2_2.pdf",sep=""),width = 16, height = 4)
#DimPlot(mammal.combined, reduction = "tsne", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/TSNE_split","_v2_2.pdf",sep=""),width = 16, height = 4)
#DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_split","_v2_2.pdf",sep=""),width = 16, height = 4)

#mammal.combined$ID1 <- Idents(mammal.combined)

#Rename to keep variable convention for the join species modelling
#mammal.anchors <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data3, human_dataA1 ), dims = 1:20, anchor.features = 4000, k.filter = 60)
#mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
#DefaultAssay(mammal.combined) <- "integrated"
#mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
#mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE)
#mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
#DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/PCA_split","_v2_3.pdf",sep=""),width = 16, height = 4)
#DimPlot(mammal.combined, reduction = "tsne", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/TSNE_split","_v2_3.pdf",sep=""),width = 16, height = 4)
#DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_split","_v2_3.pdf",sep=""),width = 16, height = 4)



#Rename to keep variable convention for the join species modelling
#mammal.anchors <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data3, human_dataA1 ), dims = 1:20, anchor.features = 4000)
#mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)

#mammal.anchors1 <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data, human_dataA1 ), dims = 1:20, anchor.features = union(rownames(mammal.combined1), intersect(intersect(intersect(Markers,rownames(rabbit_data)),rownames(cynomolgous_data) ),rownames(human_dataA1) ) ) )
#mammal.combined1 <- IntegrateData(anchorset = mammal.anchors1, dims = 1:20)
#DefaultAssay(mammal.combined) <- "integrated"
#mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
#mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE)
#mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)
#saveRDS(mammal.combined, file = paste(saveext,"mammal.combined2.rds",sep=""))




#mammal.combined <- FindClusters(mammal.combined, resolution = 0.1)
#mammal.combined$Cl1 <- Idents(mammal.combined)
#DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl1",".pdf",sep=""),width = 20, height = 8)
##Idents(mammal.combined) <- factor(c( as.character(cylabs),  as.character(RabbitBS$TimePoint) )) #c(labs2,labscy) #labs2

#mammal.combined <- FindClusters(mammal.combined, resolution = 0.2)
#mammal.combined$Cl2 <- Idents(mammal.combined)
#DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl2",".pdf",sep=""),width = 20, height = 8)

#mammal.combined <- FindClusters(mammal.combined, resolution = 0.3)
#mammal.combined$Cl3 <- Idents(mammal.combined)
#DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl3",".pdf",sep=""),width = 20, height = 8)

#mammal.combined <- FindClusters(mammal.combined, resolution = 0.4)
#mammal.combined$Cl4 <- Idents(mammal.combined)
#DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl4",".pdf",sep=""),width = 20, height = 8)

#mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)
#mammal.combined$Cl5 <- Idents(mammal.combined)
#DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl5",".pdf",sep=""),width = 20, height = 8)

#mammal.combined <- FindClusters(mammal.combined, resolution = 0.6)
#mammal.combined$Cl6 <- Idents(mammal.combined)
#DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl6",".pdf",sep=""),width = 20, height = 8)

#mammal.combined <- FindClusters(mammal.combined, resolution = 0.7)
#mammal.combined$Cl7 <- Idents(mammal.combined)
#DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl7",".pdf",sep=""),width = 20, height = 8)

#mammal.combined <- FindClusters(mammal.combined, resolution = 0.8)
#mammal.combined$Cl8 <- Idents(mammal.combined)
#DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl8",".pdf",sep=""),width = 20, height = 8)

#mammal.combined <- FindClusters(mammal.combined, resolution = 0.9)
#mammal.combined$Cl9 <- Idents(mammal.combined)
#DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl9",".pdf",sep=""),width = 20, height = 8)


#mammal.combined <- FindClusters(mammal.combined, resolution = 1.1)
#mammal.combined$Cl11 <- Idents(mammal.combined)
#DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl11",".pdf",sep=""),width = 20, height = 8)

#DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/PCA_split_Cl11",".pdf",sep=""),width = 20, height = 8)


#Idents(mammal.combined) <- mammal.combined$ID1

#DimPlot(mammal.combined, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/PCA_split_anno",".pdf",sep=""),width = 20, height = 8)

#DimPlot(mammal.combined, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/TSNE_split_anno",".pdf",sep=""),width = 20, height = 8)

#DimPlot(mammal.combined, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_split_anno",".pdf",sep=""),width = 20, height = 8)
#Idents(mammal.combined) <- factor(c( as.character(cylabs),  as.character(RabbitBS$TimePoint) )) #c(labs2,labscy) #labs2

mammal.combined2$Cl <- Idents(mammal.combined)
#uID <- as.character(Idents(mammal.combined))
#uID[which(mammal.combined$Cells=="PGC")] <- "PGC" 
#uID[which(mammal.combined$Cells=="PGCLC")] <- "PGCLC" 
#uID[which(mammal.combined$Cells=="Soma")] <- "Soma" 
#uID[which(mammal.combined$Cells=="PSC")] <- "PSC" 
#uID[which(mammal.combined$Cells=="PSC_low")] <- "PSC_low" 

#Idents(mammal.combined) <- uID

species2 <- as.character(mammal.combined2$species)
species2[which(mammal.combined2$species=="1) Rabbit")] <- "r"
species2[which(mammal.combined2$species=="2) Cynomolgous")] <- "c"
species2[which(mammal.combined2$species=="4) Human (in vitro)")] <- "h"

Idents(mammal.combined2) <- paste(species2,Idents(mammal.combined2),sep="")
DefaultAssay(mammal.combined2) <- "RNA"


library("pheatmap")
avexp  <- AverageExpression(object = mammal.combined2, slot = "data")

Markers1 <- c(panel1,panel2,panel3,panel4,panel5,panel6,panel7,panel8,panel9,panel10,panel11)

a <- avexp$RNA
a <- a[Markers1,]

orders <- c("r0","r10","r5","r8","r2","r3")
a <- a[,orders]


redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(-2, 2, length.out =40)
aa <- na.omit( log2(a+1) )

aa <- aa[which( rowSums(aa) > 0 ),]
pheatmap(aa,color =  redblue1(40),breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, scale = "row",  filename = paste(saveext, "Rabbit_HM_logged.pdf"),width=4,height=16)

mat_breaks <- seq(0, 100, length.out =40)
pheatmap( (exp(aa)-1)*100,color =  redblue1(40),breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext, "Rabbit_HM_TPM.pdf"),width=4,height=16)


#BMP2, BMP4, BMP8
#WNT3, WNT5A, WNT5B
#FGF2, FGF4, FGF5, FGF10
#NODAL, TGFB1, TGFB2, TGFB3


FeaturePlot(mammal.combined, reduction = "umap", features = "BMP2", split.by = "species", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_BMP2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "BMP4", split.by = "species", cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_BMP4.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "BMP8", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_BMP8.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "WNT3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_WNT3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "WNT5A", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_WNT5A.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "WNT5B", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_WNT5B.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "FGF2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_FGF2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "FGF4", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_FGF4.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "FGF5", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_FGF5.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "FGF10", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_FGF10.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "NODAL", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_NODAL.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "TGFB1", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TGFB1.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "TGFB2", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TGFB2.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)
FeaturePlot(mammal.combined, reduction = "umap", features = "TGFB3", split.by = "species",  cols =  c("lightgrey", "darkblue"))
ggsave(filename=paste(saveext,"/Full_TGFB3.pdf",sep=""),width = 24, height = 8,limitsize = FALSE)




orders <- c("rPGC1","rPGC2","cPGC",
            "rHyp","cHyp","hHyp","rVE","cVE","hVE", "cSYS","hSYS",
            "rICM","cICM", "hICM",
            "hEpi","cEpi",
            "rEmDisc1","cEmDisc1","hEmDisc1",
            "rEmDisc2","cEmDisc2","hEmDisc2",
            "rEmDisc3","cEmDisc3","hEmDisc3",
            "rGast","cGast","hGast",
            "cGastEnd","hGastEnd",
            "cExMes","hExMes",
            "rTb","cPTb","hTb","rPTb","hPTb")
a <- a[,orders]





mamsubset <- subset(mammal.combined, idents = c("0","6","4","8","7","PGC","PGCLC"))


#Idents(mammal.combined) <- mammal.combined$Cl6
#mamsubset <- subset(mammal.combined, idents = c("2","4","9","3","1","6","8"), invert = TRUE)
Idents(mamsubset) <- mamsubset$Cells
mamsubset <- subset(mamsubset, idents = c("Tb_CS6","ExMes_CS5","ExMes_CS6/7","Am_CS5","Am_CS6","ICM_CS3","Epi_CS3","Epi_CS4"), invert = TRUE)


mamsubset <- FindVariableFeatures(mamsubset, selection.method = "vst", nfeatures = 20000)
mamsubset <- ScaleData(mamsubset, verbose = FALSE)

mamsubset <- RunPCA(mamsubset, npcs = 30, verbose = FALSE)
mamsubset <- RunUMAP(mamsubset, reduction = "pca", dims = 1:20)
mamsubset <- RunTSNE(mamsubset, reduction = "pca", dims = 1:20)
mamsubset <- FindNeighbors(mamsubset, reduction = "pca", dims = 1:20)

DimPlot(mamsubset, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_split_subset_anno",".pdf",sep=""),width = 20, height = 8)

DimPlot(mamsubset, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_subset_anno",".pdf",sep=""),width = 20, height = 8)

DimPlot(mamsubset, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_split_subset_anno",".pdf",sep=""),width = 20, height = 8)

mamsubset <- FindClusters(mamsubset, resolution = 0.3)
mamsubset$Cl3 <- Idents(mamsubset)
DimPlot(mamsubset, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl4",".pdf",sep=""),width = 20, height = 8)

DimPlot(mamsubset, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_split_Cl4",".pdf",sep=""),width = 20, height = 8)



uID <- as.character(mamsubset$Cl3)
uID[which(mamsubset$Cl3=="2")] <- "PGC"
#uID[which(mamsubset$Cells=="PGC")] <- "PGC"
uID[which(mamsubset$Cells=="PGCLC")] <- "PGCLC"
uID[which(mamsubset$Cells=="ESC_primed")] <- "cPSC"
uID[which(mamsubset$Cells=="PSC")] <- "rPSC"
uID[which(mamsubset$Cells=="PSC_low")] <- "rPSC_low"

Idents(mamsubset) <- uID

DimPlot(mamsubset, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_new",".pdf",sep=""),width = 20, height = 8)

DimPlot(mamsubset, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_split_new",".pdf",sep=""),width = 20, height = 8)




colind <- integer( length( levels(Idents(mamsubset)) )  )
for (i in 1:length( levels(Idents(mamsubset)) ) ) {
  colind[i] <- which(cType==levels(Idents(mamsubset))[i])
}
coluse <- BaseCol[colind]

cluster_letters <- LETTERS[mamsubset$orig.ident]
cluster_letters[1:length(cluster_letters)] <- 0
cluster_letters[which(mamsubset$ID1=="E6")] <- 1
cluster_letters[which(mamsubset$ID1=="E7")] <- 2
cluster_letters[which(mamsubset$Cells=="Soma")] <- 3
cluster_letters[which(mamsubset$Cells=="PGC")] <- 3
cluster_letters[which(mamsubset$Cells=="PGCLC")] <- 4
cluster_letters[which(mamsubset$Cells=="PSC")] <- 4
cluster_letters[which(mamsubset$Cells=="PSC_low")] <- 4
cluster_letters[which(mamsubset$Cells=="ESC_primed")] <- 4
cluster_letters <- as.factor(cluster_letters)
names(cluster_letters)=rownames(mamsubset@meta.data)
mamsubset <- AddMetaData(object = mamsubset,metadata = cluster_letters,col.name = 'cell.orig')

#15,16,17,18,25

#And plot these
DimPlot(mamsubset, reduction = "pca", shape.by = "cell.orig", pt.size = 2,  label.size = 2,  split.by = "species", label = TRUE, repel = TRUE) & scale_shape_manual(values=c(15,16,17,18,25)) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_shape",".pdf",sep=""),width = 24, height = 8, useDingbats = FALSE)





Idents(mamsubset) <- colnames(mamsubset)
p <- DimPlot(mamsubset, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_split_subset_fulllabels",".pdf",sep=""),width = 100, height = 50, p, limitsize = FALSE)

Idents(mamsubset) <- colnames(mamsubset)
p <- DimPlot(mamsubset, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_subset_fulllabels",".pdf",sep=""),width = 100, height = 50, p, limitsize = FALSE)


#mamsubset <- FindClusters(mamsubset, resolution = 10.5)

#DimPlot(mamsubset, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/UMAP_split_subset_Cl",".pdf",sep=""),width = 20, height = 8)
#ggsave(filename=paste(saveext,"/DimRed/PCA_split_subset_Cl",".pdf",sep=""),width = 20, height = 8)
##DimPlot(mamsubset, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)

#wC1 <- WhichCells(mamsubset, idents = "11")
Idents(mamsubset) <- mamsubset$species
wC1 <-c("SLX.17852.i729_i515","SLX.17852.i715_i511","SLX.17852.i728_i522","SLX.17852.i712_i515")
wC2 <-c("SLX.17852.i729_i517","SLX.17852.i729_i507","SLX.17852.i729_i508","SLX.12935.i707_i506","SLX.12935.i702_i521")
wC3 <- c("SLX.12935.i704_i508","SLX.12935.i719_i505","SLX.12935.i716_i508","SLX.17852.i714_i507" )
#wC2 <- WhichCells(mamsubset, idents ="1) Rabbit")
#wC3 <- WhichCells(object = mamsubset, expression = NANOS3 > 0)
#wC4 <- intersect(wC2,wC3)

#wC3 <- subset(mammal.combined, subset = NANOS3 > 1)

Idents(mammal.combined) <- mammal.combined$Cl11
mamsubset1 <- subset(mammal.combined, idents = c("5","8","7"), invert = FALSE)
Idents(mamsubset1) <- mamsubset1$ID1
#mamsubset1 <- subset(mamsubset1, idents = c("Tb_CS5"), invert = TRUE)

mamsubset1 <- RunPCA(mamsubset1, npcs = 30, verbose = FALSE)
mamsubset1 <- RunUMAP(mamsubset1, reduction = "pca", dims = 1:20)
mamsubset1 <- RunTSNE(mamsubset1, reduction = "pca", dims = 1:20)
mamsubset1 <- FindNeighbors(mamsubset1, reduction = "pca", dims = 1:20)


#Subset on rabbit
#Idents(mamsubset) <- mamsubset$species
#subset2 <- subset(mamsubset, idents = "1) Rabbit")
Loadis <- as.data.frame(Loadings(mamsubset1, reduction = "pca")[, 1:2])
p <- ggplot(Loadis, aes(x=PC_1, y=PC_2)) + geom_point() + geom_text(label=rownames(Loadis))+ theme_classic()
ggsave(filename=paste(saveext,"/DimRed/PCA_Loadingn_Rabbit",".pdf",sep=""),width = 40, height = 40,p)


Loadis <- as.data.frame(Loadings(mamsubset, reduction = "pca")[, 1:2])
p <- ggplot(Loadis, aes(x=PC_1, y=PC_2)) + geom_point() + geom_text(label=rownames(Loadis))+ theme_classic()
ggsave(filename=paste(saveext,"/DimRed/PCA_Loadingn_EmRabbit",".pdf",sep=""),width = 40, height = 40,p)

#subset2 <- RunPCA(subset2, npcs = 30, verbose = FALSE)
#Loadis <- as.data.frame(Loadings(subset2, reduction = "pca")[, 1:2])
#p <- ggplot(Loadis, aes(x=PC_1, y=PC_2)) + geom_point() + geom_text(label=rownames(Loadis))+ theme_classic()
#ggsave(filename=paste(saveext,"/DimRed/PCA_Loadingn_Rabbit2",".pdf",sep=""),width = 40, height = 40,p)
#Idents(subset2) <- subset2$newID
#DimPlot(subset2, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/PCA_Loadingn_Rabbit3",".pdf",sep=""),width = 20, height = 8)


DimPlot(mamsubset1, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_split_subset2_anno",".pdf",sep=""),width = 20, height = 8)

DimPlot(mamsubset1, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_split_subset2_anno",".pdf",sep=""),width = 20, height = 8)

DimPlot(mamsubset1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_subset2_anno",".pdf",sep=""),width = 20, height = 8)

#mamsubset1 <- FindClusters(mamsubset1, resolution = 0.5)

DimPlot(mamsubset1, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_subset2_Cl",".pdf",sep=""),width = 20, height = 8)
DimPlot(mamsubset1, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_split_subset2_Cl",".pdf",sep=""),width = 20, height = 8)


for (i in 1:length( Markers ) ) {
#  VlnPlot(mammal.combined, features = Markers[i])
#  ggsave(filename=paste(saveext,"/Markers/Volcano_full_",Markers[i],".pdf",sep=""),split.by = "species",width = 10, height = 20,limitsize = FALSE)
  FeaturePlot(mammal.combined, features = Markers[i], split.by = "species", cols =  c("lightgrey", "darkblue"), pt.size = 1)
  ggsave(filename=paste(saveext,"/Markers/Marker_full_",Markers[i],".pdf",sep=""),width = 15, height = 5,limitsize = FALSE)
#  VlnPlot(mamsubset, features = Markers[i])
#  ggsave(filename=paste(saveext,"/Markers/Volcano_sub_",Markers[i],".pdf",sep=""),split.by = "species",width = 10, height = 20,limitsize = FALSE)
  FeaturePlot(mamsubset, features = Markers[i], split.by = "species",cols =  c("lightgrey", "darkblue"), pt.size = 1)
  ggsave(filename=paste(saveext,"/Markers/Marker_sub_",Markers[i],".pdf",sep=""),width = 15, height = 5,limitsize = FALSE)
  
  FeaturePlot(mamsubset1, features = Markers[i], split.by = "species",cols =  c("lightgrey", "darkblue"), pt.size = 1)
  ggsave(filename=paste(saveext,"/Markers/Marker_sub1_",Markers[i],".pdf",sep=""),width = 15, height = 5,limitsize = FALSE)
  
  
  FeaturePlot(mammal.combined, reduction = "pca", features = Markers[i], split.by = "species",cols =  c("lightgrey", "darkblue"), pt.size = 1)
  ggsave(filename=paste(saveext,"/Markers/Marker_pca_full_",Markers[i],".pdf",sep=""),width = 15, height = 5,limitsize = FALSE)
  #  VlnPlot(mamsubset, features = Markers[i])
  #  ggsave(filename=paste(saveext,"/Markers/Volcano_sub_",Markers[i],".pdf",sep=""),split.by = "species",width = 10, height = 20,limitsize = FALSE)
  FeaturePlot(mamsubset, reduction = "pca", features = Markers[i], split.by = "species",cols =  c("lightgrey", "darkblue"), pt.size = 1)
  ggsave(filename=paste(saveext,"/Markers/Marker_pca_sub_",Markers[i],".pdf",sep=""),width = 15, height = 5,limitsize = FALSE)
  
  FeaturePlot(mamsubset1, reduction = "pca", features = Markers[i], split.by = "species",cols =  c("lightgrey", "darkblue"), pt.size = 1)
  ggsave(filename=paste(saveext,"/Markers/Marker_pca_sub1_",Markers[i],".pdf",sep=""),width = 15, height = 5,limitsize = FALSE)
  
  
}


newID <- as.character(mammal.combined$ID1)
newID[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl11==1)] <- "Hyp"
newID[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl11==11)] <- "VE"
newID[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl11==0)] <- "EmDisc3"
newID[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl11==2)] <- "EmDisc1"
newID[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl11==10)] <- "EmDisc2"
newID[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl11==9)] <- "Epi"
newID[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl11==8)] <- "ICM"
newID[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl11==5)] <- "Tb"
newID[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl11==12)] <- "PGC"

newID[which(mammal.combined$species=="4) Human (in vitro)" & mammal.combined$Cl11==0)] <- "EmDisc3"
newID[which(mammal.combined$species=="4) Human (in vitro)" & mammal.combined$Cl11==2)] <- "EmDisc1"
newID[which(mammal.combined$species=="4) Human (in vitro)" & mammal.combined$Cl11==10)] <- "EmDisc2"

newID[which(mammal.combined$species=="2) Cynomolgous" & mammal.combined$Cl11==0)] <- "EmDisc3"
newID[which(mammal.combined$species=="2) Cynomolgous" & mammal.combined$Cl11==2)] <- "EmDisc1"
newID[which(mammal.combined$species=="2) Cynomolgous" & mammal.combined$Cl11==10)] <- "EmDisc2"


newID[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl11==3)] <- "PTb"
newID[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl11==4)] <- "PTb"
#newID[which(mammal.combined$species=="2) Cynomolgous" & mammal.combined$Cl11==3)] <- "PTb"
#newID[which(mammal.combined$species=="2) Cynomolgous" & mammal.combined$Cl11==4)] <- "PTb"
#newID[which(mammal.combined$species=="4) Human (in vitro)" & mammal.combined$Cl11==3)] <- "PTb"
#newID[which(mammal.combined$species=="4) Human (in vitro)" & mammal.combined$Cl11==4)] <- "PTb"

newID[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl11==6)] <- "Gast"
newID[which(mammal.combined$species=="2) Cynomolgous" & mammal.combined$Cl11==6)] <- "Gast"
newID[which(mammal.combined$species=="4) Human (in vitro)" & mammal.combined$Cl11==6)] <- "Gast"

newID[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl11==13)] <- "ExMes"
newID[which(mammal.combined$species=="4) Human (in vitro)" & mammal.combined$Cl11==13)] <- "ExMes"
newID[which(mammal.combined$species=="2) Cynomolgous" & mammal.combined$Cl11==13)] <- "ExMes"


newID[which(newID=="ExMes_CS5")] <- "ExMes"
newID[which(newID=="ExMes_CS6/7")] <- "ExMes"
newID[which(newID=="Tb_CS3")] <- "Tb"
newID[which(newID=="ICM_CS3")] <- "ICM"
newID[which(newID=="Epi_CS3")] <- "Epi"
#newID[which(newID=="EmDisc_CS5")] <- "eEmDisc"


#newID[which(newID=="EmDisc_CS6")] <- "lEmDisc"
newID[which(newID=="PGC_CS5")] <- "PGC"
newID[which(newID=="PGC_CS6/7")] <- "PGC"
newID[which(newID=="Tb_CS6/7")] <- "PTb"
newID[which(newID=="Tb_CS6")] <- "PTb"
newID[which(newID=="Tb_CS5")] <- "PTb"
newID[which(newID=="VE_CS4")] <- "VE"
newID[which(newID=="VE_CS5")] <- "VE"
newID[which(newID=="SYS_CS6")] <- "SYS"
newID[which(newID=="Hyp_CS3")] <- "Hyp"
newID[which(newID=="SYS_CS6/7")] <- "SYS"
newID[which(newID=="PGC_E50")] <- "PGC"
newID[which(newID=="Tb_CS4")] <- "Tb"
newID[which(newID=="Epi_CS4")] <- "Epi"

newID[which(newID=="EmDisc_CS5")] <- "GastEnd"

newID[which(newID=="EmDisc_CS6")] <- "GastEnd"
newID[which(newID=="EmDisc_CS6/7")] <- "GastEnd"

#newID[which(colnames(mammal.combined) %in% wC1)] <- "PGC1"
#newID[which(colnames(mammal.combined) %in% wC2)] <- "PGC2"
#newID[which(colnames(mammal.combined) %in% wC3)] <- "PGC3"


Idents(mammal.combined) <- as.factor(newID)

DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_subset_newanno",".pdf",sep=""),width = 20, height = 8)

DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_split_subset_newanno",".pdf",sep=""),width = 20, height = 8)

DimPlot(mammal.combined, reduction = "tsne", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_split_subset_newanno",".pdf",sep=""),width = 20, height = 8)




newID <- as.character(mamsubset$ID1)
newID[which(mamsubset$species=="1) Rabbit" & mamsubset$Cl11==1)] <- "Hyp"
newID[which(mamsubset$species=="1) Rabbit" & mamsubset$Cl11==11)] <- "VE"
newID[which(mamsubset$species=="1) Rabbit" & mamsubset$Cl11==0)] <- "EmDisc3"
newID[which(mamsubset$species=="1) Rabbit" & mamsubset$Cl11==2)] <- "EmDisc1"
newID[which(mamsubset$species=="1) Rabbit" & mamsubset$Cl11==10)] <- "EmDisc2"
newID[which(mamsubset$species=="1) Rabbit" & mamsubset$Cl11==9)] <- "Epi"
newID[which(mamsubset$species=="1) Rabbit" & mamsubset$Cl11==8)] <- "ICM"
newID[which(mamsubset$species=="1) Rabbit" & mamsubset$Cl11==5)] <- "Tb"
newID[which(mamsubset$species=="1) Rabbit" & mamsubset$Cl11==12)] <- "PGC"
newID[which(mamsubset$species=="1) Rabbit" & mamsubset$Cl11==3)] <- "Tb"
newID[which(mamsubset$species=="1) Rabbit" & mamsubset$Cl11==4)] <- "Tb"
newID[which(mamsubset$species=="1) Rabbit" & mamsubset$Cl11==6)] <- "Gast"
newID[which(mamsubset$species=="1) Rabbit" & mamsubset$Cl11==13)] <- "ExMes"



newID[which(mamsubset$species=="4) Human (in vitro)" & mamsubset$Cl11==0)] <- "EmDisc3"
newID[which(mamsubset$species=="4) Human (in vitro)" & mamsubset$Cl11==2)] <- "EmDisc1"
newID[which(mamsubset$species=="4) Human (in vitro)" & mamsubset$Cl11==10)] <- "EmDisc2"

newID[which(mamsubset$species=="2) Cynomolgous" & mamsubset$Cl11==0)] <- "EmDisc3"
newID[which(mamsubset$species=="2) Cynomolgous" & mamsubset$Cl11==2)] <- "EmDisc1"
newID[which(mamsubset$species=="2) Cynomolgous" & mamsubset$Cl11==10)] <- "EmDisc2"

newID[which(mamsubset$species=="1) Rabbit" & mamsubset$Cl11==3)] <- "PTb"
newID[which(mamsubset$species=="1) Rabbit" & mamsubset$Cl11==4)] <- "PTb"
#newID[which(mammal.combined$species=="2) Cynomolgous" & mammal.combined$Cl11==3)] <- "PTb"
#newID[which(mammal.combined$species=="2) Cynomolgous" & mammal.combined$Cl11==4)] <- "PTb"
#newID[which(mammal.combined$species=="4) Human (in vitro)" & mammal.combined$Cl11==3)] <- "PTb"
#newID[which(mammal.combined$species=="4) Human (in vitro)" & mammal.combined$Cl11==4)] <- "PTb"

newID[which(mamsubset$species=="1) Rabbit" & mamsubsetd$Cl11==6)] <- "Gast"
newID[which(mamsubset$species=="2) Cynomolgous" & mamsubset$Cl11==6)] <- "Gast"
newID[which(mamsubset$species=="4) Human (in vitro)" & mamsubset$Cl11==6)] <- "Gast"

newID[which(mamsubset$species=="1) Rabbit" & mamsubset$Cl11==13)] <- "ExMes"
newID[which(mamsubset$species=="4) Human (in vitro)" & mamsubset$Cl11==13)] <- "ExMes"
newID[which(mamsubset$species=="2) Cynomolgous" & mamsubset$Cl11==13)] <- "ExMes"



newID[which(newID=="ExMes_CS5")] <- "ExMes"
newID[which(newID=="ExMes_CS6/7")] <- "ExMes"
newID[which(newID=="Tb_CS3")] <- "Tb"
newID[which(newID=="ICM_CS3")] <- "ICM"
newID[which(newID=="Epi_CS3")] <- "Epi"
#newID[which(newID=="EmDisc_CS5")] <- "EmDisc"
#newID[which(newID=="EmDisc_CS6/7")] <- "EmDisc"
#newID[which(newID=="EmDisc_CS6")] <- "EmDisc"
newID[which(newID=="PGC_CS5")] <- "PGC"
newID[which(newID=="PGC_CS6/7")] <- "PGC"
newID[which(newID=="Tb_CS6/7")] <- "PTb"
newID[which(newID=="Tb_CS6")] <- "PTb"
newID[which(newID=="Tb_CS5")] <- "PTb"
newID[which(newID=="VE_CS4")] <- "VE"
newID[which(newID=="VE_CS5")] <- "VE"
newID[which(newID=="SYS_CS6")] <- "SYS"
newID[which(newID=="Hyp_CS3")] <- "Hyp"
newID[which(newID=="SYS_CS6/7")] <- "SYS"
newID[which(newID=="PGC_E50")] <- "PGC"
newID[which(newID=="Tb_CS4")] <- "Tb"
newID[which(newID=="Epi_CS4")] <- "Epi"
#newID[which(colnames(mamsubset) %in% wC1)] <- "PGC1"
#newID[which(colnames(mamsubset) %in% wC2)] <- "PGC2"

newID[which(newID=="EmDisc_CS5")] <- "GastEnd"

newID[which(newID=="EmDisc_CS6")] <- "GastEnd"
newID[which(newID=="EmDisc_CS6/7")] <- "GastEnd"

#newID[which(colnames(mamsubset) %in% wC1)] <- "PGC1"
#newID[which(colnames(mamsubset) %in% wC2)] <- "PGC2"
#newID[which(colnames(mamsubset) %in% wC3)] <- "PGC3"



Idents(mamsubset) <- as.factor(newID)
DimPlot(mamsubset, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_split_subset1_newanno",".pdf",sep=""),width = 20, height = 8)

DimPlot(mamsubset, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_subset1_newanno",".pdf",sep=""),width = 20, height = 8)




species2 <- mammal.combined2$species
species2[which(mammal.combined2$species=="1) Rabbit")] <- "r"
species2[which(mammal.combined2$species=="2) Cynomolgous")] <- "c"
species2[which(mammal.combined2$species=="4) Human (in vitro)")] <- "h"

Idents(mammal.combined) <- paste(species2,Idents(mammal.combined2),sep="")
DefaultAssay(mammal.combined2) <- "RNA"


library("pheatmap")
avexp  <- AverageExpression(object = mammal.combined2, slot = "data")

a <- avexp$RNA
a <- a[Markers1,]

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(-2, 2, length.out = 40)
pheatmap(log2(a+1),color =  redblue1(40),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE,  filename = paste(saveext,"Rabbit_logged.pdf",sep=""),width=10,height=16)

#pheatmap(log2(a+1),color =  redblue1(40),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE,  filename = paste(saveext,"Rabbit_logged.pdf",sep=""),width=10,height=16)


#intgenes <- rownames(avexp$integrated)
#b <- avexp$RNA[intgenes,orders]

C1 <- cor(b, method = "pearson")

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
#mat_breaks <- seq(0, 2, length.out = 20)
#pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE,  filename = "~/Desktop/HM_ara_logged.pdf",,width=10,height=16)

#mat_breaks <- seq(-1, 1, length.out = 20)
#pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, scale="row", filename = "~/Desktop/HM_ara_logged_scaled.pdf",width=10,height=16)

mat_breaks <- seq(0, 2, length.out = 20)
pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks,gaps_col=c(3,11,16,25,30,32),gaps_row=c(6,11,16,21,29), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/HM_ara_logged",".pdf",sep=""),width=10,height=16)

mat_breaks <- seq(-1, 1, length.out = 20)
pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks,gaps_col=c(3,11,16,25,30,32),gaps_row=c(6,11,16,21,29), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, scale="row", filename = paste(saveext,"/DimRed/HM_ara_logged_scaled",".pdf",sep=""),width=10,height=16)


orders <- c("rPGC1","rPGC2","cPGC",
            "rHyp","cHyp","hHyp","rVE","cVE","hVE", "cSYS","hSYS",
            "rICM","cICM", "hICM",
            "hEpi","cEpi",
            "rEmDisc1","cEmDisc1","hEmDisc1",
            "rEmDisc2","cEmDisc2","hEmDisc2",
            "rEmDisc3","cEmDisc3","hEmDisc3",
            "rGast","cGast","hGast",
            "cGastEnd","hGastEnd",
            "cExMes","hExMes",
            "rTb","cPTb","hTb","rPTb","hPTb")
a <- a[,orders]


intgenes <- rownames(avexp$integrated)
b <- avexp$RNA[intgenes,orders]

C1 <- cor(b, method = "pearson")

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
#mat_breaks <- seq(0, 2, length.out = 20)
#pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE,  filename = "~/Desktop/HM_ara_logged.pdf",,width=10,height=16)

#mat_breaks <- seq(-1, 1, length.out = 20)
#pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, scale="row", filename = "~/Desktop/HM_ara_logged_scaled.pdf",width=10,height=16)
         
mat_breaks <- seq(0, 2, length.out = 20)
pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks,gaps_col=c(3,11,16,25,30,32),gaps_row=c(6,11,16,21,29), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/HM_ara_logged",".pdf",sep=""),width=10,height=16)

mat_breaks <- seq(-1, 1, length.out = 20)
pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks,gaps_col=c(3,11,16,25,30,32),gaps_row=c(6,11,16,21,29), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, scale="row", filename = paste(saveext,"/DimRed/HM_ara_logged_scaled",".pdf",sep=""),width=10,height=16)



mat_breaks <- seq(0.5, 1, length.out = 20)
pheatmap(C1,color =  redblue1(20),breaks = mat_breaks,gaps_col=c(3,11,16,25,30,32),gaps_row=c(3,11,16,25,30,32), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/Sample_corr",".pdf",sep=""),width=16,height=16)



ronly <- c("rICM","rEmDisc1", "rEmDisc2", "rEmDisc3","rGast","rPTb","rTb","rHyp","rVE")
a <- a[,ronly]

b <- avexp$RNA[intgenes,ronly]
C1 <- cor(log2(b+1), method = "pearson")


mat_breaks <- seq(0.5, 1, length.out = 20)
pheatmap(C1,color =  redblue1(20),breaks = mat_breaks,gaps_col=c(1,4,5,7,9),gaps_row=c(1,4,5,7,9), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/Sampler_corr",".pdf",sep=""),width=16,height=16)



redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
#mat_breaks <- seq(0, 2, length.out = 20)
#pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE,  filename = "~/Desktop/HM_ara_logged.pdf",,width=10,height=16)

#mat_breaks <- seq(-1, 1, length.out = 20)
#pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE, scale="row", filename = "~/Desktop/HM_ara_logged_scaled.pdf",width=10,height=16)

mat_breaks <- seq(0, 2, length.out = 20)
pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks,gaps_col=c(1,4,5,7,9),gaps_row=c(6,10,13,18,26), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/HMr_ara_logged",".pdf",sep=""),width=4,height=6)

mat_breaks <- seq(-1, 1, length.out = 20)
pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks,gaps_col=c(1,4,5,7,9),gaps_row=c(6,10,13,18,26), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, scale="row", filename = paste(saveext,"/DimRed/HMr_ara_logged_scaled",".pdf",sep=""),width=4,height=6)


l1<-grep("rICM", Idents(mammal.combined))
l2<-grep("rEmDisc1", Idents(mammal.combined))
l3<-grep("rEmDisc2", Idents(mammal.combined))
l4<-grep("rEmDisc3", Idents(mammal.combined))
l5<-grep("rGast", Idents(mammal.combined))
l6<-grep("rPTb", Idents(mammal.combined))
l7<-grep("rTb", Idents(mammal.combined))
l8<-grep("rHyp", Idents(mammal.combined))
l9<-grep("rVE", Idents(mammal.combined))

#l10<-grep("rPGC1", Idents(mammal.combined))
#l12<-grep("rPGC3", Idents(mammal.combined))
##l11<-grep("rPGC2", Idents(mammal.combined))

lfull <- c(l1,l2,l3,l4,l5,l6,l7,l8,l9)

DN1 <- GetAssayData(mammal.combined, slot = "data")


mat_breaks <- seq(0, 2, length.out = 20)
pheatmap(DN1[Markers,lfull],color =  redblue1(20),breaks = mat_breaks,gaps_col=c(9,63,68,153,216,224,248,320,344),gaps_row=c(6,10,13,18,26), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/HMrf_ara_logged",".pdf",sep=""),width=20,height=6)

mat_breaks <- seq(-1, 1, length.out = 20)
pheatmap(DN1[Markers,lfull],color =  redblue1(20),breaks = mat_breaks,gaps_col=c(9,63,68,153,216,224,248,320,344),gaps_row=c(6,10,13,18,26), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, scale="row", filename = paste(saveext,"/DimRed/HMrf_ara_logged_scaled",".pdf",sep=""),width=20,height=6)



ronly_data <- subset(mammal.combined, idents = ronly, invert = FALSE )

markerscl <- FindAllMarkers(ronly_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use="MAST")

top100 <- markerscl %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

a <- avexp$RNA
a <- a[top100$gene,]
a <- a[,ronly]
mat_breaks <- seq(-1, 1, length.out = 20)
pheatmap(log2(a+1),color =  redblue1(20), breaks = mat_breaks,border_color = NA,gaps_col=c(1,4,5,7,9), cluster_rows=FALSE,cluster_cols=FALSE, scale="row", filename = paste(saveext,"/DimRed/HM_unbiased_ara_logged_scaled",".pdf",sep=""),width=4,height=40)

mat_breaks <- seq(0, 2, length.out = 20)
pheatmap(log2(a+1),color =  redblue1(20), breaks = mat_breaks,border_color = NA,gaps_col=c(1,4,5,7,9), cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/HM_unbiased_ara_logged",".pdf",sep=""),width=4,height=40)


write.csv(as.data.frame(markerscl), file=paste(saveext,"/Markers/Unbiased.csv",sep=""))


#DoHeatmap(mammal.combined, features = top100$gene) + NoLegend()
#ggsave(filename=paste(saveext,"/Markers/HeatMapcl_100",".pdf",sep=""),width = 5, height = 40)
#write.csv(as.data.frame(mammal.combined.markerscl), file=paste(saveext,"/Markers/Markers_cl.csv",sep=""))


#Write everythign out to csv
#write.csv(as.data.frame(Idents(object = mammal.combined)), file=paste(saveext,"EmbeddingsCl.csv",sep=""))
#write.csv(as.data.frame(Embeddings(object = mammal.combined[["umap"]])), file=paste(saveext,"Embeddings.csv",sep=""))
#write.csv(as.data.frame(mammal.combined[[]]), file=paste(saveext,"EmbeddingsKey.csv",sep=""))
#write.csv(as.data.frame(Embeddings(object = mammal.combined[["pca"]])), file=paste(saveext,"EmbeddingsPCA.csv",sep=""))

#Dis1 <- as.data.frame(as.matrix(dist(Embeddings(object = mammal.combined[["umap"]]), method = "euclidean",diag = TRUE, upper = TRUE)))
#Dis2 <- as.data.frame(as.matrix(dist(Embeddings(object = mammal.combined[["tsne"]]), method = "euclidean",diag = TRUE, upper = TRUE)))
#Dis3 <- as.data.frame(as.matrix(dist(Embeddings(object = mammal.combined[["pca"]]), method = "euclidean",diag = TRUE, upper = TRUE)))

#Dis1[Dis1 == 0] <- 1000
#Dis2[Dis2 == 0] <- 1000
#Dis3[Dis3 == 0] <- 1000

#colind1 <- integer( length( Dis1))
#colind2 <- integer( length( Dis1))
#colind3 <- integer( length( Dis1)) 

#colind4 <- integer( length( Dis1))
#colind5 <- integer( length( Dis1))
#colind6 <- integer( length( Dis1)) 

#inds1 <- which(mammal.combined$species=="2) Cynomolgous")
#inds2 <- which(mammal.combined$species=="4) Human (in vitro)")

#for (i in 1:length( Dis1 ) ) {
#  
#  colind1[i] <- which( Dis1[i,]==min(Dis1[i,]) )
#  colind2[i] <- which( Dis2[i,]==min(Dis2[i,]) )
#  colind3[i] <- which( Dis3[i,]==min(Dis3[i,]) )  #
#
#  colind4[i] <- which( Dis1[i,]==min(Dis1[i,]) )
#  colind5[i] <- which( Dis2[i,]==min(Dis2[i,]) )
#  colind6[i] <- which( Dis3[i,]==min(Dis3[i,]) )  
#}





#Idents(mammal.combined) <- Cls
#Write everythign out to csv


#Cls <- Idents(mammal.combined)
#Idents(mammal.combined) <- factor(c( as.character(cylabs),  as.character(RabbitBS$TimePoint), as.character(le3[,2]), as.character(le4[,2]), as.character(le5[,2]), as.character(le6[,2]) )) #labs3



#Epi - reds
#Tb - blues
#Hypo Green
#PGCs yellow

cType <- c(
  "cyICM_E6","cypreEpi_E7","cypreEpi_E8","cypreEpi_E9",
  "cyEmDisc_E13","cyEmDisc_E14","cyEmDisc_E16","cyEmDisc_E17",
  "cyGast_E13","cyGast_E14","cyGast_E16","cyGast_E17","cyGast_E20",
  "cyePGC","cylPGC",
  "cyHypo_E7","cyHypo_E8","cyHypo_E9","cySYS",
  "E4","E5","E6","E7",
  "cyTb_E6","cyTb_E7","cyTb_E8","cyTb_E9","cyTb_E14","cyTb_E16",
  "EPI","PE","PGC","TE")

cType <- c("cPGC","cExMes","cGastEnd","cEmDisc3","cGast","cEmDisc2","cHyp","cICM","cPTb","cEmDisc1","cEpi","cTb","cVE","cSYS","rVE","rEmDisc3","rEmDisc1","rHyp","rPTb","rTb","rEmDisc2","rGast","rICM","hPTb","hExMes","hGast","hEmDisc1","hEmDisc3","hGastEnd","hSYS","hEmDisc2","hVE","hEpi","hTb","hICM","hHyp")  
BaseCol<-c("#c51b7d","#01665e","#cb181d","#abd9e9","#cb181d","#74add1","#fd8d3c","#2d004b","#1b7837","#4575b4","#313695","#7fbc41","#feb24c","#feb24c","#fe9929","#fe9929","#2166ac","#ec7014","#00441b","#4d9221","#4393c3","#a50f15","#053061","#5aae61","#35978f","#ef3b2c","#045a8d","#3690c0","#ef3b2c","#fdd0a2","#0570b0","#fdd0a2","#023858","#b8e186","#542788","#fdae6b")

#01665e
#276419,	#4d9221,#b2182b,#053061,	#276419,"exmes",#b2182b,#4393c3,#d1e5f0,	#b2182b, #fdb863, #92c5de,#fdb863
#BaseCol <- c("#c7e9c0","#a1d99b","#74c476","#41ab5d",
#             "#ffffe5","#fff7bc","#fee391","#fec44f",
#             "#fe9929","#ec7014","#cc4c02","#993404","#662506",
#             "#252525","#000000",
#             "#41ab5d","#238b45","#006d2c","#00441b",
#             "#6a51a3","#238b45","#d94801","#2171b5",
#             "#9ecae1","#6baed6","#4292c6","#2171b5","#08519c","#08306b",
#             "#fec44f","#238b45", "#969696","#6baed6")

#"Tb_CS5","Tb_CS6","Tb_CS7", "EmDisc_CS7"
cType <- c("cPGC","cExMes","cGastEnd","cEmDisc3","cGast","cEmDisc2","cHyp","cICM","cPTb","cEmDisc1","cEpi","cTb","cVE","cSYS","rVE","rEmDisc3","rEmDisc1","rHyp","rPTb","rTb","rEmDisc2","rGast","rICM","hPTb","hExMes","hGast","hEmDisc1","hEmDisc3","hGastEnd","hSYS","hEmDisc2","hVE","hEpi","hTb","hICM","hHyp")  
BaseCol<-c("#c51b7d","#01665e","#cb181d","#abd9e9","#cb181d","#74add1","#fd8d3c","#2d004b","#1b7837","#4575b4","#313695","#7fbc41","#feb24c","#feb24c","#fe9929","#fe9929","#2166ac","#ec7014","#00441b","#4d9221","#4393c3","#a50f15","#053061","#5aae61","#35978f","#ef3b2c","#045a8d","#3690c0","#ef3b2c","#fdd0a2","#0570b0","#fdd0a2","#023858","#b8e186","#542788","#fdae6b")

#"EmDisc_CS5","EmDisc_CS6"
# "SYS_CS6"            "VE_CS5"       "Epi_CS4"     
#"VE_CS4"       "Tb_CS4"       "Epi_CS3"      "ICM_CS3"      "Hyp_CS3"      "ExMes_CS5"    "Am_CS5"       "ExMes_CS6"    "Am_CS6"      
#[46] "PGC_CS6"      "ExMes_CS7"      "Tb_CS7"       "SYS_CS7"      "Am_CS7"     


colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

cluster_letters <- LETTERS[mammal.combined$orig.ident]
cluster_letters[1:length(cluster_letters)] <- "E5"
cluster_letters[which(mammal.combined$ID3=="E4")] <- "E4"
cluster_letters[which(mammal.combined$ID3=="E5")] <- "E5"
cluster_letters[which(mammal.combined$ID3=="E6")] <- "E6"
cluster_letters[which(mammal.combined$ID3=="E7")] <- "E7"
cluster_letters <- as.factor(cluster_letters)


#cluster_letters[grep("_F", mammal.combined$ID3)] <- "E7"


#mammal.combined <- AddMetaData(object = mammal.combined,metadata = cluster_letters,col.name = 'cell.orig')
names(cluster_letters)=rownames(mammal.combined@meta.data)
mammal.combined <- AddMetaData(object = mammal.combined,metadata = cluster_letters,col.name = 'cell.orig')



DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 2, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2",".pdf",sep=""),width = 20, height = 6)


DimPlot(mammal.combined, reduction = "pca", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 2, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab2",".pdf",sep=""),width = 20, height = 6)



species2 <- mamsubset$species
species2[which(mamsubset$species=="1) Rabbit")] <- "r"
species2[which(mamsubset$species=="2) Cynomolgous")] <- "c"
species2[which(mamsubset$species=="4) Human (in vitro)")] <- "h"

Idents(mamsubset) <- paste(species2,Idents(mamsubset),sep="")


#mammal.combined$FinLab <- Idents(mammal.combined)
#mamsubset2 <- subset(mammal.combined, idents = c("2","10","0","12","6"), invert = FALSE)
#Idents(mamsubset2) <- mamsubset2$ID1
##Idents(mammal.combined) <- mammal.combined$Cl11
#mamsubset2 <- subset(mamsubset2, idents = c("Tb_CS5"), invert = TRUE)
#mamsubset2 <- RunPCA(mamsubset2, npcs = 30, verbose = FALSE)
#mamsubset2 <- RunUMAP(mamsubset2, reduction = "pca", dims = 1:20)
#mamsubset2 <- FindNeighbors(mamsubset2, reduction = "pca", dims = 1:20)
##mamsubset2 <- RunTSNE(mamsubset2, reduction = "pca", dims = 1:20)

#Idents(mamsubset2) <- mamsubse$FinLab

colind <- integer( length( levels(Idents(mamsubset)) )  )
for (i in 1:length( levels(Idents(mamsubset)) ) ) {
  colind[i] <- which(cType==levels(Idents(mamsubset))[i])
}
coluse <- BaseCol[colind]

cluster_letters <- LETTERS[mamsubset$orig.ident]
cluster_letters[1:length(cluster_letters)] <- "E5"
cluster_letters[which(mamsubset$ID3=="E4")] <- "E4"
cluster_letters[which(mamsubset$ID3=="E5")] <- "E5"
cluster_letters[which(mamsubset$ID3=="E6")] <- "E6"
cluster_letters[which(mamsubset2$ID3=="E7")] <- "E7"
cluster_letters <- as.factor(cluster_letters)


#cluster_letters[grep("_F", mammal.combined$ID3)] <- "E7"



#mammal.combined <- AddMetaData(object = mammal.combined,metadata = cluster_letters,col.name = 'cell.orig')
names(cluster_letters)=rownames(mamsubset@meta.data)
mamsubset <- AddMetaData(object = mamsubset,metadata = cluster_letters,col.name = 'cell.orig')


DimPlot(mamsubset, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 2, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/UMAP_subType_Lab2",".pdf",sep=""),width = 20, height = 6)


DimPlot(mamsubset, reduction = "pca", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 2, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/PCA_subType_Lab2",".pdf",sep=""),width = 20, height = 6)




#ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab",".pdf",sep=""),width = 20, height = 6)
##DimPlot(mammal.combined, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE, repel = TRUE)
#DimPlot(mammal.combined, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE, repel = TRUE)
#ggsave(filename=paste(saveext,"/DimRed/TSNE_Type_Lab",".pdf",sep=""),width = 20, height = 6)
#DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = .5, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
#scale_shape_manual(values = c(18,16,17,15))
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab",".pdf",sep=""),width = 20, height = 6)

#DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 1, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
#  scale_shape_manual(values = c(18,16,17,15))
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_ms=1",".pdf",sep=""),width = 20, height = 6)

#DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 2, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
#  scale_shape_manual(values = c(18,16,17,15))
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab_ms=2",".pdf",sep=""),width = 20, height = 6)

Idents(mammal.combined) <- mammal.combined$ID3

ulabs <- as.character(Idents(mammal.combined))


ulabs[which(ulabs=="ICM_CS3")] <- "ICM"
ulabs[which(ulabs=="EmDisc_CS5")] <- "EmDisc_CS5"
ulabs[which(ulabs=="EmDisc_CS6/7")] <- "EmDisc_CS6/7"

ulabs[which(ulabs=="Embryo_EPI_M_D6")] <- "Epi_D6"
ulabs[which(ulabs=="Embryo_EPI_M_D8")] <- "Epi_D8"
ulabs[which(ulabs=="Embryo_EPI_M_D10")] <- "Epi_D10"
ulabs[which(ulabs=="Embryo_EPI_M_D12")] <- "Epi_D12"
ulabs[which(ulabs=="Embryo_EPI_M_D14")] <- "Epi_D14"
ulabs[which(ulabs=="Embryo_EPI_F_D6")] <- "Epi_D6"
ulabs[which(ulabs=="Embryo_EPI_F_D8")] <- "Epi_D8"
ulabs[which(ulabs=="Embryo_EPI_F_D10")] <- "Epi_D10"
ulabs[which(ulabs=="Embryo_EPI_F_D12")] <- "Epi_D12"
ulabs[which(ulabs=="Embryo_EPI_F_D14")] <- "Epi_D14"

ulabs[which(ulabs=="Embryo_TE_M_D6")] <- "Tb_D6"
ulabs[which(ulabs=="Embryo_TE_M_D8")] <- "Tb_D8"
ulabs[which(ulabs=="Embryo_TE_M_D10")] <- "Tb_D10"
ulabs[which(ulabs=="Embryo_TE_M_D12")] <- "Tb_D12"
ulabs[which(ulabs=="Embryo_TE_M_D14")] <- "Tb_D14"
ulabs[which(ulabs=="Embryo_TE_F_D6")] <- "Tb_D6"
ulabs[which(ulabs=="Embryo_TE_F_D8")] <- "Tb_D8"
ulabs[which(ulabs=="Embryo_TE_F_D10")] <- "Tb_D10"
ulabs[which(ulabs=="Embryo_TE_F_D12")] <- "Tb_D12"
ulabs[which(ulabs=="Embryo_TE_F_D14")] <- "Tb_D14"

ulabs[which(ulabs=="Embryo_PE_M_D6")] <- "Hyp_D6"
ulabs[which(ulabs=="Embryo_PE_M_D8")] <- "Hyp_D8"
ulabs[which(ulabs=="Embryo_PE_M_D10")] <- "Hyp_D10"
ulabs[which(ulabs=="Embryo_PE_M_D12")] <- "Hyp_D12"
ulabs[which(ulabs=="Embryo_PE_M_D14")] <- "Hyp_D14"
ulabs[which(ulabs=="Embryo_PE_F_D6")] <- "Hyp_D6"
ulabs[which(ulabs=="Embryo_PE_F_D8")] <- "Hyp_D8"
ulabs[which(ulabs=="Embryo_PE_F_D10")] <- "Hyp_D10"
ulabs[which(ulabs=="Embryo_PE_F_D12")] <- "Hyp_D12"
ulabs[which(ulabs=="Embryo_PE_F_D14")] <- "Hyp_D14"

ulabs[which(ulabs=="PGC_4W_M")] <- "PGC_4W"
ulabs[which(ulabs=="PGC_4W_F")] <- "PGC_4W"

ulabs[which(ulabs=="PGC_5W_M")] <- "PGC_5W"
ulabs[which(ulabs=="PGC_5W_F")] <- "PGC_5W"
ulabs[which(ulabs=="PGC_7W_M")] <- "PGC_7W"
ulabs[which(ulabs=="PGC_7W_F")] <- "PGC_7W"
ulabs[which(ulabs=="PGC_8W_M")] <- "PGC_8W"
ulabs[which(ulabs=="PGC_8W_F")] <- "PGC_8W"

ulabs[which(ulabs=="PGC_9W_M")] <- "PGC_9W"
ulabs[which(ulabs=="PGC_9W_F")] <- "PGC_9W"

ulabs[which(ulabs=="PGC_10W_M")] <- "PGC_10-12W"
ulabs[which(ulabs=="PGC_10W_F")] <- "PGC_10-12W"
ulabs[which(ulabs=="PGC_11W_M")] <- "PGC_10-12W"
ulabs[which(ulabs=="PGC_11W_F")] <- "PGC_10-12W"
ulabs[which(ulabs=="PGC_12W_M")] <- "PGC_10-12W"
ulabs[which(ulabs=="PGC_12W_F")] <- "PGC_10-12W"

ulabs[which(ulabs=="PGC_14W_M")] <- "PGC_14-20W"
ulabs[which(ulabs=="PGC_14W_F")] <- "PGC_14-20W"
ulabs[which(ulabs=="PGC_17W_M")] <- "PGC_14-20W"
ulabs[which(ulabs=="PGC_17W_F")] <- "PGC_14-20W"
ulabs[which(ulabs=="PGC_18W_M")] <- "PGC_14-20W"
ulabs[which(ulabs=="PGC_18W_F")] <- "PGC_14-20W"

ulabs[which(ulabs=="PGC_19W_M")] <- "PGC_14-20W"
ulabs[which(ulabs=="PGC_19W_F")] <- "PGC_14-20W"
ulabs[which(ulabs=="PGC_20W_M")] <- "PGC_14-20W"
ulabs[which(ulabs=="PGC_20W_F")] <- "PGC_14-20W"

ulabs[which(ulabs=="PGC_21W_M")] <- "PGC_20W+"
ulabs[which(ulabs=="PGC_21W_F")] <- "PGC_20W+"

ulabs[which(ulabs=="PGC_23W_M")] <- "PGC_20W+"
ulabs[which(ulabs=="PGC_23W_F")] <- "PGC_20W+"

ulabs[which(ulabs=="PGC_24W_M")] <- "PGC_20W+"
ulabs[which(ulabs=="PGC_24W_F")] <- "PGC_20W+"

ulabs[which(ulabs=="PGC_25W_M")] <- "PGC_20W+"
ulabs[which(ulabs=="PGC_25W_F")] <- "PGC_20W+"

ulabs[which(ulabs=="PGC_26W_M")] <- "PGC_20W+"
ulabs[which(ulabs=="PGC_26W_F")] <- "PGC_20W+"


ulabs[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl9==1)] <- "PGC"
ulabs[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl9==10)] <- "Hyp"
ulabs[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl9==14)] <- "eEpi"
ulabs[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl9==17)] <- "mEpi"
ulabs[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl9==16)] <- "mEpi"
#ulabs[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl9==22)] <- "lEpi"
ulabs[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl9==11)] <- "Tb"
ulabs[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl9==0)] <- "Tb"
ulabs[which(mammal.combined$species=="1) Rabbit" & mammal.combined$Cl9==8)] <- "ICM/Tb"




cType <- c("PGC_CS5","PGC_CS6/7","PGC_E50",
           "ICM","Epi_CS3","EmDisc_CS5","EmDisc_CS6/7",
           "Tb_D6","Tb_D8","Tb_D10","Tb_D12","Tb_D14",
           "Hyp_D6","Hyp_D8","Hyp_D10","Hyp_D12","Hyp_D14",
           "Epi_D6","Epi_D8","Epi_D10","Epi_D12","Epi_D14",
           "Tb_CS6/7","Tb_CS3", "Hyp_CS3", "VE_CS5","SYS_CS6/7",
           "E4","E5","E6","E7",
           "PGC_4W" ,"PGC_5W","PGC_7W", "PGC_8W","PGC_9W","PGC_10-12W","PGC_14-20W","PGC_20W+",     
           "PGC","Hyp","eEpi","mEpi","lEpi","Tb","ICM/Tb")

BaseCol <- c("#e31a1c","#bd0026","#800026",
             "#c6dbef","#6baed6","#2171b5","#08306b",        
             "#9e9ac8","#807dba","#6a51a3","#54278f","#3f007d",
             "#74c476","#41ab5d","#238b45","#006d2c","#00441b",
             "#6baed6","#4292c6","#2171b5","#08519c","#08306b",             
             "#3f007d","#9e9ac8","#74c476","#238b45","#00441b",
             "#ffffb2","#fecc5c","#fd8d3c","#e31a1c",
             "#ffeda0",
             "#fed976",
             "#feb24c",
             "#fd8d3c",
             "#fc4e2a",
             "#e31a1c",
             "#bd0026",
             "#800026",
             "#800026","#fecc5c","#6baed6","#2171b5","#08306b","#3f007d","#3f007d")




Idents(mammal.combined) <- as.factor(ulabs)
mammal.combined$ulabs <- ulabs
colind <- integer( length( levels(Idents(mammal.combined)) )  )
for (i in 1:length( levels(Idents(mammal.combined)) ) ) {
  colind[i] <- which(cType==levels(Idents(mammal.combined))[i])
}
coluse <- BaseCol[colind]

mammal.combined$NewIDS <- Idents(mammal.combined)


DimPlot(mammal.combined, reduction = "umap",  rule="evenodd", cols = coluse,  pt.size = .5, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2",".pdf",sep=""),width = 20, height = 6)


DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = .5, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2",".pdf",sep=""),width = 20, height = 6)


DimPlot(mammal.combined, reduction = "pca", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = .5, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab2",".pdf",sep=""),width = 20, height = 6)


DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 1, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2_ms=1",".pdf",sep=""),width = 20, height = 6)

DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 2, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2_ms=2",".pdf",sep=""),width = 20, height = 6)


DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 2, merge.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2_ms=2_div1",".pdf",sep=""),width = 7, height = 6)

Idents(mammal.combined) <- mammal.combined$Cl2
DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 2, merge.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2_ms=2_cl2",".pdf",sep=""),width = 7, height = 6)


Idents(mammal.combined) <- mammal.combined$Cl3
DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 2, merge.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2_ms=2_cl3",".pdf",sep=""),width = 7, height = 6)

Idents(mammal.combined) <- mammal.combined$Cl4
DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 2, merge.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2_ms=2_cl4",".pdf",sep=""),width = 7, height = 6)

Idents(mammal.combined) <- mammal.combined$Cl5
DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 2, merge.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2_ms=2_cl5",".pdf",sep=""),width = 7, height = 6)

Idents(mammal.combined) <- mammal.combined$Cl6
DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 2, merge.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2_ms=2_cl6",".pdf",sep=""),width = 7, height = 6)


Idents(mammal.combined) <- mammal.combined$Cl9
DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 2, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2_ms=2_cl9",".pdf",sep=""),width = 7, height = 6)


mammal.combined <- FindClusters(mammal.combined, resolution = 1.5)
mammal.combined$Cl15 <- Idents(mammal.combined)
Idents(mammal.combined) <- mammal.combined$Cl15
DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 2, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2_ms=2_cl15",".pdf",sep=""),width = 20, height = 6)






#Idents(mammal.combined) <- mammal.combined$Cl3
#DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd",  pt.size = 2, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
#  scale_shape_manual(values = c(18,16,17,15))
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2_Cl3",".pdf",sep=""),width = 20, height = 6)

#Idents(mammal.combined) <- mammal.combined$Cl5
#DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd",  pt.size = 2, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
#  scale_shape_manual(values = c(18,16,17,15))
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2_Cl5",".pdf",sep=""),width = 13, height = 6)


#Idents(mammal.combined) <- mammal.combined$Cl6
#DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd",  pt.size = 2, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
#  scale_shape_manual(values = c(18,16,17,15))
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2_Cl6",".pdf",sep=""),width = 13, height = 6)


#Idents(mammal.combined) <- mammal.combined$Cl7
#DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd",  pt.size = 2, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
#  scale_shape_manual(values = c(18,16,17,15))
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2_Cl7",".pdf",sep=""),width = 13, height = 6)

#Idents(mammal.combined) <- mammal.combined$Cl8
#DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd",  pt.size = 2, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
#  scale_shape_manual(values = c(18,16,17,15))
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2_Cl8",".pdf",sep=""),width = 13, height = 6)

#Idents(mammal.combined) <- mammal.combined$Cl9
#DimPlot(mammal.combined, reduction = "umap", shape.by = 'cell.orig', rule="evenodd",  pt.size = 2, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
#  scale_shape_manual(values = c(18,16,17,15))
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2_Cl9",".pdf",sep=""),width = 13, height = 6)



Idents(mammal.combined) <- mammal.combined$Cl9


subs_data <- subset(mammal.combined, idents = c("1","14","16","17","10","8","0"), invert = FALSE )
Idents(subs_data) <- subs_data$NewIDS
#DimPlot(subs_data, reduction = "umap", shape.by = 'cell.orig', cols = coluse, rule="evenodd",  pt.size = 2, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
#  scale_shape_manual(values = c(18,16,17,15))
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2_subs",".pdf",sep=""),width = 13, height = 6)

#Idents(subs_data) <-subs_data$Cl7

DimPlot(subs_data, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = .5, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2","s.pdf",sep=""),width = 20, height = 6)


DimPlot(subs_data, reduction = "pca", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = .5, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab2","s.pdf",sep=""),width = 20, height = 6)


subs_data <- RunPCA(subs_data, npcs = 30, verbose = FALSE)
subs_data <- RunUMAP(subs_data, reduction = "pca", dims = 1:20)
subs_data <- RunTSNE(subs_data, reduction = "pca", dims = 1:20)
subs_data <- FindNeighbors(subs_data, reduction = "pca", dims = 1:20)
subs_data <- FindClusters(subs_data, resolution = 0.5)
Idents(subs_data) <- subs_data$NewIDS
DimPlot(subs_data, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = .5, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab2","s2.pdf",sep=""),width = 20, height = 6)


DimPlot(subs_data, reduction = "pca", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = .5, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab2","s2.pdf",sep=""),width = 20, height = 6)
#subs_data$ucelltype <- paste( subs_data$species, subs_data$Cl7, sep = "_")
#Idents(subs_data) <- subs_data$ucelltype 




AraM <- intersect(  rownames(subs_data), c("POU5F1", "GATA3", "OTX2", "EOMES", "NANOG", "SOX2", "PRDM14", "CXCR4", "LHX1", "GSC", "MIXL1", "TBXT", "FOXA2", "SOX17", "TFCP2L1", "CDX2", "PDGFRB", "NANOS3", "KLF4"))



#annotation_col = data.frame(Stage = factor(colnames(b)))
#rownames(annotation_col) <- colnames(b)

# change the color of annotation to what you want: (eg: "navy", "darkgreen")
#Var1        <- c("navy", "darkgreen")
#names(mycolors) <- annotationL #c("Exp1", "Exp2")
#pheatmap(test, annotation = annotation, annotation_colors = anno_colors)

#names(mycolors) <- colnames(b)

#anno_colors <- list(Stage = mycolors)

 "Cy_14","Cy_4","Cy_7","Cy_0","Cy_10",     
 "Cy_3","Cy_8","Cy_15","Cy_11","Cy_2","Cy_5","Cy_9","Cy_6","R_4","R_0","R_8","R_6","R_3","R_7","R_10"          
 "R_15","R_14","R_9","R_2","H_11","H_10","H_2","H_13","H_3","H_0","H_8","H_5","H_1","H_4","H_9" 
 "2) H_6"  "H_15" "H_7"  "4) Cynomologous_6"      "4) Cynomologous_9"     
 "4) Cynomologous_12"     "4) Cynomologous_5"      "4) Cynomologous_8"      "4) Cynomologous_0"      "4) Cynomologous_4"     
 "4) Cynomologous_13"     "4) Cynomologous_2"      "4) Cynomologous_7"      "4) Cynomologous_11"     "4) Cynomologous_10"    
 "4) Cynomologous_3"      "4) Cynomologous_1"      "4) Cynomologous_14"     "4) Cynomologous_17"     "4) Cynomologous_16" 

AvExp <- AverageExpression(object = subs_data,features = AraM, return.seurat = TRUE)
a <- GetAssayData(object = AvExp , slot = 'data')
colnames(a) <-c("cPGC", "cHyp", "cEmDisc", "cGast", "cEpi", "cTb", "rHyp", "rEmDisc", "rTb", "rEpi", "rGast", "rPGC", "rTb2", "hPGC", "hGast", "hTb2", "hEpi", "hHyp", "hTb", "hEmDisc")
mat_breaks <- seq(0, 8, length.out = 30)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
#pheatmap(a,color =  redblue1(20),breaks = mat_breaks, border_color = NA,annotation_col = annotation_col, annotation_colors = anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,  filename = "~/Desktop/Ara.pdf")
pheatmap(as.data.frame(a),color =  redblue1(30),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE,  filename = "~/Desktop/AraMarkers.pdf")

FortyEight <- intersect(  rownames(subs_data),  c("OOEP","TCL1A","WEE2","IL1RN","NOV","ZNF80","SPIC","ESRRB","STAT3","KLF17","NLRP9","SOX15","POU5F1","NANOG","SOX2","SFRP2","DNMT3B","T","TFAP2C","TFAP2A","HOXD3","PRDM1","PRDM14","NANOS3","GATA6","GATA4","SOX17","APOA1","TTR","APOB","HAND2","TBX4","HGF","JAM2","FABP3","GATA3","GATA2","CGB3","CGA","PDK1","PDK3","CPT1A","IDH2"))
AvExp <- AverageExpression(object = subs_data,features = FortyEight, return.seurat = TRUE)
a <- GetAssayData(object = AvExp , slot = 'data')
colnames(a) <-c("cPGC", "cHyp", "cEmDisc", "cGast", "cEpi", "cTb", "rHyp", "rEmDisc", "rTb", "rEpi", "rGast", "rPGC", "rTb2", "hPGC", "hGast", "hTb2", "hEpi", "hHyp", "hTb", "hEmDisc")
mat_breaks <- seq(0, 8, length.out = 30)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
#pheatmap(a,color =  redblue1(20),breaks = mat_breaks, border_color = NA,annotation_col = annotation_col, annotation_colors = anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,  filename = "~/Desktop/Ara.pdf")
pheatmap(as.data.frame(a),color =  redblue1(30),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE,  filename = "~/Desktop/AraMarkers48.pdf")


Idents(subs_data) <-subs_data$Cl7
subs_data.markerscl <- FindAllMarkers(subs_data, only.pos = TRUE, test.use = "MAST")

top50 <- subs_data.markerscl %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
Idents(subs_data) <- subs_data$ucelltype 

AvExp <- AverageExpression(object = subs_data,features = unique(top50$gene), return.seurat = TRUE)

a <- GetAssayData(object = AvExp , slot = 'data')
colnames(a) <-c("cPGC", "cHyp", "cEmDisc", "cGast", "cEpi", "cTb", "rHyp", "rEmDisc", "rTb", "rEpi", "rGast", "rPGC", "rTb2", "hPGC", "hGast", "hTb2", "hEpi", "hHyp", "hTb", "hEmDisc")
mat_breaks <- seq(0, 7, length.out = 30)
redblue1<-colorRampPalette(c("#00B0F0","#FFFFFF","#FF0B07"))
#pheatmap(a,color =  redblue1(20),breaks = mat_breaks, border_color = NA,annotation_col = annotation_col, annotation_colors = anno_colors, cluster_rows=FALSE,cluster_cols=FALSE,  filename = "~/Desktop/Ara.pdf")
pheatmap(as.data.frame(a),color =  redblue1(30),breaks = mat_breaks, border_color = NA, cluster_rows=TRUE,cluster_cols=TRUE,  filename = "~/Desktop/AraTop50Markers.pdf", width=5, height=35)

write.csv(as.data.frame(subs_data.markerscl), file=paste(saveext,"/Markers/Markers_cl.csv",sep=""))


#mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#top100 <- mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
#DoHeatmap(mammal.combined, features = top100$gene) + NoLegend()
#ggsave(filename=paste(saveext,"/Markers/HeatMapcl_100",".pdf",sep=""),width = 5, height = 40)
#write.csv(as.data.frame(mammal.combined.markerscl), file=paste(saveext,"/Markers/Markers_cl.csv",sep=""))

library("MAST")
Idents(subs_data) <-subs_data$Cl7
subs_data.markerscl <- FindAllMarkers(subs_data, only.pos = TRUE, test.use = "MAST")

#do.return = TRUE) + 
#  scale_shape_manual(values = c("ConGroup" = 3, "ExpGroup" = 16))

#DimPlot(mammal.combined, reduction = "tsne", split.by = "cond", label = TRUE)
#ggsave(filename=paste(saveext,"TSNE_Type_Lab",".pdf",sep=""),width = 16, height = 4)


write.csv(as.data.frame(Idents(object = mammal.combined)), file=paste(saveext,"/DimRed/EmbeddingsCl.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["tsne"]])), file=paste(saveext,"/DimRed/EmbeddingsTSNE.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["umap"]])), file=paste(saveext,"/DimRed/Embeddings.csv",sep=""))
write.csv(as.data.frame(mammal.combined[[]]), file=paste(saveext,"/DimRed/EmbeddingsKey.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["pca"]])), file=paste(saveext,"/DimRed/EmbeddingsPCA.csv",sep=""))

#mammal.combined.markerscl <- FindAllMarkers(mammal.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#top100 <- mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
#DoHeatmap(mammal.combined, features = top100$gene) + NoLegend()
#ggsave(filename=paste(saveext,"/Markers/HeatMapcl_100",".pdf",sep=""),width = 5, height = 40)
#write.csv(as.data.frame(mammal.combined.markerscl), file=paste(saveext,"/Markers/Markers_cl.csv",sep=""))


#Assign cluster based on annotated cell type (read in from a seperate file)

#Idents(mammal.combined) <- factor(c( as.character(labs2),as.character(cylabs) )) #labs3

#Idents(mammal.combined) <- factor(c( as.character(labs2),as.character(labscy) )) #labs3

#And plot these
#DimPlot(mammal.combined, reduction = "pca", label.size = 2, no.legend = "none", label = TRUE, dim.1 = 20, dim.2 = 10, repel = TRUE) # + NoLegend()
#ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Labs",".pdf",sep=""),width = 8, height = 4)
#DimPlot(mammal.combined, reduction = "umap", label.size = 2, no.legend = TRUE, label = TRUE, dim.1 = 20, dim.2 = 10, repel = TRUE) #+ NoLegend()
#ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Labs",".pdf",sep=""),width = 8, height = 4)
#DimPlot(mammal.combined, reduction = "tsne", label.size = 2, no.legend = TRUE, label = TRUE, dim.1 = 20, dim.2 = 10, repel = TRUE) #+ NoLegend()
#ggsave(filename=paste(saveext,"/DimRed/TSNE_Type_Labs",".pdf",sep=""),width = 8, height = 4)






#Identify markers based on cluster
#mammal.combined.markerscl <- FindAllMarkers(mammal.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
##top100 <- mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
#mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
#DoHeatmap(mammal.combined, features = top100$gene) + NoLegend()
#ggsave(filename=paste(saveext,"/Markers/HeatMapcl_100",".pdf",sep=""),width = 5, height = 40)
#write.csv(as.data.frame(mammal.combined.markerscl), file=paste(saveext,"/Markers/Markers_cl.csv",sep=""))

#Take top 100 of each for inspection
#top100 <- mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
#DoHeatmap(mammal.combined, features = top100$gene) + NoLegend()
#ggsave(filename=paste(saveext,"/Markers/HeatMap_100",".pdf",sep=""),width = 5, height = 40)
#TFstrue <- merge(x = as.data.frame(mammal.combined.markerscl), y = as.data.frame(TF), by="row.names", all.x=TRUE)
#write.csv(as.data.frame(TFstrue), file=paste(saveext,"/Markers/Markers.csv",sep=""))

#Specific cluster comparison ...
#cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#head(cluster5.markers, n = 5)

##Take top 20 for plotting
#top20 <- mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
#write.csv(as.data.frame(top20), file=paste(saveext,"/Markers/Top20Markers.csv",sep=""))
##top20<-read.table("/Users/christopherpenfold/Desktop/Thorsten/ALL-E25specific_Modelling20k/Top20Markers.csv",sep=",",header = T, row.names=1)
#for (i in 1:length(unique(top20$cluster))) {
#  genes<-top20$gene[which(top20$cluster==unique(top20$cluster)[i])]
#  VlnPlot(mammal.combined, features = genes)
#  ggsave(filename=paste(saveext,"/Markers/Marker_",i,".pdf",sep=""),width = 24, height = 24)
#  ggsave(filename=paste(saveext,"/Markers/Markerscatter_",i,".pdf",sep=""),width = 24, height = 24)
#  FeaturePlot(mammal.combined, features = genes, combine=TRUE)
#}






Idents(mammal.combined) <- Cls
split <- SplitObject(mammal.combined, split.by = "species")
rab <- split$rabbit

#Markers, extract out the embryonal data
sub_data <- subset(rab, idents = c("15","6","3"), invert = FALSE )
#sublab <- Idents(rab)
#sub_data$species <- "rabsub"
#sub_data <- NormalizeData(sub_data, verbose = FALSE)
#sub_data <- FindVariableFeatures(sub_data, selection.method = "vst", nfeatures = 20000)
#sub_data <- ScaleData(mammal.combined, verbose = FALSE)
#mammal.combined <- ScaleData(mammal.combined,  vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(mammal.combined), verbose = FALSE)
##sub_data <- RunUMAP(sub_data, reduction = "pca", dims = 1:20)
#sub_data <- RunPCA(sub_data, npcs = 30, verbose = FALSE)
#sub_data <- RunTSNE(sub_data, reduction = "pca", dims = 1:20)
#sub_data <- FindNeighbors(sub_data, reduction = "pca", dims = 1:20)
#sub_data <- FindClusters(sub_data, resolution = 0.5)
#Idents(sub_data)<- sublab 

#Do various plots
DimPlot(sub_data, reduction = "pca", split.by = "species", label = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_subs",".pdf",sep=""),width = 8, height = 4)
DimPlot(sub_data, reduction = "umap", split.by = "species", label = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_subs",".pdf",sep=""),width = 8, height = 4)

sub_data.markerscl <- FindAllMarkers(sub_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top100 <- sub_data.markerscl %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
DoHeatmap(sub_data.markerscl, features = top100$gene) + NoLegend()
ggsave(filename=paste(saveext,"/Markers/HeatMapcl_100Rab",".pdf",sep=""),width = 5, height = 40)
write.csv(as.data.frame(sub_data.markerscl), file=paste(saveext,"/Markers/Rab_subset_Markers_cl.csv",sep=""))

#Take top 100 of each for inspection
TFstrue <- merge(x = as.data.frame(sub_data.markerscl), y = as.data.frame(TF), by="row.names", all.x=TRUE)
write.csv(as.data.frame(TFstrue), file=paste(saveext,"/Markers/RabsubMarkers.csv",sep=""))

#Take top 20 for plotting
top20 <- sub_data.markerscl %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.csv(as.data.frame(top20), file=paste(saveext,"/Markers/RabsubTop20Markers.csv",sep=""))
#top20<-read.table("/Users/christopherpenfold/Desktop/Thorsten/ALL-E25specific_Modelling20k/Top20Markers.csv",sep=",",header = T, row.names=1)
for (i in 1:length(unique(top100$cluster))) {
  genes<-top100$gene[which(top100$cluster==unique(top100$cluster)[i])]
  VlnPlot(sub_data, features = genes)
  ggsave(filename=paste(saveext,"/Markers/RabsubMarker_",unique(top100$cluster)[i],".pdf",sep=""),width = 24, height = 100,limitsize = FALSE)
  FeaturePlot(sub_data, features = genes, combine=TRUE)
  ggsave(filename=paste(saveext,"/Markers/RabsubMarkerscatter_",unique(top100$cluster)[i],".pdf",sep=""),width = 24, height = 100,limitsize = FALSE)
}


DimPlot(rab, reduction = "pca", split.by = "species", label = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_subs",".pdf",sep=""),width = 8, height = 4)
DimPlot(rab, reduction = "umap", split.by = "species", label = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_subs",".pdf",sep=""),width = 8, height = 4)

rab.markerscl <- FindAllMarkers(rab, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top100 <- rab.markerscl %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
DoHeatmap(rab.markerscl, features = top100$gene) + NoLegend()
ggsave(filename=paste(saveext,"/Markers/HeatMapcl_100Raball",".pdf",sep=""),width = 5, height = 40)
write.csv(as.data.frame(sub_data.markerscl), file=paste(saveext,"/Markers/Rab_Markers_cl.csv",sep=""))

#Take top 100 of each for inspection
TFstrue <- merge(x = as.data.frame(rab.markerscl), y = as.data.frame(TF), by="row.names", all.x=TRUE)
write.csv(as.data.frame(TFstrue), file=paste(saveext,"/Markers/RabMarkers.csv",sep=""))

#Take top 20 for plotting
top20 <- rab.markerscl %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
top100 <- rab.markerscl %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)

write.csv(as.data.frame(top20), file=paste(saveext,"/Markers/RabTop20Markers.csv",sep=""))
#top20<-read.table("/Users/christopherpenfold/Desktop/Thorsten/ALL-E25specific_Modelling20k/Top20Markers.csv",sep=",",header = T, row.names=1)
for (i in 1:length(unique(top100$cluster))) {
  genes<-top100$gene[which(top100$cluster==unique(top100$cluster)[i])]
  VlnPlot(rab, features = genes)
  ggsave(filename=paste(saveext,"/Markers/RabMarker_",unique(top100$cluster)[i],".pdf",sep=""),width = 24, height = 100,limitsize = FALSE)
  FeaturePlot(sub_data, features = genes, combine=TRUE)
  ggsave(filename=paste(saveext,"/Markers/RabMarkerscatter_",unique(top100$cluster)[i],".pdf",sep=""),width = 24, height = 100,limitsize = FALSE)
}








#FeaturePlot(mammal.combined, features = top20$gene[which(top20$cluster==unique(top20$cluster)[1])], split.by = "species")
#ggsave(filename=paste(saveext,"UMAP_FeatureCluster1_Split",".pdf",sep=""),width = 24, height = 24)
#dev.off()

#Can do the same for violin plot representations
#VlnPlot(mammal.combined, features = mammal.combined.markerscl$gene)
#ggsave(filename=paste(saveext,"Cluster_Split",".pdf",sep=""),width = 24, height = 24)
#dev.off()


#Manually compare clusters to get markers out for specific clusters rather than 1-vs-all in previous
#mammal.combined.specificmarkers1_6 <- FindMarkers(mammal.combined, ident.1 = "1", ident.2 = "6", only.pos = TRUE)
#mammal.combined.specificmarkers4_6 <- FindMarkers(mammal.combined, ident.1 = "4", ident.2 = "6", only.pos = TRUE)
#mammal.combined.specificmarkers0_6 <- FindMarkers(mammal.combined, ident.1 = "0", ident.2 = "6", only.pos = TRUE)
#write.csv(as.data.frame(mammal.combined.specificmarkers1_6), file=paste(saveext,"Markers_cl_1_6.csv",sep=""))
#write.csv(as.data.frame(mammal.combined.specificmarkers4_6), file=paste(saveext,"Markers_cl_4_6.csv",sep=""))
#write.csv(as.data.frame(mammal.combined.specificmarkers0_6), file=paste(saveext,"Markers_cl_0_6.csv",sep=""))

#...
#mammal.combined.specificmarkers6_1 <- FindMarkers(mammal.combined, ident.1 = "6", ident.2 = "1", only.pos = TRUE)
#mammal.combined.specificmarkers6_4 <- FindMarkers(mammal.combined, ident.1 = "6", ident.2 = "4", only.pos = TRUE)
#mammal.combined.specificmarkers6_0 <- FindMarkers(mammal.combined, ident.1 = "6", ident.2 = "0", only.pos = TRUE)
#write.csv(as.data.frame(mammal.combined.specificmarkers1_6), file=paste(saveext,"Markers_cl_1_6.csv",sep=""))
#write.csv(as.data.frame(mammal.combined.specificmarkers4_6), file=paste(saveext,"Markers_cl_4_6.csv",sep=""))
#write.csv(as.data.frame(mammal.combined.specificmarkers0_6), file=paste(saveext,"Markers_cl_0_6.csv",sep=""))

#Assign cluster based on annotated cell type (read in from a seperate file)
#Cls <- Idents(mammal.combined)
#Idents(mammal.combined) <- factor(c( as.character(labs2), as.character(cylabs),  as.character(labsubs1) ,  as.character(label2) )) #c(labs2,labscy) #labs2


#And plot these

DimPlot(mammal.combined, reduction = "pca", label.size = 2, no.legend = "none", label = TRUE, dim.1 = 10, dim.2 = 10) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Type",".pdf",sep=""),width = 10, height = 4)
DimPlot(mammal.combined, reduction = "umap", label.size = 2, no.legend = TRUE, label = TRUE, dim.1 = 10, dim.2 = 10) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type",".pdf",sep=""),width = 10, height = 4)
DimPlot(mammal.combined, reduction = "tsne", label.size = 2, no.legend = TRUE, label = TRUE, dim.1 = 10, dim.2 = 10) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/TSNE_Type",".pdf",sep=""),width = 10, height = 4)

DimPlot(mammal.combined, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Lab",".pdf",sep=""),width = 20, height = 6)
DimPlot(mammal.combined, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_Type_Lab",".pdf",sep=""),width = 20, height = 6)
DimPlot(mammal.combined, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Lab",".pdf",sep=""),width = 20, height = 6)

#DimPlot(mammal.combined, reduction = "tsne", split.by = "cond", label = TRUE)
#ggsave(filename=paste(saveext,"TSNE_Type_Lab",".pdf",sep=""),width = 16, height = 4)



#Assign cluster based on annotated cell type (read in from a seperate file)

#Idents(mammal.combined) <- factor(c( as.character(labs2),as.character(cylabs) )) #labs3

#Idents(mammal.combined) <- factor(c( as.character(labs2),as.character(labscy) )) #labs3

#And plot these
DimPlot(mammal.combined, reduction = "pca", label.size = 2, no.legend = "none", label = TRUE, dim.1 = 20, dim.2 = 10, repel = TRUE) + NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_Type_Labs",".pdf",sep=""),width = 8, height = 4)
DimPlot(mammal.combined, reduction = "umap", label.size = 2, no.legend = TRUE, label = TRUE, dim.1 = 20, dim.2 = 10, repel = TRUE) + NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_Type_Labs",".pdf",sep=""),width = 8, height = 4)
DimPlot(mammal.combined, reduction = "tsne", label.size = 2, no.legend = TRUE, label = TRUE, dim.1 = 20, dim.2 = 10, repel = TRUE) + NoLegend()
ggsave(filename=paste(saveext,"/DimRed/TSNE_Type_Labs",".pdf",sep=""),width = 8, height = 4)

#And plot these
#DimPlot(mammal.combined, reduction = "pca", label.size = 2, no.legend = "none", label = TRUE, dim.1 = 20, dim.2 = 10) + NoLegend()
#ggsave(filename=paste(saveext,"PCA_Type_Labs",".pdf",sep=""),width = 8, height = 4)
#DimPlot(mammal.combined, reduction = "umap", label.size = 2, no.legend = TRUE, label = TRUE, dim.1 = 20, dim.2 = 10) + NoLegend()
#ggsave(filename=paste(saveext,"UMAP_Type_Labs",".pdf",sep=""),width = 8, height = 4)
#DimPlot(mammal.combined, reduction = "tsne", label.size = 2, no.legend = TRUE, label = TRUE, dim.1 = 20, dim.2 = 10) + NoLegend()
#ggsave(filename=paste(saveext,"TSNE_Type_Labs",".pdf",sep=""),width = 8, height = 4)



AvExp <- AverageExpression(object = mammal.combined, use.counts = TRUE)
AvExpDF <- as.data.frame(AvExp, row.names = rownames(AvExp$RNA))
AvExpDFTF <- merge(x = AvExpDF, y=as.data.frame(TF), by="row.names", all.x=TRUE)
write.csv(as.data.frame(AvExpDFTF), file=paste(saveext,"AvExp_CPM.csv",sep=""))

AvExp <- AverageExpression(object = mammal.combined, use.scale = TRUE)
AvExpDF <- as.data.frame(AvExp, row.names = rownames(AvExp$RNA))
AvExpDFTF <- merge(x = AvExpDF, y=as.data.frame(TF), by="row.names", all.x=TRUE)
write.csv(as.data.frame(AvExpDFTF), file=paste(saveext,"AvExp_Scaled.csv",sep=""))


Idents(mammal.combined) <- Cls


#Write everythign out to csv
write.csv(as.data.frame(Idents(object = mammal.combined)), file=paste(saveext,"/DimRed/EmbeddingsClr.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["tsne"]])), file=paste(saveext,"/DimRed/EmbeddingsTSNEr.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["umap"]])), file=paste(saveext,"/DimRed/Embeddingsr.csv",sep=""))
write.csv(as.data.frame(mammal.combined[[]]), file=paste(saveext,"/DimRed/EmbeddingsKey.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["pca"]])), file=paste(saveext,"/DimRed/EmbeddingsPCAr.csv",sep=""))

Idents(mammal.combined) <- factor(c( as.character(cylabs),  as.character(RabbitBS$TimePoint), as.character(le3[,3]), as.character(le4[,3]), as.character(le5[,4]), as.character(le6[,4]) )) #labs3

#Idents(mammal.combined) <- factor(c( as.character(labs2), as.character(cylabs),  as.character(labsubs1) ,  as.character(label2) )) #c(labs2,labscy) #labs2

write.csv(as.data.frame(Idents(object = mammal.combined)), file=paste(saveext,"/DimRed/EmbeddingsCl_labs.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["tsne"]])), file=paste(saveext,"/DimRed/EmbeddingsTSNE_labs.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["umap"]])), file=paste(saveext,"/DimRed/Embeddings_labs.csv",sep=""))
write.csv(as.data.frame(mammal.combined[[]]), file=paste(saveext,"/DimRed/EmbeddingsKey_labs.csv",sep=""))
write.csv(as.data.frame(Embeddings(object = mammal.combined[["pca"]])), file=paste(saveext,"/DimRed/EmbeddingsPCA_labs.csv",sep=""))


#We can now manually plot expresison in a bunch of markers
FeaturePlot(mammal.combined, features = "NANOG", split.by = "species", max.cutoff = 100, cols = c("grey", "red"))
ggsave(filename=paste(saveext,"/Markers/UMAP_NANOG_Split",".pdf",sep=""),width = 4, height = 4)

FeaturePlot(mammal.combined, features = "SP2", split.by = "species", max.cutoff = 100, cols = c("grey", "red"))
ggsave(filename=paste(saveext,"/Markers/UMAP_SP2_Split",".pdf",sep=""),width = 4, height = 4)

FeaturePlot(mammal.combined, features = "TUBB4A", split.by = "species", max.cutoff = 100, cols = c("grey", "red"))
ggsave(filename=paste(saveext,"/Markers/UMAP_TUBB4A_Split",".pdf",sep=""),width = 4, height = 4)

#Or we can plot by groups of makers
epi <- c("POU5F1","SOX2","NANOG","SALL2","SOX11","SFRP1","SFRP2","PRDM14","NODAL","TDGF1","NLRP7","DDX43")
gast <- c("EOMES","SNAI1","SNAI2","MIXL1","T","SOX17")
troph<-c("TFAP2C","CGA","GATA2","GATA3","CDX2","GCM1","HAND1","DAB2","FGFR1","KRT7")
endo<- c("TDGF1","NODAL","OTX2","RSPO3","GATA6","SOX17","PDGFRA","CER1","IGF1","GATA4","HNF4A")
EXMC<-c("TBX20","WNT5A","BMP4","HGF","HAND2","FOXF1")
amnion<-c("SOX2","NANOG","POU5F1","GCM1","TFAP2C","TFAP2A","VTCN1","GATA2","HAND1")
PGC<-c("PRDM1","TFAP2C","T","SOX17","NANOG","NANOS3","PRDM14")

VE<-c("SOX7","SOX17","MIXL1","WNT3A")

#features <- c("NFYA", "SOX2", "NANOG", "KLF4", "OTX2", "NR5A2")
#FeaturePlot(mammal.combined, features = features, split.by = "species")
#ggsave(filename=paste(saveext,"/Markers/UMAP_ZGA_Split",".pdf",sep=""),width = 4, height = 24)

FeaturePlot(mammal.combined, features = epi, split.by = "species")
ggsave(filename=paste(saveext,"/Markers/UMAP_EPI_Split",".pdf",sep=""),width = 4, height = 24)

FeaturePlot(mammal.combined, features = gast, split.by = "species")
ggsave(filename=paste(saveext,"/Markers/UMAP_GAST_Split",".pdf",sep=""),width = 4, height = 15)

FeaturePlot(mammal.combined, features = troph, split.by = "species")
ggsave(filename=paste(saveext,"/Markers/UMAP_TROPH_Split",".pdf",sep=""),width = 4, height = 24)

FeaturePlot(mammal.combined, features = EXMC, split.by = "species")
ggsave(filename=paste(saveext,"/Markers/UMAP_EXMC_Split",".pdf",sep=""),width = 4, height = 16)

FeaturePlot(mammal.combined, features = PGC, split.by = "species")
ggsave(filename=paste(saveext,"/Markers/UMAP_PGC_Split",".pdf",sep=""),width = 4, height = 24)

FeaturePlot(mammal.combined, features = endo, split.by = "species")
ggsave(filename=paste(saveext,"/Markers/UMAP_ENDO_Split",".pdf",sep=""),width = 4, height = 24)


FeaturePlot(mammal.combined, features = VE, split.by = "species")
ggsave(filename=paste(saveext,"/Markers/UMAP_VE_Split",".pdf",sep=""),width = 4, height = 24)


#Can do the same for violin plot representations
VlnPlot(mammal.combined, features = features)
ggsave(filename=paste(saveext,"/Markers/Cluster_ZGA_Split",".pdf",sep=""),width = 24, height = 24)

VlnPlot(mammal.combined, features = epi)
ggsave(filename=paste(saveext,"/Markers/Cluster_EPI_Split",".pdf",sep=""),width = 24, height = 24)

VlnPlot(mammal.combined, features = gast)
ggsave(filename=paste(saveext,"/Markers/Cluster_GAST_Split",".pdf",sep=""),width = 24, height = 24)

VlnPlot(mammal.combined, features = troph)
ggsave(filename=paste(saveext,"/Markers/Cluster_TROPH_Split",".pdf",sep=""),width = 24, height = 24)

VlnPlot(mammal.combined, features = endo)
ggsave(filename=paste(saveext,"/Markers/Cluster_ENDO_Split",".pdf",sep=""),width = 24, height = 24)

VlnPlot(mammal.combined, features = PGC)
ggsave(filename=paste(saveext,"/Markers/Cluster_PGC_Split",".pdf",sep=""),width = 24, height = 24)

#Identify markers in an unbiased way based on cell type
#mammal.combined.markers <- FindAllMarkers(mammal.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#mammal.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#Take top 100 of each for inspection
#top100 <- mammal.combined.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
#DoHeatmap(mammal.combined, features = top100$gene) + NoLegend()
#ggsave(filename=paste(saveext,"HeatMap_100",".pdf",sep=""),width = 5, height = 40)
#write.csv(as.data.frame(mammal.combined.markers), file=paste(saveext,"Markers.csv",sep=""))

#Take top 20 for plotting
#top20 <- mammal.combined.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
#write.csv(as.data.frame(top20), file=paste(saveext,"Top20Markers.csv",sep=""))
#top20<-read.table("/Users/christopherpenfold/Desktop/Thorsten/ALL-E25specific_Modelling20k/Top20Markers.csv",sep=",",header = T, row.names=1)
#for (i in 1:length(unique(top20$cluster))) {
#  genes<-top20$gene[which(top20$cluster==unique(top20$cluster)[i])]
#  VlnPlot(mammal.combined, features = genes)
#  ggsave(filename=paste(saveext,unique(top20$cluster)[i],"_marker_",i,".pdf",sep=""),width = 24, height = 24)
#  FeaturePlot(mammal.combined, features = genes, combine=TRUE)
#  ggsave(filename=paste(saveext,unique(top20$cluster)[i],"_markerscatter_",i,".pdf",sep=""),width = 24, height = 24)
#}

FeaturePlot(mammal.combined, features = c(top20$gene[which(top20$cluster==unique(top20$cluster)[1])][1],top20$gene[which(top20$cluster==unique(top20$cluster)[1])][2]), split.by = "species", max.cutoff = c(3,3),blend = TRUE,blend.threshold = .2)
ggsave(filename=paste(saveext,"/Markers/UMAP_top2markers_Cl1Split",".pdf",sep=""),width = 20, height = 5)

FeaturePlot(mammal.combined, features = c(top20$gene[which(top20$cluster==unique(top20$cluster)[2])][1],top20$gene[which(top20$cluster==unique(top20$cluster)[2])][2]), split.by = "species", max.cutoff = c(3,3),blend = TRUE,blend.threshold = .2)
ggsave(filename=paste(saveext,"/Markers/UMAP_top2markers_Cl2Split",".pdf",sep=""),width = 20, height = 5)

FeaturePlot(mammal.combined, features = c("CD52","CA12"), split.by = "species", max.cutoff = c(3,3),blend = TRUE,blend.threshold = .2)
ggsave(filename=paste(saveext,"/Markers/UMAP_top2markers_CD52_CA12_Split",".pdf",sep=""),width = 20, height = 5)

FeaturePlot(mammal.combined, features = c("ITGB6","PTGS2"), split.by = "species", max.cutoff = c(3,3),blend = TRUE,blend.threshold = .2)
ggsave(filename=paste(saveext,"/Markers/UMAP_top2markers_ITGB6_PTGS2_Split",".pdf",sep=""),width = 20, height = 5)

FeaturePlot(mammal.combined, features = c("CD47","PLBD2"), split.by = "species", max.cutoff = c(3,3),blend = TRUE,blend.threshold = .2)
ggsave(filename=paste(saveext,"/Markers/UMAP_top2markers_CD47_PLBD2_Split",".pdf",sep=""),width = 20, height = 5)

FeaturePlot(mammal.combined, features = c("CD47","ITGB6"), split.by = "species", max.cutoff = c(3,3),blend = TRUE,blend.threshold = .2)
ggsave(filename=paste(saveext,"/Markers/UMAP_top2markers_CD47_ITGB6_Split",".pdf",sep=""),width = 20, height = 5)

FeaturePlot(mammal.combined, features = c("TFAP2C","SOX17"), split.by = "species",  max.cutoff = c(3,2),blend = TRUE, blend.threshold = .1)
ggsave(filename=paste(saveext,"/Markers/UMAP_Cl1_Top2_Split",".pdf",sep=""),width = 20, height = 5)
FeaturePlot(mammal.combined, features = c("TFAP2C","SOX17"), split.by = "species",  max.cutoff = c(3,2),blend = TRUE, blend.threshold = .1)
ggsave(filename=paste(saveext,"/Markers/UMAP_Cl2_Top2_Split",".pdf",sep=""),width = 20, height = 5)

#Plot some examples of double markers, will be good to pick out
#FeaturePlot(mammal.combined, features = c("TFAP2A","SOX17"), split.by = "species", max.cutoff = c(3,2),blend = TRUE,blend.threshold = .1)
#ggsave(filename=paste(saveext,"UMAP_TFAP2A_SOX17_Split",".pdf",sep=""),width = 20, height = 5)

##FeaturePlot(mammal.combined, features = c("TFAP2A","SOX17"), split.by = "species", max.cutoff = c(3,2),blend = TRUE,blend.threshold = .1)
#ggsave(filename=paste(saveext,"UMAP_TFAP2A_SOX17_Split",".pdf",sep=""),width = 20, height = 5)

##This snippet of code will colour the plots by markers identified in another context e.g., in the clean E25A we can call markers and
##then colour those markers in another dataset. Block out for now
#top100<-read.table("/Users/christopherpenfold/Desktop/Thorsten/ALL-E25specific_Modelling20k/Markers.csv",sep=",",header = T, row.names=1) 
#mat1 <- as.matrix(GetAssayData(mammal.combined, slot = "scale.data"))
#N_1 <- matrix(0, nrow = dim(mat1)[2], ncol = length(unique(top100$cluster)))
#N_2 <- matrix(0, nrow = dim(mat1)[2], ncol = length(unique(top100$cluster)))
#N_3 <- matrix(0, nrow = dim(mat1)[2], ncol = length(unique(top100$cluster)))
#N_4 <- matrix(0, nrow = dim(mat1)[2], ncol = 1)
#for (j in 1:length(unique(top100$cluster))  ){
#  mg<-top100$gene[which(top100$cluster==unique(top100$cluster)[j])
#  for (i in 1:dim(mat1)[2]){
#    N_1[i,j] <- length(intersect(rownames(mat1)[which(mat1[,i]>1)],mg))
#    N_2[i,j] <- length((rownames(mat1)[which(mat1[,i]>1)]))      
#    #ps <- phyper(N_1[i,j]-1, N_2[i,j], dim(mat1)[1]-N_2[i,j], length(mg), lower.tail= FALSE)
#    N_3[i,j] <- length(intersect(rownames(mat1)[which(mat1[,i]>-1000)],mg))
#  }
#}

#for (j in 1:dim(N_1)[2]) {
#  vec2add <- N_1[,j] / N_3[,j] #pvalss[,6] #N_5 #1-padj
#  names(vec2add) <- colnames(x = mammal.combined)
#  mammal.combined <- AddMetaData(object = mammal.combined,metadata = vec2add, col.name = 'markerno')
#  FeaturePlot(mammal.combined, features = 'markerno', cols = c("grey", "red"))
#  ggsave(filename=paste(saveext,"/OtherMarkers/",unique(top100$cluster)[j], "markerscatter_Enrich_",j,".pdf",sep=""))
#}



Idents(mammal.combined) <- Cls



mammal.combined.markerscl <- FindAllMarkers(mammal.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top100 <- mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
DoHeatmap(mammal.combined, features = top100$gene) + NoLegend()
ggsave(filename=paste(saveext,"/Markers/HeatMapcl_100",".pdf",sep=""),width = 5, height = 40)
write.csv(as.data.frame(mammal.combined.markerscl), file=paste(saveext,"/Markers/Markers_cl.csv",sep=""))

#Take top 100 of each for inspection
TFstrue <- merge(x = as.data.frame(mammal.combined.markerscl), y = as.data.frame(TF), by="row.names", all.x=TRUE)
write.csv(as.data.frame(TFstrue), file=paste(saveext,"/Markers/Markers.csv",sep=""))

#Take top 20 for plotting
top20 <- mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.csv(as.data.frame(top20), file=paste(saveext,"/Markers/Top20Markers.csv",sep=""))
#top20<-read.table("/Users/christopherpenfold/Desktop/Thorsten/ALL-E25specific_Modelling20k/Top20Markers.csv",sep=",",header = T, row.names=1)
for (i in 1:length(unique(top20$cluster))) {
  genes<-top20$gene[which(top20$cluster==unique(top20$cluster)[i])]
  VlnPlot(mammal.combined, features = genes)
  ggsave(filename=paste(saveext,"/Markers/Marker_",unique(top20$cluster)[i],".pdf",sep=""),width = 24, height = 24)
  FeaturePlot(mammal.combined, features = genes, combine=TRUE)
  ggsave(filename=paste(saveext,"/Markers/Markerscatter_",unique(top20$cluster)[i],".pdf",sep=""),width = 24, height = 24)
  
}





split <- SplitObject(mammal.combined, split.by = "species")

rab <- split$rabbit


rab.markerscl <- FindAllMarkers(rab, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top100 <- rab.markerscl %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
DoHeatmap(rab.markerscl, features = top100$gene) + NoLegend()
ggsave(filename=paste(saveext,"/Markers/HeatMapcl_100Rab",".pdf",sep=""),width = 5, height = 40)
write.csv(as.data.frame(rab.markerscl), file=paste(saveext,"/Markers/RabMarkers_cl.csv",sep=""))

#Take top 100 of each for inspection
TFstrue <- merge(x = as.data.frame(rab.markerscl), y = as.data.frame(TF), by="row.names", all.x=TRUE)
write.csv(as.data.frame(TFstrue), file=paste(saveext,"/Markers/RabMarkers.csv",sep=""))

#Take top 20 for plotting
top20 <- rab.markerscl %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.csv(as.data.frame(top20), file=paste(saveext,"/Markers/RabTop20Markers.csv",sep=""))
#top20<-read.table("/Users/christopherpenfold/Desktop/Thorsten/ALL-E25specific_Modelling20k/Top20Markers.csv",sep=",",header = T, row.names=1)
for (i in 1:length(unique(top20$cluster))) {
  genes<-top20$gene[which(top20$cluster==unique(top20$cluster)[i])]
  VlnPlot(rab, features = genes)
  ggsave(filename=paste(saveext,"/Markers/RabMarker_",unique(top20$cluster)[i],".pdf",sep=""),width = 24, height = 24)
  FeaturePlot(rab, features = genes, combine=TRUE)
  ggsave(filename=paste(saveext,"/Markers/RabMarkerscatter_",unique(top20$cluster)[i],".pdf",sep=""),width = 24, height = 24)
  
}




Idents(mammal.combined) <- Cls
split <- SplitObject(mammal.combined, split.by = "species")
rab <- split$rabbit



#Markers, extract out the embryonal data
sub_data <- subset(rab, idents = c("6","14","3"), invert = FALSE )
sublab <- Idents(rab)
sub_data$species <- "rabsub"
#sub_data <- NormalizeData(sub_data, verbose = FALSE)
#sub_data <- FindVariableFeatures(sub_data, selection.method = "vst", nfeatures = 20000)
#sub_data <- ScaleData(mammal.combined, verbose = FALSE)
#mammal.combined <- ScaleData(mammal.combined,  vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(mammal.combined), verbose = FALSE)
sub_data <- RunPCA(sub_data, npcs = 30, verbose = FALSE)
sub_data <- RunUMAP(sub_data, reduction = "pca", dims = 1:20)
sub_data <- RunTSNE(sub_data, reduction = "pca", dims = 1:20)
sub_data <- FindNeighbors(sub_data, reduction = "pca", dims = 1:20)
sub_data <- FindClusters(sub_data, resolution = 0.5)
#Idents(sub_data)<- sublab 


#Do various plots
DimPlot(sub_data, reduction = "pca", split.by = "species", label = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_subs",".pdf",sep=""),width = 8, height = 4)
DimPlot(sub_data, reduction = "umap", split.by = "species", label = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_subs",".pdf",sep=""),width = 8, height = 4)





sub_data.markerscl <- FindAllMarkers(sub_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#mammal.combined.markerscl %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top100 <- sub_data.markerscl %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
DoHeatmap(sub_data.markerscl, features = top100$gene) + NoLegend()
ggsave(filename=paste(saveext,"/Markers/HeatMapcl_100Rab",".pdf",sep=""),width = 5, height = 40)
write.csv(as.data.frame(sub_data.markerscl), file=paste(saveext,"/Markers/Rab_subset_Markers_cl.csv",sep=""))

#Take top 100 of each for inspection
TFstrue <- merge(x = as.data.frame(rab.markerscl), y = as.data.frame(TF), by="row.names", all.x=TRUE)
write.csv(as.data.frame(TFstrue), file=paste(saveext,"/Markers/RabMarkers.csv",sep=""))

#Take top 20 for plotting
top20 <- sub_data.markerscl %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.csv(as.data.frame(top20), file=paste(saveext,"/Markers/RabsubTop20Markers.csv",sep=""))
#top20<-read.table("/Users/christopherpenfold/Desktop/Thorsten/ALL-E25specific_Modelling20k/Top20Markers.csv",sep=",",header = T, row.names=1)
for (i in 1:length(unique(top100$cluster))) {
  genes<-top100$gene[which(top100$cluster==unique(top100$cluster)[i])]
  VlnPlot(sub_data, features = genes)
  ggsave(filename=paste(saveext,"/Markers/RabsubMarker_",unique(top100$cluster)[i],".pdf",sep=""),width = 24, height = 100,limitsize = FALSE)
  FeaturePlot(sub_data, features = genes, combine=TRUE)
  ggsave(filename=paste(saveext,"/Markers/RabsubMarkerscatter_",unique(top100$cluster)[i],".pdf",sep=""),width = 24, height = 100,limitsize = FALSE)
}



Idents(mammal.combined) <- factor(c( as.character(cylabs),  as.character(RabbitBS$TimePoint), as.character(le3[,2]), as.character(le4[,2]), as.character(le5[,2]), as.character(le6[,2]) )) #labs3


split2 <- SplitObject(mammal.combined, split.by = "species")
rab <- split2$rabbit
rab4 <- subset(rab, idents = c("E4"), invert = FALSE)
rab5 <- subset(rab, idents = c("E5"), invert = FALSE)
rab6 <- subset(rab, idents = c("E6"), invert = FALSE)
rab7 <- subset(rab, idents = c("E7"), invert = FALSE)


plist <- DimPlot(rab4, reduction = "umap", split.by = "species", label = TRUE)
plist + xlim(c(-13,5)) + ylim(c(-3,13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_rab4",".pdf",sep=""),width = 8, height = 8)

plist <- DimPlot(rab5, reduction = "umap", split.by = "species", label = TRUE)
plist + xlim(c(-13,5)) + ylim(c(-3,13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_rab5",".pdf",sep=""),width = 8, height = 8)


plist <- DimPlot(rab6, reduction = "umap", split.by = "species", label = TRUE)
plist + xlim(c(-13,5)) + ylim(c(-3,13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_rab6",".pdf",sep=""),width = 8, height = 8)


plist <- DimPlot(rab7, reduction = "umap", split.by = "species", label = TRUE)
plist + xlim(c(-13,5)) + ylim(c(-3,13))
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_rab7",".pdf",sep=""),width = 8, height = 8)






Idents(mammal.combined) <- factor(c( as.character(cylabs),  as.character(RabbitBS$TimePoint), as.character(le3[,4]), as.character(le4[,4]), as.character(le5[,5]), as.character(le6[,5]) )) #labs3
split2 <- SplitObject(mammal.combined, split.by = "species")
rab <- split2$human
rab4 <- subset(rab, idents = c("D6"), invert = FALSE)
rab5 <- subset(rab, idents = c("D8"), invert = FALSE)
rab6 <- subset(rab, idents = c("D10"), invert = FALSE)
rab7 <- subset(rab, idents = c("D12"), invert = FALSE)
rab8 <- subset(rab, idents = c("D14"), invert = FALSE)


plist <- DimPlot(rab4, reduction = "umap", split.by = "species", label = TRUE)
plist + xlim(c(-10,10)) + ylim(c(-10,6))
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_hum6",".pdf",sep=""),width = 8, height = 8)

plist <- DimPlot(rab5, reduction = "umap", split.by = "species", label = TRUE)
plist + xlim(c(-10,10)) + ylim(c(-10,6))
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_hum8",".pdf",sep=""),width = 8, height = 8)

plist <- DimPlot(rab6, reduction = "umap", split.by = "species", label = TRUE)
plist + xlim(c(-10,10)) + ylim(c(-10,6))
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_hum10",".pdf",sep=""),width = 8, height = 8)

plist <- DimPlot(rab7, reduction = "umap", split.by = "species", label = TRUE)
plist + xlim(c(-10,10)) + ylim(c(-10,6))
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_hum12",".pdf",sep=""),width = 8, height = 8)

plist <- DimPlot(rab8, reduction = "umap", split.by = "species", label = TRUE)
plist + xlim(c(-10,10)) + ylim(c(-10,6))
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_hum14",".pdf",sep=""),width = 8, height = 8)




Idents(mammal.combined) <- factor(c( as.character(cylabs),  as.character(RabbitBS$TimePoint), as.character(le3[,2]), as.character(le4[,2]), as.character(le5[,2]), as.character(le6[,2]) )) #labs3
split2 <- SplitObject(mammal.combined, split.by = "species")
rab <- split2$human
rab4 <- subset(rab, idents = c("4W","5W"), invert = FALSE)
rab5 <- subset(rab, idents = c("7W","8W"), invert = FALSE)
rab6 <- subset(rab, idents = c("10W","11W"), invert = FALSE)
rab7 <- subset(rab, idents = c("12W","14W"), invert = FALSE)
rab8 <- subset(rab, idents = c("18W","19W"), invert = FALSE)
rab9 <- subset(rab, idents = c("20W","23W","24W"), invert = FALSE)

plist <- DimPlot(rab4, reduction = "umap", split.by = "species", label = TRUE)
plist + xlim(c(-10,10)) + ylim(c(4,20))
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_humpgc4",".pdf",sep=""),width = 8, height = 8)

plist <- DimPlot(rab5, reduction = "umap", split.by = "species", label = TRUE)
plist + xlim(c(-10,10)) + ylim(c(4,20))
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_humpgc5",".pdf",sep=""),width = 8, height = 8)

plist <- DimPlot(rab6, reduction = "umap", split.by = "species", label = TRUE)
plist + xlim(c(-10,10)) + ylim(c(4,20))
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_humpgc6",".pdf",sep=""),width = 8, height = 8)

plist <- DimPlot(rab7, reduction = "umap", split.by = "species", label = TRUE)
plist + xlim(c(-10,10)) + ylim(c(4,20))
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_humpgc7",".pdf",sep=""),width = 8, height = 8)

plist <- DimPlot(rab8, reduction = "umap", split.by = "species", label = TRUE)
plist + xlim(c(-10,10)) + ylim(c(4,20))
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_humpgc8",".pdf",sep=""),width = 8, height = 8)

plist <- DimPlot(rab9, reduction = "umap", split.by = "species", label = TRUE)
plist + xlim(c(-10,10))+ ylim(c(4,20))
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_humpgc9",".pdf",sep=""),width = 8, height = 8)
