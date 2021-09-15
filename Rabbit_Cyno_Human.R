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
saveext = "/Rabbit_Analysis1/"
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

oldrabbit <- readRDS(paste(saveext,"/OldAnnotations.rds",sep=""))

#Read second batch of data
rabbit_data <- readRDS("Data/Rabbit_SCE_1.rds")
rabbit_data2 <- readRDS("Data/Rabbit_SCE_2.rds")
#Human
human_dataA1 <- readRDS("Data/Human_SCE.rds")
#Cyno
cynomolgous_data <- readRDS("Data/Cyno_SCE.rds")

TF<-read.table("Data/Human_TF_MasterList_v1_02.csv",sep=",",header = F, row.names=2)

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

#First just look at second dataset
nrabbit_data3 <- rabbit_data2 #merge(rabbit_data, y = c(rabbit_data2), project = "merged")
nrabbit_data3 <- FindVariableFeatures(nrabbit_data3, selection.method = "vst", nfeatures = 3000)
nrabbit_data3 <- ScaleData(nrabbit_data3, verbose = FALSE)
nrabbit_data3 <- RunPCA(nrabbit_data3, npcs = 30, verbose = FALSE)
nrabbit_data3 <- RunUMAP(nrabbit_data3, reduction = "pca", dims = 1:20)
nrabbit_data3 <- RunTSNE(nrabbit_data3, reduction = "pca", dims = 1:20)
nrabbit_data3 <- FindNeighbors(nrabbit_data3, reduction = "pca", dims = 1:20)

DimPlot(nrabbit_data3, reduction = "tsne",  label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_nrab3",".pdf",sep=""),width = 8, height = 8)
DimPlot(nrabbit_data3, reduction = "umap",  label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_nrab3",".pdf",sep=""),width = 8, height = 8)
DimPlot(nrabbit_data3, reduction = "pca",  label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_split_nrab3",".pdf",sep=""),width = 8, height = 8)

nrabbit_data3 <- rabbit_data2 #merge(rabbit_data, y = c(rabbit_data2), project = "merged")
nrabbit_data3 <- FindVariableFeatures(nrabbit_data3, selection.method = "vst", nfeatures = 20000)
nrabbit_data3 <- ScaleData(nrabbit_data3, verbose = FALSE)
nrabbit_data3 <- RunPCA(nrabbit_data3, npcs = 30, verbose = FALSE)
nrabbit_data3 <- RunUMAP(nrabbit_data3, reduction = "pca", dims = 1:20)
nrabbit_data3 <- RunTSNE(nrabbit_data3, reduction = "pca", dims = 1:20)
nrabbit_data3 <- FindNeighbors(nrabbit_data3, reduction = "pca", dims = 1:20)

DimPlot(nrabbit_data3, reduction = "tsne",  label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_nrab4",".pdf",sep=""),width = 8, height = 8)
DimPlot(nrabbit_data3, reduction = "umap",  label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_nrab4",".pdf",sep=""),width = 8, height = 8)
DimPlot(nrabbit_data3, reduction = "pca",  label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_split_nrab4",".pdf",sep=""),width = 8, height = 8)


nrabbit_data3 <- FindClusters(nrabbit_data3, resolution = 1.5)
DimPlot(nrabbit_data3, reduction = "tsne",  label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_rab3","_cl.pdf",sep=""),width = 8, height = 8)
DimPlot(nrabbit_data3, reduction = "umap",  label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_rab3","_cl.pdf",sep=""),width = 8, height = 8)
DimPlot(nrabbit_data3, reduction = "pca",  label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_split_rab3","_cl.pdf",sep=""),width = 8, height = 8)


library("pheatmap")
avexp  <- AverageExpression(object = nrabbit_data3, slot = "data")

a <- avexp$RNA
a <- a[Markers,c("1","0","2","3")]

intgenes <- rownames(avexp$RNA)
b <- avexp$RNA[intgenes,c("1","0","2","3")]

C1 <- cor(b, method = "pearson")

redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))

mat_breaks <- seq(0, 2, length.out = 20)
pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks,gaps_row=c(6,11,16,21,29), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE,  filename = paste(saveext,"/DimRed/HM_ara_logged",".pdf",sep=""),width=8,height=16)

mat_breaks <- seq(-1, 1, length.out = 20)
pheatmap(log2(a+1),color =  redblue1(20),breaks = mat_breaks,gaps_row=c(6,11,16,21,29), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, scale="row", filename = paste(saveext,"/DimRed/HM_ara_logged_scaled",".pdf",sep=""),width=8,height=16)

mat_breaks <- seq(0.8, 0.92, length.out = 20)
pheatmap(C1[c("1","0","2","3"),c("1","0","2","3")],color =  redblue1(20), border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/Sample_corr",".pdf",sep=""),width=4,height=4)



FeaturePlot(nrabbit_data3, features = "SOX17", cols =  c("lightgrey", "darkblue"), pt.size = 5)
ggsave(filename=paste(saveext,"/Markers/Marker_full_SOX17.pdf",sep=""),width = 10, height = 10,limitsize = FALSE)
FeaturePlot(nrabbit_data3, features = "NANOS3", cols =  c("lightgrey", "darkblue"), pt.size = 5)
ggsave(filename=paste(saveext,"/Markers/Marker_full_NANOS3.pdf",sep=""),width = 10, height = 10,limitsize = FALSE)
FeaturePlot(nrabbit_data3, features = "SOX2", cols =  c("lightgrey", "darkblue"), pt.size = 5)
ggsave(filename=paste(saveext,"/Markers/Marker_full_SOX2.pdf",sep=""),width = 10, height = 10,limitsize = FALSE)

#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data, rabbit_data2, human_dataA1 ), dims = 1:20, anchor.features = 4000, k.filter = 20)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)

DimPlot(mammal.combined, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_1",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_1",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_1",".pdf",sep=""),width = 16, height = 4)

#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data, rabbit_data2, human_dataA1 ), dims = 1:20, anchor.features = 3000, k.filter = 20)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_2",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "tsne", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_2",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_2",".pdf",sep=""),width = 16, height = 4)



#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data, rabbit_data2, human_dataA1 ), dims = 1:20, anchor.features = 4000, k.filter = 40)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_3",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "tsne", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_3",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_3",".pdf",sep=""),width = 16, height = 4)

#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data, rabbit_data2, human_dataA1 ), dims = 1:20, anchor.features = 3000, k.filter = 40)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_4",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "tsne", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_4",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_4",".pdf",sep=""),width = 16, height = 4)



#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data, rabbit_data2, human_dataA1 ), dims = 1:20, anchor.features = 4000, k.filter = 60)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_5",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "tsne", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_5",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_5",".pdf",sep=""),width = 16, height = 4)

#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data, rabbit_data2, human_dataA1 ), dims = 1:20, anchor.features = 3000, k.filter = 60)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_6",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "tsne", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_6",".pdf",sep=""),width = 16, height = 4)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_6",".pdf",sep=""),width = 16, height = 4)


#Rename to keep variable convention for the join species modelling
mammal.anchors <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data, rabbit_data2, human_dataA1 ), dims = 1:20, anchor.features = 4000, k.filter = 40)
mammal.combined <- IntegrateData(anchorset = mammal.anchors, dims = 1:20)

#mammal.anchors1 <- FindIntegrationAnchors(object.list = list(cynomolgous_data, rabbit_data, human_dataA1 ), dims = 1:20, anchor.features = union(rownames(mammal.combined1), intersect(intersect(intersect(Markers,rownames(rabbit_data)),rownames(cynomolgous_data) ),rownames(human_dataA1) ) ) )
#mammal.combined1 <- IntegrateData(anchorset = mammal.anchors1, dims = 1:20)
DefaultAssay(mammal.combined) <- "integrated"
mammal.combined <- ScaleData(mammal.combined, verbose = FALSE)
mammal.combined <- RunPCA(mammal.combined, npcs = 30, verbose = FALSE)
mammal.combined <- RunUMAP(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- RunTSNE(mammal.combined, reduction = "pca", dims = 1:20)
mammal.combined <- FindNeighbors(mammal.combined, reduction = "pca", dims = 1:20)
#mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)
saveRDS(mammal.combined, file = paste(saveext,"mammal.combined.rds",sep=""))

#Do various plots
DimPlot(mammal.combined, reduction = "pca", group.by = "species", label = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA",".pdf",sep=""),width = 8, height = 4)
ElbowPlot(mammal.combined)
ggsave(filename=paste(saveext,"/DimRed/PCA_var",".pdf",sep=""),width = 8, height = 8)
VizDimLoadings(mammal.combined, dims = 1:6, nfeatures = 30, col = "blue",
               reduction = "pca", projected = FALSE, balanced = FALSE,
               ncol = NULL, combine = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_loadings",".pdf",sep=""),width = 8, height = 20)

DimPlot(mammal.combined, reduction = "tsne", group.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE",".pdf",sep=""),width = 8, height = 4)

DimPlot(mammal.combined, reduction = "umap", group.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP",".pdf",sep=""),width = 8, height = 4)

DimPlot(mammal.combined, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_split",".pdf",sep=""),width = 16, height = 4)

DimPlot(mammal.combined, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/TSNE_split",".pdf",sep=""),width = 16, height = 4)

DimPlot(mammal.combined, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split",".pdf",sep=""),width = 16, height = 4)


panel1 <- c("SALL4",
             "ESRRB",
             "KLF5",
             "KLF4",
             "DPPA5",
             "TBX3",
             "POU5F1",
             "SOX2",
             "PRDM14",
             "LIN28",
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

pannel6 <- c("DAZL",
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


Data <- GetAssayData(mammal.combined,assay = "RNA")


mammal.combined <- FindClusters(mammal.combined, resolution = 0.1)
mammal.combined$Cl1 <- Idents(mammal.combined)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl1",".pdf",sep=""),width = 20, height = 8)

mammal.combined <- FindClusters(mammal.combined, resolution = 0.2)
mammal.combined$Cl2 <- Idents(mammal.combined)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl2",".pdf",sep=""),width = 20, height = 8)

mammal.combined <- FindClusters(mammal.combined, resolution = 0.3)
mammal.combined$Cl3 <- Idents(mammal.combined)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl3",".pdf",sep=""),width = 20, height = 8)

mammal.combined <- FindClusters(mammal.combined, resolution = 0.4)
mammal.combined$Cl4 <- Idents(mammal.combined)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl4",".pdf",sep=""),width = 20, height = 8)

mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)
mammal.combined$Cl5 <- Idents(mammal.combined)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl5",".pdf",sep=""),width = 20, height = 8)

mammal.combined <- FindClusters(mammal.combined, resolution = 0.6)
mammal.combined$Cl6 <- Idents(mammal.combined)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl6",".pdf",sep=""),width = 20, height = 8)

mammal.combined <- FindClusters(mammal.combined, resolution = 0.7)
mammal.combined$Cl7 <- Idents(mammal.combined)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl7",".pdf",sep=""),width = 20, height = 8)

mammal.combined <- FindClusters(mammal.combined, resolution = 0.8)
mammal.combined$Cl8 <- Idents(mammal.combined)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl8",".pdf",sep=""),width = 20, height = 8)

mammal.combined <- FindClusters(mammal.combined, resolution = 0.9)
mammal.combined$Cl9 <- Idents(mammal.combined)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl9",".pdf",sep=""),width = 20, height = 8)


mammal.combined <- FindClusters(mammal.combined, resolution = 1.1)
mammal.combined$Cl11 <- Idents(mammal.combined)
DimPlot(mammal.combined, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl11",".pdf",sep=""),width = 20, height = 8)

DimPlot(mammal.combined, reduction = "pca", split.by = "species2", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_split_Cl11",".pdf",sep=""),width = 20, height = 8)

Idents(mammal.combined) <- mammal.combined$Cl6
mamsubset <- subset(mammal.combined, idents = c("2","4","9","3","1","6","8"), invert = TRUE)
Idents(mamsubset) <- mamsubset$ID1
mamsubset <- subset(mamsubset, idents = c("Tb_CS6","ExMes_CS5","ExMes_CS6/7"), invert = TRUE)

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

Idents(mamsubset) <- colnames(mamsubset)
p <- DimPlot(mamsubset, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/PCA_split_subset_fulllabels",".pdf",sep=""),width = 100, height = 50, p, limitsize = FALSE)

Idents(mamsubset) <- colnames(mamsubset)
p <- DimPlot(mamsubset, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_subset_fulllabels",".pdf",sep=""),width = 100, height = 50, p, limitsize = FALSE)

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
  #VlnPlot(mammal.combined, features = Markers[i])
  #ggsave(filename=paste(saveext,"/Markers/Volcano_full_",Markers[i],".pdf",sep=""),split.by = "species",width = 10, height = 20,limitsize = FALSE)
  FeaturePlot(mammal.combined, features = Markers[i], split.by = "species", cols =  c("lightgrey", "darkblue"), pt.size = 1)
  ggsave(filename=paste(saveext,"/Markers/Marker_full_",Markers[i],".pdf",sep=""),width = 15, height = 5,limitsize = FALSE)
  #VlnPlot(mamsubset, features = Markers[i])
  #ggsave(filename=paste(saveext,"/Markers/Volcano_sub_",Markers[i],".pdf",sep=""),split.by = "species",width = 10, height = 20,limitsize = FALSE)
  FeaturePlot(mamsubset, features = Markers[i], split.by = "species",cols =  c("lightgrey", "darkblue"), pt.size = 1)
  ggsave(filename=paste(saveext,"/Markers/Marker_sub_",Markers[i],".pdf",sep=""),width = 15, height = 5,limitsize = FALSE)
  
  FeaturePlot(mamsubset1, features = Markers[i], split.by = "species",cols =  c("lightgrey", "darkblue"), pt.size = 1)
  ggsave(filename=paste(saveext,"/Markers/Marker_sub1_",Markers[i],".pdf",sep=""),width = 15, height = 5,limitsize = FALSE)
  
  
  FeaturePlot(mammal.combined, reduction = "pca", features = Markers[i], split.by = "species",cols =  c("lightgrey", "darkblue"), pt.size = 1)
  ggsave(filename=paste(saveext,"/Markers/Marker_pca_full_",Markers[i],".pdf",sep=""),width = 15, height = 5,limitsize = FALSE)
  #VlnPlot(mamsubset, features = Markers[i])
  #ggsave(filename=paste(saveext,"/Markers/Volcano_sub_",Markers[i],".pdf",sep=""),split.by = "species",width = 10, height = 20,limitsize = FALSE)
  FeaturePlot(mamsubset, reduction = "pca", features = Markers[i], split.by = "species",cols =  c("lightgrey", "darkblue"), pt.size = 1)
  ggsave(filename=paste(saveext,"/Markers/Marker_pca_sub_",Markers[i],".pdf",sep=""),width = 15, height = 5,limitsize = FALSE)
  
  FeaturePlot(mamsubset1, reduction = "pca", features = Markers[i], split.by = "species",cols =  c("lightgrey", "darkblue"), pt.size = 1)
  ggsave(filename=paste(saveext,"/Markers/Marker_pca_sub1_",Markers[i],".pdf",sep=""),width = 15, height = 5,limitsize = FALSE)  
  
}


#Update cell lineages by cluster
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

Idents(mammal.combined) <- as.factor(newID)
DimPlot(mammal.combined, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_subset_newanno",".pdf",sep=""),width = 20, height = 8)
DimPlot(mammal.combined, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_split_subset_newanno",".pdf",sep=""),width = 20, height = 8)
DimPlot(mammal.combined, reduction = "tsne", split.by = "species", label = TRUE, repel = TRUE)
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

Idents(mamsubset) <- as.factor(newID)
DimPlot(mamsubset, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/PCA_split_subset1_newanno",".pdf",sep=""),width = 20, height = 8)
DimPlot(mamsubset, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_subset1_newanno",".pdf",sep=""),width = 20, height = 8)

#Split by species and cell type for cross-correlation analyses
species2 <- mammal.combined$species
species2[which(mammal.combined$species=="1) Rabbit")] <- "r"
species2[which(mammal.combined$species=="2) Cynomolgous")] <- "c"
species2[which(mammal.combined$species=="4) Human (in vitro)")] <- "h"

Idents(mammal.combined) <- paste(species2,Idents(mammal.combined),sep="")
DefaultAssay(mammal.combined) <- "RNA"

library("pheatmap")
avexp  <- AverageExpression(object = mammal.combined, slot = "data")

a <- avexp$RNA
a <- a[Markers,]

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

cType <- c("cPGC","cExMes","cGastEnd","cEmDisc3","cGast","cEmDisc2","cHyp","cICM","cPTb","cEmDisc1","cEpi","cTb","cVE","cSYS","rVE","rEmDisc3","rEmDisc1","rHyp","rPTb","rTb","rEmDisc2","rGast","rICM","hPTb","hExMes","hGast","hEmDisc1","hEmDisc3","hGastEnd","hSYS","hEmDisc2","hVE","hEpi","hTb","hICM","hHyp")  
BaseCol<-c("#c51b7d","#01665e","#cb181d","#abd9e9","#cb181d","#74add1","#fd8d3c","#2d004b","#1b7837","#4575b4","#313695","#7fbc41","#feb24c","#feb24c","#fe9929","#fe9929","#2166ac","#ec7014","#00441b","#4d9221","#4393c3","#a50f15","#053061","#5aae61","#35978f","#ef3b2c","#045a8d","#3690c0","#ef3b2c","#fdd0a2","#0570b0","#fdd0a2","#023858","#b8e186","#542788","#fdae6b")

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

names(cluster_letters)=rownames(mamsubset@meta.data)
mamsubset <- AddMetaData(object = mamsubset,metadata = cluster_letters,col.name = 'cell.orig')

DimPlot(mamsubset, reduction = "umap", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 2, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/UMAP_subType_Lab2",".pdf",sep=""),width = 20, height = 6)


DimPlot(mamsubset, reduction = "pca", shape.by = 'cell.orig', rule="evenodd", cols = coluse, pt.size = 2, split.by = "species", label = TRUE, repel = TRUE, do.return = TRUE) +
  scale_shape_manual(values = c(18,16,17,15))
ggsave(filename=paste(saveext,"/DimRed/PCA_subType_Lab2",".pdf",sep=""),width = 20, height = 6)
