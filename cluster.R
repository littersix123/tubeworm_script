library(Seurat)
library(dplyr)
library(patchwork)
library(reshape2)
library(ggplot2)

rds <- c("/dellfsqd2/ST_OCEAN/USER/zhangxianghui/05.Deep_sea/1.Riftia_achyptila/quality_control/QC_3/SY530-GC03-PL-1/03.analysis/Clustering/clustering_annotation_object.RDS",
"/dellfsqd2/ST_OCEAN/USER/zhangxianghui/05.Deep_sea/1.Riftia_achyptila/quality_control/QC_7/SY530-GC03-PL-cDNA-2-J/03.analysis/Clustering/clustering_annotation_object.RDS",
"/dellfsqd2/ST_OCEAN/USER/zhangxianghui/05.Deep_sea/1.Riftia_achyptila/quality_control/QC_3/SY530-GC03-PL-3/SY530-GC03-PL-3/03.analysis/Clustering/clustering_annotation_object.RDS",
"/dellfsqd2/ST_OCEAN/USER/zhangxianghui/05.Deep_sea/1.Riftia_achyptila/quality_control/QC_7/SY530-GC03-PL-cDNA-4-J/03.analysis/Clustering/clustering_annotation_object.RDS",
"/dellfsqd2/ST_OCEAN/USER/zhangxianghui/05.Deep_sea/1.Riftia_achyptila/quality_control/QC_3/SY530-GC04-PL-1/03.analysis/Clustering/clustering_annotation_object.RDS",
"/dellfsqd2/ST_OCEAN/USER/zhangxianghui/05.Deep_sea/1.Riftia_achyptila/quality_control/QC_7/SY530-GC04-PL-cDNA-2-J/03.analysis/Clustering/clustering_annotation_object.RDS",
"/dellfsqd2/ST_OCEAN/USER/zhangxianghui/05.Deep_sea/1.Riftia_achyptila/quality_control/QC_7/SY530-GC04-PL-cDNA-3-J/03.analysis/Clustering/clustering_annotation_object.RDS")

lib_name <- c("SY530-GC03-PL-1",
              "SY530-GC03-PL-cDNA-2-J",
              "SY530-GC03-PL-3",
              "SY530-GC03-PL-cDNA-4-J",
              "SY530-GC04-PL-1",
              "SY530-GC04-PL-cDNA-2-J",
              "SY530-GC04-PL-cDNA-3-J")

sc_list <- list()

for (i in 1:length(rds)){
                object <- readRDS(rds[i])
                name <- lib_name[i]
                object <- RenameCells(object, add.cell.id=name)
                object$stim <- name
                object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = "ONT.9786|PE-Scaf9739-0.0|PE-Scaf9739-0.1|PE-Scaf9739-0.3|ONT.9790|ONT.9783|ONT.9783|ONT.9783|ONT.5018|ONT.9780")
                object <- subset(object,percent.mt < 20)
                sc_list <- c(sc_list,object)}

sct.list <- lapply(X = sc_list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = sct.list,nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = sct.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",anchor.features = features)
combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
pdf("./ElbowPlot.pdf")
ElbowPlot(combined,ndims=30)
dev.off()
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

#UMAP
pdf("./umap.pdf",width=10,height=5)
p1 <- DimPlot(combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
dev.off()

#Different libraries umap
pdf("./DimPlot.pdf",width=10,height=5)
DimPlot(combined,reduction="umap",split.by="stim",ncol=3)
dev.off()

#Generate a list of differentially expressed genes
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)
combined <- ScaleData(combined,features=rownames(combined))

markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
topn <- markers %>% group_by(cluster)
write.table(topn[topn$p_val_adj <= 0.05,],file="RNA_GeneDiffExpFilter.xls",sep="\t",col.names = TRUE,row.names = F,quote=F)

#Heatmap
markers %>% group_by(cluster) %>%  top_n(n = 50, wt = avg_log2FC) -> top50
write.table(top50,file="top50_RNA_GeneDiffExpFilter.xls",sep="\t",col.names = TRUE,row.names = F,quote=F)
pdf("DoHeatmap.pdf")
DoHeatmap(combined, features = top50$gene) + NoLegend()
dev.off()

