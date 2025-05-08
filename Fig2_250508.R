##################
.libPaths("C:/Rprojects/covid")
print(.libPaths())
#R 4.4.1
##################

library(dplyr, lib.loc="C:/Rprojects/covid")
library(Seurat, lib.loc="C:/Rprojects/covid")
library(patchwork, lib.loc="C:/Rprojects/covid")
library(remotes, lib.loc="C:/Rprojects/covid")
library(ggplot2, lib.loc="C:/Rprojects/covid")
library(reshape2, lib.loc="C:/Rprojects/covid")
library(BiocManager, lib.loc="C:/Rprojects/covid")
library(limma, lib.loc="C:/Rprojects/covid")
library(magrittr, lib.loc="C:/Rprojects/covid")
library(scCustomize, lib.loc="C:/Rprojects/covid")
library(qs, lib.loc="C:/Rprojects/covid")
library(harmony, lib.loc="C:/Rprojects/covid")
library(DoubletFinder, lib.loc="C:/Rprojects/covid")
library(rstudioapi, lib.loc="C:/Rprojects/covid")
library(org.Hs.eg.db, lib.loc="C:/Rprojects/covid")
library(AnnotationDbi, lib.loc="C:/Rprojects/covid")
library(clusterProfiler, lib.loc="C:/Rprojects/covid")
library(tidyverse, lib.loc="C:/Rprojects/covid")
library(enrichplot, lib.loc="C:/Rprojects/covid")
library(GOSemSim, lib.loc="C:/Rprojects/covid")
library(igraph, lib.loc="C:/Rprojects/covid")
library(qgraph, lib.loc="C:/Rprojects/covid")
library(ggrepel, lib.loc="C:/Rprojects/covid")
library(viridis, lib.loc="C:/Rprojects/covid")
library(scales, lib.loc="C:/Rprojects/covid")
library(GO.db, lib.loc="C:/Rprojects/covid")
library(nnet, lib.loc="C:/Rprojects/covid")
library(pheatmap, lib.loc="C:/Rprojects/covid")


document_path <- getActiveDocumentContext()$path
document_dir <- dirname(document_path)
setwd(document_dir)


#read data
CSF <- readRDS("CSF.rds")
CSF@meta.data$spikeres <- "No"
CSF@meta.data$spikeres[CSF$CTaa.x == "CAAYGNNRLAF_CASTLGTGGSEQYF"] <- "Yes"
CSF@meta.data$spikeres[CSF$CTaa.x == "CLVGGGNKLTF_CSVGGGQDNTEAFF"] <- "Yes"

#UMAP
new_order <- c("Naive CD4+ T cells", "Effector CD4+ T cells", "Effector CD8+ T cells",
               "Natural killer cells", "B cells", "Monocytes","Myeloid Dendritic cells", "Plasmacytoid Dendritic cells")
CSF$celltypem <- factor(CSF$celltypem, levels = new_order)
Idents(object = CSF) = "celltypem"
p1 <- DimPlot(CSF, reduction = "umap", pt.size = 0.1) + theme_void()+
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "brown", "pink", "grey", "orange"))
p2 <- DimPlot(CSF, reduction = "umap", group.by = "condition1", pt.size = 0.1) + ggtitle(NULL)+ theme_void()
p1+p2

#marker genes
cd_genes <- c("CCR7", "SELL","CD4", "IL7R", "CXCR3", "CD8A", "CCL4", "GNLY", "KLRD1", "CD79A", "CD40",
              "S100A8", "FCN1", "FCGR2A", "C1QC", "LILRA4", "DERL3")
plot1 = DotPlot(object = CSF, features = cd_genes, col.max = 2)
plot2 = plot1 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
plot2

#proportion of the cells
propo <- as.data.frame.matrix(prop.table(table(Idents(CSF), CSF$sample), margin = 2))
propo$average <- rowMeans(propo[, 2:8], na.rm = TRUE)
propo <- propo[, c(1, 10, 2:9)]
colnames(propo) <- c("Patient     ", "Average    ", "Healthy Control1", "Healthy Control2", "Healthy Control3",
                     "Healthy Control4", "Healthy Control5", "Healthy Control6", "Healthy Control7", "Healthy Control8")
propo$Celltype <- rownames(propo)
propo_melt <- melt(propo, id.vars = "Celltype", variable.name = "Sample", value.name = "Value")
propo_melt$Celltype <- factor(propo_melt$Celltype, levels = new_order)
ggplot(propo_melt, aes(x = Sample, y = Value, fill = Celltype)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Proportion") +
  theme_minimal() +
  scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "brown", "pink", "grey", "orange")) +
  theme(
    axis.title = element_text(color = "black"), 
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1.4, vjust = 0, size = 11),
    panel.grid.major.x = element_blank() 
  )
#5.18x4.30

#test
table_result <- table(CSF$celltypem, CSF$condition1)
chisq.test(table_result)
cell_types <- unique(CSF$celltypem)
chisq_results <- lapply(cell_types, function(cell_type) {
  sub_table <- table(CSF$celltypem == cell_type, CSF$condition1)
  chisq.test(sub_table)
})
names(chisq_results) <- cell_types
chisq_results
p_values <- sapply(chisq_results, function(x) x$p.value)
p_adjusted <- p.adjust(p_values, method = "bonferroni")
write.csv(p_adjusted, "celltypetest.csv")
prop_result <- as.data.frame(prop.table(table_result, margin = 2))

################################################
#comparison between patient and healthy control
################################################
#DEGs
Idents(object = CSF) = "condition1"
PrepSCTFindMarkers(CSF, assay = "SCT", verbose = TRUE)
all.markers <- FindMarkers(CSF, ident.1 = "CSF", ident.2 = "healthy_CSF", logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
thre <- 0.25
all.up <- subset(all.markers, p_val_adj < 0.05 & avg_log2FC > thre)
all.down <- subset(all.markers, p_val_adj < 0.05 & -avg_log2FC > thre)

#GO enrichment (BP, UP)
all.up.degs = sort(rownames(all.up), decreasing = TRUE)
all.up.BP = enrichGO(
  all.up.degs,
  OrgDb = "org.Hs.eg.db",
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "fdr",
  qvalueCutoff = 0.01,
  minGSSize = 10,
  maxGSSize = 600,
  readable = TRUE,
  pool = FALSE
)
all.up.BP.df <- as.data.frame(all.up.BP) 
head(all.up.BP.df, n=5) -> all.up.BP.df2
all.up.BP.df2$qvalue <- -log(all.up.BP.df2$qvalue)
all.up.BP.df2$Description <- factor(all.up.BP.df2$Description, levels = all.up.BP.df2$Description[order(all.up.BP.df2$qvalue)])
ggplot(all.up.BP.df2, aes(y = Description, x = Count, fill = qvalue)) +
  geom_bar(stat = "identity") + 
  theme_classic() + 
  scale_fill_gradient(high = "Red") + 
  labs(fill='-log(q.value)') + 
  theme(axis.text.x = element_text(size = 7, color = "black"), axis.text.y = element_text(size = 11, color = "black"), axis.title.y = element_blank()) + 
  theme(legend.key.size = unit(1,"cm")) +
  theme(    axis.title = element_text(color = "black"), 
            axis.text = element_text(color = "black"))+
  scale_x_continuous(expand = c(0,0))
#3x7


################################################
#comparison between patient and healthy control according to Count
################################################
#########DEGs#########
#naive cd4
naivecd4 <- subset(CSF, celltypem == "Naive CD4+ T cells")
PrepSCTFindMarkers(naivecd4, assay = "SCT", verbose = TRUE)
naivecd4.markers <- FindMarkers(naivecd4, ident.1 = "CSF", ident.2 = "healthy_CSF", logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
naivecd4.up <- subset(naivecd4.markers, p_val_adj < 0.05 & avg_log2FC > thre)
naivecd4.down <- subset(naivecd4.markers, p_val_adj < 0.05 & -avg_log2FC > thre)
#effector cd4
effectorcd4 <- subset(CSF, celltypem == "Effector CD4+ T cells")
PrepSCTFindMarkers(effectorcd4, assay = "SCT", verbose = TRUE)
effectorcd4.markers <- FindMarkers(effectorcd4, ident.1 = "CSF", ident.2 = "healthy_CSF", logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
effectorcd4.up <- subset(effectorcd4.markers, p_val_adj < 0.05 & avg_log2FC > thre)
effectorcd4.down <- subset(effectorcd4.markers, p_val_adj < 0.05 & -avg_log2FC > thre)
#effector cd8
effectorcd8 <- subset(CSF, celltypem == "Effector CD8+ T cells")
PrepSCTFindMarkers(effectorcd8, assay = "SCT", verbose = TRUE)
VlnPlot(effectorcd8, c("GZMM", "GZMK", "GZMH", "GZMA", "PRF1", "TNFSF14"))
effectorcd8.markers <- FindMarkers(effectorcd8, ident.1 = "CSF", ident.2 = "healthy_CSF", logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
effectorcd8.up <- subset(effectorcd8.markers, p_val_adj < 0.05 & avg_log2FC > thre)
effectorcd8.down <- subset(effectorcd8.markers, p_val_adj < 0.05 & -avg_log2FC > thre)
#nk
nk <- subset(CSF, celltypem == "Natural killer cells")
PrepSCTFindMarkers(nk, assay = "SCT", verbose = TRUE)
nk.markers <- FindMarkers(nk, ident.1 = "CSF", ident.2 = "healthy_CSF", logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
nk.up <- subset(nk.markers, p_val_adj < 0.05 & avg_log2FC > thre)
nk.down <- subset(nk.markers, p_val_adj < 0.05 & -avg_log2FC > thre)
#bcell
bcell <- subset(CSF, celltypem == "B cells")
PrepSCTFindMarkers(bcell, assay = "SCT", verbose = TRUE)
bcell.markers <- FindMarkers(bcell, ident.1 = "CSF", ident.2 = "healthy_CSF", logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
bcell.up <- subset(bcell.markers, p_val_adj < 0.05 & avg_log2FC > thre)
bcell.down <- subset(bcell.markers, p_val_adj < 0.05 & -avg_log2FC > thre)
#mono
mono <- subset(CSF, celltypem == "Monocytes")
PrepSCTFindMarkers(mono, assay = "SCT", verbose = TRUE)
mono.markers <- FindMarkers(mono, ident.1 = "CSF", ident.2 = "healthy_CSF", logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
mono.up <- subset(mono.markers, p_val_adj < 0.05 & avg_log2FC > thre)
mono.down <- subset(mono.markers, p_val_adj < 0.05 & -avg_log2FC > thre)
#Myeloid Dendritic cells
MDC <- subset(CSF, celltypem == "Myeloid Dendritic cells")
PrepSCTFindMarkers(MDC, assay = "SCT", verbose = TRUE)
MDC.markers <- FindMarkers(MDC, ident.1 = "CSF", ident.2 = "healthy_CSF", logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
MDC.up <- subset(MDC.markers, p_val_adj < 0.05 & avg_log2FC > thre)
MDC.down <- subset(MDC.markers, p_val_adj < 0.05 & -avg_log2FC > thre)

naivecd4.up$celltype <- "naive CD4T"
effectorcd4.up$celltype <- "effector CD4T"
effectorcd8.up$celltype <- "effector CD8T"
nk.up$celltype <- "nk"
bcell.up$celltype <- "bcell"
mono.up$celltype <- "mono"
MDC.up$celltype <- "MDC"
degs <- rbind(naivecd4.up, effectorcd4.up, effectorcd8.up, nk.up, bcell.up, mono.up, MDC.up)

#####comparecluster function (#PDC too few cells to be compared)######
compare_result <- compareCluster(
  geneClusters = list(
    'Naive CD4+ T cells' = rownames(naivecd4.up),
    'Effector CD4+ T cells' = rownames(effectorcd4.up),
    'Effector CD8+ T cells' = rownames(effectorcd8.up),
    'Natural killer cells' = rownames(nk.up),
    'B cells' = rownames(bcell.up),
    'Monocytes' = rownames(mono.up),
    'Myeloid Dendritic cells' = rownames(MDC.up)
  ),
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "fdr",
  qvalueCutoff = 0.01,
  minGSSize = 10,
  maxGSSize = 600,
  keyType = "SYMBOL"
)

#number of GO
ccdf <- as.data.frame(compare_result)
countdf <- as.data.frame(table(ccdf$Cluster))
countdf <- rename(countdf, Celltype = Var1, Count = Freq)
new_row <- data.frame(Celltype = "Plasmacytoid Dendritic cells", Count = 0)
countdf <- rbind(countdf, new_row)
par(mfrow = c(1,1), mar = c(3, 18, 3, 3))
ggplot(countdf, aes(x = reorder(Celltype, -Count), y = Count)) +  
  geom_bar(stat = "identity", fill = "black", width = 0.5) +
  labs(x = "Celltype", y = "N of upregulated pathways", title = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.2, "cm"),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()   
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  
  theme(
    axis.line.x.bottom = element_line(color = "black")  
  )
#4x3 portlait

#emap plot using TOP50 GO
compare_result2 <- pairwise_termsim(compare_resultf)
filtered_result <- compare_resultf@compareClusterResult %>%
  group_by(Cluster) %>%
  arrange(qvalue) %>%
  slice_head(n = 50)
filtered_result2 <- compare_result2 %>% 
  filter(Description %in% filtered_result$Description)
set.seed(123)
emapplot(filtered_result2, max.overlaps=500, showCategory = 200, cex_label_category=0.8, cex_line=0.1)
#portrait 10x12


#Visualize upregulated genes in the pathways related to cytotoxicity
godf <- as.data.frame(filtered_result2)
lmc <- godf[godf$Description %in% c("leukocyte mediated cytotoxicity", "natural killer cell mediated immunity",
                                      "cell killing", "natural killer cell mediated cytotoxicity"), ]
genes.cd8 <- unique(unlist(strsplit(lmc$geneID[lmc$Cluster == "Effector CD8+ T cells"], "/")))
genes.nk <- unlist(strsplit(lmc$geneID[lmc$Cluster == "Natural killer cells"], "/"))

#CD8
effectorcd8.markers$p_val_adj[effectorcd8.markers$p_val_adj < 1e-100] <- 1e-100 # replace too small p value
effectorcd8.markers$avg_log2FC[effectorcd8.markers$avg_log2FC < -1] <- -1 # replace too small p value
effectorcd8.markers$avg_log2FC[effectorcd8.markers$avg_log2FC > 2] <- 2
effectorcd8.markers$'-log(p_val_adj)' <- -log10(effectorcd8.markers$p_val_adj)
cd8label <- subset(effectorcd8.markers, (-log(p_val_adj) > 25 | (-log(p_val_adj) > 5 & avg_log2FC > 0.5)) & rownames(effectorcd8.markers) %in% genes.cd8)
cd8label$gene_name <- rownames(cd8label)

p1 <- ggplot(effectorcd8.markers, aes(avg_log2FC, -log(p_val_adj))) +
  geom_point(aes(color = ifelse(rownames(effectorcd8.markers) %in% genes.cd8, "red", "dimgray")), size = 2/5) +
  geom_text_repel(data = cd8label, aes(label = gene_name), vjust = 1.5, size = 2, color = "red", max.overlaps = 25) +
  scale_color_identity() +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))  + 
  ggtitle("Effector CD8+ T cells") + theme(plot.title = element_text(hjust = 0.5))

#nk
nk.markers$p_val_adj[nk.markers$p_val_adj < 1e-100] <- 1e-100 # replace too small p value
nk.markers$avg_log2FC[nk.markers$avg_log2FC < -1] <- -1 # replace too small p value
nk.markers$avg_log2FC[nk.markers$avg_log2FC > 2] <- 2
nk.markers$'-log(p_val_adj)' <- -log10(nk.markers$p_val_adj)
nklabel <- subset(nk.markers, (-log(p_val_adj) > 25 | (-log(p_val_adj) > 5 & avg_log2FC > 0.5)) & rownames(nk.markers) %in% genes.nk)
nklabel$gene_name <- rownames(nklabel)

p2 <- ggplot(nk.markers, aes(avg_log2FC, -log(p_val_adj))) +
  geom_point(aes(color = ifelse(rownames(nk.markers) %in% genes.nk, "red", "dimgray")), size = 2/5) +
  geom_text_repel(data = nklabel, aes(label = gene_name), vjust = 1.5, size = 2, color = "red", max.overlaps = 25) +
  scale_color_identity() +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))  + 
  ggtitle("Natural killer cells") + theme(plot.title = element_text(hjust = 0.5))

p1+p2



##############################
#Explore subcluster of NK
##############################
#UMAP
nk <- subset(CSF, celltypem == "Natural killer cells")
DefaultAssay(object = nk) <- "integrated" 
nk <- RunPCA(nk)
DimHeatmap(nk, dims = 1:20, cells = 200, balanced = TRUE)
nk <- FindNeighbors(nk, dims = 1:10)
nk <- FindClusters(nk, resolution = 0.8)
nk <- RunUMAP(nk, dims = 1:10)
p1 = DimPlot(nk, reduction = "umap")
p2 = DimPlot(nk, reduction = "umap", group.by = "condition1")
p1+p2

#proportion
proponk <- as.data.frame.matrix(prop.table(table(Idents(nk), nk$condition1), margin = 2))
colnames(proponk) <- c("            Patient", "Healthy Control")
proponk$Cluster <- rownames(proponk)
proponk_melt <- melt(proponk, id.vars = "Cluster", variable.name = "Sample", value.name = "Value")
propo_melt$Celltype <- factor(propo_melt$Celltype, levels = new_order)
ggplot(proponk_melt, aes(x = Sample, y = Value, fill = Cluster)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Proportion") +
  theme_minimal() +
  theme(
    axis.title = element_text(color = "black"), 
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, size = 11),
    panel.grid.major.x = element_blank() 
  )
#4x3 portrait

#heatmap
DefaultAssay(object = nk) <- "SCT" 
PrepSCTFindMarkers(nk, assay = "SCT", verbose = TRUE)
nksub.markers <- FindAllMarkers(nk, only.pos = TRUE, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
nksub.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DefaultAssay(object = nk) <- "RNA" 
nk <- NormalizeData(nk, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(nk)
nk <- ScaleData(nk, features = all.genes)
DoHeatmap(nk, features = top10$gene)
DefaultAssay(object = nk) <- "SCT" 

#Heatmap of the selected genes
avg_exp <- AverageExpression(nk, features = c("FCGR3A", "GZMB", "PRF1", "IFNG", "KLRD1", "NCAM1","CXCR3", "CCR7", "IL7R", "KIT"))
mat <- avg_exp$RNA
mat_t <- t(mat)

pheatmap(mat_t,
         cluster_rows = FALSE,  
         cluster_cols = FALSE,  
         scale = "column")      
#6x3 landscape
