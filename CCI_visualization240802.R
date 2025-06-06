#visualize CCI data

library(devtools)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(dplyr)
library(Seurat)
library(patchwork)
library(remotes)
library(ggplot2)
library(reshape2)
library(BiocManager)
library("glmGamPoi")
library("limma")
library("clusterProfiler")
library("tidyverse")
library("RColorBrewer")
library("org.Hs.eg.db")
library("DOSE")
library("enrichplot")
library("ggupset")
library("pathview")
library("enrichplot")
library("ggrepel")
library(RColorBrewer)
library(rstudioapi)
library(data.table)
library(reshape2)
library(Cairo)
library(ktplots)
library(SingleCellExperiment)
library(CCPlotR)
library(igraph)

document_path <- getActiveDocumentContext()$path
document_dir <- dirname(document_path)
setwd(document_dir)

#########################################
#CCI visualization network plot
#########################################

#read single cell data
CSF = readRDS("240214_rand_csf_norm.rds")
DefaultAssay(object = CSF) <- "SCT" 

#rename the cluster(s)
CSF@meta.data$celltypem[CSF$seurat_clusters %in% c(0, 3, 4, 8, 14)] <- "Effector CD4+ T cells"
CSF@meta.data$celltypem[CSF$seurat_clusters %in% c(1)] <- "Naive CD4+ T cells"
CSF@meta.data$celltypem[CSF$seurat_clusters %in% c(2, 5, 10)] <- "Effector CD8+ T cells"
CSF@meta.data$celltypem[CSF$seurat_clusters %in% c(6)] <- "Natural killer cells"
CSF@meta.data$celltypem[CSF$seurat_clusters %in% c(7, 11)] <- "Myeloid Dendritic cells"
CSF@meta.data$celltypem[CSF$seurat_clusters %in% c(9)] <- "Monocytes"
CSF@meta.data$celltypem[CSF$seurat_clusters %in% c(12)] <- "Plasmacytoid Dendritic cells"
CSF@meta.data$celltypem[CSF$seurat_clusters %in% c(13)] <- "B cells"

new_order <- c("Plasmacytoid Dendritic cells", "Myeloid Dendritic cells", "Monocytes",
               "B cells", "Natural killer cells",
               "Effector CD8+ T cells", 
               "Effector CD4+ T cells", "Naive CD4+ T cells")
CSF$celltypem <- factor(CSF$celltypem, levels = new_order)

#split the samples
scpatient <- subset(CSF, condition1 == "CSF")
schc <- subset(CSF, condition1 == "healthy_CSF")

#read result from cellphoneDB for patient
pvals <- read.delim("cellphone/results_patient/statistical_analysis_pvalues_06_21_2024_115853.txt", check.names = FALSE)
means <- read.delim("cellphone/results_patient/statistical_analysis_means_06_21_2024_115853.txt", check.names = FALSE)
decon = fread("cellphone/results_patient/statistical_analysis_deconvoluted_06_21_2024_115853.txt", sep="\t")
res <- plot_cpdb_heatmap(pvals = pvals, cellheight = 10, cellwidth = 10, 
                  symmetrical = FALSE, return_tables = TRUE)
patient <- res$count_network
patient <- t(patient) #sender as rownames

#read result from cellphoneDB for healthy control
pvals <- read.delim("cellphone/results_hc/statistical_analysis_pvalues_06_21_2024_120907.txt", check.names = FALSE)
means <- read.delim("cellphone/results_hc/statistical_analysis_means_06_21_2024_120907.txt", check.names = FALSE)
decon = fread("cellphone/results_hc/statistical_analysis_deconvoluted_06_21_2024_120907.txt", sep="\t")
res <- plot_cpdb_heatmap(pvals = pvals, cellheight = 10, cellwidth = 10, 
                         symmetrical = FALSE, return_tables = TRUE)
hc <- res$count_network
hc <- t(hc) #sender as rownames

par(mfrow = c(1,3), mar = c(1, 2, 1, 5))
# plot for Patient
groupSizedf <- as.data.frame(table(scpatient@meta.data$celltypem))
groupSize <- as.numeric(groupSizedf[,2])
order <- match(groupSizedf$Var1, rownames(patient))
patient1s <- patient[order, order]
patient2 <- as.matrix(patient1s) 
netVisual_circle(patient2, vertex.weight = groupSize, weight.scale = TRUE, 
                 edge.weight.max = 32, edge.width.max = 6, label.edge = FALSE, 
                 title.name = "Patient")

# plot for Healthy Control
groupSizedf <- as.data.frame(table(schc@meta.data$celltypem))
groupSize <- as.numeric(groupSizedf[,2])
order <- match(groupSizedf$Var1, rownames(hc))
hc1s <- hc[order, order]
hc2 <- as.matrix(hc1s) 
netVisual_circle(hc2, vertex.weight = groupSize, weight.scale = TRUE, 
                 edge.weight.max = 32, edge.width.max = 6, label.edge = FALSE, 
                 title.name = "Healthy Control")

# plot for the Difference
dif <- patient1s - hc1s
dif[dif < 0] <- 0
dif1 <- as.matrix(dif) 
netVisual_circle(dif1, weight.scale = TRUE, label.edge = TRUE, 
                 title.name = "Difference")

#creat a bar graph for the increase according to the celltype
par(mfrow = c(1,1), mar = c(3, 18, 3, 3))
dif <- as.data.frame(dif)
dif$su <- rowSums(dif)
dif$Label <- rownames(dif)
dif <- dif %>%
  mutate(Label = factor(Label, levels = Label[order(-su)])) #order
ggplot(dif, aes(x = Label, y = su)) +
  geom_bar(stat = "identity", fill = "black", width = 0.5) +
  labs(x = "Source", y = "Increase of interactions", title = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = "black"),
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




###############
#Visualization using for loop
##############
#To visualize all cell types, use the other one.
#interests = c("Natural killer cells", "Effector CD8+ T cells", "Effector CD4+ T cells", "Monocytes", "Myeloid Dendritic cells")
interests = c("Effector CD8+ T cells")

for (intcell in interests) {
  #get patient data
  pvalp <- read.delim("cellphone/results_patient/statistical_analysis_pvalues_06_21_2024_115853.txt", check.names = FALSE)
  names(pvalp) <- gsub("\\|", "_", names(pvalp))
  pvalpl <- melt(pvalp, id.vars = names(pvalp)[1:13], measure.vars = 14:77, variable.name = "variable", value.name = "pvalue")
  pvalpl$hash <- paste0(pvalpl$id_cp_interaction, pvalpl$variable)
  intscop <- read.delim("cellphone/results_patient/statistical_analysis_interaction_scores_06_21_2024_115853.txt", check.names = FALSE)
  names(intscop) <- gsub("\\|", "_", names(intscop))
  intscopl <- melt(intscop, id.vars = names(intscop)[1:13], measure.vars = 14:77, variable.name = "variable", value.name = "intscore")
  intscopl$hash <- paste0(intscopl$id_cp_interaction, intscopl$variable)
  intscopl <- intscopl %>% dplyr::select(hash, intscore)
  mergedp <- merge(pvalpl, intscopl, by="hash")
  pvalpl1 <- mergedp %>% filter(pvalue<0.05)
  pvalpl1$source <- sapply(strsplit(as.character(pvalpl1$variable), "_"), `[`, 1)
  pvalpl1$target <- sapply(strsplit(as.character(pvalpl1$variable), "_"), `[`, 2)
  pvalpl2 <- pvalpl1 %>% filter(source %in% intcell)
  
  #get hc long
  pvalhc <- read.delim("cellphone/results_hc/statistical_analysis_pvalues_06_21_2024_120907.txt", check.names = FALSE)
  names(pvalhc) <- gsub("\\|", "_", names(pvalhc))
  pvalhcl <- melt(pvalhc, id.vars = names(pvalhc)[1:13], measure.vars = 14:77, variable.name = "variable", value.name = "pvalue")
  pvalhcl$hash <- paste0(pvalhcl$id_cp_interaction, pvalhcl$variable)
  intscohc <- read.delim("cellphone/results_hc/statistical_analysis_interaction_scores_06_21_2024_120907.txt", check.names = FALSE)
  names(intscohc) <- gsub("\\|", "_", names(intscohc))
  intscohcl <- melt(intscohc, id.vars = names(intscohc)[1:13], measure.vars = 14:77, variable.name = "variable", value.name = "intscore")
  intscohcl$hash <- paste0(intscohcl$id_cp_interaction, intscohcl$variable)
  intscohcl <- intscohcl %>% dplyr::select(hash, intscore)
  mergedhc <- merge(pvalhcl, intscohcl, by="hash")
  pvalhcl1 <- mergedhc %>% filter(pvalue<0.05)
  pvalhcl1$source <- sapply(strsplit(as.character(pvalhcl1$variable), "_"), `[`, 1)
  pvalhcl1$target <- sapply(strsplit(as.character(pvalhcl1$variable), "_"), `[`, 2)
  pvalhcl2 <- pvalhcl1 %>% filter(source %in% intcell)
  
  #found only in patient
  pvalpl3 <- pvalpl2 %>% filter(!hash %in% pvalhcl2$hash)
  pvalpl3$temp <- sapply(strsplit(as.character(pvalpl3$partner_a), ":"), `[`, 2) 
  pvalpl3$ligand <-pvalpl3$gene_a 
  pvalpl3$ligand[pvalpl3$ligand == ""] <- pvalpl3$temp[pvalpl3$ligand == ""]
  pvalpl3$temp <- sapply(strsplit(as.character(pvalpl3$partner_b), ":"), `[`, 2) 
  pvalpl3$receptor <-pvalpl3$gene_b
  pvalpl3$receptor[pvalpl3$receptor == ""] <- pvalpl3$temp[pvalpl3$receptor == ""]
  pvalpl4 <- pvalpl3 %>% filter(!source==target) #exclude autocrine
  lrp <- dplyr::select(pvalpl4, id_cp_interaction, source, target, ligand, receptor, pvalue, intscore)
  outname <- paste0("source", intcell, ".csv")
  write.csv(lrp, outname)
  
  #plot
  lrp$pvalue[lrp$pvalue == 0] <- 0.001 #minimum cut off = 0.001
  lrp$logp <- -log10(lrp$pvalue) 
  lrp$pair <-paste0(lrp$ligand, " > ", lrp$receptor)
  
  print(ggplot(lrp, aes(x = target, y = pair)) +
    geom_point(aes(size = logp, color = intscore)) +
    scale_color_gradient(low = "blue", high = "red", name = "Interaction Score") +
    scale_size(range = c(1, 3), name = "-log10(pval)") +
    labs(title = paste0("Source: ", intcell),
         x = "Target",
         y = "Ligand -> Receptor",
         size = "logp",
         color = "intscore") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),  
      axis.ticks.length = unit(0.2, "cm"),  
      legend.text = element_text(size = 6),  
      legend.title = element_text(size = 8) 
    ))
}