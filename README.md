Data can be downloaded from Google drive https://drive.google.com/drive/folders/1mLBZx9fhw4eDRdua4umsr61UsxNuHVE8?usp=sharing

## Fig2
### (B)
```R
so <- source('so_scf')
DimPlot(so, group.by='celltype')
DimPlot(so, group.by='condition1')
```
### (C)
```R
DefaultAssay(so) <- 'RNA'
features <- c('CCR7', 'SELL', 'CD4', 'IL7R', 'CXCR3', 'CD8A', 'CCL4', 'GNLY', 'KLRD1', 'CD79A', 'CD40', 'S100A9', 'FCN1', 'FCGR2A', 'C1QC', 'LILRA4', 'DERL3')
DotPlot(so,
        features = features,
        group.by = 'celltype'
       ) + theme(axis.text.x = element_text(angle = 90, hjust=1)) 
```
### (D)
```R
csv <- read.csv('goenrichment_csf_vs_publiccsf.csv')
csv %>% dplyr::filter(Cluster == 'up') -> csv
head(csv, n=5) -> csv2
csv2$qvalue <- -log(csv2$qvalue)
ggplot(csv2, aes(y = Description, x = Count, fill = qvalue)) +
  geom_bar(stat = "identity") + 
  theme_classic() + 
  scale_fill_gradient(high = "Red") + 
  labs(fill='-log(q.value)') + 
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_blank()) + 
  theme(legend.key.size = unit(1,"cm")) + 
  scale_x_continuous(expand = c(0,0))
```
### (E)
```R
csv <- read.csv('goenrichment_csf_vs_publiccsf_bp.csv')
csv %>% group_by(Cluster) %>% slice_head(n = 5) -> csv2
lists <- unique(csv2$Description)
csv3 <- csv %>% dplyr::filter(Description %in% lists)
csv4 <- transform(csv3, Cluster = factor(Cluster, levels = c("Effector CD8+ T cells", "Effector CD4+ T cells", "Naive CD4+ T cells", "Myeloid Dendritic cells")))
g <- ggplot(csv4) + 
geom_point(aes(x = Cluster, y = Description, size= GeneRatio, color = qvalue)) + 
theme_classic() + 
scale_color_gradient(low = 'red', high = 'blue') + guides(colour = guide_colourbar(reverse = TRUE)) +
theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), axis.title = element_blank())
```
### (F)
### (G)
### (H)

## Fig.3
### (B)
library(Seurat)
data <- readRDS('')
