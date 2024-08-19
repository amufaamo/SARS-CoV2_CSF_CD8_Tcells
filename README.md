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
