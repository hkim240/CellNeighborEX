### Add my localhome path to Server
myPaths <- .libPaths()
myPaths <- c("/localhome/bric/kqg803/R/library/", .libPaths())
.libPaths(myPaths)

### Load R packages
library(Seurat)
library(Matrix)

### Import data 
datadir <- system.file("extdata",'SpatialRNA/Vignette',package = 'RCTD')
counts <- read.csv(file.path(datadir,"Hippocampus_MappedDGEForR.csv"))
rownames(counts) <- counts[,1]; counts[,1] <- NULL

### Create Seurat Object
hippo <- CreateSeuratObject(counts)

### Extract final barcodes (unassigned barcodes were removed) instead of QC
list_barcodes <- read.csv(file.path(datadir,"final_barcodes_hippo.csv"))
list_barcodes <- c(list_barcodes)
list_barcodes <- list_barcodes$barcode
list_barcodes <- sapply(list_barcodes, as.character) # command for only our server(bricliinux02)
hippo <- hippo[,colnames(hippo) %in% list_barcodes]

### Normalize data
hippo <- NormalizeData(object = hippo, normalization.method = "LogNormalize", scale.factor = 10000)

### Import cell type markers obtained from mouse hippocampus scRNA-seq for 17 cell types
celltype_markers <- read.csv(file.path(datadir,"hippocampus_celltype_markers.txt"))
celltype_markers <- c(celltype_markers)
celltype_markers <- celltype_markers$gene
celltype_markers <- tolower(celltype_markers)

### Import top genes from Slide-seq (2000 genes)
top_genes <- read.csv(file.path(datadir,"hippocampus_top2000_wo_slasht.txt"))
top_genes <- c(top_genes)
top_genes <- top_genes$gene
top_genes <- tolower(top_genes)

### Get the union of cell type markers and top 2000 genes
integrated <- union(celltype_markers,top_genes) # 2035 genes

### Check whether the cell type markers exist in Slide-seq
selected <- c()
for(i in 1:length(rownames(hippo))) {
  gene <- rownames(hippo)[i]
  gene_lower <- tolower(gene)
  if(gene_lower %in% integrated){ # top_genes plus cell type markers: integrated 
    selected <- c(selected,gene) 
  }
}

### Extract genes of interest in Slide-seq
# selected from integrated: 2034 genes (one cell type markers from scRNA-seq does not exist in Slide-seq)
hippo <- hippo[rownames(hippo) %in% selected, ]
hippo <- FindVariableFeatures(hippo, selection.method = "vst", nfeatures = length(selected))
# export log-normalized data
write.table(hippo[["RNA"]][VariableFeatures(hippo)],"log_data_top_plus.txt")
# export gene_name
df <- data.frame(gene = VariableFeatures(hippo))
write.csv(df, "gene_name_top_plus.txt", row.names=F)
# export cell_id
write.table(colnames(hippo@assays$RNA@data),"cell_id_top_plus.txt")
