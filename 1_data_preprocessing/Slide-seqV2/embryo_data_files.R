### Add my localhome path to Server
myPaths <- .libPaths()
myPaths <- c("/localhome/bric/kqg803/R/library/", .libPaths())
.libPaths(myPaths)

### Load R packages
library(Seurat)
library(Matrix)

### Import mouse embryo Slide-seq count data: 51649 spots and 23124 genes
data <- read.table("data_mouse_embryo_slide-seq.txt",sep="\t", header = T, as.is = T, check.names = F, row.names=1)

### Create seurat object and perform quality control
slideSeq <- CreateSeuratObject(data)
slideSeq <- subset(slideSeq, subset = nCount_RNA > 200) # mouse embryo 42362 spots remained

### Normalize data
slideSeq <- NormalizeData(object = slideSeq, normalization.method = "LogNormalize", scale.factor = 10000)

### Import cell type markers obtained from mouse embryo scRNA-seq for 37 cell types
celltype_markers <- read.csv(file = "embryo_celltype_markers.txt", sep=' ')
celltype_markers <- c(celltype_markers)
celltype_markers <- celltype_markers$gene
celltype_markers <- tolower(celltype_markers)

### Import top genes from Slide-seq (2000 genes)
top_genes <- read.csv(file = "embryo_top2000_wo_slasht.txt", sep=' ')
top_genes <- c(top_genes)
top_genes <- top_genes$gene
top_genes <- tolower(top_genes)

### Get the union of cell type markers and top 2000 genes
integrated <- union(celltype_markers,top_genes) # 2091 genes

### Check whether the cell type markers exist in Slide-seq
selected <- c()
for(i in 1:length(rownames(data))) {
  gene <- rownames(data)[i]
  gene_lower <- tolower(gene)
  if(gene_lower %in% integrated){ # top_genes plus cell type markers: integrated 
    selected <- c(selected,gene) 
  }
}

### Extract genes of interest in Slide-seq
# selected from integrated: 2089 genes (two cell type markers from scRNA-seq does not exist in Slide-seq)
slideSeq <- slideSeq[rownames(slideSeq) %in% selected, ]
slideSeq <- FindVariableFeatures(slideSeq, selection.method = "vst", nfeatures = length(selected))
# export log-normalized data
write.table(slideSeq[["RNA"]][VariableFeatures(slideSeq)],"log_data_top_plus.txt")
# export gene_name
# write.table(rownames(slideSeq@assays$RNA@data),"gene_name_top_plus.txt") # alphabetically ordered gene names
# write.table(rownames(slideSeq),"gene_name_top_plus.txt") # alphabetically ordered gene names
df <- data.frame(gene = VariableFeatures(slideSeq))
write.csv(df, "gene_name_top_plus.txt", row.names=F)
# export cell_id
write.table(colnames(slideSeq@assays$RNA@data),"cell_id_top_plus.txt")
