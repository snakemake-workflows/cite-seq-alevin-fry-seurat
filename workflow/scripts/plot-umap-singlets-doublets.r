log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
    library(devtools)
    library(ggplot2)
    library(SingleCellExperiment)
    library(Seurat)
})

object <- readRDS(snakemake@input[[1]])

# perform PCA and generate UMAP emdeddings
object <- RunPCA(object, reduction.name = "hto.pca", reduction.key = "HPC_", verbose = F, approx=FALSE)
object <- RunUMAP(object, reduction = "hto.pca", dims = 1:9, reduction.name = "hto.umap", 
                  reduction.key = "HUMAP_", umap.method = 'umap-learn', metric='correlation', verbose = F)

pdf(file = snakemake@output[["by_xlet"]])
DimPlot(object, reduction = "hto.umap", label = F)
dev.off()

pdf(file = snakemake@output[["by_hashtag"]])
DimPlot(object, reduction = "hto.umap", label = T, group.by = "hash.ID" )
dev.off()