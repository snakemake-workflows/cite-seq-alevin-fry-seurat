log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
    library(devtools)
    library(ggplot2)
    library(SingleCellExperiment)
    library(Seurat)
})

seurat_object <- readRDS(snakemake@input[[1]])

pdf(file = snakemake@output[[1]])
VlnPlot(seurat_object, features = c("nCount_HTO"), log = TRUE)
dev.off()