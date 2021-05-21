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

pdf(file = snakemake@output[[1]])
VlnPlot(object, features = "nCount_HTO", pt.size = 0.1, log = TRUE)
dev.off()

pdf(file = snakemake@output[[1]])
VlnPlot(object, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()