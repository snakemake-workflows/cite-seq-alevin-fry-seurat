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

# remove the negatives
object <- subset(object, idents = "Negative", invert = TRUE)

saveRDS(object, snakemake@output[[1]])