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

# subsetting the data to Singlets
object <- subset(object, idents = "Singlet")

pdf(file = snakemake@output[["pdf"]])
VlnPlot(object, features = c("nCount_ADT"),  pt.size = 0.1)
dev.off()

saveRDS(object, file = snakemake@output[["rds"]])