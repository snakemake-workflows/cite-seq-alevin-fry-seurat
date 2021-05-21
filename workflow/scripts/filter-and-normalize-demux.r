log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
    library(devtools)
    library(ggplot2)
    library(SingleCellExperiment)
    library(Seurat)
})

object <- subset(object, subset = nCount_HTO < snakemake@params[["max_hto_count"]])

# HTO Normalization
DefaultAssay(object) <- "HTO"
object <- NormalizeData(object,   normalization.method = "CLR", margin = 2, verbose = F)
VariableFeatures(object) <- rownames(object[["HTO"]]@counts)
object <- ScaleData(object, assay = "HTO", verbose = F)

# demultiplex cellular barcodes based on the HTO enrichment and assign single cells back to their samples of origin
object <- HTODemux(object, assay = "HTO", positive.quantile = 0.99, verbose = F)
Idents(object) <- "HTO_classification.global"

saveRDS(object, snakemake@output[[1]])