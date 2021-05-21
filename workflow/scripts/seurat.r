log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
    library(devtools)
    library(ggplot2)
    library(SingleCellExperiment)
    library(Seurat)
})

# set the seed 
set.seed(271828)

#' Read alevin-fry quantifications into a SingleCellExperiment object
load_fry <- function(frydir, which_counts = c('S', 'A'), verbose = FALSE) {
  suppressPackageStartupMessages({
    library(rjson)
    library(Matrix)
    library(SingleCellExperiment)
  })
  
  # read in metadata
  meta_info <- fromJSON(file = file.path(frydir, "meta_info.json"))
  ng <- meta_info$num_genes
  usa_mode <- meta_info$usa_mode
  
  if (usa_mode) {
    if (length(which_counts) == 0) {
      stop("Please at least provide one status in 'U' 'S' 'A' ")
    }
    if (verbose) {
      message("processing input in USA mode, will return ", paste(which_counts, collapse = '+'))
    }
  } else if (verbose) {
    message("processing input in standard mode, will return spliced count")
  }

  # read in count matrix
  af_raw <- readMM(file = file.path(frydir, "alevin", "quants_mat.mtx"))
  # if usa mode, each gene gets 3 rows, so the actual number of genes is ng/3
  if (usa_mode) {
    if (ng %% 3 != 0) {
      stop("The number of quantified targets is not a multiple of 3")
    }
    ng <- as.integer(ng/3)
  }
  
  # read in gene name file and cell barcode file
  afg <- read.csv(file.path(frydir, "alevin", "quants_mat_cols.txt"), 
                  strip.white = TRUE, header = FALSE, nrows = ng, 
                  col.names = c("gene_ids"), row.names = 1)
  afc <- read.csv(file.path(frydir, "alevin", "quants_mat_rows.txt"), 
                  strip.white = TRUE, header = FALSE,
                  col.names = c("barcodes"), row.names = 1)

  # if in usa_mode, sum up counts in different status according to which_counts
  if (usa_mode) {
    rd <- list("S" = seq(1, ng), "U" =  seq(ng + 1, 2 * ng),
               "A" =  seq(2 * ng + 1, 3 * ng))
    o <- af_raw[, rd[[which_counts[1]]], drop = FALSE]
    for (wc in which_counts[-1]) {
      o <- o + af_raw[, rd[[wc]], drop = FALSE]
    }
  } else {
    o <- af_raw
  }
  
  # create SingleCellExperiment object
  sce <- SingleCellExperiment(list(counts = t(o)),
                              colData = afc,
                              rowData = afg
  )
  sce
}

hto_q <- load_fry(snakemake@input[["hto"]], verbose = TRUE)
adt_q <- load_fry(snakemake@input[["adt"]], verbose = TRUE)
rna_q <- load_fry(snakemake@input[["rna"]], verbose = TRUE)

common.cells <- intersect(colnames(rna_q), colnames(adt_q))
common.cells <- intersect(common.cells , colnames(hto_q))

gid_to_gname <- read.table(snakemake@input[["geneid2name"]])
rownames(rna_q) <- gid_to_gname$V2[match(rownames(rna_q), gid_to_gname$V1)]

# seurat
seurat_object <- CreateSeuratObject(counts(rna_q)[, which(colnames(rna_q) %in% common.cells)])
seurat_object[["ADT"]] <- CreateAssayObject(counts = counts(adt_q)[, which(colnames(adt_q) %in% common.cells)])
seurat_object[["HTO"]] <- CreateAssayObject(counts = counts(hto_q)[, which(colnames(hto_q) %in% common.cells)])


saveRDS(seurat_object, snakemake@output[[1]])