log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
    library(eisaR)
    library(Biostrings)
    library(BSgenome)
    library(stringr)
})

# fl is the flank length, here we set it to 
# the read length - 5 (91bp - 5bp = 86bp) 
grl <- getFeatureRanges(
  gtf = snakemake@input[["gtf"]],
  featureType = c("spliced", "intron"), 
  intronType = "separate", 
  flankLength = 86L, 
  joinOverlappingIntrons = FALSE, 
  verbose = TRUE
)

# load the genome sequence
x <- Biostrings::readDNAStringSet(snakemake@input[["fasta"]])
# fix the names
names(x) <- sapply(strsplit(names(x), " "), .subset, 1)

# add the levels and lengths
seqlevels(grl) <- seqlevels(x)
seqlengths(grl) <- seqlengths(x)

# make sure no annotations overhang the 
# references.
grl <- trim(grl)

# extract the sequences
seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = x,
  transcripts = grl
)

# get rid of sequence duplicates
seqs <- unique(seqs)
grl <- grl[names(seqs)]

df <- getTx2Gene(grl)
df[, "status"] = sapply(strsplit(df$transcript_id, "-"), function(x) if(length(x) == 2){"U"} else {"S"})
df[, "gene_id"] = sapply(strsplit(df$gene_id, "-"), function(x) x[1])

writeXStringSet(seqs, snakemake@output[["seq"]], format = "fasta")
write.table(df, snakemake@output[["t2g"]], sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)