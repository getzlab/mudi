#!/usr/bin/Rscript
suppressMessages(library(edgeR))
suppressMessages(library(argparse))

parser <- ArgumentParser("Estimate Gene-wide Dispersion estimates with EdgeR.")
parser$add_argument("-i", "--input", help="Raw counts.", required=TRUE, type="character")
parser$add_argument("-m", "--metadata", help="Metadata.", required=TRUE, type="character")
parser$add_argument("-o", "--output", help="Output file name.", required=TRUE, type="character")
args <- parser$parse_args()

message("Loading inputs")
# Load Counts
counts <- read.table(args$input, sep=',', header=TRUE, check.names=FALSE)
rownames(counts) <- counts$index
hv_genes <- counts$index

counts$index <- NULL

counts[] <- lapply(counts, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
})

# Load metadata
meta <- read.table(args$metadata, sep=',', header=TRUE)
meta <- data.frame(meta[,-1], row.names = meta[,1])

# Filter counts
counts <- counts[,colnames(counts) %in% rownames(meta)]

# Broad cell-type clusters to use
clusters <- factor(meta$Broad.cell.type)

# Create DE object
dge <- DGEList(counts, group = clusters)

# Calculate norm factors
message("Computing norm factors")
dge <- calcNormFactors(dge)

# Detection Rate
cdr <- scale(colMeans(counts > 0))

# Design matrix
de.design <- model.matrix(~ 1 + cdr + clusters)

# Estimate Dispersions
message("Estimating Dispersions")
dge <- estimateDisp(dge, de.design)

disp_df <- data.frame(hv_genes, dge$trended.dispersion, dge$tagwise.dispersion, dge$common.dispersion)

write.table(disp_df, file=args$output, sep="\t", quote=FALSE)
