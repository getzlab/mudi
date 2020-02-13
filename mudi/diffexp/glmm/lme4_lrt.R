#!/usr/bin/Rscript

suppressMessages(library("lme4"))
suppressMessages(library("lmerTest"))
suppressMessages(library("nlme"))
suppressMessages(library("argparse"))

parser <- ArgumentParser(description="lme4 Negative Binomial Mixed Effect Model for Differential Expression")
parser$add_argument("-i", "--input_counts_file", type="character",
    help="Raw counts file per gene.")
parser$add_argument("-m", "--meta_file", type="character",
    help="Metadata file with covariates.")
parser$add_argument("-d", "--dispersions_file", type="character",
    help="Estimated dispersions file.")
parser$add_argument("-f", "--model", type="character", default="Y ~ 1",
    help="GLMM model to use.")
parser$add_argument("-o", "--out_dir", type="character",
    help="Output directory.")

args <- parser$parse_args()

# --------------------------
# Load inputs
# --------------------------
meta.df <- read.table(args$meta_file, sep=',', header=TRUE)
rownames(meta.df) <- meta.df$index
meta.df$index <- NULL

# Join gene counts
meta.df$Y <- read.table(args$input_counts_file, sep=',')$V2

# Get gene name
gene <- unlist(strsplit(basename(args$input_counts_file), '\\.'))[[1]]

# Pull theta.est
disp.df <- read.table(args$dispersions_file, sep='\t')
rownames(disp.df) <- disp.df$hv_genes
theta.est <- disp.df[gene,"dge.trended.dispersion"]

# Subtype
groups <- colnames(meta.df)[grepl("^groupby_", colnames(meta.df))]

# --------------------------
# Model Fitting
# --------------------------
# Fit NULL model
null.result <- glmer(formula=args$model, data=meta.df, verbose=FALSE, family=MASS::negative.binomial(theta = theta.est))

# Fit each group
counter <- 1
lmer.results = list()

for (fm in as.vector(paste(args$model, groups, sep=" + "))){
  tryCatch(
      expr = {
          lmer.results[[groups[[counter]]]] <- glmer(formula=fm, data=meta.df, verbose=FALSE, family=MASS::negative.binomial(theta = theta.est))
          message(" * Fit model for ", groups[[counter]])
      },
      error = function(e){
          message(" * ERROR: could not fit model for  ", groups[[counter]])
          print(e)
      }
  )

  counter <- counter + 1
}

# Liklihood Ratio Tests
sum.list = vector("list", length(lmer.results))
lrt.list = vector("list", length(lmer.results))

for (i in 1:length(lmer.results)){
    # Summary
    sum.df <- summary(lmer.results[[i]])$coefficient
    colnames(sum.df) <- c('coef','stderr', 'z', 'p_val')
    sum.df <- data.frame(sum.df)
    sum.df$celltype <- names(lmer.results)[[i]]
    sum.list[[i]] <- sum.df

    # Liklihood-Ratio Test
    lrt.df <- anova(lmer.results[[i]], null.result)
    rownames(lrt.df) <- c('null', names(lmer.results)[[i]])
    lrt.list[[i]] <- lrt.df[2,]
}

# Summary dataframe
sum.results.df <- do.call(rbind, sum.list)

# Liklihood ratio test dataframe
lrt.list[[i+1]] <- lrt.df[1,]
lrt.results.df <- do.call(rbind, lrt.list)

# --------------------------
# Save Results
# --------------------------
print(paste("Saving... ", file.path(args$out_dir, paste(gene, ".", "summary.tsv", sep=""))))
write.table(sum.results.df, file=file.path(args$out_dir, paste(gene, ".", "summary.tsv", sep="")), quote=FALSE, sep='\t')

# Liklihood Ratio
print(paste("Saving... ", file.path(args$out_dir, paste(gene, ".", "lrt.tsv", sep=""))))
write.table(lrt.results.df, file=file.path(args$out_dir, paste(gene, ".", "lrt.tsv", sep="")), quote=FALSE, sep='\t')
