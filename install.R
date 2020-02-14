#!/usr/bin/Rscript

# Scran
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scran")
BiocManager::install("qvalue")
