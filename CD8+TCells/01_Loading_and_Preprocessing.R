
library(data.table)
library(tidyverse)
library(stringr)
library(Seurat)

# CD8+ T cells data: GSE153931

GSE153931 <- as.data.frame(fread("./data/GSE153931_cd8_t24_processed_data_umi_counts.txt.gz"))
genename <- GSE153931[, 1]
rownames(GSE153931) <- genename
GSE153931 <- GSE153931[, -1]

GSE153931_metadata <- fread("./data/GSE153931_cd8_t24_processed_data_annotations.txt")

# Preprocessing

GSE153931 <- as.data.frame(GSE153931)[, !(GSE153931_metadata$orig.severity %in% c("Healthy_FLU", "Healthy_RSV", "void"))]
GSE153931_metadata <- GSE153931_metadata[!(GSE153931_metadata$orig.severity %in% c("Healthy_FLU", "Healthy_RSV", "void")), ]
GSE153931_metadata$group <- "Healthy"
GSE153931_metadata$group <- ifelse(GSE153931_metadata$orig.severity == "Mild", "Mild", GSE153931_metadata$group)
GSE153931_metadata$group <- ifelse(GSE153931_metadata$orig.severity %in% c("ITU", "NITU"), "Severe", GSE153931_metadata$group)
table(GSE153931_metadata$group)

# Filtering

GSE153931_metadata <- as.data.frame(GSE153931_metadata)
rownames(GSE153931_metadata) <- GSE153931_metadata$barcode
noninfo_gene_posi <- which(apply(GSE153931, 1, function(x) sum(x == 0)) >= ncol(GSE153931) * 0.975)
GSE153931 <- GSE153931[-noninfo_gene_posi, ]
dim(GSE153931)
GSE153931_metadata <- GSE153931_metadata[colSums(GSE153931) != 0, ]
GSE153931 <- GSE153931[, colSums(GSE153931) != 0]
dim(GSE153931)

# Create a Seurat object

CD8dat <- CreateSeuratObject(counts=GSE153931, meta.data=GSE153931_metadata)
CD8dat <- NormalizeData(CD8dat)
CD8dat <- ScaleData(CD8dat, features=rownames(CD8dat))
saveRDS(CD8dat, "./result/CD8dat.rds")

