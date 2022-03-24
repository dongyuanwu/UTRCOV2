
library(data.table)
library(tidyverse)
library(stringr)
library(Seurat)

# B cells data: GSE164381
# The original expression data are sets of h5 files

filelist <- list.files("./data", pattern=".h5")

datlist <- vector(mode="list", length=length(filelist))
ids <- sapply(strsplit(filelist, "_"), function(x) x[1])
group <- NULL
cname <- NULL
for (i in 1:length(filelist)) {
    datlist[[i]] <- as.data.frame(Read10X_h5(paste0("./data/", filelist[i])))
    datlist[[i]] <- datlist[[i]][-which(rowSums(datlist[[i]]) == 0), ]
    if (grepl("_201657", filelist[i])) {
        group <- c(group, rep("COVID19", ncol(datlist[[i]])))
    } else {
        group <- c(group, rep("Healthy", ncol(datlist[[i]])))
    }
    cname <- c(cname, rep(ids[i], ncol(datlist[[i]])))
    datlist[[i]]$gene_id <- rownames(datlist[[i]])
}
GSE164381 <- reduce(datlist, full_join, by="gene_id")
genename <- GSE164381$gene_id
GSE164381 <- GSE164381[, -which(colnames(GSE164381) == "gene_id")]
rownames(GSE164381) <- genename
GSE164381[is.na(GSE164381)] <- 0

print("Whole dataset done.")


# Filtering

noninfo_gene_posi_bcell <- which(apply(GSE164381, 1, function(x) sum(x == 0)) >= ncol(GSE164381) * 0.975)
GSE164381 <- GSE164381[-noninfo_gene_posi_bcell, ]
dim(GSE164381)
group <- group[colSums(GSE164381) != 0]
cname <- cname[colSums(GSE164381) != 0]
GSE164381 <- GSE164381[, colSums(GSE164381) != 0]
dim(GSE164381)

Bcelldat_metadata <- data.frame(group=group, cname=cname)
rownames(Bcelldat_metadata) <- colnames(GSE164381)

# Create a Seurat object

Bcelldat <- CreateSeuratObject(counts=GSE164381, meta.data=Bcelldat_metadata)
Bcelldat <- NormalizeData(Bcelldat)
Bcelldat <- ScaleData(Bcelldat, features=rownames(Bcelldat))
saveRDS(Bcelldat, "./result/Bcelldat.rds")

