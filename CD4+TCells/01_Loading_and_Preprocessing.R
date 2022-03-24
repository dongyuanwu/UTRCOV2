
library(data.table)
library(tidyverse)
library(Seurat)

# CD4+ T cells data: GSE162086
# The original expression data are sets of tsv files

dir <- list.files("./data", pattern="counts.tsv")
datlist <- vector(mode="list", length=length(dir))
for (i in 1:length(dir)) {
    datlist[[i]] <- fread(paste0("./data", dir[i]))
}

metadata <- read.table("./data/GSE162086_seurat_metadata.tsv", sep="\t", header=TRUE)
dict <- metadata[!duplicated(metadata$sample), c("sample", "diagnosis")]

# Merge all datasets into one
library(plyr)
dat <- join_all(datlist, by="V1", type="inner", match="all")
rm(datlist)

# Load the information about disease severity
meta_add <- readxl::read_xlsx("./data/10X samples Demographics_ID.xlsx")
meta_add$GEO_ID...1 <- gsub("_RNA", "", meta_add$GEO_ID...1)
colnames(meta_add)[1] <- "sample"
metadata <- left_join(metadata, meta_add[, c("sample", "Hosptalized")], by="sample")
dict <- metadata[!duplicated(metadata$sample), c("sample", "Hosptalized")]

dat <- dat[rowSums(dat[, -1]) != 0, ]
dat <- as.data.frame(dat)
genes <- as.vector(dat[, 1])
dat <- dat[, -1]
rownames(dat) <- genes

temp <- strsplit(metadata$X, "-")
temp <- sapply(temp, function(x) x[1])
metadata$cells <- paste0(temp, "-", metadata$sample)

dat <- dat[, colnames(dat) %in% metadata$cells]

# Filtering
noninfo_gene_posi <- which(apply(dat, 1, function(x) sum(x == 0)) >= ncol(dat) * 0.99)
dat <- dat[-noninfo_gene_posi, ]
dim(dat)

dat <- dat[, colSums(dat) != 0]

cells <- colnames(dat)
temp <- strsplit(cells, "-")
cname <- sapply(temp, function(x) x[2])
group <- newcluster <- rep("", length(cname))
for (i in 1:length(cname)) {
    group[i] <- dict$Hosptalized[cname[i] == dict$sample]
    newcluster[i] <- metadata$new_cluster_names[cells[i] == metadata$cells]
}

metadata <- data.frame(cname=cname, group=group, newcluster=newcluster, row.names=cells)

rm(temp, noninfo_gene_posi, dir, i, dict, cname, group, newcluster, genes, cells)

# Create a Seurat object

dat <- CreateSeuratObject(counts=dat, meta.data=metadata)
dat <- NormalizeData(dat)
dat <- ScaleData(dat, features=rownames(dat))

saveRDS(dat, "./result/CD4dat.rds")

