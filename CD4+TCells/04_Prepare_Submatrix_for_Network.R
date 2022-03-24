
library(Seurat)
library(data.table)
library(tidyverse)

CD4_mild <- na.omit(fread("./result/CD4_MildvsHC.txt", sep="\t"))
CD4_severe <- na.omit(fread("./result/CD4_SeverevsHC.txt", sep="\t"))
CD4_severemild <- na.omit(fread("./result/CD4_SeverevsMild.txt", sep="\t"))
CD4 <- full_join(CD4_mild, CD4_severe, by="gene") %>% full_join(CD4_severemild, by="gene")

dat <- readRDS("./result/CD4dat.rds")
exprmat <- as.matrix(dat@assays$RNA@data)
exprmat <- exprmat[rownames(exprmat) %in% CD4$gene, ]

# Split the whole data into three sub data for the development of network

Healthy <- exprmat[, dat$group == "Healthy"]
Healthy <- Healthy[!(apply(Healthy, 1, function(x) sum(x == 0)) >= ncol(Healthy) * 0.99), ]
Mild <- exprmat[, dat$group == "No"]
Mild <- Mild[!(apply(Mild, 1, function(x) sum(x == 0)) >= ncol(Mild) * 0.99), ]
Severe <- exprmat[, dat$group == "Yes"]
Severe <- Severe[!(apply(Severe, 1, function(x) sum(x == 0)) >= ncol(Severe) * 0.99), ]

write.table(Healthy, "./CD4network/Healthy.txt", sep="\t")
write.table(Mild, "./CD4network/Mild.txt", sep="\t")
write.table(Severe, "./CD4network/Severe.txt", sep="\t")

HealthyGeneMean <- rowMeans(Healthy)
MildGeneMean <- rowMeans(Mild)
SevereGeneMean <- rowMeans(Severe)

write.table(HealthyGeneMean, "./CD4network/HealthyGeneMean.txt", sep="\t")
write.table(MildGeneMean, "./CD4network/MildGeneMean.txt", sep="\t")
write.table(SevereGeneMean, "./CD4network/SevereGeneMean.txt", sep="\t")


###############
# Generate "new" data for permutation test

# Mild vs Healthy

HealthyMild <- exprmat[, dat$group == "Healthy" | dat$group == "No"]
metadata <- dat@meta.data
metadata <- metadata[dat$group == "Healthy" | dat$group == "No", ]

(n1 <- sum(metadata$group == "Healthy"))
(n2 <- sum(metadata$group == "No"))
(p <- nrow(HealthyMild))
(n <- ncol(HealthyMild))
num.permutations <- 500
set.seed(2021)
permutation.list <- sample(1:n, n1)
for (i in 2:num.permutations) {
    permutation.list <- rbind(permutation.list, sample(1:n, n1))
}
num.perm <- nrow(permutation.list)
perm.sN <- rep(0, num.perm)

registerDoParallel(10)

foreach(i=1:num.perm) %dopar% {
    i1 <- as.vector(permutation.list[i, ])
    perm.X1 <- HealthyMild[, i1]
    perm.X2 <- HealthyMild[, -i1]
    write.table(perm.X1, paste0("./CD4network/HealthyMildPermDat/permH", i, ".txt"), quote=FALSE, sep="\t")
    write.table(perm.X2, paste0("./CD4network/HealthyMildPermDat/permM", i, ".txt"), quote=FALSE, sep="\t")
}


# Severe vs Healthy

HealthySevere <- exprmat[, dat$group == "Healthy" | dat$group == "Yes"]
metadata <- dat@meta.data
metadata <- metadata[dat$group == "Healthy" | dat$group == "Yes", ]

(n1 <- sum(metadata$group == "Healthy"))
(n2 <- sum(metadata$group == "Yes"))
(p <- nrow(HealthySevere))
(n <- ncol(HealthySevere))
num.permutations <- 500
set.seed(2021)
permutation.list <- sample(1:n, n1)
for (i in 2:num.permutations) {
    permutation.list <- rbind(permutation.list, sample(1:n, n1))
}
num.perm <- nrow(permutation.list)
perm.sN <- rep(0, num.perm)

registerDoParallel(10)

foreach(i=1:num.perm) %dopar% {
    i1 <- as.vector(permutation.list[i, ])
    perm.X1 <- HealthySevere[, i1]
    perm.X2 <- HealthySevere[, -i1]
    write.table(perm.X1, paste0("./CD4network/HealthySeverePermDat/permH", i, ".txt"), quote=FALSE, sep="\t")
    write.table(perm.X2, paste0("./CD4network/HealthySeverePermDat/permS", i, ".txt"), quote=FALSE, sep="\t")
}


# Severe vs Mild

MildSevere <- exprmat[, dat$group == "No" | dat$group == "Yes"]
metadata <- dat@meta.data
metadata <- metadata[dat$group == "No" | dat$group == "Yes", ]

(n1 <- sum(metadata$group == "No"))
(n2 <- sum(metadata$group == "Yes"))
(p <- nrow(MildSevere))
(n <- ncol(MildSevere))
num.permutations <- 500
set.seed(2021)
permutation.list <- sample(1:n, n1)
for (i in 2:num.permutations) {
    permutation.list <- rbind(permutation.list, sample(1:n, n1))
}
num.perm <- nrow(permutation.list)
perm.sN <- rep(0, num.perm)

registerDoParallel(10)

foreach(i=1:num.permutations) %dopar% {
    i1 <- as.vector(permutation.list[i, ])
    perm.X1 <- MildSevere[, i1]
    perm.X2 <- MildSevere[, -i1]
    write.table(perm.X1, paste0("./CD4network/MildSeverePermDat/permM", i, ".txt"), quote=FALSE, sep="\t")
    write.table(perm.X2, paste0("./CD4network/MildSeverePermDat/permS", i, ".txt"), quote=FALSE, sep="\t")
}
