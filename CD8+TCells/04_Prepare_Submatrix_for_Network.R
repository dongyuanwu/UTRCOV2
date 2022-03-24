
library(Seurat)
library(data.table)
library(tidyverse)

CD8_mild <- na.omit(fread("./result/CD8_MildvsHC.txt", sep="\t"))
CD8_severe <- na.omit(fread("./result/CD8_SeverevsHC.txt", sep="\t"))
CD8_severemild <- na.omit(fread("./result/CD8_SeverevsMild.txt", sep="\t"))
CD8 <- full_join(CD8_mild, CD8_severe, by="gene") %>% full_join(CD8_severemild, by="gene")

dat <- readRDS("./result/CD8dat.rds")
exprmat <- as.matrix(dat@assays$RNA@data)
exprmat <- exprmat[rownames(exprmat) %in% CD8$gene, ]

# Split the whole data into three sub data for the development of network

Healthy <- exprmat[, dat$group == "Healthy"]
Healthy <- Healthy[!(apply(Healthy, 1, function(x) sum(x == 0)) >= ncol(Healthy) * 0.99), ]
Mild <- exprmat[, dat$group == "No"]
Mild <- Mild[!(apply(Mild, 1, function(x) sum(x == 0)) >= ncol(Mild) * 0.99), ]
Severe <- exprmat[, dat$group == "Yes"]
Severe <- Severe[!(apply(Severe, 1, function(x) sum(x == 0)) >= ncol(Severe) * 0.99), ]

write.table(Healthy, "./CD8network/Healthy.txt", sep="\t")
write.table(Mild, "./CD8network/Mild.txt", sep="\t")
write.table(Severe, "./CD8network/Severe.txt", sep="\t")

HealthyGeneMean <- rowMeans(Healthy)
MildGeneMean <- rowMeans(Mild)
SevereGeneMean <- rowMeans(Severe)

write.table(HealthyGeneMean, "./CD8network/HealthyGeneMean.txt", sep="\t")
write.table(MildGeneMean, "./CD8network/MildGeneMean.txt", sep="\t")
write.table(SevereGeneMean, "./CD8network/SevereGeneMean.txt", sep="\t")


###############
# Generate "new" data for permutation test

# Mild vs Healthy

HealthyMild <- exprmat[, dat$group == "Healthy" | dat$group == "Mild"]
metadata <- dat@meta.data
metadata <- metadata[dat$group == "Healthy" | dat$group == "Mild", ]

(n1 <- sum(metadata$group == "Healthy"))
(n2 <- sum(metadata$group == "Mild"))
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
    write.table(perm.X1, paste0("./CD8network/HealthyMildPermDat/permH", i, ".txt"), quote=FALSE, sep="\t")
    write.table(perm.X2, paste0("./CD8network/HealthyMildPermDat/permM", i, ".txt"), quote=FALSE, sep="\t")
}


# Severe vs Healthy

HealthySevere <- exprmat[, dat$group == "Healthy" | dat$group == "Severe"]
metadata <- dat@meta.data
metadata <- metadata[dat$group == "Healthy" | dat$group == "Severe", ]

(n1 <- sum(metadata$group == "Healthy"))
(n2 <- sum(metadata$group == "Severe"))
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
    write.table(perm.X1, paste0("./CD8network/HealthySeverePermDat/permH", i, ".txt"), quote=FALSE, sep="\t")
    write.table(perm.X2, paste0("./CD8network/HealthySeverePermDat/permS", i, ".txt"), quote=FALSE, sep="\t")
}


# Severe vs Mild

MildSevere <- exprmat[, dat$group == "Mild" | dat$group == "Severe"]
metadata <- dat@meta.data
metadata <- metadata[dat$group == "Mild" | dat$group == "Severe", ]

(n1 <- sum(metadata$group == "Mild"))
(n2 <- sum(metadata$group == "Severe"))
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
    write.table(perm.X1, paste0("./CD8network/MildSeverePermDat/permM", i, ".txt"), quote=FALSE, sep="\t")
    write.table(perm.X2, paste0("./CD8network/MildSeverePermDat/permS", i, ".txt"), quote=FALSE, sep="\t")
}
