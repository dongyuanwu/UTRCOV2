
library(data.table)
library(tidyverse)
library(stringr)
library(Seurat)

CD8dat <- readRDS("./result/CD8dat.rds")
Idents(CD8dat) <- "group"

# Mild vs Healthy

ids_healthy <- which(CD8dat$group == "Healthy")
output <- data.frame(gene=rownames(CD8dat))
for (i in 1:100) {
    set.seed(i*i)
    ids_covid <- sample(which(CD8dat$group == "Mild"), 4700)
    ids <- ifelse(c(1:length(CD8dat$group)) %in% c(ids_covid, ids_healthy), TRUE, FALSE)
    names(ids) <- colnames(CD8dat)
    newdat <- AddMetaData(CD8dat, ids, col.name="ids")
    newdat <- subset(newdat, subset = ids)
    
    output_new <- as.data.frame(FindMarkers(newdat, ident.1="Mild", ident.2="Healthy", test.use="MAST"))
    output_new$gene <- rownames(output_new)
    output_new <- subset(output_new, select=c("gene", "p_val_adj"))
    colnames(output_new)[2] <- paste0("round", i)
    output <- output %>% full_join(output_new, by="gene")
    
    print(i)
}

save.image("./result/CD8_MildvsHC.RData")
write.table(output, "./result/CD8_MildvsHC.txt", sep="\t", quote=FALSE)


# Severe vs Healthy

ids_healthy <- which(CD8dat$group == "Healthy")
output <- data.frame(gene=rownames(CD8dat))
for (i in 1:100) {
    set.seed(i*i)
    ids_covid <- sample(which(CD8dat$group == "Severe"), 4700)
    ids <- ifelse(c(1:length(CD8dat$group)) %in% c(ids_covid, ids_healthy), TRUE, FALSE)
    names(ids) <- colnames(CD8dat)
    newdat <- AddMetaData(CD8dat, ids, col.name="ids")
    newdat <- subset(newdat, subset = ids)
    
    output_new <- as.data.frame(FindMarkers(newdat, ident.1="Severe", ident.2="Healthy", test.use="MAST"))
    output_new$gene <- rownames(output_new)
    output_new <- subset(output_new, select=c("gene", "p_val_adj"))
    colnames(output_new)[2] <- paste0("round", i)
    output <- output %>% full_join(output_new, by="gene")
    
    print(i)
}

save.image("./result/CD8_SeverevsHC.RData")
write.table(output, "./result/CD8_SeverevsHC.txt", sep="\t", quote=FALSE)


# Severe vs Mild

ids_mild <- which(CD8dat$group == "Mild")
output <- data.frame(gene=rownames(CD8dat))
for (i in 1:100) {
    set.seed(i*i)
    ids_covid <- sample(which(CD8dat$group == "Severe"), 13200)
    ids <- ifelse(c(1:length(CD8dat$group)) %in% c(ids_covid, ids_mild), TRUE, FALSE)
    names(ids) <- colnames(CD8dat)
    newdat <- AddMetaData(CD8dat, ids, col.name="ids")
    newdat <- subset(newdat, subset = ids)
    
    output_new <- as.data.frame(FindMarkers(newdat, ident.1="Severe", ident.2="Mild", test.use="MAST"))
    output_new$gene <- rownames(output_new)
    output_new <- subset(output_new, select=c("gene", "p_val_adj"))
    colnames(output_new)[2] <- paste0("round", i)
    output <- output %>% full_join(output_new, by="gene")
    
    print(i)
}

save.image("./result/CD8_SeverevsMild.RData")
write.table(output, "./result/CD8_SeverevsMild.txt", sep="\t", quote=FALSE)


