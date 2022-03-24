
library(Seurat)
library(MAST)
library(tidyverse)

dat <- readRDS("./result/CD4dat.rds")

Idents(dat) <- "group"

# Mild vs Healthy

ids_healthy <- which(dat$group == "Healthy")
output <- data.frame(gene=rownames(dat))
for (i in 1:100) {
    set.seed(i*i)
    ids_covid <- sample(which(dat$group == "No"), 4000)
    ids <- ifelse(c(1:length(dat$group)) %in% c(ids_covid, ids_healthy), TRUE, FALSE)
    names(ids) <- colnames(dat)
    newdat <- AddMetaData(dat, ids, col.name="ids")
    newdat <- subset(newdat, subset = ids)
    
    output_new <- as.data.frame(FindMarkers(newdat, ident.1="No", ident.2="Healthy", test.use="MAST"))
    output_new$gene <- rownames(output_new)
    output_new <- subset(output_new, select=c("gene", "p_val_adj"))
    colnames(output_new)[2] <- paste0("round", i)
    output <- output %>% full_join(output_new, by="gene")
    
    print(i)
}

write.table(output, "./result/CD4_MildvsHC.txt", sep="\t", quote=FALSE)
save.image("./result/CD4_MildvsHC.RData")


# Severe vs Healthy

ids_healthy <- which(dat$group == "Healthy")
output <- data.frame(gene=rownames(dat))
for (i in 1:100) {
    set.seed(i*i)
    ids_covid <- sample(which(dat$group == "Yes"), 4000)
    ids <- ifelse(c(1:length(dat$group)) %in% c(ids_covid, ids_healthy), TRUE, FALSE)
    names(ids) <- colnames(dat)
    newdat <- AddMetaData(dat, ids, col.name="ids")
    newdat <- subset(newdat, subset = ids)
    
    output_new <- as.data.frame(FindMarkers(newdat, ident.1="Yes", ident.2="Healthy", test.use="MAST"))
    output_new$gene <- rownames(output_new)
    output_new <- subset(output_new, select=c("gene", "p_val_adj"))
    colnames(output_new)[2] <- paste0("round", i)
    output <- output %>% full_join(output_new, by="gene")
    
    print(i)
}

write.table(output, "./result/CD4_SeverevsHC.txt", sep="\t", quote=FALSE)
save.image("./result/CD4_SeverevsHC.RData")


# Severe vs Mild

ids_healthy <- which(dat$group == "No")
output <- data.frame(gene=rownames(dat))
for (i in 1:100) {
    set.seed(i*i)
    ids_covid <- sample(which(dat$group == "Yes"), 23600)
    ids <- ifelse(c(1:length(dat$group)) %in% c(ids_covid, ids_healthy), TRUE, FALSE)
    names(ids) <- colnames(dat)
    newdat <- AddMetaData(dat, ids, col.name="ids")
    newdat <- subset(newdat, subset = ids)
    
    output_new <- as.data.frame(FindMarkers(newdat, ident.1="Yes", ident.2="No", test.use="MAST"))
    output_new$gene <- rownames(output_new)
    output_new <- subset(output_new, select=c("gene", "p_val_adj"))
    colnames(output_new)[2] <- paste0("round", i)
    output <- output %>% full_join(output_new, by="gene")
    
    print(i)
}

write.table(output, "./result/CD4_SeverevsMild.txt", sep="\t", quote=FALSE)
save.image("./result/CD4_SeverevsMild.RData")

