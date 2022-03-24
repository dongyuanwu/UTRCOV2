
library(doParallel)
library(foreach)
library(dna)

# Correlation Matrix Function
cor_mat <- function(hdat, cdat, flag=0) {
    
    allgene <- unique(c(cdat$V1, cdat$V2))
    if(flag == 1) {
        hdatadd <- data.frame(V1=rep(allgene, each=length(allgene)),
                              V2=rep(allgene, times=length(allgene)),
                              V3=0)
        hdat <- rbind(hdat, hdatadd)
        hdat <- subset(hdat, !duplicated(hdat[, -3]))
    }
    
    hmat <- cmat <- matrix(rep(0, length(allgene)*length(allgene)), ncol=length(allgene), nrow=length(allgene))
    
    rownames(hmat) <- allgene
    colnames(hmat) <- allgene
    rownames(cmat) <- allgene
    colnames(cmat) <- allgene
    
    hdat$V3 <- hdat$V3 / 2
    cdat$V3 <- cdat$V3 / 2
    
    diag(hmat) <- 1
    diag(cmat) <- 1
    
    for(i in 1:nrow(hmat)) {
        for(j in 1:nrow(hmat)) {
            if(i != j) {
                # healthy
                rowposi_h <- hdat$V1 %in% rownames(hmat)[i]
                colposi_h <- hdat$V2 %in% colnames(hmat)[j]
                hmat[i, j] <- hdat$V3[rowposi_h & colposi_h]
                # covid
                rowposi_c <- cdat$V1 %in% rownames(cmat)[i]
                colposi_c <- cdat$V2 %in% colnames(cmat)[j]
                cmat[i, j] <- cdat$V3[rowposi_c & colposi_c]
            }
        }
    }
    return(list(hmat, cmat))
}


############################
# Mild vs Healthy

healthy <- read.table("./CD4network/HealthyPIDC.txt", header=FALSE)
covid <- read.table("./CD4network/MildPIDC.txt", header=FALSE)

mat0 <- cor_mat(healthy, covid, flag=1)

hmat0 <- mat0[[1]]
cmat0 <- mat0[[2]]

nG <- nrow(cmat0)
di <- rep(0, nG)
for (g in 1:nG) {
    di[g] <- sum(abs(hmat0[g, ] - cmat0[g, ]))
}
di <- di / (nG - 1)

perm.di <- rep(0, nG)
count.ind <- rep(0, nG)
num.perm <- 500

cat("Starting permutation test:\n")

registerDoParallel(60)
perm.sN <- foreach(i = 1:num.perm, .combine=rbind) %dopar% {
    cat("permutation", i, "out of", num.perm,"\n")
    hperm <- read.table(paste0("./CD4network/HealthyMildPIDC/netH", i, ".txt"), header=FALSE)
    cperm <- read.table(paste0("./CD4network/HealthyMildPIDC/netM", i, ".txt"), header=FALSE)
    matperm <- cor_mat(hperm, cperm)
    perm.s1 <- matperm[[1]]
    perm.s2 <- matperm[[2]]
    for (g in 1:nG) {
        perm.di[g] <- sum(abs(perm.s1[g, ] - perm.s2[g, ]))
    }
    perm.di <- perm.di / (nG - 1)
    perm.di
}

p.value.ind <- rep(0, nG)
for (g in 1:nG) {
    p.value.ind[g] <- sum(perm.sN[, g] >= di[g]) / num.perm
}

names(p.value.ind) <- rownames(hmat0)
names(di) <- rownames(hmat0)
new("resultsIndTest", p.values=p.value.ind, d=di)


save.image("./CD4network/Permtest_HM.RData")


############################
# Severe vs Healthy

healthy <- read.table("./CD4network/HealthyPIDC.txt", header=FALSE)
covid <- read.table("./CD4network/SeverePIDC.txt", header=FALSE)

mat0 <- cor_mat(healthy, covid, flag=1)

hmat0 <- mat0[[1]]
cmat0 <- mat0[[2]]

nG <- nrow(cmat0)
di <- rep(0, nG)
for (g in 1:nG) {
    di[g] <- sum(abs(hmat0[g, ] - cmat0[g, ]))
}
di <- di / (nG - 1)

perm.di <- rep(0, nG)
count.ind <- rep(0, nG)
num.perm <- 500

cat("Starting permutation test:\n")

registerDoParallel(60)
perm.sN <- foreach(i = 1:num.perm, .combine=rbind) %dopar% {
    cat("permutation", i, "out of", num.perm,"\n")
    hperm <- read.table(paste0("./CD4network/HealthySeverePIDC/netH", i, ".txt"), header=FALSE)
    cperm <- read.table(paste0("./CD4network/HealthySeverePIDC/netS", i, ".txt"), header=FALSE)
    matperm <- cor_mat(hperm, cperm)
    perm.s1 <- matperm[[1]]
    perm.s2 <- matperm[[2]]
    for (g in 1:nG) {
        perm.di[g] <- sum(abs(perm.s1[g, ] - perm.s2[g, ]))
    }
    perm.di <- perm.di / (nG - 1)
    perm.di
}

p.value.ind <- rep(0, nG)
for (g in 1:nG) {
    p.value.ind[g] <- sum(perm.sN[, g] >= di[g]) / num.perm
}

names(p.value.ind) <- rownames(hmat0)
names(di) <- rownames(hmat0)
new("resultsIndTest", p.values=p.value.ind, d=di)


save.image("./CD4network/Permtest_HS.RData")




############################
# Severe vs Mild

mild <- read.table("./CD4network/MildPIDC.txt", header=FALSE)
severe <- read.table("./CD4network/SeverePIDC.txt", header=FALSE)

mat0 <- cor_mat(mild, severe, flag=1)

hmat0 <- mat0[[1]]
cmat0 <- mat0[[2]]

nG <- nrow(cmat0)
di <- rep(0, nG)
for (g in 1:nG) {
    di[g] <- sum(abs(hmat0[g, ] - cmat0[g, ]))
}
di <- di / (nG - 1)

perm.di <- rep(0, nG)
count.ind <- rep(0, nG)
num.perm <- 500

cat("Starting permutation test:\n")

registerDoParallel(60)
perm.sN <- foreach(i = 1:num.perm, .combine=rbind) %dopar% {
    cat("permutation", i, "out of", num.perm,"\n")
    hperm <- read.table(paste0("./CD4network/MildSeverePIDC/netM", i, ".txt"), header=FALSE)
    cperm <- read.table(paste0("./CD4network/MildSeverePIDC/netS", i, ".txt"), header=FALSE)
    matperm <- cor_mat(hperm, cperm)
    perm.s1 <- matperm[[1]]
    perm.s2 <- matperm[[2]]
    for (g in 1:nG) {
        perm.di[g] <- sum(abs(perm.s1[g, ] - perm.s2[g, ]))
    }
    perm.di <- perm.di / (nG - 1)
    perm.di
}

p.value.ind <- rep(0, nG)
for (g in 1:nG) {
    p.value.ind[g] <- sum(perm.sN[, g] >= di[g]) / num.perm
}

names(p.value.ind) <- rownames(hmat0)
names(di) <- rownames(hmat0)
new("resultsIndTest", p.values=p.value.ind, d=di)


save.image("./CD4network/Permtest_MS.RData")
