

library(igraph)

healthy <- read.table("./CD4network/HealthyPIDC.txt", header=FALSE)
mild <- read.table("./CD4network/MildPIDC.txt", header=FALSE)
severe <- read.table("./CD4network/SeverePIDC.txt", header=FALSE)
healthy_genemean <- read.table("./CD4network/HealthyGeneMean.txt")
mild_genemean <- read.table("./CD4network/MildGeneMean.txt")
severe_genemean <- read.table("./CD4network/SevereGeneMean.txt")

toDelete <- seq(1, nrow(healthy), 2) + 1
healthy1 <- healthy[-toDelete, ]
mild1 <- mild[-toDelete, ]
severe1 <- severe[-toDelete, ]

################

healthy1$V3 <- healthy1$V3 / 2
mild1$V3 <- mild1$V3 / 2
severe1$V3 <- severe1$V3 / 2
healthy1 <- subset(healthy1, V3 > 0.999)
mild1 <- subset(mild1, V3 > 0.999)
severe1 <- subset(severe1, V3 > 0.999)

pal <- colorRampPalette(c("mistyrose", "indianred3"))(100)
range1.100 <- function(x){1 + 99*(x-min(x))/(max(x)-min(x))}

# healthy network

net1 <- graph_from_data_frame(d=healthy1, directed=FALSE)
weight1 <- 0.4+0.6*(healthy1$V3-min(healthy1$V3)) / (max(healthy1$V3)-min(healthy1$V3))
E(net1)$width <- weight1 * 5

healthy_genemean <- subset(healthy_genemean, rownames(healthy_genemean) %in% healthy$V1)
genemean1 <- healthy_genemean$x
names(genemean1) <- rownames(healthy_genemean)
genemean1 <- genemean1[order(match(names(genemean1), names(V(net1))))]
genemean1 <- 0.5+0.5*(genemean1-min(genemean1)) / (max(genemean1)-min(genemean1))
V(net1)$size <- genemean1 * 8

cornum1 <- table(healthy$V1)
cornum1 <- cornum1[order(match(names(cornum1), names(V(net1))))]
V(net1)$color <- pal[round(range1.100(cornum1))]

# mild network

net2 <- graph_from_data_frame(d=mild1, directed=FALSE)
weight2 <- 0.4+0.6*(mild1$V3-min(mild1$V3)) / (max(mild1$V3)-min(mild1$V3))
E(net2)$width <- weight2 * 5

mild_genemean <- subset(mild_genemean, rownames(mild_genemean) %in% mild$V1)
genemean2 <- mild_genemean$x
names(genemean2) <- rownames(mild_genemean)
genemean2 <- genemean2[order(match(names(genemean2), names(V(net2))))]
genemean2 <- 0.5+0.5*(genemean2-min(genemean2)) / (max(genemean2)-min(genemean2))
V(net2)$size <- genemean2 * 8

cornum2 <- table(mild$V1)
cornum2 <- cornum2[order(match(names(cornum2), names(V(net2))))]
V(net2)$color <- pal[round(range1.100(cornum2))]

# severe network

net3 <- graph_from_data_frame(d=severe1, directed=FALSE)
weight3 <- 0.4+0.6*(severe1$V3-min(severe1$V3)) / (max(severe1$V3)-min(severe1$V3))
E(net3)$width <- weight3 * 5

severe_genemean <- subset(severe_genemean, rownames(severe_genemean) %in% severe$V1)
genemean3 <- severe_genemean$x
names(genemean3) <- rownames(severe_genemean)
genemean3 <- genemean3[order(match(names(genemean3), names(V(net3))))]
genemean3 <- 0.5+0.5*(genemean3-min(genemean3)) / (max(genemean3)-min(genemean3))
V(net3)$size <- genemean3 * 8

cornum3 <- table(severe$V1)
cornum3 <- cornum2[order(match(names(cornum3), names(V(net3))))]
V(net3)$color <- pal[round(range1.100(cornum3))]



################################################
# detect hub gene
################################################

cor_mat <- function(hdat, mdat, sdat, flag=0) {
    
    if(flag == 1) {
        allgene <- unique(c(mdat$V1, mdat$V2))
        hdatadd <- data.frame(V1=rep(allgene, each=length(allgene)),
                              V2=rep(allgene, times=length(allgene)),
                              V3=0)
        hdat <- rbind(hdat, hdatadd)
        hdat <- subset(hdat, !duplicated(hdat[, -3]))
    }
    
    hmat <- mmat <- smat <- matrix(rep(0, length(allgene)*length(allgene)), 
                                   ncol=length(allgene), nrow=length(allgene))
    
    rownames(hmat) <- allgene
    colnames(hmat) <- allgene
    rownames(mmat) <- allgene
    colnames(mmat) <- allgene
    rownames(smat) <- allgene
    colnames(smat) <- allgene
    
    hdat$V3 <- hdat$V3 / 2
    mdat$V3 <- mdat$V3 / 2
    sdat$V3 <- sdat$V3 / 2
    
    diag(hmat) <- 1
    diag(mmat) <- 1
    diag(smat) <- 1
    
    for(i in 1:nrow(hmat)) {
        for(j in 1:nrow(hmat)) {
            if(i != j) {
                # healthy
                rowposi_h <- hdat$V1 %in% rownames(hmat)[i]
                colposi_h <- hdat$V2 %in% colnames(hmat)[j]
                hmat[i, j] <- hdat$V3[rowposi_h & colposi_h]
                # mild
                rowposi_m <- mdat$V1 %in% rownames(mmat)[i]
                colposi_m <- mdat$V2 %in% colnames(mmat)[j]
                mmat[i, j] <- mdat$V3[rowposi_m & colposi_m]
                # severe
                rowposi_s <- sdat$V1 %in% rownames(smat)[i]
                colposi_s <- sdat$V2 %in% colnames(smat)[j]
                smat[i, j] <- sdat$V3[rowposi_s & colposi_s]
            }
        }
        print(i)
    }
    return(list(hmat, mmat, smat))
}

mat0 <- cor_mat(healthy, mild, severe, flag=1)
saveRDS(mat0, "./CD4network/mat0.rds")

hmat0 <- mat0[[1]]
mmat0 <- mat0[[2]]
smat0 <- mat0[[3]]


min.module.size <- 5
epsilon <- 0.999
healthy_module <- dna::network.modules(hmat0, m=min.module.size, epsilon=epsilon)
mild_module <- dna::network.modules(mmat0, m=min.module.size, epsilon=epsilon)
severe_module <- dna::network.modules(smat0, m=min.module.size, epsilon=epsilon)
F1 <- dna::get.modules(healthy_module)
F2 <- dna::get.modules(mild_module)
F3 <- dna::get.modules(severe_module)

healthy$V3 <- healthy$V3 / 2
mild$V3 <- mild$V3 / 2
severe$V3 <- severe$V3 / 2
healthy <- subset(healthy, V3 > 0.999)
mild <- subset(mild, V3 > 0.999)
severe <- subset(severe, V3 > 0.999)

netlist1 <- vector(mode="list", length=length(table(F1))-1)
netlist2 <- vector(mode="list", length=length(table(F2))-1)
netlist3 <- vector(mode="list", length=length(table(F3))-1)
# healthy
for(i in 1:(length(table(F1))-1)) {
    netlist1[[i]] <- healthy[healthy$V1 %in% names(F1)[F1 == i] | healthy$V2 %in% names(F1)[F1 == i], ]
}
# mild
for(i in 1:(length(table(F2))-1)) {
    netlist2[[i]] <- mild[mild$V1 %in% names(F2)[F2 == i] | mild$V2 %in% names(F2)[F2 == i], ]
}
# severe
for(i in 1:(length(table(F3))-1)) {
    netlist3[[i]] <- severe[severe$V1 %in% names(F3)[F3 == i] | severe$V2 %in% names(F3)[F3 == i], ]
}
hubgene <- function(x) {
    summ <- table(x[, 1])
    summ <- summ[order(summ, decreasing=TRUE)]
    hubname1 <- names(summ)[which((abs(summ - max(summ)) <= max(summ)*0.1) & (var(summ) != 0))]
    hubname2 <- names(summ)[which((summ == max(summ)) & (var(summ) != 0))]
    if (length(hubname2) >= length(summ)*0.3) {
        hubname <- NULL
    } else if (length(hubname1) >= length(summ)*0.3) {
        hubname <- hubname2
    } else {
        hubname <- hubname1
    }
    return(hubname)
}
hubgene1 <- unlist(sapply(netlist1, function(x) hubgene(x)))
hubgene2 <- unlist(sapply(netlist2, function(x) hubgene(x)))
hubgene3 <- unlist(sapply(netlist3, function(x) hubgene(x)))

V(net1)$color <- adjustcolor(ifelse(names(V(net1)) %in% hubgene1, "indianred2", "lightblue1"), alpha.f=0.8)
V(net2)$color <- adjustcolor(ifelse(names(V(net2)) %in% hubgene2, "indianred2", "lightblue1"), alpha.f=0.8)
V(net3)$color <- adjustcolor(ifelse(names(V(net3)) %in% hubgene3, "indianred2", "lightblue1"), alpha.f=0.8)



################################################
# final network visualization
################################################


test.layout1 <- layout_nicely(net1)
test.layout2 <- layout_nicely(net2)
test.layout3 <- layout_nicely(net3)

par(mar=c(0, 0, 0, 0))
plot(net1, vertex.label.color="black", vertex.label.font=2, vertex.label.family="sans",
     vertex.frame.color="gray50", edge.color="gray", vertex.label.cex=1.5, layout=test.layout1)
par(mar=c(0, 0, 0, 0))
plot(net2, vertex.label.color="black", vertex.label.font=2, vertex.label.family="sans",
     vertex.frame.color="gray50", edge.color="gray", vertex.label.cex=1.5, layout=test.layout2)
par(mar=c(0, 0, 0, 0))
plot(net3, vertex.label.color="black", vertex.label.font=2, vertex.label.family="sans",
     vertex.frame.color="gray50", edge.color="gray", vertex.label.cex=1.5, layout=test.layout3)