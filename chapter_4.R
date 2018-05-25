#install.packages("ade4")
library(ade4)
install.packages("vegan")
library(vegan)
install.packages("gclus")
library(gclus)
install.packages('cluster')
install.packages("RColorBrewer")
install.packages('labdsv')
install.packages("mvpart")
install.packages("MVPARTwrap")
library(cluster)
library(RColorBrewer)
library(labdsv)
library(mvpart)
library(MVPARTwrap)
#env <- read.csv("DoubesEnv.csv", row.names=1)
#dat <- data(doubs.species)
#species <- doubs[["species"]]
#spe.norm <- decostand(species, "normalize")
spe <- spe[-8,]
spe.norm <- decostand(spe,"normalize")
spe.ch <- vegdist(spe.norm,"euc")
spe.ch.single <- hclust(spe.ch, method="single")
plot(spe.ch.single)

spe.ch.complete <- hclust(spe.ch,method="complete")
plot(spe.ch.complete)

spe.ch.UPGMA <- hclust(spe.ch, method="average")
plot(spe.ch.UPGMA)

spe.ch.centroid <- hclust(spe.ch, method="centroid")
plot(spe.ch.centroid)

spe.ch.ward <- hclust(spe.ch, method="ward")
plot(spe.ch.ward)

spe.ch.single.coph <- cophenetic(spe.ch.single)
cor(spe.ch, spe.ch.single.coph)
spe.ch.comp.coph <- cophenetic(spe.ch.complete)
cor(spe.ch, spe.ch.comp.coph)
spe.ch.UPGMA.coph <- cophenetic(spe.ch.UPGMA)
cor(spe.ch, spe.ch.UPGMA.coph)
spe.ch.ward.coph <- cophenetic(spe.ch.ward)
cor(spe.ch, spe.ch.ward.coph)
cor(spe.ch, spe.ch.ward.coph, method="spearman")

spe.chwo <- reorder.hclust(spe.ch.ward, spe.ch)
dend <- as.dendrogram(spe.chwo)
heatmap(as.matrix(spe.ch),Rowv = dend, symm = TRUE, margin=c(3,3))
or <- vegemite(spe, spe.chwo)
heatmap(t(spe[rev(or$species)]),Rowv = NA, Colv = dend,
        col = c("white",brewer.pal(5,"Greens")),scale = "none",
        magin = c(4,4),ylab = "Species(weighted averages of sites)", xlab="Sites")

spe.kmeans <- kmeans(spe.norm, centers=4, nstart=100)
k <- 4
spebc.ward.g <- cutree(spe.ch.ward, k)
table(spe.kmeans$cluster, spebc.ward.g)



