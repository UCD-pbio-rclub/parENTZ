###########
# Upendra Kumar Devisetty
# upendrakumar.devisetty@googlemail.com
# B rapa parental RNA-seq 
# Modified December 12, 2015
###########

# GC pool with v1.5 Reference 
# Clustering all tissues together......
rm(list=ls())
setwd("~/Dropbox/Brassica_Parental_stuff/DE_v1.5/GC_only")

# Load packages
#source("http://bioconductor.org/biocLite.R")
library(edgeR)
library(reshape2)
library(ggplot2)
#biocLite("genefilter")
library("genefilter") 
#biocLite("hopach")
library(hopach)
#install.packages("devtools")
library(devtools)
#install_github("rafalib","ririzarr")
library(rafalib)

# Reading the couns file
counts_GC_final_merged <- read.delim("GC_merged_v1.5_mapping.tsv",row.names=NULL)
names(counts_GC_final_merged)
head(counts_GC_final_merged)
dim(counts_GC_final_merged) # [1] 41921    67
counts_GC_final_merged<-counts_GC_final_merged[counts_GC_final_merged$gene!="*",]
row.names(counts_GC_final_merged) <- counts_GC_final_merged$gene
counts_GC_final_merged <- counts_GC_final_merged[,-1]
counts_GC_final_merged[is.na(counts_GC_final_merged)] <- 0
counts_GC_final_merged <- round(counts_GC_final_merged)
head(counts_GC_final_merged)
counts_GC_final_merged <- counts_GC_final_merged[,colSums(counts_GC_final_merged) > 100000]
head(counts_GC_final_merged)[1:3]
dim(counts_GC_final_merged) # [1] 41920    64 (2 samples lost...)

#creating a samples dataframe
samples <- data.frame(
  file = names(counts_GC_final_merged),
  gt = sub("(IMB211|R500)_(SHADE|SUN)_([1-3])_(INTERNODE|LEAF|ROOT|SEEDLING|SILIQUE|APICALMERISTEM|FLORALMERISTEM).[[:print:]]+(.bam)","\\1",names(counts_GC_final_merged)),
  trt  = sub("(IMB211|R500)_(SHADE|SUN)_([1-3])_(INTERNODE|LEAF|ROOT|SEEDLING|SILIQUE|APICALMERISTEM|FLORALMERISTEM).[[:print:]]+(.bam)","\\2",names(counts_GC_final_merged)),
  rep  = sub("(IMB211|R500)_(SHADE|SUN)_([1-3])_(INTERNODE|LEAF|ROOT|SEEDLING|SILIQUE|APICALMERISTEM|FLORALMERISTEM).[[:print:]]+(.bam)","\\3",names(counts_GC_final_merged)),
  tiss  = sub("(IMB211|R500)_(SHADE|SUN)_([1-3])_(INTERNODE|LEAF|ROOT|SEEDLING|SILIQUE|APICALMERISTEM|FLORALMERISTEM).[[:print:]]+(.bam)","\\4",names(counts_GC_final_merged))
)
samples

samples$group <- paste(samples$gt, samples$trt, samples$tiss, sep = "_")
group <- samples$group
gt <- samples$gt
trt <- samples$trt
tiss <- samples$tiss
samples

# creating DGEList and normalizing (edgeR)
data_GC_int <- DGEList(counts_GC_final_merged,group=group)
#data_GC_int <- DGEList(counts_GC_final_merged,group = rep(c("IMB211","R500"),c(5,6))) # same as above...
data_GC_int$samples
data_GC_int <- data_GC_int[rowSums(cpm(data_GC_int) > 1) >= 3,] # [1] 34395    64 (7,525 genes lost because of cut-off)
dim(data_GC_int)
data_GC_int <- calcNormFactors(data_GC_int)

# normalized counts
countspermi <- cpm(data_GC_int, normalized.lib.sizes=TRUE)
head(countspermi)[,1:10]
class(countspermi)
colnames(countspermi)

# Clustering (from edx course)
clst <- countspermi
x <- t(clst)
tiss1 <- data.frame(tiss)
te1 <- cbind(x, tiss1)
head(te1)[1:10]

#distance calculation
d <- dist(te1) # takes few minutes 
head(d)
class(d)
as.matrix(d)[1:2, 1:2]

# Heirarchial clustering
hc <- hclust(d) # very quick
hc
plot(hc, labels=te1$tiss) # Not so good..
pdf(file = "cluster_Dendrogram_GC.pdf", height =7, width =12)
myplclust(hc, labels=te1$tiss, lab.col=as.numeric(te1$tiss))
dev.off()

hclusters <- cutree(hc, h=40000)
table(true=te1$tiss, cluster=hclusters) ##heriarchial clusters (may be ok?)

# kmeans clustering for only two genes
# Look at two genes for all 66 samples
plot(x[,1], x[,2])
# k-means with two genes and three clusters
km <- kmeans(x[,1:2], centers=3)
names(km)
plot(x[,1], x[,2], col=km$cluster, pch=16)
# k-means with all genes and three clusters
km <- kmeans(x, centers=3)
plot(x[,1], x[,2], col=km$cluster, pch=16)
# k-means with all genes and seven clusters
km <- kmeans(x, centers=6)
plot(x[,1], x[,2], col=km$cluster, pch=16)
table(true=te1$tiss, cluster=km$cluster) ##km clusters

# mds plots
mds <- cmdscale(dist(x))
plot(mds, col=km$cluster, pch=16)  
plot(mds, type="n")
text(mds, colnames(clst), col=km$cluster)
plot(mds, type="n")
text(mds, colnames(clst), col=as.numeric(te1$tiss))

# heatmaps
#install colour package
#install.packages("RColorBrewer")
library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# Now pick the genes with the top variance over all samples:
library("genefilter")
test1 <- gsub(".1.merged.bam", "", colnames(countspermi))
head(test1)
colnames(countspermi) <- test1
head(countspermi)
rv <- rowVars(countspermi)
head(rv)
idx <- order(-rv)[1:40]
idx

#normal heatmap
pdf(file = "normal_heatmap.pdf", height = 10, width = 7)
heatmap(countspermi[idx,],col=hmcol)
dev.off()

# heatmap2
#install.packages("gplots")
library("gplots")
cols <- palette(brewer.pal(7, "Dark2"))[samples$tiss]
length(cols)
cbind(colnames(countspermi),cols)
pdf(file = "colored_heatmap.pdf", height = 10, width = 8)
heatmap.2(countspermi[idx,],trace="none",ColSideColors=cols,col=hmcol)


heatmap.2(countspermi[idx,],col=redgreen(64), scale = "none", ColSideColors=col, 
                key = FALSE, symkey = TRUE, density.info = "none", margins=c(12,8), 
                cexCol=1, trace="none",srtCol=45)
dev.off()


# heatmap3
pdf(file = "All_tissues_no_log_scale_genefiltered_heatmap.pdf", height = 7, width = 12)
heatmap(cor(countspermi, method = "spearman"), main = "All_tissues_no_log_scale_genefiltered_heatmap")
dev.off()

# Save the normalized counts file
write.csv(countspermi, file = "clustering_all_normalized_counts.csv")

# Running heatmap and hopach for clustering of arrays (samples) with mean of the replicates....
  
# Read in the data file
countspermi2 <- read.csv("clustering_all_normalized_counts.csv", sep = ",")
head(countspermi2[1:10])
dim(countspermi2) # [1] 34395    65

# Convert to log scale
rownames(countspermi2) <- countspermi2$X 
head(countspermi2)
data <- log2(countspermi2[,-1])
head(data[1:10])
summary(data[1:10])
#rownames(data)
# heatmaps
pdf(file = "All_tissues_log_scale_genefiltered_heatmap.pdf", height = 7, width = 12)
heatmap(cor(data, method = "spearman")) #almost same as with normalzied counts...
dev.off()

# Prepare file names
data2 <- data
names(data2)
data2$ITAG <- factor(rownames(data2))
head(data2[1:10])
dim(data2)

# Melted dataframe
data3 <- melt(data2,id.var="ITAG")
head(data3)
data3$value[data3$value==-Inf] <- 0
summary(data3)
data3$variable <- gsub("_[1-3]_", "_", data3$variable) # modified on 10/16/13 to remove the replicate number and replace it with "_" to calculate mean...
head(data3)

# Calculating average mean of all samples
mean.table <- tapply(data3$value,list(data3$ITAG,data3$variable),mean)
head(mean.table)
dim(mean.table)

#writing the mean table
write.csv(file = "all_tissues_mean_table.csv", mean.table)

#Scaling and Centering of Matrix-like Objects 
scaled.table <- t(scale(t(mean.table)))
summary(scaled.table)
scaled.table <- scaled.table[complete.cases(scaled.table),]
head(scaled.table)
#mean.table <- mean.table[apply(mean.table,1,sum)!=0,]

#heatmap for clustering here..
pdf(file = "all_tissues_scaled_heatmap.pdf", height = 7, width = 12)
heatmap(cor(scaled.table, method = "spearman") )
dev.off()

#genefilter before doing heatmap..
#biocLite("genefilter")
f1 <- pOverA(0.25, log2(100))
f2 <- function(x) (IQR(x) > 0.5)
f3 <- filterfun(f1, f2)
selected <- genefilter(mean.table, f3)
sum(selected)
scaled.table.gf <- scaled.table[selected,]
nrow(scaled.table.gf)# 2300
#head(mean.table.gf)
pdf(file = "all_tissues_scaled_genefiltered_heatmap.pdf", height = 7, width = 12)
heatmap(cor(scaled.table.gf, method = "spearman")) # Pretty Good (Use this for heatmap).
dev.off()

# Hopach clustering...
# Run distance matrix first
scaled.dist <- distancematrix(scaled.table.gf, "cosangle", na.rm = TRUE)
dim(scaled.dist)

# Now run hopach for clustering of genes
scaled.dist.hopach <- hopach(scaled.table.gf, dmat = scaled.dist) # takes a while
scaled.dist.hopach$clust$k # The hopach algorithm identiﬁes 111 gene clusters
table(scaled.dist.hopach$clust$sizes) # Most of the clusters are between 1 to 5 genes.
pdf(file = "All_tissues_scaled_genefiltered_gene_cluster_hopach.pdf")
dplot(scaled.dist, scaled.dist.hopach, col=heat.colors(12), showclusters=F, ord = "final", main = "All_tissues_scaled_genefiltered_gene_cluster_hopach")
dev.off()

# Now run hopach for clustering of arrays (samples)
scaled.array.hobj <- hopach(t(scaled.table.gf), d = "euclid")
scaled.array.hobj$clust$k # 7 clusters
pdf(file = "All_tissues_scaled_genefiltered_sample_cluster_hopach.pdf", height = 5, width = 8)
dplot(distancematrix(t(scaled.table.gf), d = "euclid"), scaled.array.hobj,labels = colnames(scaled.table.gf), main = "HOPACH - distance = correlation")
dev.off() # Good

# HOPACH hierarchical clustering in MapleTree
hopach2tree(scaled.table.gf, file = "all_tissues_kidneyTree", hopach.genes = scaled.dist.hopach, hopach.arrays = scaled.array.hobj, dist.genes = scaled.dist)

# Running heatmap and hopach for clustering of arrays (samples) with all the replicates included....
  
# Read in the data file
countspermi2 <- read.csv("clustering_all_normalized_counts.csv", sep = ",")#
head(countspermi2)
dim(countspermi2)

# Convert to log scale
rownames(countspermi2) <- countspermi2$X 
head(countspermi2)
data <- log2(countspermi2[,-1])
head(data)
class(data)
str(data)
summary(data)
#rownames(data)
# heatmaps
heatmap(cor(data, method = "spearman")) #almost same as with normalzied counts...

# Prepare file names
data2 <- data
names(data2)
data2$ITAG <- factor(rownames(data2))
head(data2)
dim(data2)

# Melted dataframe
data3 <- melt(data2,id.var="ITAG")
head(data3)
data3$value[data3$value==-Inf] <- 0
summary(data3)

# Calculating average mean of all samples
mean.table <- tapply(data3$value,list(data3$ITAG,data3$variable),mean)
head(mean.table)
dim(mean.table)

#Scaling and Centering of Matrix-like Objects 
scaled.table <- t(scale(t(mean.table)))
head(scaled.table)
summary(scaled.table)
scaled.table <- scaled.table[complete.cases(scaled.table),]
mean.table <- mean.table[apply(mean.table,1,sum)!=0,]

#heatmap for clustering here..
pdf(file = "all_tissues_all_replicates_scaled_heatmap.pdf", height = 7, width = 12)
heatmap(cor(scaled.table, method = "spearman") ) # Good
dev.off()

#genefilter before doing heatmap (to remove variant genes)..
#biocLite("genefilter")
f1 <- pOverA(0.25, log2(100))
f2 <- function(x) (IQR(x) > 0.5)
f3 <- filterfun(f1, f2)
selected <- genefilter(mean.table, f3)
sum(selected)
scaled.table.gf <- scaled.table[selected,]
nrow(scaled.table.gf)
#head(mean.table.gf)
pdf(file = "all_tissues_all_replicates_scaled_genefiltered_heatmap.pdf", height = 7, width = 12)
heatmap(cor(scaled.table.gf, method = "spearman")) # Good...
dev.off()

# Trying hopach...
# Run distance matrix first
scaled.dist <- distancematrix(scaled.table.gf, "cosangle", na.rm = TRUE)
dim(scaled.dist)

# Now run hopach for clustering of arrays (samples)
scaled.array.hobj <- hopach(t(scaled.table.gf), d = "euclid")
scaled.array.hobj$clust$k # 8 clusters
pdf(file = "All_tissues_all_replicates_scaled_genefiltered_sample_cluster_hopach.pdf", height = 7, width = 12)
dplot(distancematrix(t(scaled.table.gf), d = "euclid"), scaled.array.hobj,labels = colnames(scaled.table.gf), main = "All_tissues_scaled_genefiltered_sample_cluster_hopach") # Good 
dev.off()

# Trying clara clustering (Julin's)
  
# creating an empty df for different number of clusters and metrics
#many clusters?  R says the number of clusters with the largest average silhouette width
  
results <- data.frame(k=4:30,sil.euc=NA,sil.manh=NA)
head(results)

# Loop to generate values for the metrics
for (k in results$k) {
  clara.e <- clara(scaled.table,k=k,metric="euclidean",sample=100)
  clara.m <- clara(scaled.table,k=k,metric="manhattan",sample=100)
  results$sil.euc[results$k==k] <- clara.e$silinfo$avg.width
  results$sil.manh[results$k==k] <- clara.m$silinfo$avg.width
  #print(plot(clara.e,which=2))
  #print(plot(clara.m,which=2))
  #plot.cluster(clara.e)
  #plot.cluster(clara.m)
}

results
max(results[,2])
max(results[,3])

#not so different really.  k = 5 looks good

clara.e <- clara(scaled.table,k=5,metric="euclidean",sample=100)
clara.m <- clara(scaled.table,k=5,metric="manhattan",sample=100)

plot(clara.e,which=2)
plot(clara.m,which=2)

plot.cluster <- function(clara.obj,ncol=3,title="") {
  k <- clara.obj$call$k
  print(k)
  if (k < ncol) ncol <- 4
  metric <- clara.obj$call$metric
  #medoids <- melt(clara.obj$medoids,id.var=factor(rownames(clara.obj$medoids)))
  #names(medoids) <- c("ITAG","species","value")
  #ggplot(data=medoids,aes(x=species,y=value)) + geom_line(aes(group=ITAG)) + geom_point() + facet_wrap(~ITAG,ncol=4)
  
  #make data frame which gives cluster number and (scaled) expression values for each ITAG
  clusterdata <- data.frame(clara.obj$data,ITAG=rownames(clara.obj$data),cluster=clara.obj$cluster)
  
  #use melt to reshape this
  clara.melt <- melt(clusterdata,id.vars=c("ITAG","cluster"),variable.name="species")
  clara.melt$species <- factor(clara.melt$species,levels=c("IMB211_SHADE_INTERNODE","IMB211_SHADE_LEAF","IMB211_SHADE_ROOT","IMB211_SHADE_SEEDLING",
                                                           "IMB211_SHADE_SILIQUE","IMB211_SHADE_APICALMERISTEM","IMB211_SUN_INTERNODE","IMB211_SUN_LEAF",
                                                           "IMB211_SUN_ROOT","IMB211_SUN_SEEDLING","IMB211_SUN_SILIQUE","R500_SHADE_APICALMERISTEM", 
                                                           "R500_SHADE_FLORALMERISTEM", "R500_SHADE_INTERNODE", "R500_SHADE_LEAF", "R500_SHADE_ROOT", 
                                                           "R500_SHADE_SEEDLING", "R500_SHADE_SILIQUE", "R500_SUN_APICALMERISTEM", "R500_SUN_INTERNODE",
                                                           "R500_SUN_LEAF", "R500_SUN_ROOT", "R500_SUN_SEEDLING", "R500_SUN_SILIQUE"))
  
  #sort so that plotting works
  clara.melt <- clara.melt[order(clara.melt$ITAG,clara.melt$species),]
  
  p <- ggplot(data=clara.melt,aes(x=species,y=value)) + 
    facet_wrap(~cluster,ncol) + geom_line(aes(group=ITAG),alpha=.005) +
    ggtitle(label=paste(title," k = ",k," metric = ",metric,sep=""))
  #ggsave(filename = label, plot = p, width = 5, height = 7)
  print(p)
}

pdf(file = "ALL tissues clara.e k = 5 metric = euclidean", width = 10, height=5)
plot.cluster(clara.e, ncol=5, title="ALL tissues clara.e")
dev.off()

pdf(file = "clustering/ALL tissues clara.e k = 5 metric = manhattan", width = 10, height=5)
plot.cluster(clara.m, ncol=4, title="INTERNODE clara.m")
dev.off()

# Can't cluster all tissues this way..instead use heatmap
# heatmap (hierachial clustering) based on normalized counts
countspermi_test <- countspermi
colnames(countspermi_test) <- sub("(GC|GH|F)(IMB211|R500)_(SUN|SHADE)_([1-3])_(INTERNODE|LEAF|ROOT|SEEDLING|SILIQUE|APICALMERISTEM|FLORALMERISTEM)_[[:print:]]+(.bam)","\\1\\2_\\3_\\4_\\5",colnames(countspermi_test))
pdf(file = "all_tissues_heatmap.pdf", height = 7, width = 12)
heatmap(cor(countspermi_test, method = "spearman"))
dev.off()

------------------------------------------------------------------------------------------------------------------------------------
##### PCA and SOM##########

# PCA analysis
pca <- prcomp(scaled.table.gf, scale=T)
plot(pca)
summary(pca) # The first two PC's account for 50% of the variation..
print(pca$rotation)
pdf(file = "pca_biplot.pdf", height = 7, width = 7)
biplot(pca, scale = TRUE) # good separation of samples # takes a few minutes
dev.off()

# all PC's
pcr <- print(pca$rotation)
pcr1 <- melt(pcr)
head(pcr1)

pcr1$gt <- sub("(IMB211|R500)_(SUN|SHADE)_(INTERNODE|LEAF|ROOT|SEEDLING|SILIQUE|APICALMERISTEM|FLORALMERISTEM)","\\1",pcr1[,1])
pcr1$trt <- sub("(IMB211|R500)_(SUN|SHADE)_(INTERNODE|LEAF|ROOT|SEEDLING|SILIQUE|APICALMERISTEM|FLORALMERISTEM)","\\2",pcr1[,1])
pcr1$tiss <- sub("(IMB211|R500)_(SUN|SHADE)_(INTERNODE|LEAF|ROOT|SEEDLING|SILIQUE|APICALMERISTEM|FLORALMERISTEM)","\\3",pcr1[,1])

head(pcr1)
tail(pcr1)

p <- ggplot(pcr1, aes(Var1, value)) + geom_bar(stat = "identity") + facet_grid(Var2~., scales="free") + theme(axis.text = element_text(size = 9), axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))
p <- p + guides(fill=FALSE)
p
ggsave(filename = "pcr.pdf", plot = p, width = 10, height = 20)

# To look to see which pc correspond to genotype
# For genotype
p <- ggplot(pcr1, aes(Var1, value, fill = gt)) + geom_bar(stat = "identity") + facet_grid(Var2~., scales="free") + theme(axis.text = element_text(size = 9), axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))
p <- p + guides(fill=FALSE)
p
ggsave(filename = "pcr_gt.pdf", plot = p, width = 10, height = 20)
# PC3 = genotype

# For tissue, we first should divide tissues into meristem and non meristem
head(pcr1)
test1 <- pcr1[grep("MERISTEM", pcr1$tiss),]
head(test1)
tail(test1)
test1$type <- "MERISTEM"

test2 <- pcr1[!grepl("MERISTEM", pcr1$tiss),]
head(test2)
tail(test2)
test2$type <- "NONMERISTEM"

test_comb <- rbind(test1, test2)
head(test_comb)
tail(test_comb)

p <- ggplot(test_comb, aes(Var1, value, fill = type)) + geom_bar(stat = "identity") + facet_grid(Var2~., scales="free") + theme(axis.text = element_text(size = 9), axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))
p <- p + guides(fill=FALSE)
p
ggsave(filename = "pcr_tiss_type.pdf", plot = p, width = 10, height = 20)
# PC2 = meristem vs nonmeristem

####
# all four PC's together - not doing now.....
pca_result <- data.frame(row.names=rownames(pca$rotation),PC1 = pca$rotation[,1], PC2 = pca$rotation[,2], PC3 = pca$rotation[,3], PC12 = pca$rotation[,12])
#PC1 = Tissue (Main and Lateral organs)
#PC2 = Tissue (Meristematic and Non meristematic)
#PC3 = Genotype
#PC12 = Treatment?
pca_result$samples <- rownames(pca$rotation)
head(pca_result)
pca_result_melt <- melt(pca_result)
head(pca_result_melt)
tail(pca_result_melt)
ggplot(pca_result_melt, aes(samples, value)) + geom_bar(stat = "identity") + facet_grid(variable~., scales="free") + theme(axis.text = element_text(size = 7), axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))
ggsave(filename = "pca_result_melt.pdf", plot = last_plot(), width = 4, height = 7)
pca_result_melt$gt <- sub("(IMB211|R500)_(SUN|SHADE)_(INTERNODE|LEAF|ROOT|SEEDLING|SILIQUE|APICALMERISTEM|FLORALMERISTEM)","\\1",pca_result_melt[,1])
pca_result_melt$trt <- sub("(IMB211|R500)_(SUN|SHADE)_(INTERNODE|LEAF|ROOT|SEEDLING|SILIQUE|APICALMERISTEM|FLORALMERISTEM)","\\2",pca_result_melt[,1])
pca_result_melt$tiss <- sub("(IMB211|R500)_(SUN|SHADE)_(INTERNODE|LEAF|ROOT|SEEDLING|SILIQUE|APICALMERISTEM|FLORALMERISTEM)","\\3",pca_result_melt[,1])
head(pca_result_melt)

# By genotype
p <- ggplot(pca_result_melt, aes(samples, value, fill = gt)) + geom_bar(stat = "identity") + facet_grid(variable~., scales="free") + theme(axis.text = element_text(size = 9), axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))
p <- p + guides(fill=FALSE)
p # PC3 is almost confirmed as Genotype

# By treatment (not sure?)
p1 <- ggplot(pca_result_melt, aes(samples, value, fill = trt)) + geom_bar(stat = "identity") + facet_grid(variable~., scales="free") + theme(axis.text = element_text(size = 9), axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))
p1 <- p1 + guides(fill=FALSE)
p1

# By tissue
p2 <- ggplot(pca_result_melt, aes(samples, value, fill = tiss)) + geom_bar(stat = "identity") + facet_grid(variable~., scales="free") + theme(axis.text = element_text(size = 9), axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))
p2 <- p2 + guides(fill=FALSE)
p2


# Individual PC's

# PC1
pca_1_result <- data.frame(row.names=rownames(pca$rotation),PC1 = pca$rotation[,1])
pca_1_result$samples <- rownames(pca$rotation)
head(pca_1_result)
pca_result_1_melt <- melt(pca_1_result)
head(pca_result_1_melt)
pca_result_1_melt$tis <- sub("(IMB211|R500)_(SUN|SHADE)_(INTERNODE|LEAF|ROOT|SEEDLING|SILIQUE|APICALMERISTEM|FLORALMERISTEM)","\\3",pca_result_1_melt[,1])
ggplot(pca_result_1_melt, aes(samples, value, colour = tis, fill = tis)) + geom_bar(stat = "identity") + theme(axis.text = element_text(size = 7), axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))
ggsave(filename = "pca_result_1_melt.pdf", plot = last_plot(), width = 12, height = 7)

# Removing A.M and F.M samples and replotting PC1
#te <- pca_1_result[grep("APICALMERISTEM|FLORALMERISTEM", pca_1_result$samples, invert = TRUE),]
#te_1_melt <- melt(te)
#head(te_1_melt)
#te_1_melt$tis <- sub("(GC|GH|F)(IMB211|R500)_(SUN|SHADE)_(INTERNODE|LEAF|ROOT|SEEDLING|SILIQUE)","\\4",te_1_melt[,1])
#ggplot(te_1_melt, aes(samples, value, colour = tis, fill = tis)) + geom_bar(stat = "identity")
#ggsave(filename = "te_1_melt.pdf", plot = last_plot(), width = 12, height = 7)

#INT <- te_1_melt[grep("INTERNODE", te_1_melt$tis),]
#mean(INT[,3])
#LEA <- te_1_melt[grep("LEA", te_1_melt$tis),]
#mean(LEA[,3])
#ROOT <- te_1_melt[grep("ROOT", te_1_melt$tis),]
#mean(ROOT[,3])
#SEED <- te_1_melt[grep("SEEDLING", te_1_melt$tis),]
#mean(SEED[,3])
#SIL <- te_1_melt[grep("SILIQUE", te_1_melt$tis),]
#mean(SIL[,3])
#comb <- cbind(mean(INT[,3]), mean(LEA[,3]), mean(ROOT[,3]), mean(SEED[,3]), mean(SIL[,3]))
#colnames(comb) <- c("INTERNODE", "LEAF", "ROOT", "SEEDLING", "SILIQUE")
#comb
#comb_mlet <- melt(comb)
#comb_mlet
#c1 <- ggplot(comb_mlet, aes(Var2, value, fill = Var2)) + geom_bar(stat = "identity") + guides(fill=FALSE) + xlab(NULL)
#ggsave(file = "pc1.pdf", plot = c1)

# Using tapply for the same..and also using mean of the samples irresptive of gt
# Taking the mean of tissues irrespective of gt.
comb_mlet_1 <- melt(tapply(pca_result_1_melt$value, pca_result_1_melt$tis, FUN = mean))
comb_mlet_1
c1 <- ggplot(comb_mlet_1, aes(Var1, value, fill = Var1)) + geom_bar(stat = "identity") + guides(fill=FALSE) + xlab(NULL) + theme(axis.text = element_text(size = 10), axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))
c1
ggsave(file = "pc1.pdf", plot = c1)

# PC1 explains horizontal and vertical separation of tissue specificity

# PC2
pca_2_result <- data.frame(row.names=rownames(pca$rotation),PC2 = pca$rotation[,2])
pca_2_result$samples <- rownames(pca$rotation)
head(pca_2_result)
pca_result_2_melt <- melt(pca_2_result)
pca_result_2_melt$tis <- sub("(GC|GH|F)(IMB211|R500)_(SUN|SHADE)_(INTERNODE|LEAF|ROOT|SEEDLING|SILIQUE|APICALMERISTEM|FLORALMERISTEM)","\\4",pca_result_2_melt[,1])
ggplot(pca_result_2_melt, aes(samples, value, colour = tis, fill = tis)) + geom_bar(stat = "identity") + theme(axis.text = element_text(size = 7), axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5)) 
ggsave(filename = "pca_result_2_melt.pdf", plot = last_plot(), width = 12, height = 7)

# Taking the mean of tissues irrespective of gt.
comb_mlet_2 <- melt(tapply(pca_result_2_melt$value, pca_result_2_melt$tis, FUN = mean))
comb_mlet_2
c2 <- ggplot(comb_mlet_2, aes(Var1, value, fill = Var1)) + geom_bar(stat = "identity") + guides(fill=FALSE) + xlab(NULL) + theme(axis.text = element_text(size = 15), axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))
ggsave(file = "pc2.pdf", plot = c2)

#PC3
pca_3_result <- data.frame(row.names=rownames(pca$rotation),PC3 = pca$rotation[,3])
pca_3_result$samples <- rownames(pca$rotation)
head(pca_3_result)
pca_result_3_melt <- melt(pca_3_result)
pca_result_3_melt$gt <- sub("(GC|GH|F)(IMB211|R500)_(SUN|SHADE)_(INTERNODE|LEAF|ROOT|SEEDLING|SILIQUE|APICALMERISTEM|FLORALMERISTEM)","\\2",pca_result_3_melt[,1])
ggplot(pca_result_3_melt, aes(samples, value, colour = gt, fill = gt)) + geom_bar(stat = "identity") + theme(axis.text = element_text(size = 7), axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))
ggsave(filename = "pca_result_3_melt.pdf", plot = last_plot(), width = 12, height = 7)

# Removing A.M and F.M samples and replotting PC3
te3 <- pca_3_result[grep("APICALMERISTEM|FLORALMERISTEM", pca_3_result$samples, invert = TRUE),]
te_3_melt <- melt(te3)
head(te_3_melt)
te_3_melt$tis <- sub("(GC|GH|F)(IMB211|R500)_(SUN|SHADE)_(INTERNODE|LEAF|ROOT|SEEDLING|SILIQUE)","\\4",te_3_melt[,1])
te_3_melt
comb_mlet_3 <- melt(tapply(te_3_melt$value, te_3_melt$tis, FUN = mean))
c3 <- ggplot(comb_mlet_3, aes(Var1, value, fill = Var1)) + geom_bar(stat = "identity") + guides(fill=FALSE) + xlab(NULL)
ggsave(file = "pc3.pdf", plot = c3)


#PC12
pca_12_result <- data.frame(row.names=rownames(pca$rotation),PC12 = pca$rotation[,12])
pca_12_result$samples <- rownames(pca$rotation)
head(pca_12_result)
pca_result_12_melt <- melt(pca_12_result)
pca_result_12_melt$trt <- sub("(GC|GH|F)(IMB211|R500)_(SUN|SHADE)_(INTERNODE|LEAF|ROOT|SEEDLING|SILIQUE|APICALMERISTEM|FLORALMERISTEM)","\\3",pca_result_12_melt[,1])
ggplot(pca_result_12_melt, aes(samples, value, colour = trt, fill = trt)) + geom_bar(stat = "identity") + theme(axis.text = element_text(size = 7), axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))
ggsave(filename = "pca_result_12_melt.pdf", plot = last_plot(), width = 12, height = 7)


#summary(pca)$importance[, 1:6]
#mycolors <- c("red", "green", "blue", "magenta", "black")
#plot(pca$x, pch=20, col=mycolors[sort(rep(1:5, 500))]) 
#plot(pca$x, type="n"); text(pca$x, rownames(pca$x), cex=0.8, col=mycolors[sort(rep(1:5, 500))]) 
#library(geneplotter)
#smoothScatter(pca$x)
#library(scatterplot3d)
#library(rgl); rgl.open(); offset <- 50; par3d(windowRect=c(offset, offset, 640+offset, 640+offset)); rm(offset); rgl.clear(); rgl.viewpoint(theta=45, phi=30, fov=60, zoom=1); spheres3d(pca$x[,1], pca$x[,2], pca$x[,3], radius=0.3, color=mycolors, alpha=1, shininess=20); aspect3d(1, 1, 1); axes3d(col='black'); title3d("", "", "PC1", "PC2", "PC3", col='black'); bg3d("white") 

##########

# No genotype or treatment effect here...

# Use the mean table for performing PCA and SOM #
data_R500_SUN <- mean.table[,grep("R500_SUN*", colnames(mean.table))]
head(data_R500_SUN)
tail(data_R500_SUN)
dim(data_R500_SUN)
data_R500_SUN <- data_R500_SUN[,-1] # Remove APICAL MERISTEM tissue?
data_R500_SUN <- cbind(data_R500_SUN,SP = rep("R500",34395), trt = rep("SUN", 34395))

data_R500_SHADE <- mean.table[,grep("R500_SHADE*", colnames(mean.table))]
head(data_R500_SHADE)
dim(data_R500_SHADE)
data_R500_SHADE <- data_R500_SHADE[,-c(1,2)]
data_R500_SHADE <- cbind(data_R500_SHADE,SP = rep("R500",34395), trt = rep("SHADE", 34395))


data_IMB211_SUN <- mean.table[,grep("IMB211_SUN*", colnames(mean.table))]
head(data_IMB211_SUN)
dim(data_IMB211_SUN)
data_IMB211_SUN <- cbind(data_IMB211_SUN,SP = rep("IMB211",34395), trt = rep("SUN", 34395))


data_IMB211_SHADE <- mean.table[,grep("IMB211_SHADE*", colnames(mean.table))]
head(data_IMB211_SHADE)
dim(data_IMB211_SHADE)
data_IMB211_SHADE <- data_IMB211_SHADE[,-1]
data_IMB211_SHADE <- cbind(data_IMB211_SHADE,SP = rep("IMB211",34395), trt = rep("SHADE", 34395))

data_R500_IMB211 <- rbind(data_R500_SUN, data_R500_SHADE, data_IMB211_SUN, data_IMB211_SHADE)
head(data_R500_IMB211)
tail(data_R500_IMB211)
dim(data_R500_IMB211)

#write the combined file
write.csv(data_R500_IMB211, file = "R500_IMB211.csv")

data_2 <- read.csv("R500_IMB211.csv", h = T)
head(data_2)
colnames(data_2) <- c("X", "INTERNODE", "LEAF", " ROOT", "SEEDLING", " SILIQUE", "SP", "trt")
tail(data_2)
dim(data_2)
str(data_2)
data_2$Mean <- rowMeans(data_2[,2:6])
head(data_2)
data_2$sd <- apply(data_2[,2:6],1,sd)
data_2$cv <- data_2[,10]/data_2[,9]
head(data_2)

## Select only those genes in the 5% CV
quantile(data_2$cv, probs=seq(.05,.95,.05), na.rm=T) 
data_2.data <- subset(data_2, cv > 3.65941510) # 95%
head(data_2.data)
length(data_2.data$X) # 6805
summary(data_2.data$SP)
#IMB211   R500 
#3201   3604 

## Create a matrix of the data to perform a PCA on and scale it
m.sub.data <- as.matrix(data_2.data[2:6])
head(m.sub.data)
sc.sub.data <- t(scale(t(m.sub.data)))
head(sc.sub.data)
tisdata <- as.matrix(sc.sub.data, dimnames=list(rownames(X)) )
head(tisdata)

## Perform the PCA first
tispca <- prcomp(tisdata, scale=TRUE)
#head(tispca)[1:10]
summary(tispca)
#Importance of components:
#  PC1    PC2    PC3    PC4       PC5
#Standard deviation     1.2112 1.1373 1.0877 1.0279 2.603e-15
#Proportion of Variance 0.2934 0.2587 0.2366 0.2113 0.000e+00
#Cumulative Proportion  0.2934 0.5521 0.7887 1.0000 1.000e+00

# rotation
tispca$rotation
              #PC1        PC2         PC3         PC4        PC5
#INTERNODE -0.12064390  0.1224193 -0.89810772  0.06890989 -0.3988891
#LEAF       0.29229137 -0.4822383  0.11722189  0.72648147 -0.3748275
#ROOT     -0.73979930  0.1625866  0.36524261  0.06975596 -0.5366512
#SEEDLING   0.08696526 -0.6384217  0.02142257 -0.66066485 -0.3846009
#SILIQUE   0.58749030  0.5643076  0.21401041 -0.16165443 -0.5142757

## Retrieve PCA scores
tis.pca.scores <- data.frame(tispca$x) # what is "x" here?
head(tis.pca.scores)

## Write out master data files with original data, scaled data, and PCA results
data.val <- cbind(data_2.data,sc.sub.data,tis.pca.scores)
head(data.val)

write.table(data.val, file="tissue.pca.scores.txt")
write.table(tispca$rotation, "loadings.txt")

#PCA circular plot 

#PC1 and PC2
# without theme background
t <- ggplot(data.val, aes(PC1, PC2)) + geom_point(color = "black", size = 1.0, alpha = 0.4) + theme_bw() 
t
# with theme background
t1 <- ggplot(data.val, aes(PC1, PC2)) + geom_point(color = "black", size = 1.0, alpha = 0.4) 
t1

#PC3 and PC4
tpc34 <- ggplot(data.val, aes(PC3, PC4)) + geom_point(color = "black", size = 1.0, alpha = 0.4) + theme_bw()
tpc34

ggsave(filename = "plain.PC.no.background.pdf", plot = t, height = 5, width = 5) #there are three clusters
ggsave(filename = "plain.PC.background.pdf", plot = t1, height = 5, width = 5) #there are three clusters
ggsave(filename = "plain.PC34.no.background.pdf", plot = tpc34, height = 5, width = 5) #there are three clusters


## run SOM now....
library(class)
library(MASS)
source("http://bioconductor.org/biocLite.R")
biocLite("kohonen")
library(kohonen)

## Read in data
data <- read.table("tissue.pca.scores.txt")
head(data)
names(data)
dim(data)
#colnames(data)[1:21] <- ("X", "INTERNODE.ns",  "LEAF.ns",  "ROOT.ns", "SEEDLING.ns", "SILIQUE.ns", "SP", "trt", "Mean", "sd", "cv", "INTERNODE", "LEAF", "ROOT", "SEEDLING", "SILIQUE", "PC1", "PC2", "PC3", "PC4","PC5")

## Create a matrix and scale data (data is already scaled, so nothing happens actually besides making a matrix)
m.data <- as.matrix(data[12:16])
sc.data <- t(scale(t(m.data))) # no need but you can repeat....
head(m.data)
head(sc.data)
colnames(sc.data) <- c("INTERNODE", "LEAF", "ROOT", "SEEDLING", "SILIQUE")

## Set a random seed so that SOM results are reproducible
set.seed(2)

## Perform SOM
ssom <- som(sc.data, somgrid(1,8,"hexagonal")) # try 1X8 with my data
summary(ssom)
# som map of size 1x8 with a hexagonal topology.
# Training data included; dimension is 7220 by 5
# Mean distance to the closest unit in the map: 1.031501

plot(ssom, type ="changes") # good

pdf(file = "SOM.clusters.pdf", height = 5, width = 5)
plot(ssom, type = "codes")
dev.off()

pdf(file = "SOM.counts.pdf", height = 5, width = 5)
plot(ssom, type = "counts")
dev.off()

pdf(file = "SOM.quality.pdf", height = 5, width = 5)
plot(ssom, type = "quality")
dev.off()

## Create and write-out master SOM file
data.val2 <- cbind(data, ssom$unit.classif, ssom$distances)
head(data.val2)
write.table(data.val2, file="supersom.data.c3.txt")

## Codes for the SOM nodes
codes <- ssom$codes
head(codes)
write.table(codes, file="codes.txt")

## Visualization
data <- read.table("supersom.data.c3.txt",header=TRUE)
names(data) <- c("X", "INTERNODE.ns",  "LEAF.ns",  "ROOT.ns", "SEEDLING.ns", "SILIQUE.ns", "SP", "trt", "Mean", "sd", "cv", "INTERNODE", "LEAF", "ROOT", "SEEDLING", "SILIQUE", "PC1", "PC2", "PC3", "PC4","PC5", "ssom.unit.classif", "ssom.distances")
tail(data)

## PC graphs
library(plyr)

# Circular plot for all SOM clusters projected on PC space (no tissues known here)

# For PC1 and PC2
t2 <- ggplot(data, aes(PC1, PC2))+ geom_point(size=1.0, aes(colour=factor(ssom.unit.classif))) + 
      theme_bw() + scale_colour_manual(values=c("mediumseagreen","black", "mediumorchid4", "blue", "red", "violet", "yellow", "pink")) + theme(legend.position = "none")

t2

t34 <- ggplot(data, aes(PC3, PC4))+ geom_point(size=1.0, aes(colour=factor(ssom.unit.classif))) + 
      theme_bw() + scale_colour_manual(values=c("mediumseagreen","black", "mediumorchid4", "blue", "red", "violet", "yellow", "pink")) + theme(legend.position = "none")

t34

ggsave(file = "node.PC.c3.pdf", plot = t2, height = 5, width = 5)

ggsave(file = "node.PC34.c3.pdf", plot = t34, height = 5, width = 5)


# For Species
t3 <- ggplot(data, aes(PC1, PC2))+ geom_point(alpha=0.60,aes(colour=factor(SP))) + theme_bw() + 
      scale_colour_manual(name = "genotypes", values=c("red","green"))

ggsave(file = "node.PC.c3.species.pdf", plot = t3, height = 8, width = 8)


## Boxplots for all SOM clusters for 5 tissues 
expression <- data[c(1,12:16,22)]
head(expression)
m.expression <- melt(expression, id = c("X", "ssom.unit.classif"))
head(m.expression)
tail(m.expression)
summary(m.expression)

# no jitter
bp <- ggplot(m.expression, aes(x=variable, y=value))+ 
      geom_boxplot(outlier.size=0, alpha=0.8) + 
      facet_grid(~ssom.unit.classif) +
      theme(text = element_text(size=10)) + theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

bp
ggsave(file = "boxplot_all_som_clusters_no_jitter.pdf", plot = bp, height = 7, width = 10)


# with jitter

bpj <- ggplot(m.expression, aes(x=variable, y=value))+ 
       geom_point(position="jitter",size=1.5,alpha=0.4) +
       geom_boxplot(outlier.size=0, alpha=0.8) + 
       facet_grid(~ssom.unit.classif) +
       theme(text = element_text(size=10)) + theme_bw() +
       theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

bpj

ggsave(file = "boxplot_all_som_clusters_jitter.pdf", plot = bpj, height = 7, width = 10)

# with no outline

bpjp <- ggplot(m.expression, aes(x=variable, y=value))+ 
       geom_point(position="jitter",size=1.5,alpha=0.4) + 
       facet_grid(~ssom.unit.classif) +
       theme(text = element_text(size=10)) + theme_bw() +
       theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

bpjp

ggsave(file = "point_graph_all_som_clusters_no_jitter_no_line.pdf", plot = bpjp, height = 7, width = 10)


# node 1
sub.group <- subset(data, ssom.unit.classif=="1")
expression <- sub.group[c(1,12:16)]
m.expression <- melt(expression, id=c("X")) # multiple column into one column
head(m.expression)

n1 <- ggplot(m.expression, aes(x=variable, y=value)) + 
      geom_point(position="jitter",size=1.5) + theme_bw() +
      theme(text = element_text(size=13)) + theme(axis.title.x = element_blank())
n1

ggsave(file = "node1.boxplot.pdf", plot = n1, height = 5, width = 5)


# node 2
sub.group <- subset(data, ssom.unit.classif=="2")
expression <- sub.group[c(1,12:16)]
m.expression <- melt(expression, id=c("X")) # multiple column into one column
head(m.expression)

n2 <- ggplot(m.expression, aes(x=variable, y=value)) + 
  geom_point(position="jitter",size=1.5) + theme_bw() +
  theme(text = element_text(size=13)) + theme(axis.title.x = element_blank())
n2

ggsave(file = "node2.boxplot.pdf", plot = n2, height = 5, width = 5)

# node 3
sub.group <- subset(data, ssom.unit.classif=="3")
expression <- sub.group[c(1,12:16)]
m.expression <- melt(expression, id=c("X")) # multiple column into one column
head(m.expression)

n3 <- ggplot(m.expression, aes(x=variable, y=value)) + 
  geom_point(position="jitter",size=1.5) + theme_bw() +
  theme(text = element_text(size=13)) + theme(axis.title.x = element_blank())
n3

ggsave(file = "node3.boxplot.pdf", plot = n3, height = 5, width = 5)

# node 4
sub.group <- subset(data, ssom.unit.classif=="4")
expression <- sub.group[c(1,12:16)]
m.expression <- melt(expression, id=c("X")) # multiple column into one column
head(m.expression)

n4 <- ggplot(m.expression, aes(x=variable, y=value)) + 
  geom_point(position="jitter",size=1.5) + theme_bw() +
  theme(text = element_text(size=13)) + theme(axis.title.x = element_blank())
n4

ggsave(file = "node4.boxplot.pdf", plot = n4, height = 5, width = 5)

# node 5
sub.group <- subset(data, ssom.unit.classif=="5")
expression <- sub.group[c(1,12:16)]
m.expression <- melt(expression, id=c("X")) # multiple column into one column
head(m.expression)

n5 <- ggplot(m.expression, aes(x=variable, y=value)) + 
  geom_point(position="jitter",size=1.5) + theme_bw() +
  theme(text = element_text(size=13)) + theme(axis.title.x = element_blank())
n5

ggsave(file = "node5.boxplot.pdf", plot = n5, height = 5, width = 5)

# node 6
sub.group <- subset(data, ssom.unit.classif=="6")
expression <- sub.group[c(1,12:16)]
m.expression <- melt(expression, id=c("X")) # multiple column into one column
head(m.expression)

n6 <- ggplot(m.expression, aes(x=variable, y=value)) + 
  geom_point(position="jitter",size=1.5) + theme_bw() +
  theme(text = element_text(size=13)) + theme(axis.title.x = element_blank())
n6

ggsave(file = "node6.boxplot.pdf", plot = n6, height = 5, width = 5)

# node 7
sub.group <- subset(data, ssom.unit.classif=="7")
expression <- sub.group[c(1,12:16)]
m.expression <- melt(expression, id=c("X")) # multiple column into one column
head(m.expression)

n7 <- ggplot(m.expression, aes(x=variable, y=value)) + 
  geom_point(position="jitter",size=1.5) + theme_bw() +
  theme(text = element_text(size=13)) + theme(axis.title.x = element_blank())
n7

ggsave(file = "node7.boxplot.pdf", plot = n7, height = 5, width = 5)

# node 8
sub.group <- subset(data, ssom.unit.classif=="8")
expression <- sub.group[c(1,12:16)]
m.expression <- melt(expression, id=c("X")) # multiple column into one column
head(m.expression)

n8 <- ggplot(m.expression, aes(x=variable, y=value)) + 
  geom_point(position="jitter",size=1.5) + theme_bw() +
  theme(text = element_text(size=13)) + theme(axis.title.x = element_blank())
n8

ggsave(file = "node8.boxplot.pdf", plot = n8, height = 5, width = 5)

## Line plot (Loess regression lines for all SOMs) for all SOM clusters for 5 tissues 
expression <- data[c(1,12:16,22)]
m.expression <- melt(expression, id=c("X","ssom.unit.classif"))
head(m.expression)
lp <-  ggplot(m.expression, aes(x=variable, y=value, group=ssom.unit.classif))+ theme_bw() +
      stat_smooth(method="loess",aes(colour=factor(ssom.unit.classif)),size=1) +
      scale_colour_manual(name = "ssom clusters", values=c("tomato","springgreen3","blue3","orange1","magenta3","yellow4","turquoise2","deeppink1")) +
      theme(text = element_text(size=13)) + theme(axis.title.x = element_blank())

lp

ggsave(file = "loess_regression_plot_all_som_clusters.pdf", plot = lp, height = 5, width = 8)

#==================================================================================================
# No genotype or treatment effect here (incomplete)...

data <- read.csv("all_tissues_mean_table.csv", h = T, row.names = 1)
head(data)
data <- data[,-c(1,12,13,19)]
dim(data)
data$mean <- rowMeans(data[,1:20])
data$sd <- apply(data[,1:20],1,sd)
data$cv <- data$sd/data$mean

quantile(data$cv, probs=seq(.05,.95,.05), na.rm=T) 
data_2 <- subset(data, cv > 4.69866289)
head(data_2)
dim(data_2)
nrow(data_2) # 1827

data_2_R500_SUN <- data_2[,grep("GCR500_SUN*", colnames(data_2))]
head(data_2_R500_SUN)
dim(data_2_R500_SUN)
data_2_R500_SUN <- cbind(data_2_R500_SUN,SP = rep("R500",1827), trt = rep("SUN", 1827))
#str(data_2_R500_SUN)    
colnames(data_2_R500_SUN) <- c("INTERNODE", "LEAF", "ROOT", "SEEDLING", "SILIQUE", "SP", "trt")                      

data_2_R500_SHADE <- data_2[,grep("GCR500_SHADE*", colnames(data_2))]
head(data_2_R500_SHADE)
dim(data_2_R500_SHADE)
data_2_R500_SHADE <- cbind(data_2_R500_SHADE,SP = rep("R500",1827), trt = rep("SHADE", 1827))
colnames(data_2_R500_SHADE) <- c("INTERNODE", "LEAF", "ROOT", "SEEDLING", "SILIQUE", "SP", "trt")                      


data_2_IMB211_SUN <- data_2[,grep("GCIMB211_SUN*", colnames(data_2))]
head(data_2_IMB211_SUN)
dim(data_2_IMB211_SUN)
data_2_IMB211_SUN <- cbind(data_2_IMB211_SUN,SP = rep("IMB211",1827), trt = rep("SUN", 1827))
colnames(data_2_IMB211_SUN) <- c("INTERNODE", "LEAF", "ROOT", "SEEDLING", "SILIQUE", "SP", "trt")                      


data_2_IMB211_SHADE <- data_2[,grep("GCIMB211_SHADE*", colnames(data_2))]
head(data_2_IMB211_SHADE)
dim(data_2_IMB211_SHADE)
data_2_IMB211_SHADE <- cbind(data_2_IMB211_SHADE,SP = rep("IMB211",1827), trt = rep("SHADE", 1827))
colnames(data_2_IMB211_SHADE) <- c("INTERNODE", "LEAF", "ROOT", "SEEDLING", "SILIQUE", "SP", "trt")                      


data_2_R500_IMB211 <- rbind(data_2_R500_SUN, data_2_R500_SHADE, data_2_IMB211_SUN, data_2_IMB211_SHADE)
head(data_2_R500_IMB211)
dim(data_R500_IMB211)
write.csv(file = "R500_IMB211_cv_95.csv", data_2_R500_IMB211)
