###########
# Cody Markelz
# markelz@gmail.com
# B rapa parental RNA-seq
# Modified December 14, 2015
###########
# TODO:
# subset genotype to look at shade

# load libs
library(edgeR)
library(limma)

# data directory
setwd("/Users/Cody_2/git.repos/brassica_parents/data")

##########
# growth chamber and field data import and clean up
##########
GC_counts <- read.delim("GC_merged_v1.5_mapping.tsv", row.names = NULL)
head(GC_counts)
rownames(GC_counts) <- GC_counts$gene
str(GC_counts)

# remove now row names
GC_counts <- GC_counts[,-1]

# replace all NA values with 0 
GC_counts[is.na(GC_counts)] <- 0
head(GC_counts)
tail(GC_counts)

# remove bad GC libs
colnames(GC_counts)
GC_counts <- GC_counts[,-c(13,52)]
head(GC_counts)

# remove first row of total counts
GC_counts <- GC_counts[-1,]
head(GC_counts)

##########
# greenhouse and field data import and clean up
##########

# make sure to remove bad leaf libs
# seperate field and GH data

GH_counts <- read.delim("GH_merged_v1.5_mapping.tsv", row.names = NULL)
head(GH_counts)
rownames(GH_counts) <- GH_counts$gene
GH_counts <- GH_counts[,-1]

#replace all NA values with 0 
GH_counts[is.na(GH_counts)] <- 0
head(GH_counts)
tail(GH_counts)

# remove bad GC libs
colnames(GH_counts)
GH_counts <- GH_counts[,-c(50,54,58)]
head(GH_counts)

# remove first row
GH_counts <- GH_counts[-1,]
head(GH_counts)[,1:10]

# make field dataset
colnames(GH_counts)
FIELD_counts <- GH_counts[,c(13:18,43:48)]
head(FIELD_counts)
GH_counts <- GH_counts[,-c(13:18,43:48)]
colnames(GH_counts)

###########
# clean up the colnames
###########
colnames(GC_counts)
colnames(GH_counts)
colnames(FIELD_counts)

colnames(GC_counts) <- sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)(.+)",
                             "\\1\\2\\3\\4\\5\\6\\7", colnames(GC_counts))
colnames(GH_counts) <- sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)(.+)",
                             "\\1\\2\\3\\4\\5\\6\\7", colnames(GH_counts))
colnames(FIELD_counts) <- sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)(.+)",
                             "\\1\\2\\3\\4\\5\\6\\7", colnames(FIELD_counts))

#write cleaned tables
write.table(GC_counts, file="GC_counts_counts_cleaned.csv", sep=",") 
write.table(GH_counts, file="GH_counts_counts_cleaned.csv", sep=",") 
write.table(FIELD_counts, file="FIELD_counts_counts_cleaned.csv", sep=",")

# design matrix 
GH_samples <- names(GH_counts)
GC_samples <- names(GC_counts)
FIELD_samples <- names(FIELD_counts)

GH_samples
GC_samples
FIELD_samples

#clean way to make design matrices
GC_geno <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\1", colnames(GC_counts)))
GC_geno

GC_trt <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\3", colnames(GC_counts)))
GC_trt

GC_tissue <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\7", colnames(GC_counts)))
GC_tissue

# green house 
GH_geno <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\1", colnames(GH_counts)))
GH_geno

GH_trt <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\3", colnames(GH_counts)))
GH_trt

GH_tissue <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\7", colnames(GH_counts)))
GH_tissue

# field tissue
FIELD_geno <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\1", colnames(FIELD_counts)))
FIELD_geno

FIELD_tissue <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\3", colnames(FIELD_counts)))
FIELD_tissue

#############
#############
# analyze the GC first
#############
#############
GC_geno <- relevel(GC_geno, ref = "R500")
GC_trt  <- relevel(GC_trt, ref = "SUN")
GC_tissue  <- relevel(GC_tissue, ref = "LEAF")

full_model <- model.matrix(~ GC_geno + GC_trt + GC_tissue + GC_geno:GC_trt )
full_model

# make the necessary objects
head(GC_counts)
GC_DE <- DGEList(counts = GC_counts)
GC_DE <- GC_DE[rowSums(cpm(GC_DE) > 1 ) >= 20,]
dim(GC_DE)
# [1] 27621    66

GC_DE <- calcNormFactors(GC_DE)
system.time(GC_voom <- voom(GC_DE, full_model, plot = TRUE))
system.time(GC_fit_full <-lmFit(GC_voom, full_model))
GC_fit_full <- eBayes(GC_fit_full)
toptable(GC_fit_full)
str(GC_fit_full)
head(GC_fit_full$coef)
?topTable
topTable(GC_fit_full, coef="GC_genoIMB211", n = 200)

GC_fit_full <- treat(GC_fit_full, lfc = log2(1.5))
topTreat(GC_fit_full, coef="GC_genoIMB211", n = 20)
topTreat(GC_fit_full, coef="GC_genoIMB211:GC_trtSUN", n = 20)
head(coef(GC_fit_full))


?decideTests
GC_result <- decideTests(GC_fit_full)
GC_result
head(GC_result, 50)
summary(GC_result)
?genas
?vooma

topTable(GC_fit_full, coef="GC_genoIMB211:GC_trtSHADE")

designlist <- list(
  geno <- model.matrix(~GC_geno),
  trt <- model.matrix(~GC_trt),
  tissue <- model.matrix(~GC_tissue), 
  g_trt <- model.matrix(~GC_geno*GC_trt),
  tis_geno <- model.matrix(~GC_geno*GC_tissue),
  trt_tis <- model.matrix(~GC_trt*GC_tissue),
  full <- model.matrix(~GC_geno*GC_trt*GC_trt)
)

designlist

# take a look at model selection
out <- selectModel(GC_voom,designlist, criterion="bic")

table(out$pref)
head(out$pref)
str(out)

# tissue driving most of the response
modelout <- as.data.frame(out$pref)
head(modelout)

# genotype only model
geno <- model.matrix(~GC_geno)
GC_DE <- calcNormFactors(GC_DE)
system.time(GC_geno_voom <- voom(GC_DE, geno, plot = TRUE))
system.time(GC_geno_fit <- lmFit(GC_geno_voom, geno))
GC_geno_fit <- eBayes(GC_geno_fit)
topTable(GC_geno_fit, number = 50)
# output gene lists

# treatment only model
# basically nothing. All variance is in tissue type.
trt <- model.matrix(~GC_trt)
system.time(GC_trt_voom <- voom(GC_DE, trt, plot = TRUE))
system.time(GC_trt_fit <- lmFit(GC_trt_voom, trt))
GC_trt_fit <- eBayes(GC_trt_fit)
?topTable
topTable(GC_trt_fit, number = 50)
# output gene lists

# tissue only model
tissue <- model.matrix(~GC_tissue)
system.time(GC_tissue_voom <- voom(GC_DE, tissue, plot = TRUE))
system.time(GC_tissue_fit <- lmFit(GC_tissue_voom, tissue))
GC_tissue_fit <- eBayes(GC_tissue_fit)
?topTable
topTable(GC_tissue_fit, number = 500)
# output gene lists

# interaction model
trt_geno <- model.matrix(~GC_trt:GC_geno)
system.time(GC_trt_geno_voom <- voom(GC_DE, trt_geno, plot = TRUE))
system.time(GC_trt_geno_fit <- lmFit(GC_trt_geno_voom, trt_geno))
GC_trt_geno_fit <- eBayes(GC_trt_geno_fit)
?topTable
topTable(GC_trt_geno_fit, number = 50)
# output gene lists



