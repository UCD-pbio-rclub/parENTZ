###########
# Cody Markelz
# markelz@gmail.com
# B rapa parental RNA-seq
# Modified December 12, 2015
###########
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
GC_counts <- GC_counts[,-1]

#replace all NA values with 0 
GC_counts[is.na(GC_counts)] <- 0
head(GC_counts)
tail(GC_counts)

# remove first row
GC_counts <- GC_counts[-1,]
head(GC_counts)[,1:10]


##########
# greenhouse and field data import and clean up
##########
GH_counts <- read.delim("GH_merged_v1.5_mapping.tsv", row.names = NULL)
head(GH_counts)
rownames(GH_counts) <- GH_counts$gene
GH_counts <- GH_counts[,-1]

#replace all NA values with 0 
GH_counts[is.na(GH_counts)] <- 0
head(GH_counts)
tail(GH_counts)

# remove first row
GH_counts <- GH_counts[-1,]
head(GH_counts)[,1:10]

# clean up the colnames
colnames(GC_counts)
colnames(GH_counts)

colnames(GC_counts) <- sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)(.+)",
                             "\\1\\2\\3\\4\\5\\6\\7", colnames(GC_counts))
colnames(GH_counts) <- sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)(.+)",
                             "\\1\\2\\3\\4\\5\\6\\7", colnames(GH_counts))

#write cleaned tables
write.table(GC_counts, file="GC_counts_counts_cleaned.csv", sep=",") 
write.table(GH_counts, file="GH_counts_counts_cleaned.csv", sep=",") 

# design matrix 
GH_samples <- names(GH_counts)
GC_samples <- names(GC_counts)

GH_samples
GC_samples

# make a clean way to make design matrices
GC_geno <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\1", colnames(GC_counts)))
GC_geno

GC_trt <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\3", colnames(GC_counts)))
GC_trt

GC_tissue <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\7", colnames(GC_counts)))
GC_tissue

# green house and field tissue
GH_geno <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\7", colnames(GH_counts)))
GH_geno

GH_trt <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\7", colnames(GH_counts)))
GH_trt

GH_tissue <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\7", colnames(GH_counts)))
GH_tissue

# analyze the GC first
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
topTable(GC_fit_full, coef="GC_genoIMB211", n = 2000)

GC_fit_full <- treat(GC_fit_full, lfc = 2)
topTreat(GC_fit_full, coef="GC_genoIMB211", n = 20)
topTreat(GC_fit_full, coef="GC_genoIMB211:GC_trtSUN", n = 20)
head(coef(GC_fit_full))


?decideTests
GC_result <- decideTests(GC_fit_full)
GC_result
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
  trt_tis <- model.matrix(~GC_trt*GC_tissue)
  # full <- model.matrix(~GC_geno*GC_trt*GC_trt)
)

designlist
?selectModel
out <- selectModel(GC_voom,designlist, criterion="bic")

table(out$pref)
str(out)

modelout <- out$pref
head(modelout)

# genotype only model
geno <- model.matrix(~GC_geno)
GC_DE <- calcNormFactors(GC_DE)
system.time(GC_geno_voom <- voom(GC_DE, geno, plot = TRUE))
system.time(GC_geno_fit <-lmFit(GC_geno_voom, geno))
GC_geno_fit <- eBayes(GC_geno_fit)
topTable(GC_geno_fit)

# treatment only model
trt <- model.matrix(~GC_trt)
system.time(GC_trt_voom <- voom(GC_DE, trt, plot = TRUE))
system.time(GC_trt_fit <-lmFit(GC_trt_voom, trt))
GC_trt_fit <- eBayes(GC_trt_fit)
?topTable
topTable(GC_trt_fit, number = 500)


# treatment only model
trt <- model.matrix(~GC_trt)
system.time(GC_trt_voom <- voom(GC_DE, trt, plot = TRUE))
system.time(GC_trt_fit <-lmFit(GC_trt_voom, trt))
GC_trt_fit <- eBayes(GC_trt_fit)
?topTable
topTable(GC_trt_fit, number = 500)



