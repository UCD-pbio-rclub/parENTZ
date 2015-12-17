###########
# Cody Markelz
# markelz@gmail.com
# B rapa parental RNA-seq
# Modified December 16, 2015
###########

# DE between genotypes
# DE between tissues
# GO enrichment?
# metabolic differences?

# load libs
library(edgeR)
library(limma)

# data directory
setwd("/Users/Cody_2/git.repos/brassica_parents/data")

GH_counts <- read.delim("GH_merged_v1.5_mapping.tsv", row.names = NULL)
head(GH_counts)
rownames(GH_counts) <- GH_counts$gene
GH_counts <- GH_counts[,-1]

#replace all NA values with 0 
GH_counts[is.na(GH_counts)] <- 0
head(GH_counts)
tail(GH_counts)
head(GH_counts)

# remove first row
GH_counts <- GH_counts[-1,]
head(GH_counts)[,1:10]

# make field dataset
colnames(GH_counts)
FIELD_counts <- GH_counts[,c(13:18,43:48)]
head(FIELD_counts)

###########
# clean up the colnames
###########
colnames(FIELD_counts) <- sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)(.+)",
                             "\\1\\2\\3\\4\\5\\6\\7", colnames(FIELD_counts))

FIELD_samples <- names(FIELD_counts)
FIELD_samples

# field tissue
FIELD_geno <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\1", colnames(FIELD_counts)))
FIELD_geno

FIELD_tissue <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\3", colnames(FIELD_counts)))
FIELD_tissue

FIELD_geno <- relevel(FIELD_geno, ref = "R500")
FIELD_tissue  <- relevel(FIELD_tissue, ref = "LEAF")

full_model <- model.matrix(~ FIELD_geno + FIELD_tissue + FIELD_geno*FIELD_tissue)
full_model

# make the necessary objects
head(FIELD_counts)
?DGEList
FIELD_DE <- DGEList(counts = FIELD_counts)
dim(FIELD_DE)
colnames(FIELD_DE)
cpm(FIELD_DE)

FIELD_DE <- FIELD_DE[rowSums(cpm(FIELD_DE) > 1 ) >= 5,]
dim(FIELD_DE)
# [1] 27621    66

FIELD_DE <- calcNormFactors(FIELD_DE)
system.time(FIELD_voom <- voom(FIELD_DE, full_model, plot = TRUE))
system.time(FIELD_fit_full <-lmFit(FIELD_voom, full_model))
FIELD_fit_full <- eBayes(FIELD_fit_full)

?topTable
toptable(FIELD_fit_full, n = 5000)

str(FIELD_fit_full)
head(FIELD_fit_full$coef)


