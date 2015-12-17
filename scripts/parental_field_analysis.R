###########
# Cody Markelz
# markelz@gmail.com
# B rapa parental RNA-seq
# Modified December 17, 2015
###########

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

# tissue specific data sets
FIELD_leaf <- FIELD_counts[,c(4:6,10:12)]
FIELD_fruit <- FIELD_counts[,-c(4:6,10:12)]

head(FIELD_leaf)

# field design
FIELD_geno <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\1", colnames(FIELD_counts)))
FIELD_geno

FIELD_tissue <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\3", colnames(FIELD_counts)))
FIELD_tissue

FIELD_geno <- relevel(FIELD_geno, ref = "R500")
FIELD_tissue  <- relevel(FIELD_tissue, ref = "LEAF")

full_model <- model.matrix(~ FIELD_geno + FIELD_tissue + FIELD_geno:FIELD_tissue)
full_model

# tissue models
leaf_geno <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\1", colnames(FIELD_leaf)))
leaf_geno
leaf_geno <- relevel(leaf_geno, ref = "R500")

leaf_model <- model.matrix(~ leaf_geno)
leaf_model

fruit_geno <- factor(sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)",
                       "\\1", colnames(FIELD_fruit)))
fruit_geno

fruit_model <- model.matrix(~ fruit_geno)
fruit_model

fruit_geno <- relevel(fruit_geno, ref = "R500")

# make the necessary objects
head(FIELD_counts)
FIELD_DE <- DGEList(counts = FIELD_counts)
dim(FIELD_DE)
colnames(FIELD_DE)
cpm(FIELD_DE)

FIELD_DE <- FIELD_DE[rowSums(cpm(FIELD_DE) > 1 ) >= 5,]
dim(FIELD_DE)
# [1] 25740    12

FIELD_DE <- calcNormFactors(FIELD_DE)
system.time(FIELD_voom <- voom(FIELD_DE, full_model, plot = TRUE))
system.time(FIELD_fit_full <-lmFit(FIELD_voom, full_model))
FIELD_fit_full <- eBayes(FIELD_fit_full)

?topTable
geno_eff <- toptable(FIELD_fit_full, coef = "FIELD_genoIMB211", p.value = 0.05, n = Inf)
dim(geno_eff)

tissue_eff <- toptable(FIELD_fit_full, coef = "FIELD_tissueFRUIT", p.value = 0.05, n = Inf)
dim(tissue_eff)

str(FIELD_fit_full)
head(FIELD_fit_full$coef)

#########
# tissue specific models
#########

# leaf models
head(FIELD_leaf)
leaf_de <- DGEList(counts = FIELD_leaf)
dim(leaf_de)
colnames(leaf_de)
leaf_de <- leaf_de[rowSums(cpm(leaf_de) > 1 ) >= 5,]
dim(leaf_de)
# [1] 20725     6

leaf_de <- calcNormFactors(leaf_de)
system.time(leaf_voom <- voom(leaf_de, leaf_model, plot = TRUE))
system.time(leaf_fit <-lmFit(leaf_voom, leaf_model))
leaf_fit <- eBayes(leaf_fit)
leaf_fit
leaf_genes <- topTable(leaf_fit, p.value = 0.05, n = Inf)
dim(leaf_genes)
# [1] 4435    6
head(leaf_genes)
tail(leaf_genes)

# fruit models
head(FIELD_fruit)
fruit_de <- DGEList(counts = FIELD_fruit)
dim(fruit_de)
colnames(fruit_de)
fruit_de <- fruit_de[rowSums(cpm(fruit_de) > 1 ) >= 5,]
dim(fruit_de)
# [1] 22511     6

fruit_de <- calcNormFactors(fruit_de)
system.time(fruit_voom <- voom(fruit_de, fruit_model, plot = TRUE))
system.time(fruit_fit <-lmFit(fruit_voom, fruit_model))
fruit_fit <- eBayes(fruit_fit)
fruit_fit
fruit_genes <- topTable(fruit_fit, p.value = 0.05, n = Inf)
dim(fruit_genes)
# [1] 2967    6
head(fruit_genes)

# overlap genes
merged_genes <- merge(leaf_genes, fruit_genes, by = "row.names")
head(merged_genes)
dim(merged_genes)

# write tables
setwd("/Users/Cody_2/git.repos/brassica_parents/data")
?write.table
write.table(leaf_genes, "parental_leaf_field_DE.csv", sep = ",", row.names = TRUE, col.names = TRUE)
write.table(fruit_genes, "parental_fruit_field_DE.csv", sep = ",", row.names = TRUE, col.names = TRUE)
write.table(merged_genes, "parental_merged_field_DE.csv", sep = ",", row.names = TRUE, col.names = TRUE)

#end
