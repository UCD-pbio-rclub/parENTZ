###########
# Cody Markelz
# markelz@gmail.com
# B rapa parental RNA-seq
# Modified December 4, 2015
###########

# infile 
setwd("/Users/Cody_2/git.repos/brassica_parents/data")

GC_counts <- read.delim("GC_merged_v1.5_mapping.tsv", row.names = NULL)
head(GC_counts)
rownames(GC_counts) <- GC_counts$gene
GC_counts <- GC_counts[,-1]

GH_counts <- read.delim("GH_merged_v1.5_mapping.tsv", row.names = NULL)
head(GH_counts)
rownames(GH_counts) <- GH_counts$gene
GH_counts <- GH_counts[,-1]

# clean up the colnames
colnames(GC_counts)
colnames(GH_counts)

colnames(GC_counts) <- sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)(.+)", "\\1\\2\\3\\4\\5\\6\\7", colnames(GC_counts))
colnames(GH_counts) <- sub("(\\w+)(_)(\\w+)(_)(\\d)(_)(\\w+)(.+)", "\\1\\2\\3\\4\\5\\6\\7", colnames(GH_counts))

