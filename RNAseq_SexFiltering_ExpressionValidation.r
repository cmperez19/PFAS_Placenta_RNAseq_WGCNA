library(dplyr)
library(stringr)
library(DESeq2)
library(corrplot)
library(ggplot2)
library(ggfortify) 
library(matrixStats)

# Reading unfiltered_glowing_counts.csv...
unfiltered_glowing_counts <- read.csv("unfiltered_glowing_counts.csv")

# Loading phenotype data...
load("log_PFAS_pheno.Rdata")

# Cleaning column and row names...
colnames(unfiltered_glowing_counts) <- gsub(x = colnames(unfiltered_glowing_counts), pattern = "\\.", replacement = "-")
rownames(unfiltered_glowing_counts) <- unfiltered_glowing_counts$Geneid

# Selecting female Ymasked data...
female_ymasked <- unfiltered_glowing_counts %>% select(contains("Ymasked"))
colnames(female_ymasked) <- str_extract(colnames(female_ymasked), "[^_]+")

# Selecting male YPARsMasked data...
default_alignment <- unfiltered_glowing_counts %>% select(contains("YPARsMasked"))
colnames(default_alignment) <- str_extract(colnames(default_alignment), "[^_]+")

male_yparsmasked <- default_alignment[colnames(default_alignment) %in% pheno[which(pheno$childs_sex == "Male"), "participant_id"]]

# Merging female and male data...
sexalignment <- merge(female_ymasked, male_yparsmasked, by = 'row.names')
rownames(sexalignment) <- sexalignment$Row.names
sexalignment <- sexalignment[, -1]

# Cleaning and sorting merged data...
common_ids <- intersect(colnames(sexalignment), rownames(pheno))
sexalignment <- sexalignment[, common_ids]
pheno <- pheno[common_ids, ]

# Creating DESeq dataset...
deseq_sexalignment <- DESeqDataSetFromMatrix(countData = sexalignment, colData = pheno, design = ~ 1)
# Initial number of genes:
print(nrow(deseq_sexalignment))

# Filtering low-count genes...
keep <- rowSums(fpm(deseq_sexalignment) >= 1) >= 75
deseq_sexalignment <- deseq_sexalignment[keep, ]
# Number of genes after filtering:
print(nrow(deseq_sexalignment))

sex_complement_raw_counts <- counts(deseq_sexalignment)

# Computing max expression per gene...
max <- rowMaxs(as.matrix(sex_complement_raw_counts), na.rm = TRUE)
names(max) <- rownames(sex_complement_raw_counts)
max <- max[names(sort(max, decreasing = TRUE))]

# Summarizing top expressed genes...
for (i in head(names(max))) {
  print(i)
  print(min(sex_complement_raw_counts[i, ]))
  print(max(sex_complement_raw_counts[i, ]))
  print(mean(sex_complement_raw_counts[i, ]))
  print(sd(sex_complement_raw_counts[i, ]))
}

# Sampling random genes...
set.seed(110)
sample <- sample(rownames(sex_complement_raw_counts), 10)

# Loading original RNAseq raw counts...
load("Glowing_RNAse_141QCd.RData")
GlowingPlacenta_RNAseq.rawCounts <- counts

# Checking sample overlap with original dataset...
sample %in% rownames(GlowingPlacenta_RNAseq.rawCounts)

# Identifying shared genes between datasets...
shared_genes <- rownames(GlowingPlacenta_RNAseq.rawCounts)[which(rownames(GlowingPlacenta_RNAseq.rawCounts) %in% rownames(sex_complement_raw_counts))]

# Sampling shared genes...
set.seed(110)
sample <- sample(rownames(sex_complement_raw_counts[shared_genes, ]), 10)
sample %in% rownames(GlowingPlacenta_RNAseq.rawCounts)

# Computing correlation matrix between datasets...
cor_matrix <- cor(t(GlowingPlacenta_RNAseq.rawCounts[sample, ]), t(sex_complement_raw_counts[sample, colnames(GlowingPlacenta_RNAseq.rawCounts)]))

# Plotting XIST expression...
plot(t(GlowingPlacenta_RNAseq.rawCounts["XIST", ]), t(sex_complement_raw_counts["XIST", colnames(GlowingPlacenta_RNAseq.rawCounts)]))
