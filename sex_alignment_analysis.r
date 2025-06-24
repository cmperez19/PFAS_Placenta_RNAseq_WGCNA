library(dplyr)
library(stringr)
library(DESeq2)
library(corrplot)
library(ggplot2)
library(ggfortify) 

# Load data
unfiltered_glowing_counts <- read.csv("unfiltered_glowing_counts.csv")
load("log_PFAS_pheno.Rdata")

# Fix column names and set row names
colnames(unfiltered_glowing_counts) <- gsub("\\.", "-", colnames(unfiltered_glowing_counts))
rownames(unfiltered_glowing_counts) <- unfiltered_glowing_counts$Geneid

# Separate data by masking patterns
female_ymasked <- unfiltered_glowing_counts %>% select(contains("Ymasked"))
colnames(female_ymasked) <- str_extract(colnames(female_ymasked), "[^_]+")

default_alignment <- unfiltered_glowing_counts %>% select(contains("YPARsMasked"))
colnames(default_alignment) <- str_extract(colnames(default_alignment), "[^_]+")

# Subset male samples from pheno data
male_yparsmasked <- default_alignment[, colnames(default_alignment) %in% pheno[pheno$childs_sex == "Male", "participant_id"]]

# Merge male and female matrices
sexalignment <- merge(female_ymasked, male_yparsmasked, by = 'row.names')
rownames(sexalignment) <- sexalignment$Row.names
sexalignment <- sexalignment[, -1]

# Filter pheno and alignment to match
sexalignment <- sexalignment[, order(colnames(sexalignment))]
pheno <- pheno[-c(8, 115, 141), ]
pheno <- pheno[order(rownames(pheno)), ]

# DESeq2 object
pheno_150 <- pheno

# Create DESeqDataSet
sexalignment <- sexalignment[, rownames(pheno_150)]
deseq_sexalignment <- DESeqDataSetFromMatrix(countData = sexalignment, colData = pheno_150, design = ~ 1)

# Filter for expressed genes (>=1 CPM in 75+ samples)
keep <- rowSums(fpm(deseq_sexalignment) >= 1) >= 75
deseq_sexalignment <- deseq_sexalignment[keep, ]

# Extract raw counts
sex_complement_raw_counts <- counts(deseq_sexalignment)

# Get highest expressed genes
max <- apply(sex_complement_raw_counts, 1, max, na.rm = TRUE)
max <- sort(max, decreasing = TRUE)

# Display summary for top genes
for (i in head(names(max))) {
  cat("\nGene:", i, "\n")
  cat("Min:", min(sex_complement_raw_counts[i, ]), "\n")
  cat("Max:", max(sex_complement_raw_counts[i, ]), "\n")
  cat("Mean:", mean(sex_complement_raw_counts[i, ]), "\n")
  cat("SD:", sd(sex_complement_raw_counts[i, ]), "\n")
}

# Compare with original raw counts
set.seed(110)
sample <- sample(rownames(sex_complement_raw_counts), 10)
load("Glowing_RNAse_141QCd.RData")
GlowingPlacenta_RNAseq.rawCounts <- counts

# Identify shared genes
shared_genes <- intersect(rownames(GlowingPlacenta_RNAseq.rawCounts), rownames(sex_complement_raw_counts))
set.seed(110)
sample <- sample(shared_genes, 10)

# Correlation matrix
cor_matrix <- cor(t(GlowingPlacenta_RNAseq.rawCounts[sample, ]), 
                  t(sex_complement_raw_counts[sample, colnames(GlowingPlacenta_RNAseq.rawCounts)]))

# Plot XIST expression
plot(t(GlowingPlacenta_RNAseq.rawCounts["XIST", ]), 
     t(sex_complement_raw_counts["XIST", colnames(GlowingPlacenta_RNAseq.rawCounts)]), 
     xlab = "Original", ylab = "Filtered", main = "XIST expression")
