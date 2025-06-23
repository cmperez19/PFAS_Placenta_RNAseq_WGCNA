library(dplyr)
library(stringr)
library(DESeq2)
library(corrplot)
library(ggplot2)
library(ggfortify) 
unfiltered_glowing_counts <- read.csv("unfiltered_glowing_counts.csv")
load("log_PFAS_pheno.Rdata")


colnames(unfiltered_glowing_counts) <- gsub(x = colnames(unfiltered_glowing_counts),pattern = "\\.",replacement = "-")
rownames(unfiltered_glowing_counts) <- unfiltered_glowing_counts$Geneid


female_ymasked <- unfiltered_glowing_counts %>% select(contains("Ymasked"))
colnames(female_ymasked) <- str_extract(colnames(female_ymasked), "[^_]+")


default_alignment <- unfiltered_glowing_counts %>% select(contains("YPARsMasked"))
colnames(default_alignment) <- str_extract(colnames(default_alignment), "[^_]+")

male_yparsmasked <- default_alignment[colnames(default_alignment) %in% pheno[which(pheno$childs_sex=="Male"),"participant_id"]]
#I had the phenred codes...might be deleted ): 

sexalignment <- merge(female_ymasked,male_yparsmasked, by = 'row.names')
rownames(sexalignment) = sexalignment$Row.names
sexalignment =sexalignment[,-1]


## congruency between IDs in pheno and sexalignment 
#which(colnames(sexalignment) %in% rownames(pheno) == FALSE)
#colnames(sexalignment)[12]
#G-165

#which(rownames(pheno) %in% colnames(sexalignment) == FALSE)
#rownames(pheno)[8]
#G-247
#rownames(pheno)[115]
#G-411
# rownames(pheno)[141] 
#G-316

sexalignment = sexalignment[,-12]
sexalignment = sexalignment[, order(colnames(sexalignment))]
pheno = pheno[-c(8,115,141),]
pheno = pheno[order(rownames(pheno)),]

#colnames(sexalignment) == rownames(pheno) 
#TRUE 

## create DESEQ data seq 
deseq_sexalignment = DESeqDataSetFromMatrix(countData = sexalignment, colData = pheno_150, design = ~ 1)
nrow(deseq_sexalignment) # 56852
#rerrun this section 
keep <- rowSums(fpm(deseq_sexalignment)>=1) >= 75 #1 count per million (in data with ~20M reads: cpm 1 ~ 20 counts/sample)
deseq_sexalignment <- deseq_sexalignment[keep,]
nrow(deseq_sexalignment) # 15028

sex_complement_raw_counts = counts(deseq_sexalignment)

#picking highly expressed genes from correlation plot 
max <- matrix(apply(sex_complement_raw_counts, 1, max, na.rm=TRUE))
rownames(max) = rownames(sex_complement_raw_counts)
max = max[rownames(max)[order(max[,1],decreasing = T)],1]

#MT-ND4  MT-CO1  MT-CYB    CSH1     FN1  MT-ND2 
#1400493  887247  761549  674330  567822  563444 

for( i in head(names(max))){
  print(i)
  print(min(sex_complement_raw_counts[i,]))
  print(max(sex_complement_raw_counts[i,]))
  print(mean(sex_complement_raw_counts[i,]))
  print(sd(sex_complement_raw_counts[i,]))
}

#choose random genes to compare with og data
set.seed(110)
sample = sample(rownames(sex_complement_raw_counts), 10)
#are these sampled genes in the og raw data? 
load("Glowing_RNAse_141QCd.RData")
GlowingPlacenta_RNAseq.rawCounts = counts
#dim(GlowingPlacenta_RNAseq.rawCounts)
#13831   141

sample %in% rownames(GlowingPlacenta_RNAseq.rawCounts) 
#TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE

#sum(rownames(GlowingPlacenta_RNAseq.rawCounts) %in% rownames(sex_complement_raw_counts) )
#11279
shared_genes = rownames(GlowingPlacenta_RNAseq.rawCounts)[which(rownames(GlowingPlacenta_RNAseq.rawCounts) %in% rownames(sex_complement_raw_counts) == TRUE)]


set.seed(110)
sample = sample(rownames(sex_complement_raw_counts[shared_genes,]),10)
sample %in% rownames(GlowingPlacenta_RNAseq.rawCounts)

cor_matrix <- cor(t(GlowingPlacenta_RNAseq.rawCounts[sample,]), t(sex_complement_raw_counts[sample,colnames(GlowingPlacenta_RNAseq.rawCounts)]))

# ploting XIST do it again after transformation 
plot(t(GlowingPlacenta_RNAseq.rawCounts["XIST",]), t(sex_complement_raw_counts["XIST",colnames(GlowingPlacenta_RNAseq.rawCounts)]))
#do a plot for a Y chr linked gene and one for a gene within sample 
#check highly expressed genes that are shared between both datasets 
#check how many Y chromosome genes are lost from 30 to 75 
#can do 5 count per million instead of 1 count per million 

#pca plot after transformation - unknown can be labled NA 



