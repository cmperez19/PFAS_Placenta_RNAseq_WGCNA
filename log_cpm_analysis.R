library("edgeR")
GlowingPlacenta_RNAseq.rawCounts.nonZero.141.samples_copy = t(GlowingPlacenta_RNAseq.rawCounts.nonZero.141.samples_copy)
log_cpm = cpm(GlowingPlacenta_RNAseq.rawCounts.nonZero.141.samples_copy, log = T)


library(sva)

glowing_uq_vst = data.frame(vst_data$CurrentMatrix)
row.names(glowing_uq_vst) = vst_data[["CurrentGene"]]
colnames(glowing_uq_vst) = meta_data_cp$ID


library("tidyverse")
setwd("C:/Users/CMPERE2/Downloads")
meta_data = read.csv("meta_data_cp_Sept9.csv", sep = ",", header = T)
pheno_141 = pheno %>% filter(participant_id %in% colnames(combat.adj))%>% 
  arrange(participant_id) %>% mutate(prepID = meta_data$PrepID) 
  
batch = pheno_141$prepID
## Model matrix for batch-corrections (May need to adjust model matrix to 'protect' coefficients (study specific)):
mod <- model.matrix(~1 + official_enroll_category, data=pheno_141)


## Run ComBat to remove batch effects

combat.adj = ComBat(t(log_cpm),batch = batch, mod = mod)


#PCA stuff
PCobj = prcomp(t(combat.adj), retx = T, center = T, scale. = T)
PCs = PCobj$x
ggplot(data = data.frame(PCs), aes(PC2, PC1, color = pheno_141$prepID)) + geom_point() + labs(title = "Log(CPM)")
