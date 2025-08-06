library(tidyverse)
load("Perez_Dis_PhenoData.Rdata")
load("raw_norm_combat_seq.Rdata")
meta_data = read.csv(file="meta_data_cp_Sept9.csv", header = T)
rm(glowing_uq_vst)
rm(combat.adj)

pheno = pheno %>% filter(participant_id %in% meta_data$ID)%>% 
  arrange(participant_id) %>% mutate(prepID = meta_data$PrepID) 

prep1_pre2_ID = pheno %>% filter(prepID == "prep2"| prepID == "prep3") %>% select(participant_id)


GlowingPlacenta_RNAseq.rawCounts = GlowingPlacenta_RNAseq.rawCounts.nonZero.141.samples_copy %>% select(as.vector(t(prep1_pre2_ID))) 

Lean_ID = pheno %>% filter(prepID == "prep2"| prepID == "prep3") %>% select(participant_id, official_enroll_category) %>% filter(official_enroll_category == "Lean")
Lean_ID = as.vector(t(Lean_ID$participant_id))
Overweight_ID = pheno %>% filter(prepID == "prep2"| prepID == "prep3") %>% select(participant_id, official_enroll_category) %>% filter(official_enroll_category == "Overweight")
Overweight_ID = as.vector(t(Overweight_ID$participant_id))

GlowingPlacenta_RNAseq.rawCounts$X = row.names(GlowingPlacenta_RNAseq.rawCounts)
GlowingPlacenta_RNAseq.rawCounts = GlowingPlacenta_RNAseq.rawCounts %>% relocate(X)

Lean_ID_index = c()
for (x in 1:length(Lean_ID)){ 
  Lean_ID_index[x] = which(colnames(GlowingPlacenta_RNAseq.rawCounts) == Lean_ID[x])
}

Overweight_ID_index = c()
for (x in 1:length(Overweight_ID)){ 
  Overweight_ID_index[x] = which(colnames(GlowingPlacenta_RNAseq.rawCounts) == Overweight_ID[x])
}
rownames(GlowingPlacenta_RNAseq.rawCounts) = NULL
#write.table(GlowingPlacenta_RNAseq.rawCounts, file = "GlowingPlacenta_RNAseq.rawCounts-93-samples.txt", sep = "\t", row.names = F, col.names = T)
library(TRAPR)
Sample <- TRAPR.Data.ReadExpressionTable("GlowingPlacenta_RNAseq.rawCounts-93-samples.txt", sep = "\t", Exp1 = Lean_ID_index, Exp2 = Overweight_ID_index, Tag = c('Lean', 'Overweight'))

nSample <- TRAPR.Normalize(Sample, Method = "UpperQuartile")

vst_data = TRAPR.Transformation.VSN(nSample)



library(sva)

glowing_uq_vst = data.frame(vst_data$CurrentMatrix)
row.names(glowing_uq_vst) = vst_data[["CurrentGene"]]
colnames(glowing_uq_vst) = vst_data[["CurrentSample"]]

pheno = pheno %>% filter(prepID == "prep2"| prepID == "prep3")
batch = pheno$prepID
## Model matrix for batch-corrections (May need to adjust model matrix to 'protect' coefficients (study specific)):
mod <- model.matrix(~1 + official_enroll_category , data=pheno)

## Run ComBat to remove batch effects

combat.adj = ComBat(glowing_uq_vst,batch = batch, mod = mod)
PCobj = prcomp(t(combat.adj), retx = T, center = T, scale. = T)
PCs = PCobj$x
ggplot(data = data.frame(PCs), aes(PC2, PC1, color = pheno$prepID)) + geom_point() + labs(title = "UQ +VST +Combat for prep2 and prep3 (n=93)")

PCobj = prcomp(t(glowing_uq_vst), retx = T, center = T, scale. = T)
PCs = PCobj$x
ggplot(data = data.frame(PCs), aes(PC2, PC1, color = pheno$prepID)) + geom_point() + labs(title = "UQ +VST for prep2 and prep3 (n=93)")
