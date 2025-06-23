library(tidyverse)
load("sexcomplement_raw_counts_pheno150.RData")
#Trapr requires a tab delimited text file where the first column is the sample list and the rows are genes
sex_complement_raw_counts$Genes = rownames(sex_complement_raw_counts)
sex_complement_raw_counts= sex_complement_raw_counts%>% relocate(Genes) 
rownames(sex_complement_raw_counts) = NULL
write.table(sex_complement_raw_counts, file = "sex_complement_raw_counts.txt", sep = "\t", row.names = F, col.names = T)

#exp 1 and exp 2 can be used to label the dichotomized samples in Trapr
pheno_150 = pheno_150 %>% filter(participant_id %in% colnames(sex_complement_raw_counts))
Lean_ID = pheno_150 %>% select(participant_id, official_enroll_category) %>% filter(official_enroll_category == "Lean")
Lean_ID = as.vector(t(Lean_ID$participant_id))
Overweight_ID = pheno_150 %>% select(participant_id, official_enroll_category) %>% filter(official_enroll_category == "Overweight")
Overweight_ID = as.vector(t(Overweight_ID$participant_id))

#have to label the index of the columns that are lean or overweight
Lean_ID_index = c()
for (x in 1:length(Lean_ID)){ 
  Lean_ID_index[x] = which(colnames(sex_complement_raw_counts) == Lean_ID[x])
}

Overweight_ID_index = c()
for (x in 1:length(Overweight_ID)){ 
  Overweight_ID_index[x] = which(colnames(sex_complement_raw_counts) == Overweight_ID[x])
}

#Upper Quartile and VSN 
library(TRAPR)
Sample <- TRAPR.Data.ReadExpressionTable("sex_complement_raw_counts.txt", sep = "\t", 
                                         Exp1 = Lean_ID_index, Exp2 = Overweight_ID_index, Tag = c('Lean', 'Overweight'))
nSample <- TRAPR.Normalize(Sample, Method = "UpperQuartile")
vst_data = TRAPR.Transformation.VSN(nSample)

library(sva)
#creating a dataframe for UQ and VST data 
glowing_uq_vst = data.frame(vst_data$CurrentMatrix)
row.names(glowing_uq_vst) = vst_data[["CurrentGene"]]
colnames(glowing_uq_vst) = vst_data[["CurrentSample"]]

#prep data is not in pheno_150 data for all 153 samples 
metaDataGlowing <- read.csv("~/WGCNA_PFAS/metaDataGlowing-for-TME.csv")
colnames(metaDataGlowing)[4] = "participant_id"
#adds NAs to IDs that do not have prepID 
pheno_150 = merge(pheno_150, metaDataGlowing[,c("participant_id","PrepID")], by = "participant_id", all.x = TRUE)

## Model matrix for PrepID-corrections (May need to adjust model matrix to 'protect' coefficients (study specific)):
mod <- model.matrix(~1 + official_enroll_category , data=pheno_150)

## Run ComBat to remove PrepID effects
pheno_150$PrepID[is.na(pheno_150$PrepID)] <- 0
batch = pheno_150$PrepID
sex_comp_combat.adj = ComBat(glowing_uq_vst, batch= batch, mod = mod)

PCobj = prcomp(t(sex_comp_combat.adj), retx = T, center = T, scale. = T)
PCs = PCobj$x
#ggplot(data = data.frame(PCs), aes(PC2, PC1, color = pheno_150$PrepID)) + geom_point() + labs(title = "UQ +VST +Combat")
a = ggplot(data = data.frame(PCs), aes(PC1, PC2, color = pheno_150$PrepID, label=pheno_150$participant_id)) + 
  geom_point() +
  scale_colour_discrete(labels = c("NA", "Prep 1", "Prep 2", "Prep 3"))+
  #geom_text(aes(label=ifelse(PC1>200, as.character(pheno_150$participant_id),''),hjust=0,vjust=0)) + 
  labs(title = "UQ +VST +Combat", color = "Batch ID") + 
  ylab("PC2 (10.148%)") + 
  xlab("PC1 (30.132%)")



b = ggplot(data = data.frame(PCs), aes(PC1, PC2, color = pheno_150$childs_sex, label=pheno_150$participant_id)) + 
  geom_point() +
 # geom_text(aes(label=ifelse(PC1>200, as.character(pheno_150$participant_id),''),hjust=0.5,vjust=0.5)) + 
  labs(color = "Sex of Neonate") + 
  ylab("PC2 (10.148%)") + 
  xlab("PC1 (30.132%)")

library(gridExtra)
grid.arrange(a, b) 

#PCobj = prcomp(t(glowing_uq_vst), retx = T, center = T, scale. = T)
#PCs = PCobj$x
#ggplot(data = data.frame(PCs), aes(PC2, PC1, color = pheno_150$PrepID)) + geom_point() + labs(title = "UQ +VST")

#saving variable sex_comp_combat.adj
#save(sex_comp_combat.adj, file ="sex_comp_combat.adj.Rdata")

colnamessexcomp = colnames(sex_comp_combat.adj)
sex_comp_combat.adj = data.frame(sex_comp_combat.adj)
colnames(sex_comp_combat.adj) = colnamessexcomp
load("combat.adj.Rdata")
plot(t(combat.adj["XIST",]), t(sex_comp_combat.adj["XIST",colnames(combat.adj)]), xlab = "Original UQ + VST + Combat RNAseq",
     ylab = "Sex Complement UQ + VST + Combat RNAseq", main = "XIST")

plot(t(combat.adj["DDX3Y",]), t(sex_comp_combat.adj["DDX3Y",colnames(combat.adj)]), xlab = "Original UQ + VST + Combat RNAseq",
     ylab = "Sex Complement UQ + VST + Combat RNAseq", main = "DDX3Y")


# Load necessary library
library(ggplot2)
colnamesdefault = colnames(combat.adj)
combat.adj = data.frame(combat.adj)
colnames(combat.adj) = colnamesdefault
# Extract the relevant data
original_data <- t(combat.adj["EPPK1", ])
sex_comp_data <- t(sex_comp_combat.adj["EPPK1", colnames(combat.adj)])

# Create a data frame for ggplot2
data <- data.frame(
  Original = original_data,
  SexComplement = sex_comp_data, 
  Sex = pheno_150[colnames(combat.adj),"childs_sex"]
)
colnames(data) = c("OriginalDefault", "SexComplement", "Sex")
# Create the plot
ggplot(data, aes(x = OriginalDefault, y = SexComplement)) +
  geom_point(aes(color = Sex)) +
  labs(
    x = "Original UQ + VST + Combat RNAseq",
    y = "Sex Complement UQ + VST + Combat RNAseq",
    title = "EPPK1"
  ) +
  theme_minimal()


plot(t(combat.adj["ZFY",]), t(sex_comp_combat.adj["ZFY",colnames(combat.adj)]), xlab = "Original UQ + VST + Combat RNAseq",
     ylab = "Sex Complement UQ + VST + Combat RNAseq", main = "ZFY")




#genes SRY, UTX, are not found in both datasets 
sex_genes = t(sex_comp_combat.adj[c("DDX3X", "DDX3Y","USP9X","USP9Y", "ZFX", "ZFY"),colnames(combat.adj)])
sex_genes = data.frame(sex_genes)
sex_genes$participant_id = rownames(sex_genes)
rownames(pheno_150) = pheno_150$participant_id
sex_genes = merge(sex_genes,pheno_150[colnames(combat.adj),c("childs_sex","participant_id")],by="participant_id")


library(patchwork)

# Prepare the individual plots
DDX3xy <- sex_genes %>% 
  select(-participant_id) %>%
  tidyr::pivot_longer(starts_with("DD"), names_to = "genes", values_to = "expression") %>%
  ggplot(aes(genes, expression)) + 
  geom_jitter(aes(color=childs_sex)) + 
  labs(x = "", y = "Expression (UQ + VST + Combat)") + 
  theme(legend.position = "bottom")

ZFxy <- sex_genes %>% 
 select(-participant_id) %>%
  tidyr::pivot_longer(starts_with("ZF"), names_to = "genes", values_to = "expression") %>%
  ggplot(aes(genes, expression)) + 
  geom_jitter(aes(color=childs_sex)) + 
  labs(x = "", y = "") + 
  theme(legend.position = "none")

USP9xy <- sex_genes %>% 
  select(-participant_id) %>%
  tidyr::pivot_longer(starts_with("USP9"), names_to = "genes", values_to = "expression") %>%
  ggplot(aes(genes, expression)) + 
  geom_jitter(aes(color=childs_sex)) + 
  labs(x = "", y = "") + 
  theme(legend.position = "none")

# Extract the legend
legend <- cowplot::get_legend(DDX3xy)

# Combine the plots
combined_plot <- (DDX3xy + theme(legend.position = "none") | ZFxy |
                    USP9xy) + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')

# Display the combined plot with the legend
combined_plot + plot_annotation(title = "Sex Chromosome Complement Reference Genome (N=141)") + plot_layout(nrow = 1, heights = c(10, 1))

#DDX3xy = t(combat.adj[c("DDX3X", "DDX3Y"),])
#DDX3xy = data.frame(DDX3xy)
#DDX3xy$participant_id = rownames(DDX3xy)
#rownames(pheno_150) = pheno_150$participant_id
#DDX3xy = merge(DDX3xy,pheno_150[colnames(combat.adj),c("childs_sex","participant_id")],by="participant_id")
sex_genes = t(combat.adj[c("DDX3X", "DDX3Y","USP9X","USP9Y", "ZFX", "ZFY"),])
sex_genes = data.frame(sex_genes)
sex_genes$participant_id = rownames(sex_genes)
rownames(pheno_150) = pheno_150$participant_id
sex_genes = merge(sex_genes,pheno_150[colnames(combat.adj),c("childs_sex","participant_id")],by="participant_id")

DDX3xy <- sex_genes %>% 
  select(-participant_id) %>%
  tidyr::pivot_longer(starts_with("DD"), names_to = "genes", values_to = "expression") %>%
  ggplot(aes(genes, expression)) + 
  geom_jitter(aes(color=childs_sex)) + 
  labs(x = "", y = "Expression (UQ + VST + Combat)") + 
  theme(legend.position = "bottom")



USP9xy <- sex_genes %>% 
  select(-participant_id) %>%
  tidyr::pivot_longer(starts_with("USP9"), names_to = "genes", values_to = "expression") %>%
  ggplot(aes(genes, expression)) + 
  geom_jitter(aes(color=childs_sex)) + 
  labs(x = "", y = "") + 
  theme(legend.position = "none")

# Extract the legend
legend <- cowplot::get_legend(DDX3xy)

# Combine the plots
combined_plot <- (DDX3xy + theme(legend.position = "none") | 
                    USP9xy) + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')

# Display the combined plot with the legend
combined_plot + plot_annotation(title = "Default Reference Genome (N = 141)") + plot_layout(nrow = 1, heights = c(10, 1))

#rdata file. pheno_150 should include preps 
save(pheno_150, sex_comp_combat.adj, file = "sex_comp_combat_adj_pheno.RData")
