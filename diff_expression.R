library(tidyverse)
library(gprofiler2)
library(dplyr)
load("combat.adj.Rdata")
load("eigengnes_spearman_9.Rdata")
load("log_PFAS_pheno.Rdata")

#("G-155","G-157") are not in the pfas data 
combat.adj = t(combat.adj)
combat.adj = combat.adj[!row.names(combat.adj) %in% c("G-155","G-157"),]


pc <- prcomp(combat.adj)
pc_data = pc$x
pc_data = pc_data[pheno_139$participant_id[which(pheno_139$childs_sex == "Male")],"PC1"]

PFASp = c("L-PFOA","PFHxS","PFNA","PFDA","s-PFOS")
log_pfas = log_pfas[pheno_139$participant_id[which(pheno_139$childs_sex == "Male")],PFASp]
pheno_139 = pheno[rownames(combat.adj),]
pheno_139 = pheno_139[pheno_139$participant_id[which(pheno_139$childs_sex == "Male")], ]
DEG = data.frame()

  for(PFAS in colnames(log_pfas)){
    for (gene in colnames(combat.adj)) {
  # Fit linear regression model
  model <- lm(combat.adj[pheno_139$participant_id[which(pheno_139$childs_sex == "Male")],gene] ~ log_pfas[,PFAS] + pheno_139$mom_age_at_birth 
              + pheno_139$mom_education + pheno_139$gest_age_in_weeks_edd 
              + pheno_139$MomBMI_36wks +  pheno_139$MomBMI_10wks 
              + pheno_139$official_enroll_category+ pc_data)
  regression_results <- summary(model)
  p_value <- regression_results$coefficients[2,4]
  estimate <- coef(model)[2]
  std_err <- regression_results$coefficients[2,2]
  df <- data.frame(p.value = p_value, estimate = estimate, gene_name = gene, std.err = std_err, PFAS = PFAS)
 DEG <- rbind(DEG, df)
    } 
  }

#q value for multiple testing corrections
#DEG2 <- DEG %>% group_split(.$PFAS) 
#x = p.adjust(DEG2[[1]]$p.value, method = "fdr", n = 13831)

DEG <- DEG %>%
  group_by(PFAS) %>%
  mutate(adjusted_p_value = p.adjust(p.value, method = "fdr", n = n()))


#sum(DEG[which(DEG$PFAS == "L-PFOA"), "adjusted_p_value"] == x)




# Create the MA plot

threshold = 0.05

DEG_PFNA <- ggplot(DEG[which(DEG$PFAS == "PFNA"),], aes(x = estimate, y = -log10(p.value))) +
  geom_point(data = subset(DEG[which(DEG$PFAS == "PFNA"),], adjusted_p_value <= threshold), aes(color = "Significant"), size = 3) +
  geom_point(data = subset(DEG[which(DEG$PFAS == "PFNA"),], adjusted_p_value > threshold), aes(color = "Not Significant"), size = 3) +
  scale_color_manual(values = c("Significant" = "pink", "Not Significant" = "grey")) +
  labs(title = "PFNA Volcano Plot (Female)",
       x = "Coefficient Estimate",
       y = "-log10(P-Value)",
       color = "FDR q-Value < 0.05") +
  theme_minimal()



DEG_PFNA + geom_label_repel(data = subset(DEG[which(DEG$PFAS == "PFNA"),], adjusted_p_value<= threshold), aes(label = gene_name),
                 box.padding   = 0.5, 
                 segment.color = 'grey50') 


#VDR SNP frequency 

for (snp_col in colnames(VDR_SNPs)[-1]) {
  cat("Frequency table for", snp_col, ":\n")
  allele_counts <- table(unlist(strsplit(VDR_SNPs[[snp_col]], "")))
  allele_frequencies <- allele_counts / sum(allele_counts)
  df_freq <- data.frame(
    Allele = names(allele_counts),
    Frequency = allele_frequencies
  )
  print(df_freq)
  cat("\n")
}



#snp and gene expression 
Glowing_VDR_SNPs_20230922 <- read.delim("~/Downloads/Glowing_VDR_SNPs_20230922.txt")
VDR_SNPs = Glowing_VDR_SNPs_20230922 %>% filter(.$participant_id %in% rownames(combat.adj))
summary(lm(combat.adj[VDR_SNPs$participant_id,"VDR"] ~ VDR_SNPs[,2] + VDR_SNPs[,3] + VDR_SNPs[,4] + VDR_SNPs[,5]))

#eigengene and VDR snps 
summary(lm(eigengenes[VDR_SNPs$participant_id, "MEdarkred"] ~ factor(VDR_SNPs[,2]) + factor(VDR_SNPs[,3]) + factor(VDR_SNPs[,4]) + factor(VDR_SNPs[,5])))

#summary(lm(eigengenes[VDR_SNPs$participant_id, "MEdarkred"] ~ factor(VDR_SNPs[,5]) + factor(VDR_SNPs[,5])*log_pfas[VDR_SNPs$participant_id,4]))


for( i in colnames(log_pfas)){
  for (j in 2:5){
    print(j)
    print(i)
  x = summary(lm(eigengenes[VDR_SNPs$participant_id, "MEdarkred"] ~  log_pfas[VDR_SNPs$participant_id,i] * factor(VDR_SNPs[,j])))
  print(x)
}}


#BMKR 
library(bkmr)
#mixture<-as.matrix(log_pfas[complete.cases(pheno_139$MomBMI_36wks),])

VDR_SNPs$rs7975232.Apa1 = ifelse(VDR_SNPs$rs7975232.Apa1 == "AA", 0, ifelse(VDR_SNPs$rs7975232.Apa1 == "AC", 1, 2))
VDR_SNPs$rs1544410_Bsm1 = ifelse(VDR_SNPs$rs1544410_Bsm1 == "CC", 0, ifelse(VDR_SNPs$rs1544410_Bsm1 == "CT", 1, 2))
VDR_SNPs$rs2228570.Fok1 = ifelse(VDR_SNPs$rs2228570.Fok1 == "AA", 0, ifelse(VDR_SNPs$rs2228570.Fok1 == "GA", 1, 2))
VDR_SNPs$rs731236.Taq1 = ifelse(VDR_SNPs$rs731236.Taq1 == "AA", 0, ifelse(VDR_SNPs$rs731236.Taq1 == "AG", 1, 2))
mixture = as.matrix(VDR_SNPs[,-1])
rownames(mixture) <- VDR_SNPs$participant_id


y <- eigengenes[VDR_SNPs$participant_id,"MEdarkred"]
# model <- lm(combat.adj[,gene] ~ log_pfas[,PFAS] + pheno_139$mom_age_at_birth 
#+ pheno_139$mom_education + pheno_139$gest_age_in_weeks_edd 
#+ pheno_139$MomBMI_36wks +  pheno_139$childs_sex + pheno_139$MomBMI_10wks 
#+ pheno_139$official_enroll_category+ pc_data)
#G-
covariates<- pheno_139[VDR_SNPs$participant_id,c("mom_age_at_birth", "mom_education","gest_age_in_weeks_edd","MomBMI_36wks", "childs_sex","MomBMI_10wks",
                                   "official_enroll_category")]
#G-257 has to be removed because it has an NA in mombmi_36wks
#covariates need to be numeric
pc = data.frame(pc$x)
covariates$PC1 = pc[VDR_SNPs$participant_id,]
covariates$official_enroll_category = ifelse(covariates[, 7] == "Lean", 0, 1)
covariates$childs_sex = ifelse(covariates$childs_sex == "Female", 0, 1)

covariates$mom_education = as.numeric(gsub("[^0-9]", "", covariates$mom_education))

mixture <- mixture[complete.cases(covariates$MomBMI_36wks),]
y <- y[-73]
covariates <- covariates[complete.cases(covariates$MomBMI_36wks), ]



set.seed(10)
knots100  <- fields::cover.design(mixture, nd = 50)$design
temp <-  kmbayes(y=y, Z=mixture, X=covariates, iter=1000, verbose=FALSE, 
                 varsel=TRUE, knots=knots100)


#cor plot PFAS and VDR expression 
cor_pfas = cor(cbind(log_pfas,combat.adj[,"VDR"]), method = "spearman")

cor_melt <- melt(cor_pfas)

# Create heat map
ggplot(cor_melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + geom_text(aes(label = round(value, 2)), color = "black", size = 4)

#bar chart of module freq colors with top REVIGO term 
modules_table = data.frame(table(gene_color$color))
modules_table = modules_table[c(8, 1:7, 9:22),]
modules_table = modules_table %>% slice(2:nrow(modules_table)) %>% arrange(desc(Freq)) %>% rbind(slice(modules_table,1), .)
modules_table = modules_table[c(2:22, 1),]

modules_table$REVIGO_term = c("regulation of protein stability", "snRNA 3'-end processing", "regulation of localization", "vesicle-mediated transport", "regulation of vascular endothelial cell proliferation",
                "regulation of locomotion", "cell cycle G2/M phase transition", "regulation of Ras protein signal transduction", "peptidyl-amino acid modification", "transferrin transport", "cellular response to stimulus",
                "maintenance of cell number", "ncRNA metabolic process", "regulation of cell adhesion", "regulation of immune system process", "cellular response to chemical stimulus", "regulation of cell adhesion", 
                "negative regulation of serotonin secretion", "regulation of cell killing", "regulation of transepithelial transport", "translational initiation", " ")

ggplot(modules_table, aes(Freq, Var1)) + scale_y_discrete(limits = rev(modules_table$Var1)) + geom_bar(stat = "identity", fill = modules_table$Var1) +
  geom_text(aes(label= REVIGO_term),size = 6, hjust = -0.01)  + xlim(c(0,7000)) + theme_minimal(base_size = 20) + xlab("Number of Genes") + ylab("Module")


#hub genes 
x = chooseTopHubInEachModule(datExpr = combat.adj, colorh = gene_color$color, power = 9)
summary(lm(combat.adj[VDR_SNPs$participant_id,"UTY"] ~ factor(VDR_SNPs[,3])))
summary(lm(combat.adj[VDR_SNPs$participant_id,"UTY"] ~ factor(VDR_SNPs[,5])))



#KEGG and GO gene set enrichment analysis saved to excel sheets 
library(gprofiler2)
library(openxlsx)
wb <- createWorkbook()
for (i in unique(gene_color$color)[-1]) {
  y <- gene_color[which(gene_color$color == i), "gene_name"]
  x <- gost(query = y, organism = "hsapiens", sources = "KEGG")
  addWorksheet(wb, paste(i))
  writeData(wb, paste(i), x[["result"]] )
}





saveWorkbook(wb, "KEGG_Terms_Modules.xlsx", overwrite = TRUE)


nGenes = ncol(combat.adj)
nSamples = nrow(combat.adj)

#Recalculating MEs with label colors
MEs0 = eigengenes
MEs = orderMEs(MEs0)
datTraits = pheno_139[,c("gest_age_in_weeks_edd",  "birth_wt_kg", "birth_len_cm")]
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


# Assuming you have previously calculated moduleTraitCor and moduleTraitPvalue

# Identify indices where p-values are not significant
non_sig_indices <- moduleTraitPvalue > 0.05  # Adjust the significance level as needed

# Update moduleTraitCor based on non-significant p-values
moduleTraitCor[non_sig_indices] <- NA  # Replace non-significant values with NA

# Generating text matrix with significant values and empty strings for non-significant ones
textMatrix <- ifelse(non_sig_indices, "", paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = ""))

# Reshape textMatrix to match dimensions of moduleTraitCor
dim(textMatrix) <- dim(moduleTraitCor)

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

#Displaying the correlation values in a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = c("Gestational Age", "Birth Weight (kg)", "Birth Length (cm)"),
               yLabels = gsub("ME", "", names(MEs)),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


summary(lm(pheno_139$birth_wt_kg ~ eigengenes$MEtan  +pheno_139$mom_age_at_birth + pheno_139$mom_education + pheno_139$gest_age_in_weeks_edd + pheno_139$MomBMI_36wks +  pheno_139$childs_sex + pheno_139$MomBMI_10wks + pheno_139$official_enroll_category+ pc_data))
