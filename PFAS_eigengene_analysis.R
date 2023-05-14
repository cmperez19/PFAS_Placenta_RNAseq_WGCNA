load("~/Github/Placenta_RNAseq_WGCNA/raw_norm_combat_seq.RData")
load("test_poweers.RData")
colnames(combat.adj) = colnames(glowing_uq_vst)
library(WGCNA)
library(gprofiler2)
library(dplyr)
combat.adj = t(combat.adj)
module_Eigengenes = moduleEigengenes(expr = combat.adj, colors = colorDynamicHybridTOM[[5]], impute = TRUE, nPC = 1, excludeGrey = TRUE)
module_Eigengenes = module_Eigengenes$eigengenes
module_Eigengenes = module_Eigengenes[which(row.names(module_Eigengenes) %in% pfas_df_cp$Name), ]
PFASp = c("L-PFOA","PFHxS","PFNA","PFDA","s-PFOS","sum1", "sum2" )
pheno = pheno[row.names(module_Eigengenes),]
eigengene_PFAS = data.frame()
#139 profiles have RNA-seq data and PFAS 
#add covariates to model from proposal (GA, maternal age, maternal weight gain, maternal education, PC1 from RNA-seq)
#factor analysis between PFAS chemicals
#can explore hub genes and TF that are annotated within them
#maternal age, maternal education, gestational age, maternal BMI, gestational weight gain, and ancestry probability
for (i in 1:ncol(module_Eigengenes)) {
  for (x in PFASp){
    lm = summary(lm(module_Eigengenes[,i] ~ pfas_df_cp[rownames(module_Eigengenes),x] + pheno$mom_age_at_birth + pheno$mom_education + pheno$gest_age_in_weeks_edd + pheno$MomBMI_10wks + pheno$Prob_Caucasian + pheno$child_sex))$coefficients[2,4] 
    df = data.frame(p.value = lm, module = colnames(module_Eigengenes)[i], PFAS = x)
    eigengene_PFAS = rbind(eigengene_PFAS, df)
    #eigengene_color_PFAS[paste0(colnames(module_Eigengenes_DynamicHybridTOM_subset)[i],x)] = summary(lm(module_Eigengenes_DynamicHybridTOM_subset[,i] ~ pheno[,x]))$coefficients[2,4] 
  }
}

eigengene_PFAS = eigengene_PFAS %>% filter(.$p.value <= 0.05) %>% mutate(module = gsub("ME", "", .$module))


#KEGG Analysis
PFAS_GO = list()
df_colorDynamicHybridTOM = data.frame(colorDynamicHybridTOM[[5]], colnames(combat.adj))
for (i in unique(eigengene_PFAS$module)){
  x = df_colorDynamicHybridTOM %>% filter(.$colorDynamicHybridTOM..5.. == i) %>% pull(colnames.combat.adj.)
  colnames(x) = NULL
  print(x)
  PFAS_GO[i] = gost(x, organism = "hsapiens", sources = "GO")
}

PFAS_KEGG = list()
for (i in unique(eigengene_PFAS$module)){
  x = df_colorDynamicHybridTOM %>% filter(.$colorDynamicHybridTOM..5.. == i) %>% pull(colnames.combat.adj.)
  colnames(x) = NULL
  PFAS_KEGG[i] = gost(x, organism = "hsapiens", sources = "KEGG")
}

save(module_Eigengenes,combat.adj,pfas_df_cp,pheno,pfas.p, file="april24_input.Rda")


