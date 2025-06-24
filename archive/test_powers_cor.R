knitr::opts_chunk$set(echo = TRUE)
load("/Users/cynthiaperez/Library/Group Containers/UBF8T346G9.Office/Outlook/Outlook 15 Profiles/Main Profile/Files/S0/1/Attachments/raw_norm_combat_seq[2763521].RData")
load("~/GitHub/Planet_Epigenetic_GA_Clocks/Perez_Dis_PhenoData.RData")

colnames(combat.adj) = colnames(glowing_uq_vst)
library(WGCNA)
library(gprofiler2)
library(dplyr)
combat.adj = t(combat.adj)

powers=seq(from=2,to=30,by=1) 
chosen_powers = data.frame()
pickSoftThreshold_fit = list()

for( Rsquare in c(0.80, 0.90)){
  for (cors in c("Pearson", "Spearman", "Biweight")){
    
    if (cors == "Pearson"){
     cortype = "cor"
     options = list(use = "p")
  
    } else if (cors == "Spearman"){
      cortype = "cor"
      options = list(method = "spearman")
    } else {
      cortype = "bicor"
      options = list(use = "p")
    }
    
    #pick soft threshold 
    sft = pickSoftThreshold(combat.adj,RsquaredCut = Rsquare, dataIsExpr = T,  powerVector=powers, networkType ="signed", corFnc = cortype, corOptions = options, moreNetworkConcepts=T) 
    pickSoftThreshold_fit[[paste0("power_", Rsquare, "_", cors)]] = sft[["fitIndices"]]
    v = c(Rsquare, cors, sft[["powerEstimate"]])
    chosen_powers = rbind(chosen_powers,v )
  }
}
colnames(chosen_powers) = c("R.sq", "cor_name", "power")

#loop through network building 


k = list()
colorDynamicHybridTOM  = list()
for ( i in 1:nrow(chosen_powers)) {
  powers = as.numeric(chosen_powers$power[i])
if (chosen_powers$cor_name[i] == "Pearson"){
  ADJ1 = abs(cor(combat.adj,method = "pearson"))^powers
} else if (chosen_powers$cor_name[i] == "Spearman"){ 
  ADJ1 = abs(cor(combat.adj,method = "spearman"))^powers 
} else {
  ADJ1 = abs(bicor(combat.adj))^powers
}
  k[[paste0(chosen_powers$power[i],"_",chosen_powers$cor_name[i])]]=softConnectivity(datE=combat.adj,power= powers)
  dissTOM=TOMdist(ADJ1)
  hierTOM = hclust(as.dist(dissTOM),method="average")
  colorDynamicHybridTOM[[paste0(chosen_powers$power[i],"_",chosen_powers$cor_name[i])]] = labels2colors(cutreeDynamic(hierTOM, distM= dissTOM , cutHeight = 0.99,
                                                      deepSplit=2, pamRespectsDendro = FALSE))
} 

for ( i in 1:6){
  
  hist(k[[i]], ylab = "K Connections", main = paste("R square:", chosen_powers$R.sq[i]," Correlation Test:", chosen_powers$cor_name[i], " Power:", chosen_powers$power[i] ))
  
}
table(colorDynamicHybridTOM[[3]])

#eigengenessss
module_Eigengenes = moduleEigengenes(expr = combat.adj, colors = colorDynamicHybridTOM[[5]], impute = TRUE, nPC = 1, excludeGrey = TRUE)

#hub genes 
top_hub2 = chooseTopHubInEachModule(
  combat.adj, 
  colorDynamicHybridTOM[[5]], 
  omitColors = "grey", 
  power = 2, 
  type = "signed")


#PFAS and eigenenes 

pheno_153 = pheno 
pheno_153 = pheno_153 %>% arrange(participant_id)
load("~/GitHub/Placenta_RNAseq_WGCNA/PFAS_Placenta_RNAseq.RData")

module_Eigengenes_DynamicHybridTOM_subset = module_Eigengenes$eigengenes %>% rownames_to_column("ID") %>% filter(ID %in% pheno$ID) 
rownames(module_Eigengenes_DynamicHybridTOM_subset) = module_Eigengenes_DynamicHybridTOM_subset$ID
module_Eigengenes_DynamicHybridTOM_subset = module_Eigengenes_DynamicHybridTOM_subset[,-1]

pheno = pheno %>% filter(ID %in% rownames(module_Eigengenes_DynamicHybridTOM_subset)) %>% mutate(sum = rowSums(.[c("PFBSp","PFOSp","PFNAp","PFOAp","PFHxSp")]))
eigengene_PFAS = data.frame()
PFASp = c("PFBSp","PFOSp","PFNAp","PFOAp","PFHxSp", "sum")
for (i in 1:ncol(module_Eigengenes_DynamicHybridTOM_subset)) {
  for (x in PFASp) {
    lm = summary(lm(module_Eigengenes_DynamicHybridTOM_subset[,i] ~ pheno[,x]))$coefficients[2,4] 
   df = data.frame(p.value = lm, module = colnames(module_Eigengenes_DynamicHybridTOM_subset)[i], PFAS = x)
   eigengene_PFAS = rbind(eigengene_PFAS, df)
    #eigengene_color_PFAS[paste0(colnames(module_Eigengenes_DynamicHybridTOM_subset)[i],x)] = summary(lm(module_Eigengenes_DynamicHybridTOM_subset[,i] ~ pheno[,x]))$coefficients[2,4] 
  }
}

eigengene_PFAS = eigengene_PFAS %>% filter(.$p.value <= 0.05) %>% mutate(module = gsub("ME", "", .$module))



pheno_153 = pheno_153 %>%  filter(participant_id %in% rownames(module_Eigengenes$eigengenes))
eigengene_bw = data.frame()
for (x in 1:ncol(module_Eigengenes$eigengenes)) {
  lm = summary(lm(module_Eigengenes$eigengenes[,x] ~ pheno_153$birth_wt_kg + pheno_153$mom_age_at_birth + pheno_153$mom_education + pheno_153$gest_age_in_weeks_edd + pheno_153$MomBMI_36wks + pheno_153$mom_race + pheno_153$MomBMI_10wks))$coefficients[2,4] 
  df = data.frame(p.value = lm, module = colnames(module_Eigengenes_DynamicHybridTOM_subset)[x])
  eigengene_bw = rbind(eigengene_bw, df)
  }

eigengene_bw = eigengene_bw %>% filter(.$p.value <= 0.05)%>% mutate(module = gsub("ME", "", .$module))


#KEGG Analysis
PFAS_KEGG = list()
df_colorDynamicHybridTOM = data.frame(colorDynamicHybridTOM[[5]], colnames(combat.adj))
for (i in unique(eigengene_PFAS$module)){
  x = df_colorDynamicHybridTOM %>% filter(.$colorDynamicHybridTOM..5.. == i) %>% pull(colnames.combat.adj.)
  colnames(x) = NULL
  PFAS_KEGG[i] = gost(x, organism = "hsapiens", sources = "KEGG")
}

PFAS_bw = list() 
for (i in eigengene_bw$module){
  x = df_colorDynamicHybridTOM %>% filter(.$colorDynamicHybridTOM..5.. == i) %>% pull(colnames.combat.adj.)
  colnames(x) = NULL
  PFAS_bw[i] = gost(x, organism = "hsapiens", sources = "GO", significant = T)
}

TOM = TOMsimilarityFromExpr(combat.adj, power = 9)
modules = "lightgreen"
  # Select module probes
probes = colnames(combat.adj)
inModule = is.finite(match(colorDynamicHybridTOM[[5]], modules))
modProbes = probes[inModule]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               nodeAttr = colorDynamicHybridTOM[[5]][inModule])
