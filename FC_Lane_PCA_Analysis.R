setwd("~/GitHub/Placenta_RNAseq_WGCN")
load("raw_norm_combat_seq.RData")
Glowing_Placenta_RNAseq_FCandLane <- read.csv("Glowing_Placenta_RNAseq_FCandLane.csv") 
nrow(Glowing_Placenta_RNAseq_FCandLane)
meta_data = read.csv("metaDataGlowing-for-TME.csv", header = T)
nrow(meta_data)
#Compare two vectors of different length in R
keep_IDs = which(meta_data$ID %in% Glowing_Placenta_RNAseq_FCandLane$Sample.ID == TRUE)
Glowing_Placenta_RNAseq_FCandLane2 = Glowing_Placenta_RNAseq_FCandLane[keep_IDs,]

#not all samples have "G-" in Glowing_Placenta_RNAseq_FCandLane$Sample.ID
grepl("G-", Glowing_Placenta_RNAseq_FCandLane$Sample.ID)
newIDs = which(grepl("G-", Glowing_Placenta_RNAseq_FCandLane$Sample.ID) ==  FALSE)
#sub("G", "G-",Glowing_Placenta_RNAseq_FCandLane$Sample.ID[newIDs])
for ( i in newIDs){
  Glowing_Placenta_RNAseq_FCandLane[i, "Sample.ID"] = sub("G", "G-", Glowing_Placenta_RNAseq_FCandLane[i, "Sample.ID"])
}

#G-130 has a large space 
# Glowing_Placenta_RNAseq_FCandLane$Sample.ID[38] == "G-130"
# [1] FALSE
# Glowing_Placenta_RNAseq_FCandLane$Sample.ID[38]
# [1] "G-130 "
Glowing_Placenta_RNAseq_FCandLane$Sample.ID[38] ="G-130"

colnames(Glowing_Placenta_RNAseq_FCandLane)[7] = "ID"
library("tidyverse")
meta_data = dplyr::left_join(meta_data, Glowing_Placenta_RNAseq_FCandLane[,c("FC","Lane","ID")], by = "ID")

meta_data$FC = as.factor(meta_data$FC)
meta_data$Lane = as.factor(meta_data$Lane)

#Conduct PC for glowing_uq_vst 
PCobj = prcomp(t(glowing_uq_vst), retx = T, center = T, scale. = T)
PCs = PCobj$x
#rownames(PCs) == meta_data$GlowRNAseqID #TRUE
#rownames(PCs) = meta_data$ID
meta_data["PC1_uq_vst"] = PCs[,1]
meta_data["PC1_uq_vst >= 25"] = (meta_data$PC1_uq_vst >= 25)
meta_data["PC2_uq_vst"] = PCs[,2]

x1 = ggplot(data = meta_data, aes(PC2_uq_vst, PC1_uq_vst ,colour = `PrepID`)) + geom_point(show.legend = F) + geom_text(data = subset(meta_data, `PC1_uq_vst >= 25` == TRUE & (`PrepID` == "prep2" | `PrepID` == "prep3" )),aes( label= ID), hjust=0,vjust=0 , show.legend = FALSE) + 
   labs(title = "By Batch (UQ and VST)", x = "PC2", y = "PC1") + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme_light()

#Conduct PC for combat
PCobj = prcomp(t(combat.adj), retx = T, center = T, scale. = T)
PCs = PCobj$x
meta_data["PC1"] = PCs[,1]
meta_data["PC2"] = PCs[,2]

x2= ggplot(data = meta_data, aes(PC2, PC1 ,colour = `PrepID`)) + geom_point() + geom_text(data = subset(meta_data, `PC1` >= 25 & (`PrepID` == "prep2" | `PrepID` == "prep3" )),aes( label= ID), hjust=0,vjust=0, show.legend = FALSE) + 
  labs(title = "By Batch (Combat)", x = "PC2", y = "PC1") + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))  + theme_light()

grid.arrange(x1, x2, ncol = 2)

ggplot(data = meta_data, aes(PC2_uq_vst, PC1_uq_vst ,colour = `PrepID`)) + geom_point() + geom_text(data = subset(meta_data, `PC1_uq_vst >= 25` == TRUE & (`PrepID` == "prep2" | `PrepID` == "prep3" )),aes( label= ID), hjust=0,vjust=0 , show.legend = FALSE) + 
  labs(title = "By Batch (UQ and VST)", x = "PC2", y = "PC1") + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme(legend.position = "bottom") + theme_light()

