setwd("~/GitHub/Placenta_RNAseq_WGCNA")
load("raw_norm_combat_seq.RData")
Glowing_Placenta_RNAseq_FCandLane <- read.csv("Glowing_Placenta_RNAseq_FCandLane.csv") 
nrow(Glowing_Placenta_RNAseq_FCandLane)
meta_data = read.csv("metaDataGlowing-for-TME.csv", header = T)
nrow(meta_data)
#Compare two vectors of different length in R
#only 118 IDs from meta data are represented in Glowing_Placenta_RNAseq_FCandLane 
keep_IDs = which(meta_data$ID %in% Glowing_Placenta_RNAseq_FCandLane$Sample.ID == TRUE)


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

#PCA plot for uq and vst transformed data 
ggplot(data = meta_data, aes(PC2_uq_vst, PC1_uq_vst ,colour = `PrepID`)) + geom_point() + geom_text(data = subset(meta_data, `PC1_uq_vst >= 25` == TRUE & (`PrepID` == "prep2" | `PrepID` == "prep3" )),aes( label= ID), hjust=0,vjust=0 , show.legend = FALSE) + 
labs(title = "By Batch (UQ and VST)", x = "PC2", y = "PC1") + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme(legend.position = "bottom") + theme_light()

#Conduct PC for combat
PCobj = prcomp(t(combat.adj), retx = T, center = T, scale. = T)
PCs = PCobj$x
meta_data["PC1_combat"] = PCs[,1]
meta_data["PC2_combat"] = PCs[,2]

#Raw count sum/depth 
meta_data_cp = read.csv(file = "meta_data_cp_Sept8.csv", sep = ",",header = T)
raw_data = t(GlowingPlacenta_RNAseq.rawCounts.nonZero.141.samples_copy)
raw_data = cbind(raw_data, rep(0, nrow(raw_data)))

for(i in 1:nrow(raw_data)){
  raw_data[i, ncol(raw_data)] = sum(raw_data[i, ])
}

raw_data = data.frame(raw_data)
raw_data$ID = row.names(raw_data)
colnames(raw_data)[13832] = "raw_counts"

meta_data_cp = dplyr::left_join(meta_data_cp, raw_data[,c("ID", "raw_counts")], by = "ID")
#ggplot(data = subset(meta_data, PrepID == "prep2" & PrepID == "prep3"), aes(x = raw_counts)) + geom_histogram(color="black", fill="white")

ggplot(data = meta_data_cp, aes(x = raw_counts, color = PrepID, fill = PrepID)) + geom_histogram( alpha=0.5, position="identity") + facet_wrap(~ PrepID)

ggplot(data = meta_data_cp, aes(y = raw_counts,x = PrepID, color = PrepID)) + geom_boxplot()

ggplot(data = meta_data_cp, aes(x = 1:141, y = raw_counts, color = PrepID)) + geom_point() + scale_x_continuous(breaks = c(1, 25, 50, 75,100, 125, 141 ), labels = c("G-114", "G-148" ,"G-188" ,"G-257" ,"G-312" ,"G-382", "G-419")  ) + labs(x = "IDs (ordered by row number")

