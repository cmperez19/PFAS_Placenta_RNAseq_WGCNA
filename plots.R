
library(ggplot2)
library(FactoMineR)

# Assuming 'gene_data' is your data frame of gene expression data
# Perform PCA
pca_result <- PCA(combat.adj[pheno_139$participant_id,gene_color[which(gene_color$color == "darkred"),"gene_name"]], scale.unit = TRUE, ncp = 2, graph = FALSE)

# Extract loadings for the first principal component (PC1)
loadings <- as.data.frame(pca_result$var$coord[,1])

# Add gene names as a new column
loadings$gene <- rownames(loadings)

# Rename columns for clarity
colnames(loadings) <- c("loading", "gene")

# Create a ggplot of the loadings
ggplot(loadings, aes(x = gene, y = loading)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Gene") +
  ylab("Loading on PC1 (Eigengene)") +
  ggtitle("Loadings of Dark Red Module Genes on the First Principal Component")


loadings$abs_loadings <- abs(loadings$loading)

# Optionally, you can sort the loadings to see the genes with the highest contributions
sort(loadings$abs_loadings, decreasing = TRUE)

# Set a threshold for significant contribution
# For example, consider top 10% as significant
threshold <- quantile(abs_loadings, 0.9)

# Identify genes with significant contribution
significant_genes <- names(abs_loadings[abs_loadings >= threshold])

# Print the significant genes
print(significant_genes)


color_palette <- colorRampPalette(c("green", "black", "red"))(100)
sex_data <- factor(pheno_139$childs_sex)
annotation_df <- data.frame(Sex = sex_data)
rownames(annotation_df) <- pheno_139$participant_id
sex_colors <- list(Sex= c(Male = "blue", Female = "pink"))
gene_data <- t(combat.adj[pheno_139$participant_id,gene_color[which(gene_color$color == "darkred"),"gene_name"]])


# Generate the heatmap
pheatmap(gene_data, 
         annotation_col = annotation_df, 
         annotation_colors = sex_colors,
         scale = "row", 
         color = color_palette, 
         clustering_method = "complete",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = TRUE,
         show_colnames = TRUE)


#modules that interacted with PFAS what is their hub gene - is it sexually dimorphic?
ME_sex = c("black", "blue", "brown", "green", "greenyellow","grey60", "lightyellow","magenta","midnightblue","pink","red","turquoise")

hub_genes = chooseTopHubInEachModule(datExpr = combat.adj, colorh = gene_color$color, power = 9)[ME_sex]

for( i in hub_genes){
  x =  t.test(combat.adj[which(pheno_139[,"childs_sex"] == "Male"), i], combat.adj[which(pheno_139[,"childs_sex"] == "Female"), i] )
  if(x[["p.value"]] < 0.05){
    print(i)
    print(x) 
  }
}

for ( i in ME_sex){
  color = gene_color %>% filter(color == i) %>% pull(gene_name)
  TFN = TRNModel %>% filter(.$targetGene %in% color) 
  sorted_df = sort(table(TFN$TF), decreasing = TRUE)
  for( j in names(sorted_df[1])){
    x =  t.test(combat.adj[which(pheno_139[,"childs_sex"] == "Male"), j], combat.adj[which(pheno_139[,"childs_sex"] == "Female"), j] )
    if(x[["p.value"]] < 0.05){
      print(i)
      print(j)
      print(x) 
    }
  }
  
}
#is the top TF of this module also sexually dimorphic? 



library(ggplot2)
library(sjPlot)
library(sjmisc)
df = pheno_139 %>% dplyr::select(mom_age_at_birth, mom_education, gest_age_in_weeks_edd, MomBMI_36wks, childs_sex, MomBMI_10wks, official_enroll_category) %>% rename(., sex = childs_sex) %>%
  mutate(MEdarkred = eigengenes$MEdarkred, pc = pc_data,  PFNA = log_pfas$PFNA) 
                                                                                                                       
fit <- lm(MEdarkred ~ PFNA + mom_age_at_birth + mom_education + gest_age_in_weeks_edd + MomBMI_36wks + sex + MomBMI_10wks + official_enroll_category + PFNA:sex, data = df)

plot_model(fit, type = "int")
