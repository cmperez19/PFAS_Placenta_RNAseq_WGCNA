library(WGCNA)
library(gprofiler2)
library(dplyr)
library(biomaRt)

load("sex_comp_combat_adj_pheno.RData")

# Filter out low-quality samples
exclude_ids <- c("G-321", "G-413", "G-174")
pheno_150 <- pheno_150[!pheno_150$participant_id %in% exclude_ids, ]
rownames(pheno_150) <- pheno_150$participant_id
sex_comp_combat.adj <- sex_comp_combat.adj[, !colnames(sex_comp_combat.adj) %in% exclude_ids]

# Remove sex chromosome genes to explore network without sex genes 
#combat.adj <- t(sex_comp_combat.adj)
#gene_ids <- colnames(combat.adj)
#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#gene_info <- getBM(attributes = c('hgnc_symbol', 'chromosome_name'),
                  # filters = 'hgnc_symbol',
                   #values = gene_ids,
                   #mart = ensembl)
#sex_genes <- gene_info %>% 
 # filter(chromosome_name %in% c("X", "Y")) %>% 
  #pull(hgnc_symbol)
#combat.adj <- combat.adj[, !colnames(combat.adj) %in% sex_genes] 

# Construct network using Spearman correlation and beta = 5
beta <- 5
adjacency <- abs(cor(combat.adj, method = "spearman"))^beta
dissTOM <- TOMdist(adjacency)
hierTOM <- hclust(as.dist(dissTOM), method = "average")
colors <- labels2colors(cutreeDynamic(hierTOM, distM = dissTOM, 
                                       cutHeight = 0.99, deepSplit = 2, 
                                       pamRespectsDendro = FALSE))

# Plot connectivity histogram
connectivity <- softConnectivity(datE = combat.adj, power = beta)
hist(connectivity, xlab = "K Connections", ylab = "Frequency",
     main = "Connectivity Histogram\n(R^2 = 0.80, Spearman, Power = 5)")

# Save module assignments and eigengenes
gene_color_mrna <- data.frame(gene_name = colnames(combat.adj), color = colors)
ME_list <- moduleEigengenes(combat.adj, colors = colors, impute = TRUE, nPC = 1, excludeGrey = TRUE)
module_Eigengenes_mrna <- ME_list$eigengenes

save(gene_color_mrna, module_Eigengenes_mrna, sex_comp_combat.adj, pheno_150,
     file = "bmi_mrna_ME.Rdata")

# Export network for a specific module
target_module <- "lightgreen"
TOM <- TOMsimilarityFromExpr(combat.adj, power = beta)
in_module <- colors == target_module
mod_probes <- colnames(combat.adj)[in_module]
mod_TOM <- TOM[in_module, in_module]
dimnames(mod_TOM) <- list(mod_probes, mod_probes)

exportNetworkToCytoscape(mod_TOM,
                         edgeFile = paste0("CytoscapeInput-edges-", target_module, ".txt"),
                         nodeFile = paste0("CytoscapeInput-nodes-", target_module, ".txt"),
                         weighted = TRUE,
                         threshold = 0.02,
                         nodeNames = mod_probes,
                         nodeAttr = colors[in_module])


