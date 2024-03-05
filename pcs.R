#scree plot of dark red
darkred = gene_color[which(gene_color$color == "darkred"),"gene_name"]

pcs_darkred = prcomp(combat.adj[,darkred])
variance_explained = pcs_darkred$sdev^2
scree_data <- data.frame(
  PC = 1:length(variance_explained),
  Variance_Explained = variance_explained
)

ggplot(scree_data, aes(x = PC, y = Variance_Explained)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  labs(
    x = "Principal Component",
    y = "Variance Explained",
    title = "Scree Plot Dark Red"
  ) +
  theme_minimal()


library(gprofiler2)

darkred_GO = gost(query= darkred, organism = "hsapiens", multi_query = TRUE, sources = "GO", user_thresholdap = 0.05)

#gene files 
for (color in unique(gene_color$color)){ 
  genes = gene_color[which(gene_color$color == color),"gene_name"]
  write.table(genes, file = paste(color, ".txt"), quote = F, sep = "/t", col.names = F, row.names = F)
}

#PC 1 and PC2 plot sex for darkred 


library(ggfortify) 
darkred.pca.plot <- autoplot(pcs_darkred, 
                          data = pheno_139, label = TRUE, 
                          colour = 'childs_sex')



library(tidyverse)
library(readxl)

# Set your downloads folder
downloads_folder <- "/Users/cynthiaperez/Downloads"

# Get a list of files in the downloads folder with the specified pattern
file_list <- list.files(downloads_folder, pattern = " _EnrichResults_2023-12-27.csv", full.names = TRUE)

# Create an Excel workbook
excel_wb <- createWorkbook()

# Loop through each file
for (file_path in file_list) {
  
  # Read the CSV file into a data frame
  df <- read_csv(file_path)
  
  # Filter the data based on the condition (e.g., significant p-value)
  df_filtered <- df %>%
    filter(FisherAdjustP <= 0.05) %>% dplyr::select(TF, NTargets, SigGenes, FisherTest, FisherAdjustP) # Adjust this condition based on your requirement
  
  # Extract the file name for use as a sheet name
  sheet_name <- tools::file_path_sans_ext(basename(file_path)) %>%
                    str_remove(" _EnrichResults_2023-12-27")
  # Add the filtered data to the Excel workbook
  
  addWorksheet(excel_wb, sheetName = sheet_name)
  writeData(excel_wb, sheet = sheet_name, df_filtered)
}

# Save the Excel workbook
saveWorkbook(excel_wb, "/Users/cynthiaperez/WGCNA_PFAS/TF_enrichment.xlsx")

