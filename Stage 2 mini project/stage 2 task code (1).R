#Step 3-6 (Heatmap and Generation of upregulated and downregulated data)

# Install dplyr
install.packages("dplyr")

#load library
library(dplyr)

# Load the dataset from the URL 
url <- "https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv"
glioblastoma_data <- read.csv(url, row.names = 1)

# View the first few rows of the dataset
head(glioblastoma_data)

# Remove non-numeric columns and scale the data
data_for_clustering <- glioblastoma_data[, sapply(glioblastoma_data, is.numeric)]
scaled_data <- t(scale(t(data_for_clustering)))  # Mean-center and scale by rows (genes)
# Heatmap using a diverging color palette (blue-white-red)
heatmap.2(scaled_data,
          Rowv = TRUE,
          Colv = TRUE,
          dendrogram = "both",
          #col = bluered(100),  # Diverging color palette
          col = colorRampPalette(c("blue", "white", "red"))(100),  # Blue to white to red
          trace = "none",
          margins = c(10, 10),
          key = TRUE,
          keysize = 1.5,
          density.info = "none",
          cexCol = 0.7,
          main = "Heatmap with Diverging Palette")


# Heatmap using a sequential color palette (light blue to dark blue)
x11()
heatmap.2(scaled_data,
          Rowv = TRUE,
          Colv = TRUE,
          dendrogram = "both",
          col = colorRampPalette(c("lightblue", "blue", "darkblue"))(100),  # Sequential palette
          trace = "none",
          margins = c(10, 10),
          key = TRUE,
          keysize = 1.5,
          density.info = "none",
          cexCol = 0.7,
          main = "Heatmap with Sequential Palette")

# Cluster only rows (genes)
x11()
heatmap.2(scaled_data,
          Rowv = TRUE,         # Cluster rows (genes)
          Colv = FALSE,        # Do not cluster columns (samples)
          dendrogram = "row",  # Show dendrogram for rows only
          col = bluered(100),  # Diverging color palette
          trace = "none",
          margins = c(10, 10),
          key = TRUE,
          cexCol = 0.7,
          density.info = "none",
          main = "Clustered Genes (Rows) Only")


# Cluster only columns (samples)
x11()
heatmap.2(scaled_data,
          Rowv = FALSE,        # Do not cluster rows (genes)
          Colv = TRUE,         # Cluster columns (samples)
          dendrogram = "column", # Show dendrogram for columns only
          col = bluered(100),  # Diverging color palette
          trace = "none",
          margins = c(10, 10),
          key = TRUE,
          cexCol = 0.7,
          density.info = "none",
          main = "Clustered Samples (Columns) Only")

# Cluster both rows (genes) and columns (samples)
x11()
heatmap.2(scaled_data,
          Rowv = TRUE,         # Cluster rows (genes)
          Colv = TRUE,         # Cluster columns (samples)
          dendrogram = "both", # Show dendrogram for both rows and columns
          col = bluered(100),  # Diverging color palette
          trace = "none",
          margins = c(10, 10),
          key = TRUE,
          cexCol = 0.7,
          density.info = "none",
          main = "Clustered Genes and Samples Together")




# For step 7-8(Functional enrichment analysis and visualisation)
# Data from:https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv" was queried on ShinyGo and top 5 biological pathways recorded
#New project was created and enrichment data imported into environment, then the following code was followed


# Create a data frame for top 5 pathways based on provided data
data <- data.frame(
  pathway = c("Animal organ development", "Cell differentiation", "Cellular developmental process", 
              "Regulation of developmental process", "Response to organic substance"),
  gene_count = c(196, 209, 211, 145, 173), # Number of genes associated with each pathway
  p_value = c(1.96e-21, 1.62e-17, 2.06e-17, 2.06e-17, 2.64e-15) # FDR p-values for each pathway
)

# View the data
print(data)

# Load necessary libraries
library(ggplot2)

# Create a new column for the negative log10 of the p-value
data$log_p_value <- -log10(data$p_value)

# Create the lollipop plot
ggplot(data, aes(x=reorder(pathway, gene_count), y=gene_count)) +
  geom_segment(aes(x=pathway, xend=pathway, y=0, yend=gene_count), color="skyblue") + # Lollipop stem
  geom_point(aes(size=log_p_value), color="blue", alpha=0.8) + # Lollipop point scaled by significance
  coord_flip() + # Flip the coordinates for better readability
  labs(title="Top 5 Pathways by Gene Count and Significance",
       x="Pathway", 
       y="Number of Genes",
       size="-log10(p-value)") + # Label for the point size
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5))
