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

# Normalize the data by counts per million (CPM)
# sum across each column (sample) and divide by 1 million
cpm_data <- apply(glioblastoma_data, 2, function(x) (x / sum(x)) * 1e6)

# View the CPM-normalized data
head(cpm_data)

# Log2-transform the normalized data with a pseudocount of 1 (to avoid issues with zero)
log2_transformed_data <- log2(cpm_data + 1)

# View the log2-transformed data
head(log2_transformed_data)

#Let us check the names of each sample
colnames(log2_transformed_data) #Note that we have 10 Samples

#Each of the samples is labelled 1-10
log2_transformed_data[, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)]

#To calculate fold changes, Data needs to have two conditions.
#split the data into two groups (Treatment and Control)
treatment_samples <- log2_transformed_data[, c(1, 2, 3, 4, 5)]  # Treatment group
control_samples <- log2_transformed_data[, c(6, 7, 8, 9, 10)]    # Control group

# Calculate the mean expression for control and treatment groups
mean_control <- rowMeans(control_samples)
mean_treatment <- rowMeans(treatment_samples)

# Calculate the fold change (log2 fold change: treatment - control)
fold_change <- mean_treatment - mean_control

# View the fold changes
head(fold_change)

# Perform t-tests for each gene between control and treatment
p_values <- apply(log2_transformed_data, 1, function(x) {
  t.test(x[1:2], x[3:4])$p.value  # Adjust the indices based on your actual control/treatment samples
})

# View the p-values
head(p_values)

# Combine the fold change and p-values into a single data frame to identify differentially expressed gene
results <- data.frame(Gene = rownames(log2_transformed_data), 
                      FoldChange = fold_change, 
                      PValue = p_values)

# View the results
head(results)
results


#Filter for upregulated and downregulated genes based on fold change and p-value cutoffs:

#Upregulated gene

# Set fold change and p-value cutoffs
fold_change_cutoff <- 1.5  
p_value_cutoff <- 0.05

# Filter for upregulated genes
upregulated_genes <- results %>%
  filter(FoldChange > fold_change_cutoff & PValue < p_value_cutoff)

# View the upregulated genes
head(upregulated_genes)

#Downregulated gene
# Filter for downregulated genes
downregulated_genes <- results %>%
  filter(FoldChange < -fold_change_cutoff & PValue < p_value_cutoff)

# View the downregulated genes
head(downregulated_genes)

# Save the results to CSV files
write.csv(as.data.frame(upregulated_genes), "upregulated_genes.csv")
write.csv(as.data.frame(downregulated_genes), "downregulated_genes.csv")





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
