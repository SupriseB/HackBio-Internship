if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
BiocManager::install("sesame")
#edgeR and lima for Differential Analysis
BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("EDASeq", force = TRUE)
install.packages('gplots')
BiocManager::install("sesameData", force = TRUE)
installed.packages()["sesameData", ]
install.packages("httr2")
BiocManager::install("SummarizedExperiment", force = TRUE)
install.packages("pryr")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")  # Human gene annotation database
BiocManager::install("enrichplot")
install.packages("pheatmap")
BiocManager::install("impute")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
BiocManager::install("pathview")

library(BiocManager)
library(TCGAbiolinks)
library(edgeR)
library(limma)
library(EDASeq)
library(gplots)
library(sesameData) 
library(SummarizedExperiment)
library(dplyr)
library(pryr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(pheatmap)
library(sesame)
library(impute)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(pathview)

#########################     Building Query      ##########################

query <- GDCquery(project = "TCGA-LGG")

getProjectSummary("TCGA-LGG")
?GDCquery



#Download and process data
LGGQ <- GDCquery(
  project = "TCGA-LGG",  
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value"
)

View(LGGQ)



# Check memory usage
mem_used()  # Current memory usage
gc()


GDCdownload(LGGQ)
LGG.data <- GDCprepare(LGGQ)
str(LGG.data)  
head(LGG.data)
colnames(LGG.data)
row.names(LGG.data)


#Explore metadata information
LGG.data$tumor_descriptor
LGG.data$barcode
LGG.data$tissue_type
LGG.data$paper_IDH.specific.DNA.Methylation.Cluster
LGG.data$paper_Random.Forest.Sturm.Cluster
LGG.data$paper_Supervised.DNA.Methylation.Cluster
LGG.data$paper_IDH.status
LGG.data$paper_MGMT.promoter.status
LGG.data$paper_Pan.Glioma.DNA.Methylation.Cluster  #distinct groups of gliomas characterized by unique patterns of DNA methylation

# Create metadata for our use
sMetaLGG <- data.frame(
  barcode = LGG.data$barcode,
  IDH.specific.DNA.Methylation.Cluster = LGG.data$paper_IDH.specific.DNA.Methylation.Cluster,
  IDH_status = LGG.data$paper_IDH.status,
  stringsAsFactors = FALSE  # Prevents automatic conversion of strings to factors
)

head(sMetaLGG)


#######################           Preprocessing     ##########################


#Select unstranded dataset !!!
LGG.raw.data <- assay(LGG.data)
dim(LGG.raw.data)

# Step 1: Check available IDH status
unique(LGG.data$paper_IDH.status)

# Step 2: Separate WT and mutant samples
# Adjust based on the actual labels in your data
WT_samples <- LGG.data$barcode[LGG.data$paper_IDH.status == "WT"]
mutant_samples <- LGG.data$barcode[LGG.data$paper_IDH.status == "Mutant"]

str(WT_samples)
str(mutant_samples)

# Step 4: Subset the raw count data for selected samples
LGG_samples <- c(WT_samples, mutant_samples)
### LGG.raw.data.selected <- LGG.raw.data[, selected_samples]


# Check the structure and dimensions of the raw data matrix
str(LGG_samples)
dim(LGG_samples)  # Should have genes as rows and samples as columns
head(LGG_samples)


###############         NORMALIZATION & FILTRATION    ##################

# Normalization using limma
# Perform between-array normalization
# Quantile normalization ensures that the distribution of methylation values across samples is the same,
# which is essential when comparing methylation levels across multiple samples.
LGG.norm.data <- normalizeBetweenArrays(LGG.raw.data, method = "quantile")

# Check the result
boxplot(LGG.norm.data, main="Post-normalization Boxplot", las=2)

dim(LGG.norm.data)

# Filter out probes with low variance
varianceThreshold <- 0.01  # You can adjust this threshold based on your dataset
LGG.filtered.data <- LGG.norm.data[apply(LGG.norm.data, 1, var) > varianceThreshold, ]

# Filter probes with many NA values (if any)
naThreshold <- 0.2  # Filter out probes with more than 20% missing values
LGG.filtered.data <- LGG.filtered.data[rowMeans(is.na(LGG.filtered.data)) < naThreshold, ]

# Optional: Impute missing values if required
LGG.filtered.data <- impute.knn(LGG.filtered.data)$data

dim(LGG.filtered.data)
head(LGG.filtered.data)


sMetaLGG_clean <- na.omit(sMetaLGG)

write.csv(sMetaLGG_clean, file = "/home/berkay/HackBio/stg4/smetaLGG_clean.csv", row.names = FALSE)
write.csv(results, file = "/home/berkay/HackBio/stg4/DMA_results.csv", row.names = TRUE)
write.csv(signif_results, file = "/home/berkay/HackBio/stg4/signif_results.csv", row.names = TRUE)

# Convert sMetaLGG# Convert your filtered data to a CSV file
# write.csv(sMetaLGG, file = "LGG_filtered_data.csv", row.names = TRUE)
# getwd()

setwd("home/berkay/HackBio/stg4")
file.exists("/home/berkay/HackBio/stg4")

# Compress the CSV file into a ZIP archive
zip(zipfile = "LGG_filtered_data.zip", files = "LGG_filtered_data.csv")



#####################         DMA               ######################
# Extract the sample barcodes from LGG.filtered.data
sample_barcodes <- colnames(LGG.filtered.data)

# Create the group variable by matching the barcodes
group <- ifelse(LGG.data$paper_IDH.status[match(sample_barcodes, LGG.data$barcode)] == "WT", "WT", "Mutant")

# Ensure the group has the correct length (should be 534)
length(group)  # Should be 534
design <- model.matrix(~0 + group)  # No intercept, since we are comparing WT vs. Mutant
colnames(design) <- c("WT", "Mutant")



# Check the available IDH status in LGG.data
unique(LGG.data$paper_IDH.status)

# Check the dimensions of the design matrix
dim(design)  # Should be (534, 2)

# Find samples with missing IDH status
missing_IDH_samples <- sample_barcodes[is.na(LGG.data$paper_IDH.status[match(sample_barcodes, LGG.data$barcode)])]

# Check how many samples are missing IDH status
length(missing_IDH_samples)

# Get the indices of samples that have non-missing IDH status
valid_samples <- sample_barcodes[!is.na(LGG.data$paper_IDH.status[match(sample_barcodes, LGG.data$barcode)])]

# Subset LGG.filtered.data to keep only the valid samples
LGG.filtered.data <- LGG.filtered.data[, valid_samples]

# Check the new dimensions of LGG.filtered.data
dim(LGG.filtered.data)  # Should be (probes, 513)


#############             Fit linear model            ###############


# Step 2: Fit linear model using limma
fit <- lmFit(LGG.filtered.data, design)

# Step 3: Define the contrast (Mutant vs WT)
contrast.matrix <- makeContrasts(Mutant_vs_WT = Mutant - WT, levels = design)

# Apply the contrast matrix to the fit
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

####################      Apply thresholds to identify significant methylation changes ######

# Step 4: Extract results
results <- topTable(fit2, adjust.method = "fdr", number = Inf, sort.by = "P")
head(results)

# Step 5: Apply thresholds to identify significant methylation changes
# Typically, you apply thresholds for both p-value and logFC
signif_results <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5, ]

# View significant results
head(signif_results)
dim(signif_results)

rm(LGG.norm.data)
rm(LGGQ)
rm(LGG.data)
rm(LGG.raw.data)
gc()


############      Volcano plot    #######

volcano_data <- data.frame(
  logFC = results$logFC,
  negLogPval = -log10(results$P.Value)
)

ggplot(volcano_data, aes(x = logFC, y = negLogPval)) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Methylation",
       x = "log Fold Change", y = "-log10(p-value)") +
  geom_vline(xintercept = c(-0.5, 0.5), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = "dashed")


gc()


# Assuming 'results' is your data frame from the differential analysis
# Check if logFC and P.Value columns exist
head(results)

# Define cutoffs
logFC_cutoff <- 0.4  # Adjust if needed
pvalue_cutoff <- 0.05  # Adjust if needed

# Create a new column for significance based on logFC and p-value
results$significant <- ifelse(results$P.Value < pvalue_cutoff & abs(results$logFC) > logFC_cutoff, "Significant", "Not significant")

# Check the new 'significant' column
table(results$significant)

# Create the volcano plot with customized colors
library(ggplot2)
ggplot(results, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = significant), size = 1.5) +  # Color by significance
  scale_color_manual(values = c("Significant" = "black", "Not significant" = "gray70")) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "red") +  # Vertical lines for logFC threshold
  geom_hline(yintercept = -log10(pvalue_cutoff), linetype = "dashed", color = "blue") +  # Horizontal line for p-value threshold
  labs(title = "Volcano Plot of Differential Methylation", x = "log Fold Change", y = "-log10(p-value)") +
  theme_minimal()







##################      HEATMAP  ##################



# Step 8: Heatmap of top differentially methylated probes
colnames(top_data)  # List column names in top_data
group <- factor(group)


# !!! Reassign group based on LGG data's IDH status
group <- LGG.data$paper_IDH.status
length(group)  # Should match ncol(top_data)
ncol(top_data)

pheatmap(top_data, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = FALSE,
         show_colnames = FALSE)


# Select top 50 probes
top_50_probes <- top_data[1:50, ]

# Check dimensions to ensure only 50 probes are selected
dim(top_50_probes)

# Check the length of the group factor
length(group)

annotation_df <- data.frame(group = group)
row.names(annotation_df) <- colnames(top_50_probes)

# Check if column names of top_50_probes match row names of annotation_df
all(colnames(top_50_probes) == row.names(annotation_df))



print(group)


# Align group with the columns of top_50_probes
group <- group[match(colnames(top_50_probes), names(group))]

# Recreate annotation_df
annotation_df <- data.frame(group = group)
row.names(annotation_df) <- colnames(top_50_probes)



pheatmap(top_50_probes, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE,  # Display probe names
         show_colnames = FALSE, 
         annotation_col = annotation_df)  # Use the correct annotation


rm(LGG.data)

# Check for NAs
sum(is.na(group))  # Should return 0 if no missing values

# Ensure group length matches the number of columns in top_50_probes
group <- group[1:ncol(top_50_probes)]

# Check the unique values in group
unique(group)

# Define colors for the groups
annotation_colors <- list(group = c("WT" = "blue", "Mutant" = "red"))

# Re-run pheatmap with annotation colors
pheatmap(top_50_probes, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         annotation_col = annotation_df, 
         annotation_colors = annotation_colors)
    


######################      Functional Enrichment Analysis      #################

probe_annotations <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Select probes from your signif_results (replace "signif_results" with your data)
top_probes <- rownames(signif_results)
print(top_probes)
mapped_genes <- probe_annotations[probe_annotations$Name %in% top_probes, "UCSC_RefGene_Name"]
print(mapped_genes)

# Extract unique genes from the mapped probes
genes <- unique(unlist(strsplit(as.character(mapped_genes), ";")))
print(genes)

library(clusterProfiler)
library(org.Hs.eg.db)  # Human gene annotation database

# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform GO enrichment analysis
go_results <- enrichGO(gene_entrez$ENTREZID, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID", 
                       ont = "BP",  # Biological Process
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05)

# View and plot results
head(go_results)
barplot(go_results)
cnetplot(go_results)
dotplot(go_results)



# Perform KEGG enrichment analysis
kegg_results <- enrichKEGG(gene_entrez$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
head(kegg_results)
head(gene_entrez)
gene_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
length(gene_entrez$ENTREZID)

kegg_results <- enrichKEGG(gene_entrez$ENTREZID, organism = 'hsa', pvalueCutoff = 0.1)  # Use 0.1 instead of 0.05 


sum(is.na(gene_entrez$ENTREZID))  # Check how many IDs failed to map

# Convert KEGG results to a dataframe
kegg_df <- as.data.frame(kegg_results)


# Extracting results manually if as.data.frame fails
kegg_df <- kegg_results@result

head(kegg_df)
summary(kegg_df)




# View and plot results
head(kegg_df)
dotplot(kegg_results)
cnetplot(kegg_results)
barplot(kegg_results)









