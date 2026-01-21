library(ggtext)
library(ComplexHeatmap) #<-Bioconductor 
library(circlize)

#In This script wi will try to extract on txt the most interesting features from the PCA analysis
#PLease use all the elements created in the environment from 1_Aorta_firstpart.R

#You can extract the PCA 1 and 2 values for each sample, and see how the PCA creates unsupervized clusters using K=2
#with this line:
Groups_table <- merged_data[, c(1, 2, 3, 65, 67)]

############################################################################################
# Perform hierarchical clustering based on all features
dist_matrix <- dist(pca_scores[, c("PC1", "PC2")])  # You can use the PCs 1 and 2, or other if interesting
hc <- hclust(dist_matrix, method = "ward.D2")  # Hierarchical clustering using Ward's method

# Plot the dendrogram based on PCA
plot(hc)

#with this you can see the top proteins (in the part 1 of the script were 200, but it can 
#be modified) visualized on the biplot and choosen as most important features based on
#the PCA
View(top_proteins_K)

###################################################################################### 
# GO-BP ANALYSIS of the top 200 proteins based on PCA1 and 2 most important features
# and the Derregulated Proteins fom volcano Plot
######################################################################################

# Convert protein accessions to gene symbols using org.Hs.eg.db (200 PCA most important features)
protein_to_gene <- AnnotationDbi::select(org.Hs.eg.db,
                                         keys = top_proteins_K,
                                         columns = c("SYMBOL", "UNIPROT"),
                                         keytype = "UNIPROT")

# View the table with Protein Accession and corresponding Gene Symbols
print(protein_to_gene)

###########################################################################
# Join both results: Proteins based on PCA and derregulated proteins from 
#ProteomeDiscoverer volcano plot
###########################################################################
# Check if protein accessions from protein_to_gene are found in column 5 of x2025_03_14_codice_2024_08_Vascolari_2_FINALE_Sp_upanddownproteins_fc2_Proteins

# Extract the protein accessions from column 5 of the x2025_03_14_codice_2024_08_Vascolari_2_FINALE_Sp_upanddownproteins_fc2_Proteins
regulated_proteins <- x2025_03_14_codice_2024_08_Vascolari_2_FINALE_Sp_upanddownproteins_fc2_Proteins[, 5]

# Step 2: Create a new column in protein_to_gene indicating if PCA 200 most important features are also in protein is in regulated_proteins
protein_to_gene$Regulated <- protein_to_gene$UNIPROT %in% regulated_proteins

# Step 3: Add the abundance ratio (column 32) and adj. p-value (column 33) to protein_to_gene
# Create a data frame with protein accessions, abundance ratio (Fold Change), and adj.p-value
abundance_and_pvalue <- x2025_03_14_codice_2024_08_Vascolari_2_FINALE_Sp_upanddownproteins_fc2_Proteins[, c(5, 32, 33)]
colnames(abundance_and_pvalue) <- c("UNIPROT", "Abundance_Ratio_Aneurism_vs_healthy", "Adj_p_value")

# Merge the abundance ratio and adj.p-value with protein_to_gene based on the protein accession (UNIPROT)
protein_to_gene <- merge(protein_to_gene, abundance_and_pvalue, by = "UNIPROT", all.x = TRUE)

# View the updated protein_to_gene table with the top 200 PCA features. If a feature
#is present also in de derregulated proteins list from volcano plot
#you will see "TRUE" value on the column named "Regulated", otherwise it will be "FALSE"
View(protein_to_gene)

#Take only Proteins that are on the Derregulated List and the one on the PCA most important features
regulated_proteins <- protein_to_gene %>%
  filter(Regulated == TRUE) %>%
  pull(UNIPROT)

regulated_data <- scaled_proteins_clean_filtered[ , colnames(scaled_proteins_clean_filtered) %in% regulated_proteins]

common_samples <- intersect(rownames(regulated_data), merged_data$SampleName)
print(common_samples)  # This will show the common sample names, just to check there are no errors

# Subset merged_data based on common sample names and order it according to regulated_data's row names
merged_data_ordered <- merged_data[merged_data$SampleName %in% common_samples, ]
merged_data_ordered <- merged_data_ordered[match(rownames(regulated_data), merged_data_ordered$SampleName), ]

# Step 1: Reshape the data from wide to long format
regulated_data_long <- regulated_data %>%
  as.data.frame() %>%
  rownames_to_column("SampleName") %>%
  tidyr::pivot_longer(cols = -SampleName, 
                      names_to = "Protein", 
                      values_to = "Abundance")
# Step 2: Merge with the condition information
regulated_data_long <- regulated_data_long %>%
  left_join(merged_data_ordered[, c("SampleName", "condition")], by = "SampleName")

# Step 3: Create the heatmap plot
ggplot(regulated_data_long, aes(x = Protein, y = SampleName, fill = Abundance)) +
  geom_tile() +  # Create heatmap tiles
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Color scale for abundance
  theme_minimal() +  # Clean theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),  # Rotate protein names for better visibility
        axis.text.y = element_text(size = 6),  # Reduce sample label size
        panel.grid = element_blank()) +  # Remove gridlines
  labs(title = "Heatmap Proteins present in Differential Expressed and top 25% Biplot feat",
       x = "Protein",
       y = "Sample",
       fill = "Abundance") +
  facet_wrap(~condition, scales = "free_y")  # Facet by condition for better separation between sample groups


##################################
##Total Treemap of proteins present in both 25% proteins from PCA and volcanoplot derregulated. They cluster better 
#our samples
#################################

# Step 2.1: Calculate distance matrices for both samples and proteins
#Info about the measures choosen:
#Manhattan distance (dist(..., method = "manhattan")): The sum of absolute differences (also known as L1 norm).
#Canberra distance (dist(..., method = "canberra")): A weighted version of the Manhattan distance.
#"ward.D": Minimizes the total within-cluster variance.


# For proteins (distance between proteins)
dist_samples<- dist(regulated_data, method = "canberra") #euclidean, manhattan, maximum, canberra  # Default: Euclidean distance between proteins
hc_samples <- hclust(dist_samples, method = "ward.D2") #ward.D, ward.D2, single, complete, average, centroid, median, mcquitty

# Plot dendrogram for Complete linkage with Euclidean distance
plot(hc_samples, main = "Ward.D2 with canberra Distance", xlab = "", sub = "", cex = 0.9)

protein_to_gene_map <- setNames(protein_to_gene$SYMBOL, protein_to_gene$UNIPROT)

##########################################################
colnames(regulated_data) <- protein_to_gene_map[colnames(regulated_data)]
# Check the changes
head(colnames(regulated_data))
# Step 5.1: Create the heatmap with dendrograms
col_fun = circlize::colorRamp2(c(0, 1000), c("white", "red"))
heatmap <- ComplexHeatmap::Heatmap(regulated_data,
                                   col = col_fun,
                                   name = "Hierarchical Clustering: Ward.D2; Canberra Distance",  # Title for the heatmap
                                   cluster_rows = hc_samples,  # Clustering for samples
                                   cluster_columns = hc_proteins,  # Clustering for proteins
                                   show_row_dend = TRUE,  # Show dendrogram for rows (samples)
                                   show_column_dend = TRUE,  # Show dendrogram for columns (proteins)
                                   row_title = "Samples",  # Label for rows
                                   column_title = "Protein Abundance"  # Label for columns
)
# Draw the heatmap
ComplexHeatmap::draw(heatmap)
#########################################################

# Step 2: Replace UNIPROT accessions with SYMBOLs in the column names of regulated_data
colnames(regulated_data) <- protein_to_gene_map[colnames(regulated_data)]

# Step 3: Handle any missing SYMBOLs (if some UNIPROT accessions do not have a SYMBOL)
# Optionally, you can replace NAs with a default label like "Unknown"
colnames(regulated_data)[is.na(colnames(regulated_data))] <- "Unknown"

# Step 4: Create the heatmap with the SYMBOLs on the columns
heatmap <- ComplexHeatmap::Heatmap(regulated_data,
                                   name = "Protein Abundance, Hierarchical Clustering Ward.D2 with Canberra Distance",  # Title for the heatmap
                                   cluster_rows = hc_samples,  # Clustering for samples (assumed already done)
                                   cluster_columns = hc_proteins,  # Clustering for proteins (now SYMBOLS)
                                   show_row_dend = TRUE,  # Show dendrogram for rows (samples)
                                   show_column_dend = TRUE,  # Show dendrogram for columns (proteins)
                                   row_title = "Samples",  # Label for rows
                                   column_title = "Proteins",  # Label for columns (proteins, now SYMBOLS)
                                   column_names_gp = gpar(fontsize = 8)  # Adjust font size for column names (SYMBOLS)
)

# Step 5: Plot the heatmap
ComplexHeatmap::draw(heatmap)

write.csv(scaled_proteins, file = "scaled_proteins_accessionNumber.txt", row.names = TRUE)
write.csv(merged_data, file = "PCAScore_cleanproteins_byfile.txt", row.names = TRUE)
write.csv(protein_to_gene, file = "Top25percentProteinsfromPCA_regulatedonVolcanoPLot.txt", row.names = TRUE)

abundance_and_pvalue_2 <- x2025_03_14_codice_2024_08_Vascolari_2_FINALE_Sp_upanddownproteins_fc2_Proteins[, c(5, 26, 32, 33)]
colnames(abundance_and_pvalue_2) <- c("UNIPROT", "GeneSymbol", "Abundance_Ratio_Aneurism_vs_healthy", "Adj_p_value")
write.csv(abundance_and_pvalue_2, file = "Derregulatedproteinsfromvolcanoplot_logFC2_pvalue005.txt", row.names = FALSE)

# Filter the data to keep only rows where Adj_p_value >= 0.05
abundance_and_pvalue_3 <- abundance_and_pvalue_2 %>%
  filter(Adj_p_value <= 0.05)

write.csv(abundance_and_pvalue_3, file = "Derregulatedproteinsfromvolcanoplot_logFC2_adjpvalue005.txt", row.names = FALSE)
