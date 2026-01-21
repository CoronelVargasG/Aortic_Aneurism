library(ggplot2)
library(factoextra)    
library(dplyr)
library(tibble)
library(scales)
library(ggrepel)
library(org.Hs.eg.db)

#Import all the tables on https://github.com/CoronelVargasG/Aortic_Aneurism and be sure 
#to import correctly the columns names (heading)

All_proteins <- x2025_03_14_codice_2024_08_Vascolari_2_FINALE_Sp_proteinsforR_Proteins[, c(3, 5, 6, 7, 8:69)]
scaled_proteins <- All_proteins[, c(-2, -3, -4)]
row.names(scaled_proteins) <- scaled_proteins[, 1]
scaled_proteins <- scaled_proteins[, -1]

#the columns names are in the correct order of the samples on the matrix "x2025_03_14_codice_2024_08_Vascolari_2_FINALE_Sp_proteinsforR_Proteins"
#If you want to rehuse this script with other data from proteome discoverer, be sure to name correctly each sample
scaled_proteins_column_names <- c("10364_Aortic_Aneurism 1",
                                  "10364_Aortic_Aneurism 2",
                                  "10367_Aortic_Aneurism 1",
                                  "10367_Aortic_Aneurism 2",
                                  "VD2100020_Aortic_Aneurism 1",
                                  "VD2100020_Aortic_Aneurism 2",
                                  "VD2100023_Aortic_Aneurism 1",
                                  "VD2100023_Aortic_Aneurism 2",
                                  "VD2100042_Aortic_Aneurism 1",
                                  "VD2100042_Aortic_Aneurism 2",
                                  "VD2100089_Aortic_Aneurism 1",
                                  "VD2100089_Aortic_Aneurism 2",
                                  "VD2100127_Aortic_Aneurism 1",
                                  "VD2100127_Aortic_Aneurism 2",
                                  "VD2100143_Aortic_Aneurism 1",
                                  "VD2100143_Aortic_Aneurism 2",
                                  "VD2200046_Aortic_Aneurism 1",
                                  "VD2200046_Aortic_Aneurism 2",
                                  "VD2300002_Aortic_Aneurism 1",
                                  "VD2300002_Aortic_Aneurism 2",
                                  "VD2300004_Aortic_Aneurism 1",
                                  "VD2300004_Aortic_Aneurism 2",
                                  "VD2300023_Aortic_Aneurism 1",
                                  "VD2300023_Aortic_Aneurism 2",
                                  "VD2300050_Aortic_Aneurism 1",
                                  "VD2300050_Aortic_Aneurism 2",
                                  "VD2300066_Aortic_Aneurism 1",
                                  "VD2300066_Aortic_Aneurism 2",
                                  "VD2300083_Aortic_Aneurism 1",
                                  "VD2300083_Aortic_Aneurism 2",
                                  "VD2300087_Aortic_Aneurism 1",
                                  "VD2300087_Aortic_Aneurism 2",
                                  "VD2300117_Aortic_Aneurism 1",
                                  "VD2300117_Aortic_Aneurism 2",
                                  "VD2400050_Aortic_Aneurism 1",
                                  "VD2400050_Aortic_Aneurism 2",
                                  "VD2400066_Aortic_Aneurism 1",
                                  "VD2400066_Aortic_Aneurism 2",
                                  "FA1_Aortic_Aneurism 1",
                                  "FA1_Aortic_Aneurism 2",
                                  "FA3_Aortic_Aneurism 1",
                                  "FA3_Aortic_Aneurism 2",
                                  "FA5_Aortic_Aneurism 1",
                                  "FA5_Aortic_Aneurism 2",
                                  "FA6_Aortic_Aneurism 1",
                                  "FA7_Aortic_Aneurism 1",
                                  "FA7_Aortic_Aneurism 2",
                                  "FA6_Aortic_Aneurism 2",
                                  "SCALA_50-----_Healthy 1",
                                  "SCALA_50-----_Healthy 2",
                                  "245AIT_Healthy 1",
                                  "245AIT_Healthy 2",
                                  "284AIT_Healthy 1",
                                  "284AIT_Healthy 2",
                                  "230AIT_Healthy 1",
                                  "230AIT_Healthy 2",
                                  "161FCP/FA_Healthy 1",
                                  "161FCP/FA_Healthy 2",
                                  "279TCFA/HR_Healthy 1",
                                  "279TCFA/HR_Healthy 2",
                                  "282HR_Healthy 1",
                                  "282HR_Healthy 2"
)

experimental_design <- x2024_08_02_Vascolari_amples[c(-19, -20, -48),c(3, 6, 7, 8) ]
rownames(experimental_design) <- NULL
new_row <- data.frame(File = "F65",
                      SampleName = "FA6",
                      condition = "Aortic_Aneurism",
                      instrumentareplica = "2")
experimental_design <- rbind(experimental_design[1:47, ], new_row, experimental_design[48:nrow(experimental_design), ])
rownames(experimental_design) <- NULL

colnames(scaled_proteins) <- scaled_proteins_column_names
experimental_design$SampleName <- scaled_proteins_column_names
View(scaled_proteins)
View(experimental_design)

##########DATA CLEANING: AT LEAST 3 detections in at least one group#############
# 1. Handle NaN values (imputation or removal)
scaled_proteins[scaled_proteins == ""] <- NA
# 3. Force numeric conversion and matrix format
# ----------------------------------------------------
scaled_proteins_clean <- apply(scaled_proteins, 2, as.numeric)
rownames(scaled_proteins_clean) <- rownames(scaled_proteins)
scaled_proteins_clean[is.na(scaled_proteins_clean)] <- 0  # Impute with zeros or any other method
#Invert the columns and rows
scaled_proteins_clean <- t(scaled_proteins_clean)
# Split row names into two groups: Aortic_Aneurism and Healthy
group_aortic <- grepl("_Aortic_Aneurism", rownames(scaled_proteins_clean))
group_healthy <- grepl("_Healthy", rownames(scaled_proteins_clean))

# Define a function to check if a column meets the criteria AT LEAST 3 detections in at least one group
check_non_zero <- function(col, group1, group2) {
  # Extract values for each group
  group1_values <- col[group1]
  group2_values <- col[group2]
  
  # Check if there are at least 3 non-zero values in either group
  sum(group1_values != 0) >= 3 | sum(group2_values != 0) >= 3
}

# Apply the check to each column
cols_to_keep <- apply(scaled_proteins_clean, 2, check_non_zero, group1 = group_aortic, group2 = group_healthy)

# Subset the matrix by the columns that meet the condition, while keeping the original row order
scaled_proteins_clean_filtered <- scaled_proteins_clean[, cols_to_keep]

# Display the filtered matrix
View(scaled_proteins_clean_filtered)

# Check the dimensions of the original and filtered matrices
dim_original <- dim(scaled_proteins_clean)
dim_filtered <- dim(scaled_proteins_clean_filtered)

###########ANALYSIS#################

# 2. Perform PCA
pca_result <- prcomp(scaled_proteins_clean_filtered, center = TRUE, scale. = FALSE)

# 3. Get PCA scores (principal components)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleName <- rownames(pca_scores)

# 4. Merge PCA scores with experimental design data to get conditions
merged_data <- merge(pca_scores, experimental_design, by = "SampleName")

# 5. Calculate the proportion of variance explained by each PC
explained_variance <- (pca_result$sdev^2) / sum(pca_result$sdev^2)  # Eigenvalue / Total Variance
explained_variance_percent <- round(explained_variance * 100, 1)  # Convert to percentage

# Print explained variance for each PC
cat("Explained variance by PC1:", explained_variance_percent[1], "%\n")
cat("Explained variance by PC2:", explained_variance_percent[2], "%\n")

# Define the level for stat_ellipse
ellipse_level <- 0.90

###########################PLOTS################
# ----------------------------------------------------
# Elbow method using k-means
# ----------------------------------------------------
set.seed(123)

wss <- sapply(1:50, function(k) {
  kmeans(pca_scores[, 1:50], centers = k, nstart = 25)$tot.withinss
})

# ----------------------------------------------------
# Elbow plot
# ----------------------------------------------------
elbow_df <- data.frame(
  k = 1:50,
  wss = wss
)
#kmeans elbowplot
ggplot(elbow_df, aes(x = k, y = wss)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2) +
  labs(
    title = "Elbow Plot (k-means on PCA scores)",
    x = "Number of clusters (k)",
    y = "Total within-cluster sum of squares"
  ) +
  theme_minimal()


# Plot the elbow plot based on PC components
factoextra::fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 100)) + 
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Elbow Plot for PCA 3105 proteins")

# 2) PCA Plot with Labels, Color by Condition, and Ellipses Around Groups

# Get PCA scores and merge with experimental design
pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleName <- rownames(pca_scores)
merged_data <- merge(pca_scores, experimental_design, by = "SampleName")

# Calculate explained variance percentages
explained_variance <- (pca_result$sdev^2) / sum(pca_result$sdev^2)  # Eigenvalue / Total Variance
explained_variance_percent <- round(explained_variance * 100, 1)  # Convert to percentage

# Plot PCA with labels and ellipses around groups
ggplot2::ggplot(merged_data, aes(x = PC1, y = PC2, color = condition, label = SampleName)) + 
  ggplot2::geom_point(size = 3) +
  ggplot2::stat_ellipse(level = ellipse_level) +
  ggrepel::geom_text_repel(size = 3) + 
  ggplot2::labs(title = paste("PCA Plot with Ellipses level", ellipse_level, "and Sample Labels"),
                x = paste("PC1 (", explained_variance_percent[1], "%)", sep = ""),
                y = paste("PC2 (", explained_variance_percent[2], "%)", sep = "")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom")

########################################################################################
#3D PCA
########################################################################################
library(plotly)

pca_plot_df <- as.data.frame(pca_result$x)
pca_plot_df$Group <- ifelse(
  grepl("_Aortic_Aneurism", rownames(pca_plot_df)),
  "Aortic Aneurism",
  "Healthy")
# Percent variance
var_explained <- round(100 * summary(pca_result)$importance[2, ], 1)

plot_ly(
  pca_plot_df,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~Group,
  colors = c("Healthy" = "#1B9E77", "Aortic Aneurism" = "#D95F02"),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 4)
) %>%
  layout(
    title = "3D PCA of Scaled Proteins",
    scene = list(
      xaxis = list(title = paste0("PC1 (", var_explained[1], "%)")),
      yaxis = list(title = paste0("PC2 (", var_explained[2], "%)")),
      zaxis = list(title = paste0("PC3 (", var_explained[3], "%)"))
    )
  )


########################################################################################
# 3) First Biplot Showing the Top 10 Most Important Proteins Based on the Condition
# Define a scaling factor for the vectors (you can adjust this value)
scale_factor <- 50000  # You can increase or decrease this depending on the visibility
########################################################################################

# Identify the most important proteins (around 20% of total proteins) for each condition 
#(based on absolute value of loadings)
loadings <- pca_result$rotation[, 1:2]  # Loadings for PC1 and PC2
importance <- apply(loadings, 1, function(x) sum(abs(x)))  # Sum of absolute loadings
top_proteins <- names(sort(importance, decreasing = TRUE))[1:798] #we have a total of 3191 proteins, so 25% is around 638

# Create biplot showing the top 10 most important proteins
top_loadings <- loadings[top_proteins, ]
top_loadings_df <- as.data.frame(top_loadings)

# 4) Second Biplot Showing the Top 25% Most Important Proteins Based on the Elbow Plot's K Value

# Assume K is chosen based on the elbow plot, for example, K = 3
K <- 2

#############OPTIONAL if you want to look for a specific protein
feature_to_circle <- "P17677"
p17677_coords <- data.frame(PC1 = pca_result$rotation[feature_to_circle, "PC1"] * scale_factor, PC2 = pca_result$rotation[feature_to_circle, "PC2"] * scale_factor)
#############

loadings_K <- pca_result$rotation[, 1:K]  # Loadings for first K components
importance_K <- apply(loadings_K, 1, function(x) sum(abs(x)))  # Sum of absolute loadings
top_proteins_K <- names(sort(importance_K, decreasing = TRUE))[1:798]
# Add the cluster information to the merged_data
kmeans_result <- kmeans(pca_scores[, c("PC1", "PC2")], centers = K)
merged_data$cluster <- as.factor(kmeans_result$cluster)

# Create biplot showing the top 500 most important proteins based on K components
top_loadings_K <- loadings_K[top_proteins_K, ]
top_loadings_K_df <- as.data.frame(top_loadings_K)
num_top_proteins <- nrow(top_loadings_K_df)

# Apply the scaling factor to the loadings to make the protein vectors more visible
scaled_top_loadings_K_df <- top_loadings_K_df * scale_factor


# Create the second biplot with the adjusted protein vectors
ggplot2::ggplot(merged_data, aes(x = PC1, y = PC2, color = condition)) + 
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(aes(label = SampleName), size = 5, hjust = 0.5, vjust = 0.5) + 
  ggplot2::geom_segment(data = as.data.frame(scaled_top_loadings_K_df), aes(x = 0, y = 0, xend = PC1, yend = PC2),
                        arrow = grid::arrow(type = "open", length = grid::unit(0.2, "inches")), color = "black") +
  ggplot2::geom_text(data = as.data.frame(scaled_top_loadings_K_df), aes(x = PC1 * 1.1, y = PC2 * 1.1, label = rownames(scaled_top_loadings_K_df)),
                     size = 3, color = "black") +
  ggplot2::stat_ellipse(level = ellipse_level) +
  ggplot2::labs(title = paste("PCA Biplot with Top", num_top_proteins, "Important Proteins (K =", K, ", Scale Factor =", scale_factor, ")"),
                x = paste("PC1 (", explained_variance_percent[1], "%)", sep = ""),
                y = paste("PC2 (", explained_variance_percent[2], "%)", sep = "")) +
  
  # 🔴 Red circle around P17677
  geom_point(
    data = p17677_coords,
    aes(x = PC1, y = PC2),
    shape = 21,
    size = 7,
    stroke = 2,
    color = "red",
    fill = NA,
    inherit.aes = FALSE
  ) +
  
  # Label P17677
  geom_text(
    data = p17677_coords,
    aes(x = PC1 * 1.15, y = PC2 * 1.15, label = feature_to_circle),
    color = "red",
    fontface = "bold",
    size = 4,
    inherit.aes = FALSE
  ) +
  ggplot2::theme_minimal()

