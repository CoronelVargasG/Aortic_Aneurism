library(org.Hs.eg.db)

#######GOBP - of proteins derregulated as well as important in PCA (25% of total)################
#Step0 take all the identified genes and their adj p-value
all_genes_identified <- x2025_03_14_codice_2024_08_Vascolari_2_FINALE_Sp_proteinsforR_Proteins[, c(5, 7)]
# Remove rows with any NA values
all_genes_identified <- na.omit(all_genes_identified)
# Remove rows where the first column contains the character ";"
all_genes_identified <- all_genes_identified[!grepl(";", all_genes_identified[[1]]), ]
all_symbols_identified <- all_genes_identified$Gene.Symbol
all_entrez_ids <- clusterProfiler::bitr(all_symbols_identified, fromType = "SYMBOL", toType = "ENTREZID", org.Hs.eg.db)
# Merge with the original table to add Entrez IDs
GOBP_all_entrez_ids <- merge(all_entrez_ids, all_genes_identified, by.x = "SYMBOL", by.y = "Gene.Symbol")
all_gene_list <- GOBP_all_entrez_ids$Abundance.Ratio.Adj..P.Value...Aortic_Aneurism.....Healthy.
names(all_gene_list) <- GOBP_all_entrez_ids$ENTREZID
View(all_gene_list)

#################################################################################################################
#######Genes of interest: Proteins present in Volcano Plot as Derregulated with Adj. p-val <or= to 0.05
#######as well as present in the 25& most significative features of the PCA
# Read the tables (adjust sep if needed, e.g. sep = "\t")
top25 <- read.table(
  "Top25percentProteinsfromPCA_regulatedonVolcanoPLot.txt",
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE
)
top25_regulated <- subset(top25, Regulated == TRUE, row.names=FALSE) #keep only proteins present in Derregulatedproteinsfromvolcanoplot_logFC2_adjpvalue005.txt
top25_regulated <- top25_regulated[,-1]
row.names(top25_regulated) <- NULL
# Step 1: Filter the data based on adjusted p-value < 0.05
significant_data <- top25_regulated[
  top25_regulated$Adj_p_value <= 0.05, 
  c("SYMBOL", "Adj_p_value")
]

# Step 2: Remove duplicates (if any) and ensure unique gene symbols
significant_data <- unique(significant_data)

# Step 3: Sort genes by adjusted p-value (lower p-value means higher significance)
significant_data <- significant_data[order(significant_data$Adj_p_value), ]
row.names(significant_data)<-NULL
# Step 4: Convert GeneSymbols to Entrez IDs using org.Hs.eg.db
gene_entrez_ids <- clusterProfiler::bitr(significant_data$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Check if conversion was successful
head(gene_entrez_ids)

# Step 5: Create a ranked gene list (using p-value as the ranking criterion)
# Create a named vector where names are Entrez IDs and values are the corresponding adjusted p-values
gene_list <- significant_data$Adj_p_value
names(gene_list) <- gene_entrez_ids$ENTREZID

# Step 6: Perform GO enrichment analysis (Biological Process by default), without considering p-value
#but keeping as universe all identified proteins
go_result <- clusterProfiler::enrichGO(
  gene = names(gene_list),  # Provide the Entrez IDs
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  universe = names(all_gene_list), #background genes. If missing, the all genes listed in the database (eg TERM2GENE table) will be used as background
  minGSSize	= 1, #minimal size of genes annotated by Ontology term for testing.
  ont = "BP",  # "BP" for Biological Process; you can also use "CC" for Cellular Component, "MF" for Molecular Function
  pAdjustMethod = "none",  # one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
  qvalueCutoff = 1,  # Adjust the q-value threshold as needed
  readable = TRUE  # Convert Entrez IDs back to gene symbols
)

#Step 7: Perform GSEA for all proteins detected If you want the adjusted p-values to influence the analysis
# 1. Remove NAs
geneList2 <- gene_list[!is.na(names(gene_list))]

# 2. Convert to ranking metric
geneList2 <- -log10(geneList2)

# 3. Remove infinite values (adj p = 0)
geneList2 <- geneList2[is.finite(geneList2)]

# 4. Remove duplicated Entrez IDs (KEEP the most significant)
geneList2 <- geneList2[!duplicated(names(geneList2))]

# 5. Sort decreasing (MANDATORY)
geneList2 <- sort(geneList2, decreasing = TRUE)

#run analysis
gsea_result <- clusterProfiler::gseGO(
  geneList = geneList2,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 2,
  pAdjustMethod = "none",
  verbose = TRUE,  # Adjust the q-value threshold as needed
  readable = TRUE
)

gsea_df <- gsea_result@result
gsea_df$core_enrichment_ENTREZ <- gsea_df$core_enrichment

gsea_df$core_enrichment_SYMBOL <- clusterProfiler::setReadable(
  gsea_result,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)@result$core_enrichment



# Step 7: View the GO analysis results
head(go_result)
head(gsea_result)
head(gsea_df)

# If you want to visualize the GO enrichment results, you can plot a barplot
barplot(go_result, showCategory = 10)

# Alternatively, a dotplot can be used
clusterProfiler::dotplot(go_result, showCategory = 10)
clusterProfiler::dotplot(gsea_result, showCategory = 10)

# Step 1: Convert the GO result to a data frame
go_result_df <- as.data.frame(go_result)
go_result_df_forrevigo <- go_result_df[, c(1, 8)]
# Step 2: Export the data frame to a CSV file
write.csv(go_result_df, "GO_enrichment_results.csv", row.names = FALSE)
write.table(go_result_df_forrevigo, "GO_enrichment_results_forrevigo.txt", row.names = FALSE, sep = "\t")

# Alternatively, you can export it as a tab-delimited text file
write.table(go_result_df, "GO_enrichment_results.txt", sep = "\t", row.names = FALSE)

#write the table of GSEA analysis based on adj p-value of the Top25percentProteinsfromPCA_regulatedonVolcanoPLot
gsea_export_txt <- as.data.frame(gsea_df)
write.table(gsea_export_txt, "gsea_export_txt_based_on_adj_pvalue_Top25percentProteinsfromPCA_regulatedonVolcanoPLot.txt", row.names = FALSE, sep = "\t")
