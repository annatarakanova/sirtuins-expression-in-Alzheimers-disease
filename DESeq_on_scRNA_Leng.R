# =================== Check the content of SingleCellExperiment object ====

library(SingleCellExperiment)
sce <- readRDS("~/short_term_project/EntorhinalCortex_Leng_sce_ann.rds")

sce

assayNames(sce)

colnames(colData(sce))
df <- head(as.data.frame(colData(sce)))

colnames(rowData(sce))

table(sce$celltype)

# There is no healthy samples, all of them are Alzheimer's disease on 3 stages (Braak stages - 0, 2, 6)
unique(colData(sce)$BraakStage)
table(sce$BraakStage)
# Condition: 0, 2, 6

table(sce$donor_id, sce$BraakStage)

table(sce$celltype, sce$BraakStage)

sce

colnames(reducedDim(sce))


# ====================== SIRT expression in different cell types (tSNE) ====

# Load necessary libraries 
library(ggplot2)
library(patchwork) # for combining plots
library(viridis) # for better color scales 

# Get the tSNE coordinates from your SingleCellExperiment object
tsne_coords <- reducedDim(sce, "TSNE")
colnames(tsne_coords) <- colnames(reducedDim(sce, "TSNE"))
tsne_df <- as.data.frame(tsne_coords)

# Add cell metadata
tsne_df$celltype <- colData(sce)$celltype
tsne_df$stage <- colData(sce)$BraakStage 

# Firstly, create a reference tSNE colored by cell type
celltype_tsne <- ggplot(tsne_df, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(aes(color = celltype), size = 0.1, alpha = 0.6) + # alpha - semi-transparent (0.6)
  scale_color_brewer(palette = "Set3", name = "Cell Type") + # a color-blind friendly palette, name the legend "Cell Type"
  labs(title = "Cell type clusters") + # Adds the plot title
  theme_minimal() + # Uses a clean, minimal theme without background grids
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(), # Removes axis labels and ticks
    legend.key.size = unit(0.5, "cm")  # Increase legend key size
  ) + 
  guides(color = guide_legend(
    override.aes = list(size = 3, alpha = 1)  # Larger, opaque points in legend
  ))

print(celltype_tsne)

# Secondly, create a reference tSNE colored by stage
stage_tsne <- ggplot(tsne_df, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(aes(color = stage), size = 0.1, alpha = 1) + # alpha - semi-transparent (0.6)
  scale_color_brewer(palette = "RdYlBu", name = "Braak stage") + # a color-blind friendly palette, name the legend "Cell Type"
  labs(title = "Clusters by Braak stages") + # Adds the plot title
  theme_minimal() + # Uses a clean, minimal theme without background grids
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(), # Removes axis labels and ticks
    legend.key.size = unit(0.5, "cm")  # Increase legend key size
  ) + 
  guides(color = guide_legend(
    override.aes = list(size = 3, alpha = 1)  # Larger, opaque points in legend
  ))

print(stage_tsne)




# ====================== SIRT expression in different cell types (histograms) ====

# Histograms (Kate's code)
library(dplyr)
library(tidyr)

sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")

available_sirt_genes <- sirt_genes[sirt_genes %in% rownames(sce)]

expr_matrix <- assay(sce, 'logcounts')

sirt_expr <- as.matrix(expr_matrix[available_sirt_genes, ])

metadata <- data.frame(
  CellID = colnames(sce),
  CellType = sce$celltype,        
  Stage = sce$BraakStage    
)

plot_data <- as.data.frame(t(sirt_expr)) # транспонирование
plot_data$CellID <- rownames(plot_data)

plot_data <- plot_data %>%
  left_join(metadata, by = "CellID") %>%
  pivot_longer(
    cols = all_of(available_sirt_genes),
    names_to = "Gene",
    values_to = "Expression"
  )

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Stage), !is.na(Expression))

specific_celltypes <- c("neuron", 
                        "oligodendrocyte", 
                        "oligodendrocyte precursor cell", 
                        "astrocyte", 
                        "central nervous system macrophage")  

# Plot histograms for each SIRT across 5 cell types (dividing into AD and Control)
p <- plot_data %>%
  filter(CellType %in% specific_celltypes) %>%
  ggplot(aes(x = Expression)) +
  geom_density(aes(y = after_stat(scaled), color = CellType), 
               linewidth = 0.8, alpha = 0.8) +
  facet_grid(Gene ~ Stage, scales = "free") +
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage of cells") +
  scale_x_continuous(name = "Expression level (logcounts)") +
  labs(
    color = "Cell type") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(p)



# Plot histograms for each cell type across SIRT1-7 (dividing into AD and Control)
celltype_plots <- list()

for (cell_type in specific_celltypes) {
  p <- plot_data %>%
    filter(CellType == cell_type) %>%
    ggplot(aes(x = Expression)) +
    geom_density(aes(y = after_stat(scaled), color = Gene), 
                 linewidth = 0.8, alpha = 0.8) +
    facet_grid(. ~ Stage, scales = "free") +
    scale_y_continuous(labels = scales::percent_format(), name = "Percentage of cells") +
    scale_x_continuous(name = "Expression level (logcounts)") +
    scale_color_brewer(palette = "Set1", name = "SIRT Gene") +
    labs(title = NULL) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.margin = margin(5.5, 40, 5.5, 5.5)  # Add right margin for label
    )
  
  # Add cell type name as annotation on the right
  celltype_plots[[cell_type]] <- p + 
    annotate("text", x = Inf, y = Inf, label = cell_type, 
             hjust = 1.1, vjust = 1.5, size = 5, fontface = "bold")
}

combined_plot <- wrap_plots(celltype_plots, ncol = 1) +
  plot_annotation(title = "SIRT expression distribution across cell types") +
  plot_layout(guides = "collect") &  # Collect all guides into one legend
  theme(legend.position = "bottom")   # Place the legend at the bottom

print(combined_plot)

# Extract each of the plot separately 
celltype_plots[["neuron"]]
celltype_plots[["oligodendrocyte"]]
celltype_plots[["oligodendrocyte precursor cell"]]
celltype_plots[["astrocyte"]]
celltype_plots[["central nervous system macrophage"]]



# My old histograms 
# Define threshold for low expressed genes by checking the distribution of expression values for one SIRT gene
expr_values <- assay(sce, "logcounts")["SIRT1", ]
summary(expr_values) # key statistics 
hist(expr_values[expr_values > 0], breaks = 50, main = "SIRT1 Expression Distribution")

expr_values <- assay(sce, "logcounts")["SIRT2", ]
summary(expr_values) # key statistics 
hist(expr_values[expr_values > 0], breaks = 50, main = "SIRT2 Expression Distribution")

expr_values <- assay(sce, "logcounts")["SIRT3", ]
summary(expr_values) # key statistics 
hist(expr_values[expr_values > 0], breaks = 50, main = "SIRT3 Expression Distribution")

expr_values <- assay(sce, "logcounts")["SIRT4", ]
summary(expr_values) # key statistics 
hist(expr_values[expr_values > 0], breaks = 50, main = "SIRT4 Expression Distribution")

expr_values <- assay(sce, "logcounts")["SIRT5", ]
summary(expr_values) # key statistics 
hist(expr_values[expr_values > 0], breaks = 50, main = "SIRT5 Expression Distribution")

expr_values <- assay(sce, "logcounts")["SIRT6", ]
summary(expr_values) # key statistics 
hist(expr_values[expr_values > 0], breaks = 50, main = "SIRT6 Expression Distribution")

expr_values <- assay(sce, "logcounts")["SIRT7", ]
summary(expr_values) # key statistics 
hist(expr_values[expr_values > 0], breaks = 50, main = "SIRT7 Expression Distribution")

# Based on expression distribution low_expr_threshold <- 0.1



# Add SIRT expression values (using logcounts for visualization)
logcounts_matrix <- assay(sce, "logcounts")

# Create individual UMAP plots for each SIRT gene
sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")
sirt_tsne_plots <- list()

for (sirt_gene in sirt_genes) {
  if (sirt_gene %in% rownames(logcounts_matrix)) {
    # Add expression for this SIRT gene
    tsne_df$expression <- logcounts_matrix[sirt_gene, ]
    
    # Define threshold for low expression
    low_expr_threshold <- 0.1
    
    # Create the tSNE plot
    p <- ggplot(tsne_df, aes(x = tSNE_1, y = tSNE_2)) +
      geom_point(aes(color = ifelse(expression < low_expr_threshold, NA, expression)), # Colors points based on expression, but sets low-expression cells to NA so they get the na.value color
                 size = 0.1, alpha = 0.6) +
      scale_color_viridis_c( # Applies a color scale
        name = "Expression",
        option = "plasma", # palette name 
        na.value = "grey80",  # Color for values below threshold
        limits = c(low_expr_threshold, max(tsne_df$expression, na.rm = TRUE)) # limits = c(low_expr_threshold, max(...)) - 
      ) +
      labs(
        title = paste("SIRT Expression:", sirt_gene),
        subtitle = paste("Grey: expression <", low_expr_threshold, "| Colors: expression >", low_expr_threshold)
      ) +
      theme_minimal() +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold")
      )
    
    sirt_tsne_plots[[sirt_gene]] <- p # Stores each plot in a list for later use
  } else {
    cat("SIRT gene", sirt_gene, "not found in the dataset\n")
  }
}

# Arrange all plots in a grid
combined_tsne <- wrap_plots(sirt_tsne_plots, ncol = 3) +
  plot_annotation(title = "SIRT gene expression across cell types (tSNE)")

print(combined_tsne)




# ================== SIRT expression in different cell types (heatmap) ====

# 1. Heatmap with average SIRT expression 
# Calculate mean expression of SIRTs in each cell type
cell_types <- unique(colData(sce)$celltype)

# Create an empty matrix to store average expression
expr_matrix <- matrix(NA, nrow = length(cell_types), ncol = length(sirt_genes))
rownames(expr_matrix) <- cell_types
colnames(expr_matrix) <- sirt_genes

for (cell_type in cell_types) {
  cell_indices <- which(colData(sce)$celltype == cell_type)
  
  for (sirt_gene in sirt_genes) {
    if (sirt_gene %in% rownames(logcounts_matrix)) {
      expr_values <- logcounts_matrix[sirt_gene, cell_indices]
      if (length(expr_values) > 0) {
        expr_matrix[cell_type, sirt_gene] <- mean(expr_values, na.rm = TRUE)
      }
    }
  }
}

# Convert to data frame for plotting
expr_df <- as.data.frame(expr_matrix)
expr_df$cell_type <- rownames(expr_df)

# Reshape to long format for ggplot
library(reshape2)
expr_long <- melt(expr_df, id.vars = "cell_type", 
                  variable.name = "SIRT", 
                  value.name = "expression")

# Create the heatmap
heatmap_plot <- ggplot(expr_long, aes(x = SIRT, y = reorder(cell_type, expression), fill = expression)) +
  geom_tile(color = "white", linewidth = 0.5) +
  # Add text with expression values
  geom_text(aes(label = sprintf("%.2f", expression)), color = "white", size = 3, fontface = "bold") +
  scale_fill_viridis_c(name = "Mean\nExpression", option = "viridis") +
  labs(title = "Average SIRT expression across cell types",
       x = "SIRT genes",
       y = "Cell types") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(heatmap_plot)


# 2. Heatmap with percentage of cells with non-zero SIRT expression 
# Define the specific order for cell types (order the same as in the 1st heatmap for more convenient comparison)
cell_type_order <- c("endothelial cell", "vascular associated smooth muscle cell", 
                     "pericyte", "central nervous system macrophage", 
                     "leukocyte", "fibroblast", "astrocyte", "neuron", 
                     "oligodendrocyte precursor cell", "oligodendrocyte")

# Create an empty matrix to store percentage expression
pct_matrix <- matrix(NA, nrow = length(cell_types), ncol = length(sirt_genes))
rownames(pct_matrix) <- cell_types
colnames(pct_matrix) <- sirt_genes

for (cell_type in cell_types) {
  cell_indices <- which(colData(sce)$celltype == cell_type)
  
  for (sirt_gene in sirt_genes) {
    if (sirt_gene %in% rownames(logcounts_matrix)) {
      expr_values <- logcounts_matrix[sirt_gene, cell_indices]
      if (length(expr_values) > 0) {
        # Calculate percentage of cells with expression > 0
        pct_matrix[cell_type, sirt_gene] <- mean(expr_values > 0) * 100 
      }
    }
  }
}

# Convert to data frame for plotting
pct_df <- as.data.frame(pct_matrix)
pct_df$cell_type <- rownames(pct_df)

# Reshape to long format for ggplot
library(reshape2)
pct_long <- melt(pct_df, id.vars = "cell_type", 
                 variable.name = "SIRT", 
                 value.name = "percentage")

# Convert cell_type to factor with specified order
pct_long$cell_type <- factor(pct_long$cell_type, levels = cell_type_order)

# Create the heatmap
heatmap_plot <- ggplot(pct_long, aes(x = SIRT, y = cell_type, fill = percentage)) +
  geom_tile(color = "white", linewidth = 0.5) +
  # Add text with percentage values
  geom_text(aes(label = sprintf("%.1f%%", percentage)), color = "white", size = 3, fontface = "bold") +
  scale_fill_viridis_c(name = "% Cells\nExpressing", option = "plasma") +
  labs(title = "Percentage of cells expressing SIRT genes",
       x = "SIRT genes",
       y = "Cell types") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(heatmap_plot)







# Summary table with the number of cells in each cell type with SIRTs expression > 0
# Create summary table
summary_table <- data.frame()

for (cell_type in unique(colData(sce)$celltype)) {
  # Get cells of this type
  cell_indices <- which(colData(sce)$celltype == cell_type)
  total_cells <- length(cell_indices)
  
  for (sirt_gene in sirt_genes) {
    if (sirt_gene %in% rownames(assay(sce, "logcounts"))) {
      expr_values <- assay(sce, "logcounts")[sirt_gene, cell_indices]
      
      # Calculate statistics
      cells_expressing <- sum(expr_values > 0)
      pct_expressing <- (cells_expressing / total_cells) * 100
      avg_expression <- mean(expr_values, na.rm = TRUE)
      avg_expression_expressing <- ifelse(cells_expressing > 0, 
                                          mean(expr_values[expr_values > 0], na.rm = TRUE), 
                                          0)
      
      summary_table <- rbind(summary_table, data.frame(
        cell_type = cell_type,
        sirt_gene = sirt_gene,
        total_cells = total_cells,
        cells_expressing = cells_expressing,
        pct_expressing = pct_expressing,
        avg_expression_all = avg_expression,
        avg_expression_expressing = avg_expression_expressing
      ))
    }
  }
}

# Print the table
print(summary_table)



table(sce$celltype, sce$BraakStage)




# =========================== Pseudobulk DE analysis (using edgeR) ====

library(edgeR)
library(limma)
library(ggplot2)
library(dplyr)

cat("=== ANALYSIS 1: Comparing BraakStage 0 vs 2 ===\n")

# Define key columns
donor_column <- "donor_id"
stage_column <- "BraakStage"
celltype_column <- "celltype"

# Get SIRT genes
sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")

# Create empty list to store results for 0 vs 2 comparison
edgeR_results_list_0vs2 <- list()

# Define cell types to analyze
cell_types_to_analyze <- c("neuron", "oligodendrocyte", "astrocyte", 
                           "oligodendrocyte precursor cell", 
                           "central nervous system macrophage")

# ANALYSIS 1: 0 vs 2 comparison
for (cell_type_of_interest in cell_types_to_analyze) {
  
  cat("\n\n=== Processing cell type:", cell_type_of_interest, "(0 vs 2) ===\n")
  
  # A. Subset the sce object to the specific cell type
  sce_subset <- sce[, sce[[celltype_column]] == cell_type_of_interest]
  cat("Number of cells in subset:", ncol(sce_subset), "\n")
  
  # B. Create Pseudobulk matrix
  count_matrix <- assay(sce_subset, "counts")
  donor_ids <- colData(sce_subset)[[donor_column]]
  unique_donors <- unique(donor_ids)
  
  # Initialize pseudobulk matrix
  pb_matrix <- matrix(0, nrow = nrow(count_matrix), ncol = length(unique_donors))
  rownames(pb_matrix) <- rownames(count_matrix)
  colnames(pb_matrix) <- unique_donors
  
  # Sum counts for each donor
  for (donor in unique_donors) {
    donor_cells <- which(donor_ids == donor)
    if (length(donor_cells) > 0) {
      pb_matrix[, donor] <- rowSums(count_matrix[, donor_cells, drop = FALSE])
    }
  }
  
  # C. Prepare metadata for 0 vs 2 comparison
  metadata <- unique(colData(sce_subset)[, c(donor_column, stage_column)])
  metadata <- metadata[match(colnames(pb_matrix), metadata[[donor_column]]), ]
  
  # Keep only samples with BraakStage 0 or 2
  metadata <- metadata[metadata[[stage_column]] %in% c(0, 2), ]
  pb_matrix <- pb_matrix[, colnames(pb_matrix) %in% metadata[[donor_column]], drop = FALSE]
  
  # Check if we have samples in both groups
  if (nrow(metadata) == 0 || length(unique(metadata[[stage_column]])) < 2) {
    cat("Skipping - not enough samples in both BraakStage 0 and 2\n")
    next
  }
  
  rownames(metadata) <- metadata[[donor_column]]
  metadata[[stage_column]] <- factor(metadata[[stage_column]], levels = c(0, 2))
  
  cat("Samples distribution (0 vs 2):\n")
  print(table(metadata[[stage_column]]))
  
  # D. edgeR Analysis for 0 vs 2
  # Create DGEList object
  dge <- DGEList(counts = pb_matrix, group = metadata[[stage_column]])
  
  # Pre-filter lowly expressed genes
  keep <- filterByExpr(dge, group = metadata[[stage_column]])
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  cat("Number of genes after filtering:", nrow(dge), "\n")
  
  # Normalize
  dge <- calcNormFactors(dge)
  
  # Create design matrix
  design <- model.matrix(~ metadata[[stage_column]])
  colnames(design) <- gsub("metadata\\[\\[stage_column\\]\\]", "", colnames(design))
  
  # Estimate dispersions
  dge <- estimateDisp(dge, design)
  
  # Fit model and test (BraakStage2 vs BraakStage0)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = 2)
  
  # Get results
  res <- topTags(qlf, n = Inf, sort.by = "PValue")
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Store results
  edgeR_results_list_0vs2[[cell_type_of_interest]] <- res_df
  
  # Print SIRT results
  sirt_results <- res_df[res_df$gene %in% sirt_genes, ]
  cat("\nSIRT gene results for 0 vs 2:\n")
  print(sirt_results)
}

cat("\n\n=== ANALYSIS 1 COMPLETE: 0 vs 2 comparison saved in 'edgeR_results_list_0vs2' ===\n")

# Description of the result edgeR table:
# 1. logFC - log2 fold change(how much expression changes in AD vs NC)
# 2. logCPM - log2 counts per million (average expression level)
# 3. F - F-statistic (effect size relative to variability)
# 4. PValue - raw p-value (probability of seeing this effect by chance if no real difference)
# 5. FDR - false discovery rate (adjusted p-value for multiple testing)



cat("\n\n\n=== ANALYSIS 2: Comparing BraakStage 0 vs 6 ===\n")

# Create empty list to store results for 0 vs 6 comparison
edgeR_results_list_0vs6 <- list()

# ANALYSIS 2: 0 vs 6 comparison
for (cell_type_of_interest in cell_types_to_analyze) {
  
  cat("\n\n=== Processing cell type:", cell_type_of_interest, "(0 vs 6) ===\n")
  
  # A. Subset the sce object to the specific cell type
  sce_subset <- sce[, sce[[celltype_column]] == cell_type_of_interest]
  cat("Number of cells in subset:", ncol(sce_subset), "\n")
  
  # B. Create Pseudobulk matrix
  count_matrix <- assay(sce_subset, "counts")
  donor_ids <- colData(sce_subset)[[donor_column]]
  unique_donors <- unique(donor_ids)
  
  # Initialize pseudobulk matrix
  pb_matrix <- matrix(0, nrow = nrow(count_matrix), ncol = length(unique_donors))
  rownames(pb_matrix) <- rownames(count_matrix)
  colnames(pb_matrix) <- unique_donors
  
  # Sum counts for each donor
  for (donor in unique_donors) {
    donor_cells <- which(donor_ids == donor)
    if (length(donor_cells) > 0) {
      pb_matrix[, donor] <- rowSums(count_matrix[, donor_cells, drop = FALSE])
    }
  }
  
  # C. Prepare metadata for 0 vs 6 comparison
  metadata <- unique(colData(sce_subset)[, c(donor_column, stage_column)])
  metadata <- metadata[match(colnames(pb_matrix), metadata[[donor_column]]), ]
  
  # Keep only samples with BraakStage 0 or 6
  metadata <- metadata[metadata[[stage_column]] %in% c(0, 6), ]
  pb_matrix <- pb_matrix[, colnames(pb_matrix) %in% metadata[[donor_column]], drop = FALSE]
  
  # Check if we have samples in both groups
  if (nrow(metadata) == 0 || length(unique(metadata[[stage_column]])) < 2) {
    cat("Skipping - not enough samples in both BraakStage 0 and 6\n")
    next
  }
  
  rownames(metadata) <- metadata[[donor_column]]
  metadata[[stage_column]] <- factor(metadata[[stage_column]], levels = c(0, 6))
  
  cat("Samples distribution (0 vs 6):\n")
  print(table(metadata[[stage_column]]))
  
  # D. edgeR Analysis for 0 vs 6
  # Create DGEList object
  dge <- DGEList(counts = pb_matrix, group = metadata[[stage_column]])
  
  # Pre-filter lowly expressed genes
  keep <- filterByExpr(dge, group = metadata[[stage_column]])
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  cat("Number of genes after filtering:", nrow(dge), "\n")
  
  # Normalize
  dge <- calcNormFactors(dge)
  
  # Create design matrix
  design <- model.matrix(~ metadata[[stage_column]])
  colnames(design) <- gsub("metadata\\[\\[stage_column\\]\\]", "", colnames(design))
  
  # Estimate dispersions
  dge <- estimateDisp(dge, design)
  
  # Fit model and test (BraakStage6 vs BraakStage0)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = 2)
  
  # Get results
  res <- topTags(qlf, n = Inf, sort.by = "PValue")
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Store results
  edgeR_results_list_0vs6[[cell_type_of_interest]] <- res_df
  
  # Print SIRT results
  sirt_results <- res_df[res_df$gene %in% sirt_genes, ]
  cat("\nSIRT gene results for 0 vs 6:\n")
  print(sirt_results)
}

cat("\n\n=== ANALYSIS 2 COMPLETE: 0 vs 6 comparison saved in 'edgeR_results_list_0vs6' ===\n")




# =========================== edgeR results visualisation ====
# Load necessary library for plotting
library(ggplot2)
# install.packages("ggrepel")
library(ggrepel)  # For better label placement

# I. 0 vs 2 Braak stages ====

# 1. For neuron cells

# Get the results for neurons
neuron_results <- edgeR_results_list_0vs2[["neuron"]]

# Convert to data frame for easier handling
neuron_df <- as.data.frame(neuron_results)
neuron_df$gene <- rownames(neuron_df)

# Create a column to highlight SIRT genes
neuron_df$gene_type <- "other"
neuron_df$gene_type[neuron_df$gene %in% c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")] <- "SIRT"

# Create significance threshold column with padj <= 0.1
neuron_df$significant <- neuron_df$FDR <= 0.1

# Create the volcano plot
volcano_plot <- ggplot(neuron_df, aes(x = logFC, y = -log10(PValue))) +
  # Plot all points
  geom_point(aes(color = gene_type, alpha = significant), size = 1.5) +
  
  # Highlight significant points (padj <= 0.1)
  geom_point(data = subset(neuron_df, significant == TRUE), 
             color = "red", size = 2.5) +
  
  # Add gene labels for significant points (padj <= 0.1)
  geom_text_repel(data = subset(neuron_df, significant == TRUE),
                  aes(label = gene),
                  color = "red",
                  size = 3,
                  max.overlaps = 20,  # Increase if you want to see more labels
                  box.padding = 0.5,
                  segment.color = "red",
                  segment.alpha = 0.5) +
  
  # Highlight SIRT genes with labels
  geom_point(data = subset(neuron_df, gene_type == "SIRT"), 
             color = "blue", size = 2.5) +
  geom_text_repel(data = subset(neuron_df, gene_type == "SIRT"), 
                  aes(label = gene), 
                  color = "blue", 
                  size = 3,
                  box.padding = 0.3,
                  segment.color = "blue",
                  segment.alpha = 0.3) +
  
  # Add significance thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgreen", alpha = 0.5) +
  
  # Customize colors and theme
  scale_color_manual(values = c("SIRT" = "blue", "other" = "gray60")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
  
  # Labels and title
  labs(title = "Differential Expression in Neurons: 0 vs 2 stages",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change",
       y = "-log10(p-value)") +
  
  theme_minimal() +
  theme(legend.position = "right")

print(volcano_plot)



# Extract SIRT genes from the neuron results
sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")
sirt_results <- neuron_df[neuron_df$gene %in% sirt_genes, ]

# Check if we found all SIRT genes
cat("SIRT genes found:", nrow(sirt_results), "\n")
print(sirt_results$gene)

# If some SIRT genes are missing, they might have been filtered out during DESeq2 analysis
# missing_sirts <- setdiff(sirt_genes, sirt_results$gene)
# if (length(missing_sirts) > 0) {
# cat("Missing SIRT genes (filtered out due to low expression):", missing_sirts, "\n")
# }

# Calculate standard error from p-values
sirt_results$lfcSE <- abs(sirt_results$logFC) / qnorm(1 - sirt_results$PValue/2)

# Calculate 95% confidence intervals
sirt_results$CI_lower <- sirt_results$logFC - 1.96 * sirt_results$lfcSE
sirt_results$CI_upper <- sirt_results$logFC + 1.96 * sirt_results$lfcSE

# Create forest plot
forest_plot <- ggplot(sirt_results, aes(x = reorder(gene, logFC), y = logFC)) +
  geom_point(aes(color = FDR < 0.1), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, color = FDR < 0.1), 
                width = 0.2, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     name = "Significant (FDR < 0.1)") +
  labs(title = "SIRT expression in neurons",
       x = "SIRT genes",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)



# 2. For oligodendrocytes

# Get the results for oligodendrocytes
oligodendrocyte_results <- edgeR_results_list_0vs2[["oligodendrocyte"]]

# Convert to data frame for easier handling
oligodendrocyte_df <- as.data.frame(oligodendrocyte_results)
oligodendrocyte_df$gene <- rownames(oligodendrocyte_df)

# Create a column to highlight SIRT genes
oligodendrocyte_df$gene_type <- "other"
oligodendrocyte_df$gene_type[oligodendrocyte_df$gene %in% c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")] <- "SIRT"

# Create significance threshold column with padj <= 0.1
oligodendrocyte_df$significant <- oligodendrocyte_df$FDR <= 0.1

# Create the volcano plot
volcano_plot <- ggplot(oligodendrocyte_df, aes(x = logFC, y = -log10(PValue))) +
  # Plot all points
  geom_point(aes(color = gene_type, alpha = significant), size = 1.5) +
  
  # Highlight significant points (padj <= 0.1)
  geom_point(data = subset(oligodendrocyte_df, significant == TRUE), 
             color = "red", size = 2.5) +
  
  # Add gene labels for significant points (padj <= 0.1)
  geom_text_repel(data = subset(oligodendrocyte_df, significant == TRUE),
                  aes(label = gene),
                  color = "red",
                  size = 3,
                  max.overlaps = 20,  # Increase if you want to see more labels
                  box.padding = 0.5,
                  segment.color = "red",
                  segment.alpha = 0.5) +
  
  # Highlight SIRT genes with labels
  geom_point(data = subset(oligodendrocyte_df, gene_type == "SIRT"), 
             color = "blue", size = 2.5) +
  geom_text_repel(data = subset(oligodendrocyte_df, gene_type == "SIRT"), 
                  aes(label = gene), 
                  color = "blue", 
                  size = 3,
                  box.padding = 0.3,
                  segment.color = "blue",
                  segment.alpha = 0.3) +
  
  # Add significance thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgreen", alpha = 0.5) +
  
  # Customize colors and theme
  scale_color_manual(values = c("SIRT" = "blue", "other" = "gray60")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
  
  # Labels and title
  labs(title = "Differential Expression in Oligodendrocytes: 0 vs 2 stages",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change",
       y = "-log10(p-value)") +
  
  theme_minimal() +
  theme(legend.position = "right")

print(volcano_plot)



# Extract SIRT genes from the oligodendrocytes results
sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")
sirt_results <- oligodendrocyte_df[oligodendrocyte_df$gene %in% sirt_genes, ]

# Check if we found all SIRT genes
cat("SIRT genes found:", nrow(sirt_results), "\n")
print(sirt_results$gene)

# Calculate standard error from p-values
sirt_results$lfcSE <- abs(sirt_results$logFC) / qnorm(1 - sirt_results$PValue/2)

# Calculate 95% confidence intervals
sirt_results$CI_lower <- sirt_results$logFC - 1.96 * sirt_results$lfcSE
sirt_results$CI_upper <- sirt_results$logFC + 1.96 * sirt_results$lfcSE

# Create forest plot
forest_plot <- ggplot(sirt_results, aes(x = reorder(gene, logFC), y = logFC)) +
  geom_point(aes(color = FDR < 0.1), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, color = FDR < 0.1), 
                width = 0.2, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     name = "Significant (FDR < 0.1)") +
  labs(title = "SIRT expression in oligodendrocytes",
       x = "SIRT genes",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)



# 3. For astrocytes

# Get the results for astrocytes
astrocyte_results <- edgeR_results_list_0vs2[["astrocyte"]]

# Convert to data frame for easier handling
astrocyte_df <- as.data.frame(astrocyte_results)
astrocyte_df$gene <- rownames(astrocyte_df)

# Create a column to highlight SIRT genes
astrocyte_df$gene_type <- "other"
astrocyte_df$gene_type[astrocyte_df$gene %in% c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")] <- "SIRT"

# Create significance threshold column with padj <= 0.1
astrocyte_df$significant <- astrocyte_df$FDR <= 0.1

# Create the volcano plot
volcano_plot <- ggplot(astrocyte_df, aes(x = logFC, y = -log10(PValue))) +
  # Plot all points
  geom_point(aes(color = gene_type, alpha = significant), size = 1.5) +
  
  # Highlight significant points (padj <= 0.1)
  geom_point(data = subset(astrocyte_df, significant == TRUE), 
             color = "red", size = 2.5) +
  
  # Add gene labels for significant points (padj <= 0.1)
  geom_text_repel(data = subset(astrocyte_df, significant == TRUE),
                  aes(label = gene),
                  color = "red",
                  size = 3,
                  max.overlaps = 20,  # Increase if you want to see more labels
                  box.padding = 0.5,
                  segment.color = "red",
                  segment.alpha = 0.5) +
  
  # Highlight SIRT genes with labels
  geom_point(data = subset(astrocyte_df, gene_type == "SIRT"), 
             color = "blue", size = 2.5) +
  geom_text_repel(data = subset(astrocyte_df, gene_type == "SIRT"), 
                  aes(label = gene), 
                  color = "blue", 
                  size = 3,
                  box.padding = 0.3,
                  segment.color = "blue",
                  segment.alpha = 0.3) +
  
  # Add significance thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgreen", alpha = 0.5) +
  
  # Customize colors and theme
  scale_color_manual(values = c("SIRT" = "blue", "other" = "gray60")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
  
  # Labels and title
  labs(title = "Differential Expression in Astrocytes: 0 vs 2 stages",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change",
       y = "-log10(p-value)") +
  
  theme_minimal() +
  theme(legend.position = "right")

print(volcano_plot)



# Extract SIRT genes from the astrocytes results
sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")
sirt_results <- astrocyte_df[astrocyte_df$gene %in% sirt_genes, ]

# Check if we found all SIRT genes
cat("SIRT genes found:", nrow(sirt_results), "\n")
print(sirt_results$gene)

# Calculate standard error from p-values
sirt_results$lfcSE <- abs(sirt_results$logFC) / qnorm(1 - sirt_results$PValue/2)

# Calculate 95% confidence intervals
sirt_results$CI_lower <- sirt_results$logFC - 1.96 * sirt_results$lfcSE
sirt_results$CI_upper <- sirt_results$logFC + 1.96 * sirt_results$lfcSE

# Create forest plot
forest_plot <- ggplot(sirt_results, aes(x = reorder(gene, logFC), y = logFC)) +
  geom_point(aes(color = FDR < 0.1), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, color = FDR < 0.1), 
                width = 0.2, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     name = "Significant (FDR < 0.1)") +
  labs(title = "SIRT expression in astrocytes",
       x = "SIRT genes",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)



# 4. For oligodendrocytes precursor cells

# Get the results for oligodendrocytes precursor cells
oligodendrocytePrecursor_results <- edgeR_results_list_0vs2[["oligodendrocyte precursor cell"]]

# Convert to data frame for easier handling
oligodendrocytePrecursor_df <- as.data.frame(oligodendrocytePrecursor_results)
oligodendrocytePrecursor_df$gene <- rownames(oligodendrocytePrecursor_df)

# Create a column to highlight SIRT genes
oligodendrocytePrecursor_df$gene_type <- "other"
oligodendrocytePrecursor_df$gene_type[oligodendrocytePrecursor_df$gene %in% c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")] <- "SIRT"

# Create significance threshold column with padj <= 0.1
oligodendrocytePrecursor_df$significant <- oligodendrocytePrecursor_df$FDR <= 0.1

# Create the volcano plot
volcano_plot <- ggplot(oligodendrocytePrecursor_df, aes(x = logFC, y = -log10(PValue))) +
  # Plot all points
  geom_point(aes(color = gene_type, alpha = significant), size = 1.5) +
  
  # Highlight significant points (padj <= 0.1)
  geom_point(data = subset(oligodendrocytePrecursor_df, significant == TRUE), 
             color = "red", size = 2.5) +
  
  # Add gene labels for significant points (padj <= 0.1)
  geom_text_repel(data = subset(oligodendrocytePrecursor_df, significant == TRUE),
                  aes(label = gene),
                  color = "red",
                  size = 3,
                  max.overlaps = 20,  # Increase if you want to see more labels
                  box.padding = 0.5,
                  segment.color = "red",
                  segment.alpha = 0.5) +
  
  # Highlight SIRT genes with labels
  geom_point(data = subset(oligodendrocytePrecursor_df, gene_type == "SIRT"), 
             color = "blue", size = 2.5) +
  geom_text_repel(data = subset(oligodendrocytePrecursor_df, gene_type == "SIRT"), 
                  aes(label = gene), 
                  color = "blue", 
                  size = 3,
                  box.padding = 0.3,
                  segment.color = "blue",
                  segment.alpha = 0.3) +
  
  # Add significance thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgreen", alpha = 0.5) +
  
  # Customize colors and theme
  scale_color_manual(values = c("SIRT" = "blue", "other" = "gray60")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
  
  # Labels and title
  labs(title = "Differential Expression in oligodendrocytes precursor cells: 0 vs 2 stages",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change",
       y = "-log10(p-value)") +
  
  theme_minimal() +
  theme(legend.position = "right")

print(volcano_plot)



# Extract SIRT genes from the oligodendrocyte precursor cells results
sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")
sirt_results <- oligodendrocytePrecursor_df[oligodendrocytePrecursor_df$gene %in% sirt_genes, ]

# Check if we found all SIRT genes
cat("SIRT genes found:", nrow(sirt_results), "\n")
print(sirt_results$gene)

# Calculate standard error from p-values
sirt_results$lfcSE <- abs(sirt_results$logFC) / qnorm(1 - sirt_results$PValue/2)

# Calculate 95% confidence intervals
sirt_results$CI_lower <- sirt_results$logFC - 1.96 * sirt_results$lfcSE
sirt_results$CI_upper <- sirt_results$logFC + 1.96 * sirt_results$lfcSE

# Create forest plot
forest_plot <- ggplot(sirt_results, aes(x = reorder(gene, logFC), y = logFC)) +
  geom_point(aes(color = FDR < 0.1), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, color = FDR < 0.1), 
                width = 0.2, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     name = "Significant (FDR < 0.1)") +
  labs(title = "SIRT expression in oligodendrocytes precursor cells",
       x = "SIRT genes",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)



# 5. For central nervous system macrophage

# Get the results for central nervous system macrophage
CenNervSysMacr_results <- edgeR_results_list_0vs2[["central nervous system macrophage"]]

# Convert to data frame for easier handling
CenNervSysMacr_df <- as.data.frame(CenNervSysMacr_results)
CenNervSysMacr_df$gene <- rownames(CenNervSysMacr_df)

# Create a column to highlight SIRT genes
CenNervSysMacr_df$gene_type <- "other"
CenNervSysMacr_df$gene_type[CenNervSysMacr_df$gene %in% c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")] <- "SIRT"

# Create significance threshold column with padj <= 0.1
CenNervSysMacr_df$significant <- CenNervSysMacr_df$FDR <= 0.1

# Create the volcano plot
volcano_plot <- ggplot(CenNervSysMacr_df, aes(x = logFC, y = -log10(PValue))) +
  # Plot all points
  geom_point(aes(color = gene_type, alpha = significant), size = 1.5) +
  
  # Highlight significant points (padj <= 0.1)
  geom_point(data = subset(CenNervSysMacr_df, significant == TRUE), 
             color = "red", size = 2.5) +
  
  # Add gene labels for significant points (padj <= 0.1)
  geom_text_repel(data = subset(CenNervSysMacr_df, significant == TRUE),
                  aes(label = gene),
                  color = "red",
                  size = 3,
                  max.overlaps = 20,  # Increase if you want to see more labels
                  box.padding = 0.5,
                  segment.color = "red",
                  segment.alpha = 0.5) +
  
  # Highlight SIRT genes with labels
  geom_point(data = subset(CenNervSysMacr_df, gene_type == "SIRT"), 
             color = "blue", size = 2.5) +
  geom_text_repel(data = subset(CenNervSysMacr_df, gene_type == "SIRT"), 
                  aes(label = gene), 
                  color = "blue", 
                  size = 3,
                  box.padding = 0.3,
                  segment.color = "blue",
                  segment.alpha = 0.3) +
  
  # Add significance thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgreen", alpha = 0.5) +
  
  # Customize colors and theme
  scale_color_manual(values = c("SIRT" = "blue", "other" = "gray60")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
  
  # Labels and title
  labs(title = "Differential Expression in Central Nervous System Macrophages: 0 vs 2 stages",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change",
       y = "-log10(p-value)") +
  
  theme_minimal() +
  theme(legend.position = "right")

print(volcano_plot)



# Extract SIRT genes from the Central Nervous System Macrophages results
sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")
sirt_results <- CenNervSysMacr_df[CenNervSysMacr_df$gene %in% sirt_genes, ]

# Check if we found all SIRT genes
cat("SIRT genes found:", nrow(sirt_results), "\n")
print(sirt_results$gene)

# Calculate standard error from p-values
sirt_results$lfcSE <- abs(sirt_results$logFC) / qnorm(1 - sirt_results$PValue/2)

# Calculate 95% confidence intervals
sirt_results$CI_lower <- sirt_results$logFC - 1.96 * sirt_results$lfcSE
sirt_results$CI_upper <- sirt_results$logFC + 1.96 * sirt_results$lfcSE

# Create forest plot
forest_plot <- ggplot(sirt_results, aes(x = reorder(gene, logFC), y = logFC)) +
  geom_point(aes(color = FDR < 0.1), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, color = FDR < 0.1), 
                width = 0.2, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     name = "Significant (FDR < 0.1)") +
  labs(title = "SIRT expression in central nervous system macrophages",
       x = "SIRT genes",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)



# II. 0 vs 6 Braak stages ====

# 1. For neuron cells

# Get the results for neurons
neuron_results <- edgeR_results_list_0vs6[["neuron"]]

# Convert to data frame for easier handling
neuron_df <- as.data.frame(neuron_results)
neuron_df$gene <- rownames(neuron_df)

# Create a column to highlight SIRT genes
neuron_df$gene_type <- "other"
neuron_df$gene_type[neuron_df$gene %in% c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")] <- "SIRT"

# Create significance threshold column with padj <= 0.1
neuron_df$significant <- neuron_df$FDR <= 0.1

# Create the volcano plot
volcano_plot <- ggplot(neuron_df, aes(x = logFC, y = -log10(PValue))) +
  # Plot all points
  geom_point(aes(color = gene_type, alpha = significant), size = 1.5) +
  
  # Highlight significant points (padj <= 0.1)
  geom_point(data = subset(neuron_df, significant == TRUE), 
             color = "red", size = 2.5) +
  
  # Add gene labels for significant points (padj <= 0.1)
  geom_text_repel(data = subset(neuron_df, significant == TRUE),
                  aes(label = gene),
                  color = "red",
                  size = 3,
                  max.overlaps = 20,  # Increase if you want to see more labels
                  box.padding = 0.5,
                  segment.color = "red",
                  segment.alpha = 0.5) +
  
  # Highlight SIRT genes with labels
  geom_point(data = subset(neuron_df, gene_type == "SIRT"), 
             color = "blue", size = 2.5) +
  geom_text_repel(data = subset(neuron_df, gene_type == "SIRT"), 
                  aes(label = gene), 
                  color = "blue", 
                  size = 3,
                  box.padding = 0.3,
                  segment.color = "blue",
                  segment.alpha = 0.3) +
  
  # Add significance thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgreen", alpha = 0.5) +
  
  # Customize colors and theme
  scale_color_manual(values = c("SIRT" = "blue", "other" = "gray60")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
  
  # Labels and title
  labs(title = "Differential Expression in Neurons: 0 vs 6 stages",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change",
       y = "-log10(p-value)") +
  
  theme_minimal() +
  theme(legend.position = "right")

print(volcano_plot)



# Extract SIRT genes from the neuron results
sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")
sirt_results <- neuron_df[neuron_df$gene %in% sirt_genes, ]

# Check if we found all SIRT genes
cat("SIRT genes found:", nrow(sirt_results), "\n")
print(sirt_results$gene)

# If some SIRT genes are missing, they might have been filtered out during DESeq2 analysis
# missing_sirts <- setdiff(sirt_genes, sirt_results$gene)
# if (length(missing_sirts) > 0) {
# cat("Missing SIRT genes (filtered out due to low expression):", missing_sirts, "\n")
# }

# Calculate standard error from p-values
sirt_results$lfcSE <- abs(sirt_results$logFC) / qnorm(1 - sirt_results$PValue/2)

# Calculate 95% confidence intervals
sirt_results$CI_lower <- sirt_results$logFC - 1.96 * sirt_results$lfcSE
sirt_results$CI_upper <- sirt_results$logFC + 1.96 * sirt_results$lfcSE

# Create forest plot
forest_plot <- ggplot(sirt_results, aes(x = reorder(gene, logFC), y = logFC)) +
  geom_point(aes(color = FDR < 0.1), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, color = FDR < 0.1), 
                width = 0.2, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     name = "Significant (FDR < 0.1)") +
  labs(title = "SIRT expression in neurons",
       x = "SIRT genes",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)



# 2. For oligodendrocytes

# Get the results for oligodendrocytes
oligodendrocyte_results <- edgeR_results_list_0vs6[["oligodendrocyte"]]

# Convert to data frame for easier handling
oligodendrocyte_df <- as.data.frame(oligodendrocyte_results)
oligodendrocyte_df$gene <- rownames(oligodendrocyte_df)

# Create a column to highlight SIRT genes
oligodendrocyte_df$gene_type <- "other"
oligodendrocyte_df$gene_type[oligodendrocyte_df$gene %in% c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")] <- "SIRT"

# Create significance threshold column with padj <= 0.1
oligodendrocyte_df$significant <- oligodendrocyte_df$FDR <= 0.1

# Create the volcano plot
volcano_plot <- ggplot(oligodendrocyte_df, aes(x = logFC, y = -log10(PValue))) +
  # Plot all points
  geom_point(aes(color = gene_type, alpha = significant), size = 1.5) +
  
  # Highlight significant points (padj <= 0.1)
  geom_point(data = subset(oligodendrocyte_df, significant == TRUE), 
             color = "red", size = 2.5) +
  
  # Add gene labels for significant points (padj <= 0.1)
  geom_text_repel(data = subset(oligodendrocyte_df, significant == TRUE),
                  aes(label = gene),
                  color = "red",
                  size = 3,
                  max.overlaps = 20,  # Increase if you want to see more labels
                  box.padding = 0.5,
                  segment.color = "red",
                  segment.alpha = 0.5) +
  
  # Highlight SIRT genes with labels
  geom_point(data = subset(oligodendrocyte_df, gene_type == "SIRT"), 
             color = "blue", size = 2.5) +
  geom_text_repel(data = subset(oligodendrocyte_df, gene_type == "SIRT"), 
                  aes(label = gene), 
                  color = "blue", 
                  size = 3,
                  box.padding = 0.3,
                  segment.color = "blue",
                  segment.alpha = 0.3) +
  
  # Add significance thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgreen", alpha = 0.5) +
  
  # Customize colors and theme
  scale_color_manual(values = c("SIRT" = "blue", "other" = "gray60")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
  
  # Labels and title
  labs(title = "Differential Expression in Oligodendrocytes: 0 vs 6 stages",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change",
       y = "-log10(p-value)") +
  
  theme_minimal() +
  theme(legend.position = "right")

print(volcano_plot)



# Extract SIRT genes from the oligodendrocytes results
sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")
sirt_results <- oligodendrocyte_df[oligodendrocyte_df$gene %in% sirt_genes, ]

# Check if we found all SIRT genes
cat("SIRT genes found:", nrow(sirt_results), "\n")
print(sirt_results$gene)

# Calculate standard error from p-values
sirt_results$lfcSE <- abs(sirt_results$logFC) / qnorm(1 - sirt_results$PValue/2)

# Calculate 95% confidence intervals
sirt_results$CI_lower <- sirt_results$logFC - 1.96 * sirt_results$lfcSE
sirt_results$CI_upper <- sirt_results$logFC + 1.96 * sirt_results$lfcSE

# Create forest plot
forest_plot <- ggplot(sirt_results, aes(x = reorder(gene, logFC), y = logFC)) +
  geom_point(aes(color = FDR < 0.1), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, color = FDR < 0.1), 
                width = 0.2, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     name = "Significant (FDR < 0.1)") +
  labs(title = "SIRT expression in oligodendrocytes",
       x = "SIRT genes",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)



# 3. For astrocytes

# Get the results for astrocytes
astrocyte_results <- edgeR_results_list_0vs6[["astrocyte"]]

# Convert to data frame for easier handling
astrocyte_df <- as.data.frame(astrocyte_results)
astrocyte_df$gene <- rownames(astrocyte_df)

# Create a column to highlight SIRT genes
astrocyte_df$gene_type <- "other"
astrocyte_df$gene_type[astrocyte_df$gene %in% c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")] <- "SIRT"

# Create significance threshold column with padj <= 0.1
astrocyte_df$significant <- astrocyte_df$FDR <= 0.1

# Create the volcano plot
volcano_plot <- ggplot(astrocyte_df, aes(x = logFC, y = -log10(PValue))) +
  # Plot all points
  geom_point(aes(color = gene_type, alpha = significant), size = 1.5) +
  
  # Highlight significant points (padj <= 0.1)
  geom_point(data = subset(astrocyte_df, significant == TRUE), 
             color = "red", size = 2.5) +
  
  # Add gene labels for significant points (padj <= 0.1)
  geom_text_repel(data = subset(astrocyte_df, significant == TRUE),
                  aes(label = gene),
                  color = "red",
                  size = 3,
                  max.overlaps = 20,  # Increase if you want to see more labels
                  box.padding = 0.5,
                  segment.color = "red",
                  segment.alpha = 0.5) +
  
  # Highlight SIRT genes with labels
  geom_point(data = subset(astrocyte_df, gene_type == "SIRT"), 
             color = "blue", size = 2.5) +
  geom_text_repel(data = subset(astrocyte_df, gene_type == "SIRT"), 
                  aes(label = gene), 
                  color = "blue", 
                  size = 3,
                  box.padding = 0.3,
                  segment.color = "blue",
                  segment.alpha = 0.3) +
  
  # Add significance thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgreen", alpha = 0.5) +
  
  # Customize colors and theme
  scale_color_manual(values = c("SIRT" = "blue", "other" = "gray60")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
  
  # Labels and title
  labs(title = "Differential Expression in Astrocytes: 0 vs 6 stages",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change",
       y = "-log10(p-value)") +
  
  theme_minimal() +
  theme(legend.position = "right")

print(volcano_plot)



# Extract SIRT genes from the astrocytes results
sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")
sirt_results <- astrocyte_df[astrocyte_df$gene %in% sirt_genes, ]

# Check if we found all SIRT genes
cat("SIRT genes found:", nrow(sirt_results), "\n")
print(sirt_results$gene)

# Calculate standard error from p-values
sirt_results$lfcSE <- abs(sirt_results$logFC) / qnorm(1 - sirt_results$PValue/2)

# Calculate 95% confidence intervals
sirt_results$CI_lower <- sirt_results$logFC - 1.96 * sirt_results$lfcSE
sirt_results$CI_upper <- sirt_results$logFC + 1.96 * sirt_results$lfcSE

# Create forest plot
forest_plot <- ggplot(sirt_results, aes(x = reorder(gene, logFC), y = logFC)) +
  geom_point(aes(color = FDR < 0.1), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, color = FDR < 0.1), 
                width = 0.2, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     name = "Significant (FDR < 0.1)") +
  labs(title = "SIRT expression in astrocytes",
       x = "SIRT genes",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)



# 4. For oligodendrocytes precursor cells

# Get the results for oligodendrocytes precursor cells
oligodendrocytePrecursor_results <- edgeR_results_list_0vs6[["oligodendrocyte precursor cell"]]

# Convert to data frame for easier handling
oligodendrocytePrecursor_df <- as.data.frame(oligodendrocytePrecursor_results)
oligodendrocytePrecursor_df$gene <- rownames(oligodendrocytePrecursor_df)

# Create a column to highlight SIRT genes
oligodendrocytePrecursor_df$gene_type <- "other"
oligodendrocytePrecursor_df$gene_type[oligodendrocytePrecursor_df$gene %in% c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")] <- "SIRT"

# Create significance threshold column with padj <= 0.1
oligodendrocytePrecursor_df$significant <- oligodendrocytePrecursor_df$FDR <= 0.1

# Create the volcano plot
volcano_plot <- ggplot(oligodendrocytePrecursor_df, aes(x = logFC, y = -log10(PValue))) +
  # Plot all points
  geom_point(aes(color = gene_type, alpha = significant), size = 1.5) +
  
  # Highlight significant points (padj <= 0.1)
  geom_point(data = subset(oligodendrocytePrecursor_df, significant == TRUE), 
             color = "red", size = 2.5) +
  
  # Add gene labels for significant points (padj <= 0.1)
  geom_text_repel(data = subset(oligodendrocytePrecursor_df, significant == TRUE),
                  aes(label = gene),
                  color = "red",
                  size = 3,
                  max.overlaps = 20,  # Increase if you want to see more labels
                  box.padding = 0.5,
                  segment.color = "red",
                  segment.alpha = 0.5) +
  
  # Highlight SIRT genes with labels
  geom_point(data = subset(oligodendrocytePrecursor_df, gene_type == "SIRT"), 
             color = "blue", size = 2.5) +
  geom_text_repel(data = subset(oligodendrocytePrecursor_df, gene_type == "SIRT"), 
                  aes(label = gene), 
                  color = "blue", 
                  size = 3,
                  box.padding = 0.3,
                  segment.color = "blue",
                  segment.alpha = 0.3) +
  
  # Add significance thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgreen", alpha = 0.5) +
  
  # Customize colors and theme
  scale_color_manual(values = c("SIRT" = "blue", "other" = "gray60")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
  
  # Labels and title
  labs(title = "Differential Expression in oligodendrocytes precursor cells: 0 vs 6 stages",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change",
       y = "-log10(p-value)") +
  
  theme_minimal() +
  theme(legend.position = "right")

print(volcano_plot)



# Extract SIRT genes from the oligodendrocyte precursor cells results
sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")
sirt_results <- oligodendrocytePrecursor_df[oligodendrocytePrecursor_df$gene %in% sirt_genes, ]

# Check if we found all SIRT genes
cat("SIRT genes found:", nrow(sirt_results), "\n")
print(sirt_results$gene)

# Calculate standard error from p-values
sirt_results$lfcSE <- abs(sirt_results$logFC) / qnorm(1 - sirt_results$PValue/2)

# Calculate 95% confidence intervals
sirt_results$CI_lower <- sirt_results$logFC - 1.96 * sirt_results$lfcSE
sirt_results$CI_upper <- sirt_results$logFC + 1.96 * sirt_results$lfcSE

# Create forest plot
forest_plot <- ggplot(sirt_results, aes(x = reorder(gene, logFC), y = logFC)) +
  geom_point(aes(color = FDR < 0.1), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, color = FDR < 0.1), 
                width = 0.2, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     name = "Significant (FDR < 0.1)") +
  labs(title = "SIRT expression in oligodendrocytes precursor cells",
       x = "SIRT genes",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)



# 5. For central nervous system macrophage

# Get the results for central nervous system macrophage
CenNervSysMacr_results <- edgeR_results_list_0vs6[["central nervous system macrophage"]]

# Convert to data frame for easier handling
CenNervSysMacr_df <- as.data.frame(CenNervSysMacr_results)
CenNervSysMacr_df$gene <- rownames(CenNervSysMacr_df)

# Create a column to highlight SIRT genes
CenNervSysMacr_df$gene_type <- "other"
CenNervSysMacr_df$gene_type[CenNervSysMacr_df$gene %in% c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")] <- "SIRT"

# Create significance threshold column with padj <= 0.1
CenNervSysMacr_df$significant <- CenNervSysMacr_df$FDR <= 0.1

# Create the volcano plot
volcano_plot <- ggplot(CenNervSysMacr_df, aes(x = logFC, y = -log10(PValue))) +
  # Plot all points
  geom_point(aes(color = gene_type, alpha = significant), size = 1.5) +
  
  # Highlight significant points (padj <= 0.1)
  geom_point(data = subset(CenNervSysMacr_df, significant == TRUE), 
             color = "red", size = 2.5) +
  
  # Add gene labels for significant points (padj <= 0.1)
  geom_text_repel(data = subset(CenNervSysMacr_df, significant == TRUE),
                  aes(label = gene),
                  color = "red",
                  size = 3,
                  max.overlaps = 20,  # Increase if you want to see more labels
                  box.padding = 0.5,
                  segment.color = "red",
                  segment.alpha = 0.5) +
  
  # Highlight SIRT genes with labels
  geom_point(data = subset(CenNervSysMacr_df, gene_type == "SIRT"), 
             color = "blue", size = 2.5) +
  geom_text_repel(data = subset(CenNervSysMacr_df, gene_type == "SIRT"), 
                  aes(label = gene), 
                  color = "blue", 
                  size = 3,
                  box.padding = 0.3,
                  segment.color = "blue",
                  segment.alpha = 0.3) +
  
  # Add significance thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgreen", alpha = 0.5) +
  
  # Customize colors and theme
  scale_color_manual(values = c("SIRT" = "blue", "other" = "gray60")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
  
  # Labels and title
  labs(title = "Differential Expression in Central Nervous System Macrophages: 0 vs 6 stages",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change",
       y = "-log10(p-value)") +
  
  theme_minimal() +
  theme(legend.position = "right")

print(volcano_plot)



# Extract SIRT genes from the Central Nervous System Macrophages results
sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")
sirt_results <- CenNervSysMacr_df[CenNervSysMacr_df$gene %in% sirt_genes, ]

# Check if we found all SIRT genes
cat("SIRT genes found:", nrow(sirt_results), "\n")
print(sirt_results$gene)

# Calculate standard error from p-values
sirt_results$lfcSE <- abs(sirt_results$logFC) / qnorm(1 - sirt_results$PValue/2)

# Calculate 95% confidence intervals
sirt_results$CI_lower <- sirt_results$logFC - 1.96 * sirt_results$lfcSE
sirt_results$CI_upper <- sirt_results$logFC + 1.96 * sirt_results$lfcSE

# Create forest plot
forest_plot <- ggplot(sirt_results, aes(x = reorder(gene, logFC), y = logFC)) +
  geom_point(aes(color = FDR < 0.1), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, color = FDR < 0.1), 
                width = 0.2, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     name = "Significant (FDR < 0.1)") +
  labs(title = "SIRT expression in central nervous system macrophages",
       x = "SIRT genes",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)




# =========================== SIRT gene expression changes across cell types ====
# Create a combined forest plot

# 0 vs 2 stages 

# LFC of each SIRT per cell types
all_sirt_results <- data.frame()

for (cell_type in names(edgeR_results_list_0vs2)) {
  res_df <- edgeR_results_list_0vs2[[cell_type]]
  sirt_df <- res_df[res_df$gene %in% sirt_genes, ]
  
  if (nrow(sirt_df) > 0) {
    sirt_df$cell_type <- cell_type
    
    # Calculate standard error from p-values (CORRECT METHOD)
    sirt_df$lfcSE <- abs(sirt_df$logFC) / qnorm(1 - sirt_df$PValue/2)
    
    # Calculate 95% confidence intervals
    sirt_df$CI_lower <- sirt_df$logFC - 1.96 * sirt_df$lfcSE
    sirt_df$CI_upper <- sirt_df$logFC + 1.96 * sirt_df$lfcSE
    sirt_df$padj <- p.adjust(sirt_df$FDR, method = "BH")
    
    all_sirt_results <- rbind(all_sirt_results, sirt_df)
  }
}

# Now create your forest plot
forest_plot <- ggplot(all_sirt_results, aes(x = gene, y = logFC, color = cell_type)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(title = "SIRT expression across cell types, 0 vs 2 stages",
       x = "SIRT genes",
       y = "Log2 fold change (95% CI)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)



# LFC of SIRTs per cell type
all_sirt_results <- data.frame()

for (cell_type in names(edgeR_results_list_0vs2)) {
  res_df <- edgeR_results_list_0vs2[[cell_type]]
  sirt_df <- res_df[res_df$gene %in% sirt_genes, ]
  
  if (nrow(sirt_df) > 0) {
    sirt_df$cell_type <- cell_type
    
    # Calculate standard error from p-values (CORRECT METHOD)
    sirt_df$lfcSE <- abs(sirt_df$logFC) / qnorm(1 - sirt_df$PValue/2)
    
    # Calculate 95% confidence intervals
    sirt_df$CI_lower <- sirt_df$logFC - 1.96 * sirt_df$lfcSE
    sirt_df$CI_upper <- sirt_df$logFC + 1.96 * sirt_df$lfcSE
    sirt_df$padj <- p.adjust(sirt_df$FDR, method = "BH")
    
    all_sirt_results <- rbind(all_sirt_results, sirt_df)
  }
}

# Create a color palette for SIRT genes
sirt_colors <- c("SIRT1" = "#E41A1C", "SIRT2" = "#377EB8", "SIRT3" = "#4DAF4A",
                 "SIRT4" = "#984EA3", "SIRT5" = "#FF7F00", "SIRT6" = "#A65628",
                 "SIRT7" = "#F781BF")

# Create the forest plot with cell types on X-axis and SIRT genes in legend
forest_plot <- ggplot(all_sirt_results, aes(x = cell_type, y = logFC, color = gene)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = sirt_colors, name = "SIRT Gene") +
  labs(title = "SIRT expression across cell types, 0 vs 2 stages",
       x = "Cell types",
       y = "Log2 fold change (95% CI)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

print(forest_plot)




# Create individual forest plots for each cell type with their own legends
# and save the plots 
for (cell_type in unique(all_sirt_results$cell_type)) {
  cell_data <- all_sirt_results[all_sirt_results$cell_type == cell_type, ]
  
  p <- ggplot(cell_data, aes(x = gene, y = logFC, color = gene)) +
    geom_point(size = 1.8) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                  width = 0.08, size = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
    scale_color_manual(values = sirt_colors, name = NULL) +
    scale_y_continuous(limits = c(-2, 2)) +
    scale_x_discrete(expand = expansion(mult = 0.1)) +
    labs(title = cell_type, y = "Log2FC") +
    theme_classic() +  # Classic theme has white background
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 6),
      axis.title.y = element_text(size = 7, face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
      legend.position = "right",
      legend.text = element_text(size = 5.5),
      legend.key.size = unit(0.25, "cm"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add black frame
      plot.margin = margin(1, 1, 1, 1)
    )
  
  print(p)
  ggsave(paste0("narrow_sirt_", gsub(" ", "_", tolower(cell_type)), ".png"), 
         p, width = 2, height = 3, dpi = 300, bg = "white")
}




# 0 vs 6 stages 

# LFC of each SIRT per cell types
all_sirt_results <- data.frame()

for (cell_type in names(edgeR_results_list_0vs6)) {
  res_df <- edgeR_results_list_0vs6[[cell_type]]
  sirt_df <- res_df[res_df$gene %in% sirt_genes, ]
  
  if (nrow(sirt_df) > 0) {
    sirt_df$cell_type <- cell_type
    
    # Calculate standard error from p-values (CORRECT METHOD)
    sirt_df$lfcSE <- abs(sirt_df$logFC) / qnorm(1 - sirt_df$PValue/2)
    
    # Calculate 95% confidence intervals
    sirt_df$CI_lower <- sirt_df$logFC - 1.96 * sirt_df$lfcSE
    sirt_df$CI_upper <- sirt_df$logFC + 1.96 * sirt_df$lfcSE
    sirt_df$padj <- p.adjust(sirt_df$FDR, method = "BH")
    
    all_sirt_results <- rbind(all_sirt_results, sirt_df)
  }
}

# Now create your forest plot
forest_plot <- ggplot(all_sirt_results, aes(x = gene, y = logFC, color = cell_type)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(title = "SIRT expression across cell types, 0 vs 6 stages",
       x = "SIRT genes",
       y = "Log2 fold change (95% CI)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)



# LFC of SIRTs per cell type
all_sirt_results <- data.frame()

for (cell_type in names(edgeR_results_list_0vs6)) {
  res_df <- edgeR_results_list_0vs6[[cell_type]]
  sirt_df <- res_df[res_df$gene %in% sirt_genes, ]
  
  if (nrow(sirt_df) > 0) {
    sirt_df$cell_type <- cell_type
    
    # Calculate standard error from p-values (CORRECT METHOD)
    sirt_df$lfcSE <- abs(sirt_df$logFC) / qnorm(1 - sirt_df$PValue/2)
    
    # Calculate 95% confidence intervals
    sirt_df$CI_lower <- sirt_df$logFC - 1.96 * sirt_df$lfcSE
    sirt_df$CI_upper <- sirt_df$logFC + 1.96 * sirt_df$lfcSE
    sirt_df$padj <- p.adjust(sirt_df$FDR, method = "BH")
    
    all_sirt_results <- rbind(all_sirt_results, sirt_df)
  }
}

# Create a color palette for SIRT genes
sirt_colors <- c("SIRT1" = "#E41A1C", "SIRT2" = "#377EB8", "SIRT3" = "#4DAF4A",
                 "SIRT4" = "#984EA3", "SIRT5" = "#FF7F00", "SIRT6" = "#A65628",
                 "SIRT7" = "#F781BF")

# Create the forest plot with cell types on X-axis and SIRT genes in legend
forest_plot <- ggplot(all_sirt_results, aes(x = cell_type, y = logFC, color = gene)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = sirt_colors, name = "SIRT Gene") +
  labs(title = "SIRT expression across cell types, 0 vs 6 stages",
       x = "Cell types",
       y = "Log2 fold change (95% CI)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

print(forest_plot)




# Create individual forest plots for each cell type with their own legends
# and save the plots 
for (cell_type in unique(all_sirt_results$cell_type)) {
  cell_data <- all_sirt_results[all_sirt_results$cell_type == cell_type, ]
  
  p <- ggplot(cell_data, aes(x = gene, y = logFC, color = gene)) +
    geom_point(size = 1.8) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                  width = 0.08, size = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
    scale_color_manual(values = sirt_colors, name = NULL) +
    scale_y_continuous(limits = c(-2, 2)) +
    scale_x_discrete(expand = expansion(mult = 0.1)) +
    labs(title = cell_type, y = "Log2FC") +
    theme_classic() +  # Classic theme has white background
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 6),
      axis.title.y = element_text(size = 7, face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
      legend.position = "right",
      legend.text = element_text(size = 5.5),
      legend.key.size = unit(0.25, "cm"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add black frame
      plot.margin = margin(1, 1, 1, 1)
    )
  
  print(p)
  ggsave(paste0("narrow_sirt_", gsub(" ", "_", tolower(cell_type)), ".png"), 
         p, width = 2, height = 3, dpi = 300, bg = "white")
}
