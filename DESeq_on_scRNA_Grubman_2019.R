# =================== Check the content of SingleCellExperiment object ====

library(SingleCellExperiment)
sce <- readRDS("~/short_term_project/Grubman_2019_sce_ann.rds")

# Check class
class(sce)

# Quick summary
sce

# Available assays (usually "counts", "logcounts", maybe "normcounts")
assayNames(sce)

# Cell metadata (per-cell annotations)
colnames(colData(sce))
df <- head(as.data.frame(colData(sce)))   # peek at first 6 rows by method "head" 

# Gene metadata
colnames(rowData(sce))
head(rowData(sce))

# Number of cells per cell type
table(colData(sce)$celltype) # table() - counts how many times each unique value appears

# Check unique values in "condition" column of colData (cell metadata)
unique(colData(sce)$sample_status)

# If "condition" (AD vs Control) is present
if("sample_status" %in% colnames(colData(sce))) {
  table(colData(sce)$celltype, colData(sce)$sample_status) # rows - cell types, columns - sample_status
}

# Sample info (sometimes stored as "sample", "donor", or "batch")
colnames(colData(sce))
unique(colData(sce)$donor_ID) # 6 donors

# Check UMAP (stores dimensionaly reduction coordinates)
head(reducedDim(sce))
colnames(reducedDim(sce)) # UMAP1 - X-axis in UMAP plot, UMAP2 - Y-axis in UMAP plot



# ====================== SIRT expression in different cell types (UMAP) ====

# Load necessary libraries 
library(ggplot2)
library(patchwork) # for combining plots
library(viridis) # for better color scales 

# Get the UMAP coordinates from your SingleCellExperiment object
umap_coords <- reducedDim(sce, 'UMAPHMNY')
colnames(umap_coords) <- colnames(reducedDim(sce, 'UMAPHMNY'))
umap_df <- as.data.frame(umap_coords)

# Add cell metadata
umap_df$celltype <- colData(sce)$celltype
umap_df$condition <- colData(sce)$sample_status

# Firstly, create a reference UMAP colored by cell type
celltype_umap <- ggplot(umap_df, aes(x = UMAPHMNY_1, y = UMAPHMNY_2)) +
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

print(celltype_umap)


sample_status <- colData(sce)$sample_status

# Secondly, create a reference UMAP colored by condition
condition_umap <- ggplot(umap_df, aes(x = UMAPHMNY_1, y = UMAPHMNY_2)) +
  geom_point(aes(color = sample_status), size = 0.1, alpha = 1) + # alpha - semi-transparent (0.6)
  scale_color_brewer(palette = "RdYlBu", name = "Condition") + # a color-blind friendly palette, name the legend "Cell Type"
  labs(title = "Condition clusters") + # Adds the plot title
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

print(condition_umap)









# ====================== SIRT expression in different cell types (histograms) ====

# Check what SIRT genes are actually in your dataset
available_genes <- rownames(assay(sce, 'logcounts'))
sirt_genes_in_data <- available_genes[grepl("SIRT", available_genes, ignore.case = TRUE)]

cat("SIRT genes found in dataset:\n")
print(sirt_genes_in_data)
# "SIRT5" "SIRT3" "SIRT7" "SIRT6" "SIRT2"

# Histograms (Kate's code)
library(dplyr)
library(tidyr)

sirt_genes <- c("SIRT2", "SIRT3", "SIRT5", "SIRT6", "SIRT7")

available_sirt_genes <- sirt_genes[sirt_genes %in% rownames(sce)]

expr_matrix <- assay(sce, 'logcounts')

sirt_expr <- as.matrix(expr_matrix[available_sirt_genes, ])

metadata <- data.frame(
  CellID = colnames(sce),
  CellType = sce$celltype,        
  Condition = sce$sample_status  
)

plot_data <- as.data.frame(t(sirt_expr)) # транспонирование
plot_data$CellID <- rownames(plot_data)


# Join first
plot_data <- plot_data %>%
  left_join(metadata, by = "CellID")

# Then pivot
plot_data <- plot_data %>%
  pivot_longer(
    cols = all_of(available_sirt_genes),
    names_to = "Gene",
    values_to = "Expression"
  )

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression))

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression)) %>%
  mutate(Condition = factor(Condition, levels = c("Healthy", "Alzheimer's disease")))

specific_celltypes <- c("neuron", 
                        "oligodendrocyte", 
                        "oligodendrocyte precursor cell", 
                        "astrocyte")  

# Plot histograms for each SIRT across 5 cell types (dividing into AD and Control)
p <- plot_data %>%
  filter(CellType %in% specific_celltypes) %>%
  ggplot(aes(x = Expression)) +
  geom_density(aes(y = after_stat(scaled), color = CellType), 
               linewidth = 0.8, alpha = 0.8) +
  facet_grid(Gene ~ Condition, scales = "free") +
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



# Plot histohrams for each cell type across SIRT1-7 (dividing into AD and Control)
celltype_plots <- list()

for (cell_type in specific_celltypes) {
  p <- plot_data %>%
    filter(CellType == cell_type) %>%
    ggplot(aes(x = Expression)) +
    geom_density(aes(y = after_stat(scaled), color = Gene), 
                 linewidth = 0.8, alpha = 0.8) +
    facet_grid(. ~ Condition, scales = "free") +
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



# My old histograms 
# Define threshold for low expressed genes by checking the distribution of expression values for one SIRT gene
# expr_values <- assay(sce, "logcounts")["SIRT1", ]
# summary(expr_values) # key statistics 
# hist(expr_values[expr_values > 0], breaks = 50, main = "SIRT1 Expression Distribution")
expr_values <- assay(sce, "logcounts")["SIRT2", ]
summary(expr_values) # key statistics 
hist(expr_values[expr_values > 0], breaks = 50, main = "SIRT2 Expression Distribution")

expr_values <- assay(sce, "logcounts")["SIRT3", ]
summary(expr_values) # key statistics 
hist(expr_values[expr_values > 0], breaks = 50, main = "SIRT3 Expression Distribution")

# expr_values <- assay(sce, "logcounts")["SIRT4", ]
# summary(expr_values) # key statistics 
# hist(expr_values[expr_values > 0], breaks = 50, main = "SIRT4 Expression Distribution")

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
sirt_genes <- c("SIRT2", "SIRT3", "SIRT5", "SIRT6", "SIRT7")
sirt_umap_plots <- list()

for (sirt_gene in sirt_genes) {
  if (sirt_gene %in% rownames(logcounts_matrix)) {
    # Add expression for this SIRT gene
    umap_df$expression <- logcounts_matrix[sirt_gene, ]
    
    # Define threshold for low expression
    low_expr_threshold <- 0.1
    
    # Create the UMAP plot
    p <- ggplot(umap_df, aes(x = UMAPHMNY_1, y = UMAPHMNY_2)) +
      geom_point(aes(color = ifelse(expression < low_expr_threshold, NA, expression)), # Colors points based on expression, but sets low-expression cells to NA so they get the na.value color
                 size = 0.1, alpha = 0.6) +
      scale_color_viridis_c( # Applies a color scale
        name = "Expression",
        option = "plasma", # palette name 
        na.value = "grey80",  # Color for values below threshold
        limits = c(low_expr_threshold, max(umap_df$expression, na.rm = TRUE)) # limits = c(low_expr_threshold, max(...)) - 
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
    
    sirt_umap_plots[[sirt_gene]] <- p # Stores each plot in a list for later use
  } else {
    cat("SIRT gene", sirt_gene, "not found in the dataset\n")
  }
}

# Arrange all plots in a grid
combined_umap <- wrap_plots(sirt_umap_plots, ncol = 3) +
  plot_annotation(title = "SIRT gene expression across cell types (UMAPHMNY)")

print(combined_umap)





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
      expressing_cells <- expr_values
      if (length(expressing_cells) > 0) {
        expr_matrix[cell_type, sirt_gene] <- mean(expressing_cells, na.rm = TRUE)
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
cell_type_order <- c("vascular associated smooth muscle cell", "pericyte", "central nervous system macrophage", 
                     "fibroblast", "neuron", "astrocyte", "endothelial cell", "oligodendrocyte precursor cell", 
                     "oligodendrocyte", "leukocyte")

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








# =========================== I. Pseudobulk DE analysis (using DESeq2) ====

# Load necessary library
library(DESeq2)

# Define the key columns based on your exploration
donor_column <- "donor_ID" 
condition_column <- "sample" # AD and ct (control)

# Define the cell types we want to analyze (focusing on the major ones)
cell_types_to_analyze <- c("neuron", "oligodendrocyte", "astrocyte", 
                           "oligodendrocyte precursor cell", 
                           "central nervous system macrophage")

# Create an empty list to store results for each cell type
deseq_results_list <- list()

# Loop through each cell type of interest
for (cell_type_of_interest in cell_types_to_analyze) {
  
  cat("\n\n", rep("=", 50), "\n")
  cat("Processing cell type:", cell_type_of_interest, "\n")
  cat(rep("=", 50), "\n\n")
  
  # A. Subset the SCE object to the specific cell type
  sce_subset <- sce[, sce$celltype == cell_type_of_interest]
  cat("Number of cells in subset:", ncol(sce_subset), "\n")
  
  # B. Check if we have enough cells and donors for analysis
  donor_table <- table(colData(sce_subset)[[donor_column]])
  cat("Number of unique donors:", length(donor_table), "\n")
  
  condition_table <- table(colData(sce_subset)[[condition_column]])
  cat("Cells per condition:\n")
  print(condition_table)
  
  # C. Create Pseudobulk Matrix using base R
  # Get the raw counts matrix
  count_matrix <- assay(sce_subset, "counts")
  
  # Get the donor IDs for each cell
  donor_ids <- colData(sce_subset)[[donor_column]]
  
  # Get unique donors
  unique_donors <- unique(donor_ids)
  
  # Initialize an empty matrix for pseudobulk counts
  pb_matrix <- matrix(0, nrow = nrow(count_matrix), ncol = length(unique_donors))
  rownames(pb_matrix) <- rownames(count_matrix)
  colnames(pb_matrix) <- unique_donors
  
  # Sum counts for each donor
  for (donor in unique_donors) {
    donor_cells <- which(donor_ids == donor)
    if (length(donor_cells) > 1) {
      pb_matrix[, donor] <- rowSums(count_matrix[, donor_cells, drop = FALSE])
    } else if (length(donor_cells) == 1) {
      pb_matrix[, donor] <- count_matrix[, donor_cells]
    }
  }
  
  cat("Pseudobulk matrix dimensions:", dim(pb_matrix), "(genes x donors)\n")
  
  # D. Prepare the metadata for DESeq2
  metadata <- unique(colData(sce_subset)[, c(donor_column, condition_column)])
  metadata <- metadata[match(colnames(pb_matrix), metadata[[donor_column]]), ]
  rownames(metadata) <- metadata[[donor_column]]
  
  # Ensure the condition is a factor and set the reference level to "ct" (Control)
  metadata[[condition_column]] <- factor(metadata[[condition_column]], levels = c("ct", "AD"))
  print("Metadata for this cell type:")
  print(metadata)
  
  # E. Create the DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(
    countData = pb_matrix,
    colData = metadata,
    design = as.formula(paste("~", condition_column)) # Design: ~ condition
  )
  
  # F. Pre-filter lowly expressed genes
  # Keep genes with a count of at least 10 in at least 3 donors
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  cat("Number of genes after filtering:", nrow(dds), "\n")
  
  # G. Run DESeq2 differential expression analysis
  dds <- DESeq(dds)
  
  # H. Extract results. contrast: 'AD' vs base level 'ct'
  res <- results(dds, contrast = c(condition_column, "AD", "ct"), alpha = 0.05)
  res_sorted <- res[order(res$padj), ] # Sort by adjusted p-value
  
  # I. Save the full results for this cell type
  deseq_results_list[[cell_type_of_interest]] <- res_sorted
  
  # J. Print summary and SIRT results
  cat("\nDESeq2 results summary:\n")
  print(summary(res))
  
  sirt_genes <- c("SIRT2", "SIRT3", "SIRT5", "SIRT6", "SIRT7")
  sirt_results <- res_sorted[rownames(res_sorted) %in% sirt_genes, ]
  
  cat("\nResults for SIRT genes:\n")
  if (nrow(sirt_results) > 0) {
    print(sirt_results)
  } else {
    cat("No SIRT genes found in the results for this cell type (may have been filtered out due to low expression).\n")
  }
  
  # Optional: Save the results for this cell type to a CSV file
  # write.csv(as.data.frame(res_sorted), file = paste0("deseq_results_", gsub(" ", "_", cell_type_of_interest), ".csv"))
}

# Save the entire list of results for later inspection
# saveRDS(deseq_results_list, file = "deseq_results_list_all_celltypes.rds")

# cat("\n\nAnalysis complete! Results saved to CSV files and 'deseq_results_list_all_celltypes.rds'\n")





# =========================== DESeq results visualization ====

# Load necessary library for plotting
library(ggplot2)
# install.packages("ggrepel")
library(ggrepel)  # For better label placement

# 1. For neuron cells =====

# Get the results for neurons
neuron_results <- deseq_results_list[["neuron"]]

# Convert to data frame for easier handling
neuron_df <- as.data.frame(neuron_results)
neuron_df$gene <- rownames(neuron_df)

# Create a column to highlight SIRT genes
neuron_df$gene_type <- "other"
neuron_df$gene_type[neuron_df$gene %in% c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")] <- "SIRT"

# Create significance threshold column with padj <= 0.1
neuron_df$significant <- neuron_df$padj <= 0.1

# Create the volcano plot
volcano_plot <- ggplot(neuron_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  # Plot all points
  geom_point(aes(color = gene_type), size = 1.5, alpha = 0.6) +
  
  # Highlight significant points (padj <= 0.1)
  geom_point(data = subset(neuron_df, significant == TRUE), 
             color = "red", size = 2.5) +
  
  # Add gene labels for significant points (padj <= 0.1)
  geom_text_repel(data = subset(neuron_df, significant == TRUE),
                  aes(label = gene),
                  color = "red",
                  size = 3,
                  max.overlaps = 20,
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
  
  # Customize colors - only show gene_type in legend
  scale_color_manual(values = c("SIRT" = "blue", "other" = "gray60"),
                     name = "Gene Type") +
  
  # Labels and title
  labs(title = "Differential Expression in Neurons: AD vs Control",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change (AD vs Control)",
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

# Create the forest plot
sirt_forest <- ggplot(sirt_results, aes(x = log2FoldChange, y = reorder(gene, log2FoldChange))) +
  # Add vertical line at 0 (no effect)
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  
  # Add points for the effect sizes
  geom_point(aes(color = padj < 0.1), size = 3) +
  
  # Add error bars (95% confidence intervals)
  geom_errorbarh(aes(xmin = log2FoldChange - 1.96 * lfcSE, 
                     xmax = log2FoldChange + 1.96 * lfcSE,
                     color = padj < 0.1), 
                 height = 0.2) +
  
  # Custom colors for significance
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     labels = c("TRUE" = "padj < 0.1", "FALSE" = "padj ≥ 0.1"),
                     name = "Significance") +
  
  # Labels and theme
  labs(title = "SIRT gene expression changes in Neurons",
       subtitle = "AD vs Control with 95% confidence intervals",
       x = "Log2 Fold Change with 95% CI",
       y = "SIRT Genes") +
  
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

print(sirt_forest)




# 2. For oligodendrocytes =====

# Get the results for oligodendrocytes
oligodendrocyte_results <- deseq_results_list[["oligodendrocyte"]]

# Convert to data frame for easier handling
oligodendrocyte_df <- as.data.frame(oligodendrocyte_results)
oligodendrocyte_df$gene <- rownames(oligodendrocyte_df)

# Create a column to highlight SIRT genes
oligodendrocyte_df$gene_type <- "other"
oligodendrocyte_df$gene_type[oligodendrocyte_df$gene %in% c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")] <- "SIRT"

# Create significance threshold column with padj <= 0.1
oligodendrocyte_df$significant <- oligodendrocyte_df$padj <= 0.1

# Create the volcano plot
volcano_plot <- ggplot(oligodendrocyte_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  # Plot all points
  geom_point(aes(color = gene_type), size = 1.5, alpha = 0.6) +
  
  # Highlight significant points (padj <= 0.1)
  geom_point(data = subset(oligodendrocyte_df, significant == TRUE), 
             color = "red", size = 2.5) +
  
  # Add gene labels for significant points (padj <= 0.1)
  geom_text_repel(data = subset(oligodendrocyte_df, significant == TRUE),
                  aes(label = gene),
                  color = "red",
                  size = 3,
                  max.overlaps = 20,
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
  
  # Customize colors - only show gene_type in legend
  scale_color_manual(values = c("SIRT" = "blue", "other" = "gray60"),
                     name = "Gene Type") +
  
  # Labels and title
  labs(title = "Differential expression in oligodendrocytes: AD vs Control",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change (AD vs Control)",
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

# If some SIRT genes are missing, they might have been filtered out during DESeq2 analysis
# missing_sirts <- setdiff(sirt_genes, sirt_results$gene)
# if (length(missing_sirts) > 0) {
# cat("Missing SIRT genes (filtered out due to low expression):", missing_sirts, "\n")
# }

# Create the forest plot
sirt_forest <- ggplot(sirt_results, aes(x = log2FoldChange, y = reorder(gene, log2FoldChange))) +
  # Add vertical line at 0 (no effect)
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  
  # Add points for the effect sizes
  geom_point(aes(color = padj < 0.1), size = 3) +
  
  # Add error bars (95% confidence intervals)
  geom_errorbarh(aes(xmin = log2FoldChange - 1.96 * lfcSE, 
                     xmax = log2FoldChange + 1.96 * lfcSE,
                     color = padj < 0.1), 
                 height = 0.2) +
  
  # Custom colors for significance
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     labels = c("TRUE" = "padj < 0.1", "FALSE" = "padj ≥ 0.1"),
                     name = "Significance") +
  
  # Labels and theme
  labs(title = "SIRT gene expression changes in Oligodendrocytes",
       subtitle = "AD vs Control with 95% confidence intervals",
       x = "Log2 Fold Change with 95% CI",
       y = "SIRT Genes") +
  
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

print(sirt_forest)




# 3. For astrocytes =====

# Get the results for astrocytes
astrocyte_results <- deseq_results_list[["astrocyte"]]

# Convert to data frame for easier handling
astrocyte_df <- as.data.frame(astrocyte_results)
astrocyte_df$gene <- rownames(astrocyte_df)

# Create a column to highlight SIRT genes
astrocyte_df$gene_type <- "other"
astrocyte_df$gene_type[astrocyte_df$gene %in% c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")] <- "SIRT"

# Create significance threshold column with padj <= 0.1
astrocyte_df$significant <- astrocyte_df$padj <= 0.1

# Create the volcano plot
volcano_plot <- ggplot(astrocyte_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  # Plot all points
  geom_point(aes(color = gene_type), size = 1.5, alpha = 0.6) +
  
  # Highlight significant points (padj <= 0.1)
  geom_point(data = subset(astrocyte_df, significant == TRUE), 
             color = "red", size = 2.5) +
  
  # Add gene labels for significant points (padj <= 0.1)
  geom_text_repel(data = subset(astrocyte_df, significant == TRUE),
                  aes(label = gene),
                  color = "red",
                  size = 3,
                  max.overlaps = 20,
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
  
  # Customize colors - only show gene_type in legend
  scale_color_manual(values = c("SIRT" = "blue", "other" = "gray60"),
                     name = "Gene Type") +
  
  # Labels and title
  labs(title = "Differential expression in astrocytes: AD vs Control",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change (AD vs Control)",
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

# If some SIRT genes are missing, they might have been filtered out during DESeq2 analysis
# missing_sirts <- setdiff(sirt_genes, sirt_results$gene)
# if (length(missing_sirts) > 0) {
# cat("Missing SIRT genes (filtered out due to low expression):", missing_sirts, "\n")
# }

# Create the forest plot
sirt_forest <- ggplot(sirt_results, aes(x = log2FoldChange, y = reorder(gene, log2FoldChange))) +
  # Add vertical line at 0 (no effect)
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  
  # Add points for the effect sizes
  geom_point(aes(color = padj < 0.1), size = 3) +
  
  # Add error bars (95% confidence intervals)
  geom_errorbarh(aes(xmin = log2FoldChange - 1.96 * lfcSE, 
                     xmax = log2FoldChange + 1.96 * lfcSE,
                     color = padj < 0.1), 
                 height = 0.2) +
  
  # Custom colors for significance
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     labels = c("TRUE" = "padj < 0.1", "FALSE" = "padj ≥ 0.1"),
                     name = "Significance") +
  
  # Labels and theme
  labs(title = "SIRT gene expression changes in Astrocytes",
       subtitle = "AD vs Control with 95% confidence intervals",
       x = "Log2 Fold Change with 95% CI",
       y = "SIRT Genes") +
  
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

print(sirt_forest)




# 4. For oligodendrocytes precursor cells =====

# Get the results for oligodendrocytes precursor cells
oligodendrocytePrecursor_results <- deseq_results_list[["oligodendrocyte precursor cell"]]

# Convert to data frame for easier handling
oligodendrocytePrecursor_df <- as.data.frame(oligodendrocytePrecursor_results)
oligodendrocytePrecursor_df$gene <- rownames(oligodendrocytePrecursor_df)

# Create a column to highlight SIRT genes
oligodendrocytePrecursor_df$gene_type <- "other"
oligodendrocytePrecursor_df$gene_type[oligodendrocytePrecursor_df$gene %in% c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")] <- "SIRT"

# Create significance threshold column with padj <= 0.1
oligodendrocytePrecursor_df$significant <- oligodendrocytePrecursor_df$padj <= 0.1

# Create the volcano plot
volcano_plot <- ggplot(oligodendrocytePrecursor_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  # Plot all points
  geom_point(aes(color = gene_type), size = 1.5, alpha = 0.6) +
  
  # Highlight significant points (padj <= 0.1)
  geom_point(data = subset(oligodendrocytePrecursor_df, significant == TRUE), 
             color = "red", size = 2.5) +
  
  # Add gene labels for significant points (padj <= 0.1)
  geom_text_repel(data = subset(oligodendrocytePrecursor_df, significant == TRUE),
                  aes(label = gene),
                  color = "red",
                  size = 3,
                  max.overlaps = 20,
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
  
  # Customize colors - only show gene_type in legend
  scale_color_manual(values = c("SIRT" = "blue", "other" = "gray60"),
                     name = "Gene Type") +
  
  # Labels and title
  labs(title = "Differential expression in oligodendrocytes precursor cells: AD vs Control",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change (AD vs Control)",
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

# If some SIRT genes are missing, they might have been filtered out during DESeq2 analysis
# missing_sirts <- setdiff(sirt_genes, sirt_results$gene)
# if (length(missing_sirts) > 0) {
# cat("Missing SIRT genes (filtered out due to low expression):", missing_sirts, "\n")
# }

# Create the forest plot
sirt_forest <- ggplot(sirt_results, aes(x = log2FoldChange, y = reorder(gene, log2FoldChange))) +
  # Add vertical line at 0 (no effect)
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  
  # Add points for the effect sizes
  geom_point(aes(color = padj < 0.1), size = 3) +
  
  # Add error bars (95% confidence intervals)
  geom_errorbarh(aes(xmin = log2FoldChange - 1.96 * lfcSE, 
                     xmax = log2FoldChange + 1.96 * lfcSE,
                     color = padj < 0.1), 
                 height = 0.2) +
  
  # Custom colors for significance
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     labels = c("TRUE" = "padj < 0.1", "FALSE" = "padj ≥ 0.1"),
                     name = "Significance") +
  
  # Labels and theme
  labs(title = "SIRT gene expression changes in Oligodendrocyte precursor cells",
       subtitle = "AD vs Control with 95% confidence intervals",
       x = "Log2 Fold Change with 95% CI",
       y = "SIRT Genes") +
  
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

print(sirt_forest)








# =========================== SIRT gene expression changes across cell types ====

# Create a combined data frame with all cell types
combined_sirt_data <- data.frame()

# Loop through each cell type and extract SIRT results
for (cell_type in names(deseq_results_list)) {
  # Get results for this cell type
  cell_results <- as.data.frame(deseq_results_list[[cell_type]])
  cell_results$gene <- rownames(cell_results)
  
  # Extract SIRT genes
  cell_sirt <- cell_results[cell_results$gene %in% c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7"), ]
  
  # Add cell type information
  if (nrow(cell_sirt) > 0) {
    cell_sirt$cell_type <- cell_type
    combined_sirt_data <- rbind(combined_sirt_data, cell_sirt)
  }
}

# Create a color palette for SIRT genes
sirt_colors <- c("SIRT1" = "#E41A1C", "SIRT2" = "#377EB8", "SIRT3" = "#4DAF4A",
                 "SIRT4" = "#984EA3", "SIRT5" = "#FF7F00", "SIRT6" = "#A65628",
                 "SIRT7" = "#F781BF")

# Create the combined forest plot
combined_forest <- ggplot(combined_sirt_data, aes(x = log2FoldChange, y = reorder(gene, log2FoldChange))) +
  # Add vertical line at 0 (no effect)
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  
  # Add points and error bars colored by SIRT gene
  geom_point(aes(color = gene), size = 2) +
  geom_errorbarh(aes(xmin = log2FoldChange - 1.96 * lfcSE, 
                     xmax = log2FoldChange + 1.96 * lfcSE,
                     color = gene), 
                 height = 0.15, alpha = 0.7) +
  
  # Facet by cell type
  facet_grid(. ~ cell_type, scales = "free_x") +
  
  # Use the SIRT color palette
  scale_color_manual(values = sirt_colors, name = "SIRT Gene") +
  
  # Labels and theme
  labs(title = "SIRT Gene Expression Changes Across Cell Types",
       subtitle = "AD vs Control with 95% Confidence Intervals",
       x = "Log2 Fold Change with 95% CI",
       y = "SIRT Genes") +
  
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_rect(fill = "lightgray"),
        axis.text.x = element_text(angle = 45, hjust = 1))

print(combined_forest)







# =========================== II. Pseudobulk DE analysis (using edgeR) ====

library(edgeR)
library(limma)
library(ggplot2)
library(dplyr)

colnames(colData(sce))


# Define key columns
donor_column <- "sample_ID"
condition_column <- "sample_status"
celltype_column <- "celltype"

# Get SIRT genes
sirt_genes <- c("SIRT2", "SIRT3", "SIRT5", "SIRT6", "SIRT7")

# Create empty list to store results
edgeR_results_list <- list()

# Define cell types to analyze (focus on major ones)
cell_types_to_analyze <- c("neuron", "oligodendrocyte", "astrocyte", 
                           "oligodendrocyte precursor cell", 
                           "central nervous system macrophage")

# Pseudobulk creation and edgeR analysis
for (cell_type_of_interest in cell_types_to_analyze) {
  
  cat("\n\n=== Processing cell type:", cell_type_of_interest, "===\n")
  
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
  
  # Sum counts for each donor to create a pseudobulk matrix
  for (donor in unique_donors) {
    donor_cells <- which(donor_ids == donor)
    if (length(donor_cells) > 0) {
      pb_matrix[, donor] <- rowSums(count_matrix[, donor_cells, drop = FALSE])
    }
  }
  
  # C. Prepare metadata
  metadata <- unique(colData(sce_subset)[, c(donor_column, condition_column)])
  metadata <- metadata[match(colnames(pb_matrix), metadata[[donor_column]]), ]
  rownames(metadata) <- metadata[[donor_column]]
  
  # Ensure condition is factor with correct reference level
  metadata[[condition_column]] <- factor(metadata[[condition_column]], levels = c("Healthy", "Alzheimer's disease"))
  print("Metadata for this cell type:")
  print(metadata)
  
  # D. edgeR Analysis
  # Create DGEList object
  dge <- DGEList(counts = pb_matrix, group = metadata[[condition_column]])
  
  # Pre-filter lowly expressed genes (keep genes with at least 10 counts in at least 3 samples)
  keep <- filterByExpr(dge, group = metadata[[condition_column]])
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  cat("Number of genes after filtering:", nrow(dge), "\n")
  
  # Normalize (TMM normalization)
  dge <- calcNormFactors(dge)
  
  # Create design matrix
  design <- model.matrix(~ metadata[[condition_column]])
  colnames(design) <- gsub("metadata\\[\\[condition_column\\]\\]", "", colnames(design))
  
  # Estimate dispersions
  dge <- estimateDisp(dge, design)
  
  # Fit model and test
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = 2)  # Test the condition coefficient
  
  # Get results
  res <- topTags(qlf, n = Inf, sort.by = "PValue")
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Store results
  edgeR_results_list[[cell_type_of_interest]] <- res_df
  
  # Print SIRT results
  sirt_results <- res_df[res_df$gene %in% sirt_genes, ]
  cat("\nSIRT gene results for", cell_type_of_interest, ":\n")
  print(sirt_results)
}

# Description of the result edgeR table:
# 1. logFC - log2 fold change(how much expression changes in AD vs NC)
# 2. logCPM - log2 counts per million (average expression level)
# 3. F - F-statistic (effect size relative to variability)
# 4. PValue - raw p-value (probability of seeing this effect by chance if no real difference)
# 5. FDR - false discovery rate (adjusted p-value for multiple testing)



# =========================== edgeR results visualisation ====
# Load necessary library for plotting
library(ggplot2)
# install.packages("ggrepel")
library(ggrepel)  # For better label placement

# 1. For neuron cells =====

# Get the results for neurons
neuron_results <- edgeR_results_list[["neuron"]]

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
  labs(title = "Differential Expression in Neurons: AD vs Control",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change (AD vs Control)",
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



# 2. For oligodendrocytes =====

# Get the results for oligodendrocytes
oligodendrocyte_results <- edgeR_results_list[["oligodendrocyte"]]

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
  labs(title = "Differential Expression in Oligodendrocytes: AD vs Control",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change (AD vs Control)",
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






# 3. For astrocytes =====

# Get the results for astrocytes
astrocyte_results <- edgeR_results_list[["astrocyte"]]

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
  labs(title = "Differential Expression in Astrocytes: AD vs Control",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change (AD vs Control)",
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





# 4. For oligodendrocytes precursor cells =====

# Get the results for oligodendrocytes precursor cells
oligodendrocytePrecursor_results <- edgeR_results_list[["oligodendrocyte precursor cell"]]

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
  labs(title = "Differential Expression in oligodendrocytes precursor cells: AD vs Control",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change (AD vs Control)",
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




# 5. For central nervous system macrophage =====

# Get the results for central nervous system macrophage
CenNervSysMacr_results <- edgeR_results_list[["central nervous system macrophage"]]

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
  labs(title = "Differential Expression in Central Nervous System Macrophages: AD vs Control",
       subtitle = "Blue points: SIRT genes | Red points: padj ≤ 0.1",
       x = "log2 Fold Change (AD vs Control)",
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

# LFC of each SIRT per cell types
all_sirt_results <- data.frame()

for (cell_type in names(edgeR_results_list)) {
  res_df <- edgeR_results_list[[cell_type]]
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
  labs(title = "SIRT expression across cell types",
       x = "SIRT genes",
       y = "Log2 fold change (95% CI)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(forest_plot)



# LFC of SIRTs per cell type
all_sirt_results <- data.frame()

for (cell_type in names(edgeR_results_list)) {
  res_df <- edgeR_results_list[[cell_type]]
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
  labs(title = "SIRT expression across cell types",
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
  # ggsave(paste0("narrow_sirt_", gsub(" ", "_", tolower(cell_type)), ".png"), 
         # p, width = 2, height = 3, dpi = 300, bg = "white")
}


