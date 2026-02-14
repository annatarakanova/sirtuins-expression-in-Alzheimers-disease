#### Analysis of scRNA-seq dataset from middle temporal gyrus (Alzheimer's disease vs Control) ####

# Download lemur package from Bioconductor 
# BiocManager::install("lemur")
library(lemur)

# Read SingleCellExperiment file 
library(SingleCellExperiment)
sce <- readRDS("~/short_term_project/GSE188545_sce_ann.rds")

# explore sce object
sce

colnames(colData(sce))

unique(colData(sce)$condition)

unique(colData(sce)$sample_id)

condition <- colData(sce)$condition
sample_id <- colData(sce)$sample_id

# fit basic LEMUR model (learn the model)
fit <- lemur(sce, 
             design = ~ condition, # design matrix 
             n_emb = 15,  # Start with 30 latent dimensions
             test_fraction = 0.5)  # Hold out 50% of cells for validation

fit

# Check the description of the package and possible functions 
library(help = "lemur") # description 
ls("package:lemur") # functions 

# align corresponding cells using sample_id 
fit <- align_by_grouping(fit, grouping = sample_id) # This ensures that corresponding cells are close to each other in the fit$embedding

# Samples are correlated with conditions => we don't put both sample_id and condition
# in design matrix. Use 'align_by_grouping' function to remove the sample-to-sample 
# technical variation while preserving the true biological condition effects.
# This tells LEMUR: "Please align cells from different samples so that biologically 
# similar cells are close together, regardless of which sample they came from."
# It aligns the coordinate systems so that the same cell types from different 
# samples have similar coordinates in the latent space.

fit

# UMAP of fit$embedding
umap <- uwot::umap(t(fit$embedding))

# UMAP colored by sample_id
as_tibble(fit$colData) |>
  mutate(umap = umap) |>
  ggplot(aes(x = umap[,1], y = umap[,2])) +
  geom_point(aes(color = sample_id), size = 0.5) +
  facet_wrap(vars(condition)) + coord_fixed()

# UMAP colored by condition
as_tibble(fit$colData) |>
  mutate(umap = umap) |>
  ggplot(aes(x = umap[,1], y = umap[,2])) +
  geom_point(aes(color = condition), size = 0.5) +
  facet_wrap(vars(condition)) + coord_fixed()

# UMAP colored by cell type
as_tibble(fit$colData) |>
  mutate(umap = umap) |>
  ggplot(aes(x = umap[,1], y = umap[,2])) +
  geom_point(aes(color = celltype), size = 0.5) +
  facet_wrap(vars(condition)) + coord_fixed()


# Test differential expression - predict the effect of the Alzheimer's disease for each cell and each gene
# New "DE" slot in assays contains the predicted logarithmic fold changes between the two conditions specified in contrast
fit <- test_de(fit, 
               contrast = cond(condition = "AD") - cond(condition = "HC"))

fit # new slot "DE" is added into the assays 

# Make a dataframe with predicted logarithmic fold changes between the two conditions 
df_de <- as.data.frame(assay(fit, "DE"))

# Differential expression (log fold changes) of SIRT4 on the UMAP
df <- tibble(umap = umap) |>
  mutate(de = assay(fit, "DE")["SIRT4", ])

ggplot(df, aes(x = umap[,1], y = umap[,2])) +
  geom_point(aes(color = de)) +
  scale_color_gradient2(low = "#FFD800", high= "#0056B9") + coord_fixed()

ggplot(df, aes(x = de)) + geom_histogram(bins = 100)

# Search for cell neighborhoods that show consistent differential expression
# The group_by argument determines how the pseudobulk samples are formed.

# Group by both sample and condition (most conservative)
nei <- find_de_neighborhoods(fit, 
                             group_by = vars(sample_id, condition))


