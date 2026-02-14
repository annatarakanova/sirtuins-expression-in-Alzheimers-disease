# sirtuins-expression-in-Alzheimers-disease
The dysregulation of sirtuins across neural cell types and brain regions in the context of Alzheimer’s disease.

## The data
snRNA-seq data from 3 human brain regions (Prefrontal cortex, Middle temporal gyrus and Entorhinal cortex) of healthy patients and patients with Alzheimer's disease.

## The goal
To investigate the dysregulation of sirtuins across five brain cell types of three brain regions in the context of Alzheimer’s disease.

## The workflow
1. Analyze snRNA-seq dataset and SIRT expression
2. Subset the dataset by cell types
3. For each cell type create a pseudobulk profile
4. Pre-filter lowly expressed genes
5. Find DEGs between AD and control groups (edgeR)
6. Analyze SIRT expression in each cell type

## The results
Pseudobulk analysis revealed no statistically significant changes in the expression of sirtuin family genes in Alzheimer's disease. However, there was a tendency to increase
SIRT4 expression in neurons and oligodendrocytes. This result requires additional verification using more advanced methods (LEMUR) capable of accounting for heterogeneity within a single cell
type.
