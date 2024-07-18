# scooter

## Exploratory single-cell data analysis at the sample level

scxplor is designed to help you explore your single-cell RNA-seq data in a simple and time-efficient way. It summarizes your cell type annotations by providing tools for compositional data analysis as well as tools for gene expression on the sample and cell type level.
- **Supervised analysis:** visualize your samples with box plots comparing groups or PCA colored by group
- **Unsupervised analysis:** cluster your samples based on their similarity

Why unsupervised analysis?
- **Quality Control**: Identify outliers, biases, and potential technical artifacts.
- **Dimensionality Reduction**: Condense 1000s of dimensions into a few highly interpretable features.
- **Detection of Biological Variability**: reveal important insights into population heterogeneity, developmental trajectories, and responses to stimuli or disease.

### Installation

``` r
# install.packages("remotes")
remotes::install_github("carmonalab/scooter")
```

<br>

## Summarize your scRNA-seq data

A list of annotated Seurat objects can be summarized into a list of **scoot objects** using the `scoot` function. Compositional cell type distribution and aggregated transcriptomic profile (pseudobulk) are returned for each sample.

``` r
obj.list <- SplitObject(obj, split.by = "Sample")

scoot_object_list <- scoot(obj.list)

scoot_summary <- merge_scoot_objects(scoot_object_list)
```

### scoot object content

The scoot object summarize the cell type annotation and contain the following slots:

- Seurat object metadata (dataframe): `metadata`
- Cell type composition for each layer of cell type prediction: `composition`. Including:
  - Cell counts
  - Frequency
  - CLR (Centred log ratio)-transformed counts (useful for downstream analyses such as PCA/[Logratio analysis](https://doi.org/10.1146/annurev-statistics-042720-124436) )
- Aggregated profile of predicted cell types: `aggregated_profile`. Including:
  - Aggregated expression per cell type.
  - Mean of [UCell](https://github.com/carmonalab/UCell) scores per cell type, if additional signatures are provided, for example from [SignatuR](https://github.com/carmonalab/SignatuR).
