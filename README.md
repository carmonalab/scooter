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

#########
BELOW TBD
#########
<br>

# Summarized cell annotation

`Run.HiTME` will return the Seurat object or list of them with new metadata indicating cell type annotation.

Annotated Seurat objects can be summarized into **HiT objects** using `get.HiTObject` function. For this function the grouping variable `group.by` resulting from `Run.HiTME` annotation or additional annotations need to be indicated. Compositional cell type distribution and aggregated transcriptomic profile (pseudobulk) are returned for each sample.

``` r
HiT_summary <- get.HiTObject(annotated.obj ,
                            group.by = list("layer1" = "scGate_multi",
                                            "layer2" = "functional.cluster"))
```

Alternatively, HiT summarizing object can be obtained directly using `Run.HiTME` with parameters `return.Seurat = FALSE`.

``` r
HiT_summary <- Run.HiTME(object = obj,
                        scGate.model = models.TME,
                        ref.maps = ref.maps,
                        return.Seurat = FALSE)
```

## Hit Object content

The Hit object summarize the cell type annotation and contain the following slots:

1.  Seurat object metadata (dataframe): `metadata`

2.  Cell type predictions for each cell in the data set (list): `predictions`

3.  Cell type composition for each layer of cell type prediction: `composition`. Including:


    3.1. Cell counts

    3.2. Frequency

    3.3. CLR (Centred log ratio)-transformed counts (useful for downstream analyses such as PCA/[Logratio analysis](https://doi.org/10.1146/annurev-statistics-042720-124436) )


4.  Aggregated profile of predicted cell types: `aggregated_profile`. Including:

    4.1. Average and aggregated expression per cell type of all genes in the dataset and a subset of them.

    4.2. Mean of [UCell](https://github.com/carmonalab/UCell) scores per cell type, if additional signatures are provided, for example from [SignatuR](https://github.com/carmonalab/SignatuR).
