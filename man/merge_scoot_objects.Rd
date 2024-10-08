% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{merge_scoot_objects}
\alias{merge_scoot_objects}
\title{Merge scoot objects}
\usage{
merge_scoot_objects(
  scoot_object = NULL,
  group_by = NULL,
  metadata_vars = NULL,
  pseudobulk_matrix = "list",
  verbose = FALSE
)
}
\arguments{
\item{scoot_object}{List of scoot objects}

\item{group_by}{If only merging for certain layers of annotation is intended, layers names can be indicated here as vector. Otherwise all layers present in all scoot object will be merged.}

\item{metadata_vars}{Variables to keep as metadata. (Default: NULL, keeping unique metadata columns per sample, dropping single-cell metadata)}

\item{pseudobulk_matrix}{Paramater to determine whther obtain the pseudobulk matrix as a single matrix (\code{"unique"}), or as one matrix for each cell type in the layer (\code{"list"})}

\item{verbose}{Whether to show optional messages or not}

\item{ncores}{The number of cores to use, by default, all available cores - 2.}

\item{bparam}{A \code{BiocParallel::bpparam()} object that tells how to parallelize. If provided, it overrides the \code{ncores} parameter.}

\item{progressbar}{Whether to show a progressbar or not}
}
\value{
Merged scootObject
}
\description{
Merge scoot objects
}
