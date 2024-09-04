#' Get a scoot object summarizing cell type classification and aggregated profiles.
#'
#' @param object A Seurat object or a list of Seurat objects.
#' @param ann_layer_cols List or vector with one or multiple Seurat object metadata columns with cell type annotations (e.g. layer 1 broad cell type annotation (CD4, CD8, ...) and layer 2 with high-resolution subtypes (CD4.treg, CD4.tfh, CD8.naive, CD8.exhausted, ...))
#' @param min_cells Minimum number of cells for any annotated cell type in layer_1.
#' @param name_additional_signatures Names of additional signatures as found in object metadata to take into account.
#' @param useNA logical whether to return aggregated profile for NA (undefined) cell types.
#' @param clr_zero_impute_perc Parameter for internal \link{get_celltype_composition}.
#' @param layer_links Parameter for internal \link{get_celltype_composition}
#' @param assay Parameter for internal \link{get_aggregated_profile}.
#' @param layer_1_link Column of metadata linking layer_1 prediction (e.g. scGate ~ CellOntology_ID) in order to perform subsetting for second layer classification.
#'
#' @importFrom methods setClass new
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallelly availableCores
#'
#' @import SeuratObject
#'
#' @return scoot object summarizing cell type classification and aggregated profiles.
#' @export scoot
#'


scoot <- function(object,
                  ann_layer_cols = list("layer_1" = "scGate_multi",
                                        "layer_2" = "functional.cluster"
                  ),
                  min_cells = 10,
                  name_additional_signatures = c("IFN_UCell", "HeatShock_UCell",
                                                 "cellCycle.G1S_UCell", "cellCycle.G2M_UCell"),
                  useNA = FALSE,
                  clr_zero_impute_perc = 1,
                  layer_links = c("scGate_multi" = "functional.cluster"),
                  assay = "RNA",
                  layer_1_link = "CellOntology_ID",

                  ncores = parallelly::availableCores() - 2,
                  bparam = NULL,
                  progressbar = TRUE) {

  if (is.null(ann_layer_cols)) {
    stop("Please provide at least one annotation column")
  }

  if (!is.list(ann_layer_cols)) {
    ann_layer_cols <- as.list(ann_layer_cols)
  }

  # Run extended if object got scGate info
  # Extract values from misc slot from object
  if (any(ann_layer_cols == "scGate_multi")) {
    if (is.null(name_additional_signatures)) {
      name_additional_signatures <- c("IFN_UCell", "HeatShock_UCell",
                                      "cellCycle.G1S_UCell", "cellCycle.G2M_UCell")
    }
  }

  if (!is.numeric(min_cells) &
      !length(min_cells) == 1) {
    stop("Please provide an number for min_cells")
  }

  if (min_cells < 1) {
    message("Please provide an number larger than 1 for min_cells. min_cells was changed to 1")
    min_cells <- 1
  }

  args <- list(ann_layer_cols,
               min_cells,
               name_additional_signatures,
               useNA,
               clr_zero_impute_perc,
               layer_links,
               assay,
               layer_1_link)

  # if object is a single Seurat object, turn into a list
  if (!is.list(object)) {
    if (all(sapply(table(object[[ann_layer_cols[[1]]]]), function(x) {x < min_cells}))) {
      message("Not enough annotated cells in object")
      return(NULL)
    }

    scoot_list <- scoot_helper(object,
                               ann_layer_cols,
                               min_cells,
                               name_additional_signatures,
                               useNA,
                               clr_zero_impute_perc,
                               layer_links,
                               assay,
                               layer_1_link)
  } else {
    if (is.null(names(object))) {
      names(object) <- seq_along(object)
    }
    nams <- names(object)

    # set parallelization parameters
    param <- set_parallel_params(ncores = ncores,
                                 bparam = bparam,
                                 progressbar = progressbar)

    scoot_list <- BiocParallel::bplapply(X = names(object),
                                         BPPARAM = param,
                                         function(x) {
                                           if (all(sapply(table(object[[x]][[ann_layer_cols[[1]]]]), function(x) {x < min_cells}))) {
                                             message(paste0("Not enough annotated cells in object ", x, ". Returning NULL."))
                                             return(NULL)
                                           }

                                           do.call(scoot_helper,
                                                   c(object[[x]], args)
                                           )
                                         })
    names(scoot_list) <- nams
  }

  return(scoot_list)
}

scoot_helper <- function(object,
                         ann_layer_cols,
                         min_cells,
                         name_additional_signatures,
                         useNA,
                         clr_zero_impute_perc,
                         layer_links,
                         assay,
                         layer_1_link) {

  if (is.null(object)) {
    stop("Please provide a Seurat object")
  } else if (!inherits(object, "Seurat")) {
    stop("Not Seurat object included, cannot be processed.\n")
  }


  if (!any(ann_layer_cols %in% names(object@meta.data))) {
    stop("ann_layer_cols ", paste(ann_layer_cols, collapse = ", "), " not found in metadata")
  } else if (!all(ann_layer_cols %in% names(object@meta.data))) {
    ann_layer_cols <- ann_layer_cols[ann_layer_cols %in% names(object@meta.data)]
    message("Only found ", paste(ann_layer_cols, collapse = ", ") , " as grouping variable for scoot object.")
  }

  for (v in seq_along(ann_layer_cols)) {
    if (is.null(names(ann_layer_cols)[[v]]) || is.na(names(ann_layer_cols)[[v]])) {
      names(ann_layer_cols)[[v]] <- paste0("layer", v)
    }
  }


  # Get unique metadata columns
  # Get metadata column names which have all the same value.
  # Each object is from a specific sample, so each object has e.g. a metadata column "Sample" which contains all the same value (e.g. "sample1" for sample 1, "sample2" for sample 2, etc.)
  # Other metadata columns, e.g. scGate_multi can be dropped, as the merged scoot object is a summary per sample, so single-cell metadata columns don't make sense.
  umc <- apply(object@meta.data, 2, function(y) length(unique(y))) == 1
  umc_names <- names(umc)[umc == TRUE]
  metadata <- object@meta.data[1, umc_names]


  # Compute proportions
  # message("\nComputing cell type composition...\n")

  comp_prop <- get_celltype_composition(object,
                                        ann_layer_cols = ann_layer_cols,
                                        min_cells_composition = min_cells,
                                        useNA = useNA,
                                        clr_zero_impute_perc = clr_zero_impute_perc,
                                        layer_1_link = layer_1_link,
                                        layer_links = layer_links)
  # Compute avg expression
  # message("\nComputing aggregated profile...\n")

  avg_expr <- get_aggregated_profile(object,
                                     ann_layer_cols = ann_layer_cols,
                                     min_cells_aggregated = min_cells,
                                     assay = assay,
                                     useNA = useNA)

  aggr_sig <- get_aggregated_signature(object,
                                       ann_layer_cols = ann_layer_cols,
                                       min_cells_aggregated = min_cells,
                                       name_additional_signatures = name_additional_signatures,
                                       useNA = useNA)

  scoot <- methods::new("scoot",
                        metadata = metadata,
                        aggregated_profile = list("pseudobulk" = avg_expr,
                                                  "signatures" = aggr_sig),
                        composition = comp_prop
  )
  return(scoot)
}


#' Calculate cell type composition or frequencies
#'
#' @param object A seurat object
#' @param ann_layer_cols The Seurat object metadata column(s) containing celltype annotations (provide as character vector, containing the metadata column name(s))
#' @param split_by A Seurat object metadata column to split by (e.g. sample names)
#' @param min_cells_composition Set a minimum threshold for number of cells to calculate relative abundance (e.g. less than 10 cells -> no relative abundance will be calculated)
#' @param useNA Whether to include not annotated cells or not (labelled as "NA" in the ann_layer_cols). Can be defined separately for each ann_layer_cols (provide single boolean or vector of booleans)
#' @param layer_links Named vector indicating the relation of multiple layer classification, default is \code{c("scGate_multi", "functional.cluster")}. Name of the vector element ought to be layer_1 and element layer_2. This parameter is used to compute relative compositional data of layer_2 within layer_1 classes.
#' @param clr_zero_impute_perc To calculate the clr-transformed relative abundance ("clr_freq"), zero values are not allowed and need to be imputed (e.g. by adding a pseudo cell count). Instead of adding a pseudo cell count of flat +1, here a pseudo cell count of +1% of the total cell count will be added to all cell types, to better take into consideration the relative abundance ratios (e.g. adding +1 cell to a total cell count of 10 cells would have a different, i.e. much larger effect, than adding +1 to 1000 cells).
#' @param layer_1_link Column of metadata linking layer_1 prediction (e.g. scGate ~ CellOntology_ID) in order to perform subsetting for second layer classification.

#' @importFrom Hotelling clr
#' @importFrom dplyr group_by summarize filter ungroup mutate select left_join n coalesce bind_rows across all_of
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom rrapply rrapply
#'
#' @return Cell type compositions as a list of data.frames containing cell counts, relative abundance (freq) and clr-transformed freq (freq_clr), respectively.
#' @export get_celltype_composition
#'
#' @examples
#' devtools::install_github('satijalab/seurat-data')
#' library(SeuratData)
#' options(timeout = max(300, getOption("timeout")))
#' InstallData("panc8")
#' data("panc8")
#' panc8 = UpdateSeuratObject(object = panc8)
#' # Calculate overall composition
#' celltype_compositions.overall <- get_celltype_composition(object = panc8, ann_layer_cols = "celltype")
#'
#' # Calculate sample-wise composition

#' celltype_compositions.sample_wise <- get_celltype_composition(object = panc8, ann_layer_cols = "celltype", split_by = "orig.ident")
#'
get_celltype_composition <- function(object = NULL,
                                     ann_layer_cols = NULL,
                                     split_by = NULL,
                                     min_cells_composition = 10,
                                     useNA = FALSE,
                                     layer_links = c("scGate_multi" = "functional.cluster"),
                                     clr_zero_impute_perc = 1,
                                     layer_1_link = "CellOntology_ID") {

  if (is.null(object)) {
    stop("Please provide a Seurat object or metadata as dataframe")
  }

  # input can be a Seurat object or a dataframe containing its meta.data
  # convert object to metadata if seurat object is provided
  if (inherits(object, "Seurat")) {
    meta.data <- object@meta.data
    if (is.null(meta.data)) {
      stop("No metadata found in this Seurat object")
    }
  } else {
    stop("Not Seurat object or dataframe included, cannot be processed.\n")
  }

  if (is.null(ann_layer_cols)) {
    stop("Please specificy an ann_layer_cols variable")
  }

  # Assess wheter split_by variable is in metadata
  if (!is.null(split_by) &&
      !split_by %in% names(meta.data)) {
    stop("split_by variable not found in meta.data!\n")
  }

  if (length(useNA) != 1 &&
      length(useNA) != length(ann_layer_cols)) {
    stop("useNA variable must be of length 1 or the same length as ann_layer_cols")
  }

  # convert ann_layer_cols to list if not
  if (!is.list(ann_layer_cols)) {
    ann_layer_cols <- as.list(ann_layer_cols)
  }


  # Rename ann_layer_cols if not indicated
  for (v in seq_along(ann_layer_cols)) {
    if (is.null(names(ann_layer_cols)[[v]]) ||
        is.na(names(ann_layer_cols)[[v]])) {
      names(ann_layer_cols)[[v]] <- ann_layer_cols[[v]]
    }

  }



  # Keep only grouping variables in metadata
  if (!any(ann_layer_cols %in% names(meta.data))) {
    stop("ann_layer_cols ", paste(ann_layer_cols, collapse = ", "), " not found in metadata")
  } else if (!all(ann_layer_cols %in% names(meta.data))) {
    group_present <- ann_layer_cols %in% names(meta.data)
    # accommodate useNA vector
    if (length(useNA) == length(ann_layer_cols)) {
      useNA <- useNA[group_present]

    }
    # keep only grouping variables found in metadata
    ann_layer_cols <- ann_layer_cols[group_present]
    message("Only found ", paste(ann_layer_cols, collapse = ", ") , " as grouping variables.")
  }

  # Factorize cell type to keep all in compositional data
  meta.data <- meta.data %>%
    dplyr::mutate_at(unlist(ann_layer_cols), as.factor)

  # evaluate which layers are included in ann_layer_cols and layer_links
  for (l in seq_along(layer_links)) {
    if (!all(c(names(layer_links[l]), layer_links[l]) %in%
             ann_layer_cols)) {
      message("Not computing relative proportion values of layer_2 within layer_1 for ",
              paste0(names(layer_links[l]), " -> ", layer_links[l]))
    }
  }


  # Rep useNa parameters according to ann_layer_cols variable
  if (length(useNA) == 1) {
    useNA <- rep(useNA, length(ann_layer_cols))
  }

  # list to store compositions
  celltype_compositions <- list()


  for (i in seq_along(ann_layer_cols)) {

    suppressMessages(
      {
        ctable <- compositional_data(data = meta.data,
                                     split_by = split_by,
                                     ann_layer_col_1 = ann_layer_cols[[i]],
                                     useNA = useNA[i],
                                     clr_zero_impute_perc = clr_zero_impute_perc)

        # stop if very few cells
        if (sum(ctable$cell_counts) < min_cells_composition) {
          # Return empty data.frame
          ctable <- ctable[0,]
        }
        # get proportion relative to layer_1 types
        if (ann_layer_cols[[i]] %in% layer_links) {
          # stop if very few cells
          if (sum(ctable$cell_counts) < min_cells_composition) {
            # Return empty list
            ctable <- list()
          } else {
            lay <- layer_links[which(layer_links == ann_layer_cols[[i]])]
            if (lay %in% names(object@misc$layer2_param)) {
              meta_split <- split(meta.data,
                                  meta.data[[names(lay)]])
              # switch list names to cell ontology ID
              cellonto_dic <- lapply(meta_split, function(x) {
                nam <- unique(x[[names(lay)]])
                val <- unique(x[[layer_1_link]])
                names(val) <- nam
                return(val)
              }) %>%
                unname() %>%
                unlist()

              levs <- object@misc$layer2_param[[lay]]$levels2_per_levels1
              names(levs) <- names(cellonto_dic[match(names(levs), cellonto_dic)])

              # If a celltype was not detected, drop it
              levs <- levs[!is.na(names(levs))]

              # Filter only layer_1 cells with representation in layer_2
              meta_split <- meta_split[names(levs)]

              # add factor to group by layer_1
              for (f in names(meta_split)) {
                meta_split[[f]]$functional.cluster <- factor(meta_split[[f]]$functional.cluster,
                                                             levels = levs[[f]])
              }

              ctable_split <- lapply(meta_split,
                                     compositional_data,
                                     split_by = split_by,
                                     ann_layer_col_1 = ann_layer_cols[[i]],
                                     useNA = useNA[i])

              # If not enough cells, return empty dataframe
              ctable_split <- lapply(ctable_split,
                                     function(x) {
                                       if (sum(x$cell_counts) < min_cells_composition) {
                                         x[0,]
                                       } else {x}
                                     })

              ctable <- c(ctable_split)
            }
          }
        }
      })

    ## Append
    celltype_compositions[[names(ann_layer_cols)[i]]] <- ctable
  }

  return(celltype_compositions)

}



#' Compute aggregated gene expression
#'
#' Function to compute aggregated expression (pseudobulk, i.e. sum counts per ident), and average expression by indicated cell type or grouping variable.
#'
#'
#' @param object A seurat object or a list of seurat objects
#' @param ann_layer_cols The Seurat object metadata column(s) containing celltype annotations
#' @param min_cells_aggregated Minimum number of cells per sample and cell type to aggregate data for pseudobulk or aggregated signatures. Aggregated data for cell types with less than that value will not be returned in aggregated data. Default value is 10.
#' @param assay Assay to retrieve information. By default "RNA".
#' @param ... Extra parameters for internal Seurat functions: AverageExpression, AggregateExpression, FindVariableFeatures

#' @importFrom Seurat AverageExpression AggregateExpression FindVariableFeatures
#' @importFrom dplyr filter
#' @return Average and aggregated expression as a list of matrices for all genes and indicated gene lists filtering.
#' @export get_aggregated_profile


get_aggregated_profile <- function(object,
                                   ann_layer_cols = NULL,
                                   min_cells_aggregated = 10,
                                   assay = "RNA",
                                   ...) {

  if (is.null(object)) {
    stop("Please provide a Seurat object")
    if (!inherits(object, "Seurat")) {
      stop("Please provide a Seurat object")
    }
  }

  if (is.null(ann_layer_cols)) {
    stop("Please specificy an ann_layer_cols variable")
  }

  # convert ann_layer_cols to list if not
  if (!is.list(ann_layer_cols)) {
    ann_layer_cols <- as.list(ann_layer_cols)
  }

  # Rename ann_layer_cols if not indicated
  for (v in seq_along(ann_layer_cols)) {
    if (is.null(names(ann_layer_cols)[[v]]) || is.na(names(ann_layer_cols)[[v]])) {
      names(ann_layer_cols)[[v]] <- ann_layer_cols[[v]]
    }
  }

  if (!any(ann_layer_cols %in% names(object@meta.data))) {
    stop("ann_layer_cols ", paste(ann_layer_cols, collapse = ", "), " not found in metadata")
  } else if (!all(ann_layer_cols %in% names(object@meta.data))) {
    ann_layer_cols <- ann_layer_cols[ann_layer_cols %in% names(object@meta.data)]
    message("Only found ", paste(ann_layer_cols, collapse = ", ") , " as grouping variables.")
  }


  avg_exp <- list()

  # loop over different grouping
  for (i in names(ann_layer_cols)) {

    # compute pseudobulk
    suppressWarnings({

      # Calculate total (sample) pseudobulks
      if (i == names(ann_layer_cols)[1]) {

        # Calculate pseudobulk for ALL cells in sample
        avg_exp[[i]] <- object@assays[["RNA"]]["counts"]
        row_names <- row.names(avg_exp[[i]])
        avg_exp[[i]] <- Matrix::Matrix(rowSums(avg_exp[[i]]))
        row.names(avg_exp[[i]]) <- row_names
        colnames(avg_exp[[i]]) <- "all"


        # Calculate pseudobulk for annotated cells in sample
        # Subset to remove not annotated cells
        obj <- object[, which(!is.na(object[[ann_layer_cols[[i]]]]))]

        # Calculate pseudobulk of all annotated cells
        mat <- obj@assays[["RNA"]]["counts"]
        row_names <- row.names(mat)
        mat <- Matrix::Matrix(rowSums(mat))
        row.names(mat) <- row_names
        colnames(mat) <- "all.annotated_only"

        avg_exp[[i]] <- cbind(avg_exp[[i]], mat)


        # Calculate pseudobulk for not annotated cells in sample
        # Subset to remove annotated cells
        obj <- object[, which(is.na(object[[ann_layer_cols[[i]]]]))]

        # Calculate pseudobulk of all annotated cells
        mat <- obj@assays[["RNA"]]["counts"]
        row_names <- row.names(mat)
        mat <- Matrix::Matrix(rowSums(mat))
        row.names(mat) <- row_names
        colnames(mat) <- "all.not_annotated_only"

        avg_exp[[i]] <- cbind(avg_exp[[i]], mat)
      }

      # remove from aggregated data cell with less than min_cells_aggregated
      cnts <- compositional_data(object@meta.data,
                                 ann_layer_col_1 = ann_layer_cols[[i]],
                                 only_counts = TRUE)

      # Remove cell types with less than min_cells_aggregated
      cnts <- cnts[!is.na(cnts[["celltype"]]),]
      cnts <- cnts[cnts[["cell_counts"]] > min_cells_aggregated, 1]
      keep <- unlist(cnts)

      if (length(keep) == 0) {
        avg_exp[[i]] <- NULL
      } else {
        object@meta.data$fltr <- object@meta.data[[ann_layer_cols[[i]]]]
        object <- object[, as.character(object$fltr) %in% keep]
        object <- object[, !is.na(object$fltr)]

        if (length(unique(object@meta.data[[ann_layer_cols[[i]]]])) >= 2) {
          suppressWarnings(
            suppressMessages(
              mat <-
                Seurat::AggregateExpression(object,
                                            group.by = ann_layer_cols[[i]],
                                            assays = assay,
                                            verbose = FALSE,
                                            ...)[[assay]]
            )
          )
          avg_exp[[i]] <- cbind(avg_exp[[i]], mat)
        } else {
          # Handle case if there is only one cell type
          col_name <- as.character(unique(object@meta.data[[ann_layer_cols[[i]]]]))
          avg_exp[[i]] <- object@assays[["RNA"]]["counts"]
          row_names <- row.names(avg_exp[[i]])
          avg_exp[[i]] <- Matrix::Matrix(rowSums(avg_exp[[i]]))
          row.names(avg_exp[[i]]) <- row_names
          colnames(avg_exp[[i]]) <- col_name
        }
      }
    })
  }
  return(avg_exp)
}


#' Compute aggregated additional signatures by cell type
#'
#' Function to compute aggregated signatures of predicted cell types.
#'
#'
#' @param object A seurat object or metadata data frame.
#' @param ann_layer_cols The Seurat object metadata column(s) containing celltype annotations (idents).
#' @param min_cells_aggregated Minimum number of cells per sample and cell type to aggregate data for pseudobulk or aggregated signatures. Aggregated data for cell types with less than that value will not be returned in aggregated data. Default value is 10.
#' @param name_additional_signatures Names of additional signatures to compute the aggregation per cell type.
#' @param fun Function to aggregate the signature, e.g. mean or sum.
#' @param useNA logical whether to return aggregated signatures for NA (undefined) cell types, default is FALSE.

#' @importFrom dplyr group_by summarize_at filter
#' @return Aggregated signature score for each indicated cell type grouping Results is NULL of not additional signatures are indicated or present in metadata.
#' @export get_aggregated_signature


get_aggregated_signature <- function(object,
                                     ann_layer_cols = NULL,
                                     min_cells_aggregated = 10,
                                     name_additional_signatures = NULL,
                                     fun = mean,
                                     useNA = FALSE) {

  # input can be a Seurat object or a dataframe containing its meta.data
  # convert object to metadata if seurat object is provided
  if (inherits(object, "Seurat")) {
    meta.data <- object@meta.data
    if (is.null(meta.data)) {
      stop("No metadata found in this Seurat object")
    }
  } else if (inherits(object, "data.frame")) {
    meta.data <- object
  } else {
    stop("Not Seurat object or dataframe included, cannot be processed.\n")
  }


  if (is.null(ann_layer_cols)) {
    stop("Please specificy an ann_layer_cols variable")
  }

  # convert ann_layer_cols to list if not
  if (!is.list(ann_layer_cols)) {
    ann_layer_cols <- as.list(ann_layer_cols)
  }

  # Rename ann_layer_cols if not indicated
  for (v in seq_along(ann_layer_cols)) {
    if (is.null(names(ann_layer_cols)[[v]]) ||
        is.na(names(ann_layer_cols)[[v]])) {
      names(ann_layer_cols)[[v]] <- ann_layer_cols[[v]]
    }
  }

  if (!any(ann_layer_cols %in% names(meta.data))) {
    stop("ann_layer_cols ", paste(ann_layer_cols, collapse = ", "), " not found in metadata")
  } else if (!all(ann_layer_cols %in% names(meta.data))) {
    ann_layer_cols <- ann_layer_cols[ann_layer_cols %in% names(meta.data)]
    message("Only found ", paste(ann_layer_cols, collapse = ", ") , " as grouping variable.")
  }

  if (is.null(name_additional_signatures)) {
    aggr_sig <- NULL
  } else {

    if (!any(name_additional_signatures %in% colnames(meta.data))) {
      aggr_sig <- NULL
      message("No signatures columns found in this object's metadata")
    } else {
      add_sig_cols <- name_additional_signatures[name_additional_signatures %in%
                                                   colnames(meta.data)]

      if (length(add_sig_cols) < length(name_additional_signatures)) {
        not_found <- name_additional_signatures[!name_additional_signatures %in%
                                                  colnames(meta.data)]
        message(paste("The following additional signature columns were not found
                    in the object's metadata colnames: ", not_found))
      }

      aggr_sig <- list()

      for (e in names(ann_layer_cols)) {
        # remove from aggregated data cell with less than min_cells_aggregated
        cnts <- compositional_data(meta.data,
                                   ann_layer_col_1 = ann_layer_cols[[e]],
                                   only_counts = TRUE,
                                   useNA = useNA)

        keep <- cnts[cnts[["cell_counts"]] > min_cells_aggregated, 1] %>%
          unlist()


        aggr_sig[[e]] <- meta.data %>%
          dplyr::filter(.data[[ann_layer_cols[[e]]]] %in% keep) %>%
          dplyr::group_by(.data[[ann_layer_cols[[e]]]]) %>%
          dplyr::summarize_at(add_sig_cols, fun, na.rm = TRUE)

        # filter out NA if useNA=F
        if (!useNA) {
          aggr_sig[[e]] <- aggr_sig[[e]] %>%
            dplyr::filter(!is.na(.data[[ann_layer_cols[[e]]]]))
        }

        colnames(aggr_sig[[e]])[1] <- "celltype"
      }
    }
  }
  return(aggr_sig)
}



#' Merge scoot objects
#'
#' @param scoot_object List of scoot objects
#' @param group_by If only merging for certain layers of annotation is intended, layers names can be indicated here as vector. Otherwise all layers present in all scoot object will be merged.
#' @param metadata_vars Variables to keep as metadata. (Default: NULL, keeping unique metadata columns per sample, dropping single-cell metadata)
#' @param pseudobulk_matrix Paramater to determine whther obtain the pseudobulk matrix as a single matrix (\code{"unique"}), or as one matrix for each cell type in the layer (\code{"list"})
#' @param ncores The number of cores to use, by default, all available cores - 2.
#' @param bparam A \code{BiocParallel::bpparam()} object that tells how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param progressbar Whether to show a progressbar or not
#' @param verbose Whether to show optional messages or not

#' @importFrom dplyr mutate mutate_if filter %>% mutate_all full_join
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom data.table rbindlist
#' @importFrom methods slot

#' @return Merged scootObject
#' @export merge_scoot_objects
#'

merge_scoot_objects <- function(scoot_object = NULL,
                                group_by = NULL,
                                metadata_vars = NULL,
                                pseudobulk_matrix = "list",
                                verbose = FALSE) {

  if (is.null(scoot_object) ||
      !is.list(scoot_object)) {
    stop("Please provide a list of scoot object containing more than one sample")
    if (suppressWarnings(!all(lapply(scoot_object, function(x) {inherits(x, "scoot")})))) {
      stop("Not all components of the list are scoot objects.")
    }
  }
  scoot_object[sapply(scoot_object, is.null)] <- NULL

  # fetch which layers are included in all scoot objects
  layers_in_scoot <- lapply(scoot_object, function(x) {
    a <- names(x@composition)
    b <- names(x@aggregated_profile$pseudobulk)
    c <- names(x@aggregated_profile$signatures)
    u <- unique(c(a,b,c))
    return(u)
  })

  # if group_by is not NULL retrieve the layers present in scoot object
  if (is.null(group_by)) {
    group_by <- unique(unlist(layers_in_scoot))
  }


  if (suppressWarnings(!all(lapply(layers_in_scoot, function(x) {any(group_by %in% x)})))) {
    stop("None of the supplied scoot object contain at least one of these layers: ", paste(group_by, collapse = ", "),
         ".\nPlease make sure to indicate the name of the layer.")
  } else {
    present <- sapply(group_by,
                      function(char) {
                        unlist(lapply(layers_in_scoot,
                                      function(df) {char %in% df}
                        ))
                      }
    )

    present_sum <- colSums(present)

    for (l in names(group_by)) {
      message("** ", l, " present in ", present_sum[[l]], " / ", length(scoot_object), " scoot objects.")
    }
  }

  if (!is.null(metadata_vars)) {
    if (verbose) {message("\n#### Metadata ####")}
    if (suppressWarnings(!all(lapply(scoot_object, function(x) {any(metadata_vars %in% names(x@metadata))})))) {
      message("Not all supplied scoot object contain ", paste(metadata_vars, collapse = ", "),
              "metadata elements in their metadata")
    }
    if (verbose) {
      in_md <- sapply(metadata_vars,
                      function(char) {
                        unlist(lapply(scoot_object,
                                      function(df) {char %in% colnames(df@metadata)}
                        ))
                      }
      )
      in_md_sum <- colSums(in_md)

      for (l in metadata_vars) {
        message("** ", l, " present in ",
                in_md_sum[[l]], " / ", length(scoot_object),
                " scoot objects.")
      }
    }
  }

  # Join data from all scoot objects

  # Metadata
  if (is.null(metadata_vars)) {
    # Per scootObject from input, get metadata column names which have all the same value.
    # E.g. each scootObject is from a specific sample, so each scootObject has metadata column "Sample" which contains all the same value (e.g. "sample1" for sample 1, "sample2" for sample 2, etc.)
    # Other metadata columns, e.g. scGate_multi can be dropped, as the merged scootObject is a summary per sample, so single-cell metadata columns don't make sense.
    all_names <- list()
    for (x in names(scoot_object)) {
      all_names[[x]] <- names(scoot_object[[x]]@metadata)
    }
    all_names_in_all <- Reduce(intersect, all_names)
  }
  metadata <- lapply(names(scoot_object),
                     function(x) {scoot_object[[x]]@metadata[1, all_names_in_all] %>%
                         mutate(scoot_sample = x)
                     }
  )
  metadata <- data.table::rbindlist(metadata, use.names=TRUE, fill=TRUE)

  # Remove columns that contain identical values across all rows
  metadata <- metadata[, sapply(metadata, function(x) length(unique(x)) > 1)]

  comp_prop <- list()
  avg_expr <- list()
  aggr_sig <- list()

  for (gb in group_by) {

    layer_present <- row.names(present)[present[,gb]]

    # Composition
    message("Merging compositions of " , gb, "...")

    type <- "composition"

    # assess that all scoot objects in the list contain the composition data for that layer
    is_df_check <- scoot_object[layer_present] %>%
      lapply(methods::slot, name = type) %>%
      lapply("[[", gb) %>%
      lapply(is.data.frame) %>%
      unlist()
    is_list_check <- scoot_object[layer_present] %>%
      lapply(methods::slot, name = type) %>%
      lapply("[[", gb) %>%
      lapply(function(x) {is.list(x) &!inherits(x, "data.frame")}) %>%
      unlist()

    if (all(is_df_check)) {
      df <- lapply(X = layer_present,
                   function(x) {
                     scoot_object[[x]]@composition[[gb]] %>%
                       mutate(scoot_sample = x)
                   })
      df <- data.table::rbindlist(df, use.names=TRUE, fill=TRUE)
      comp_prop[[gb]] <- df

    } else if (all(is_list_check)) {
      gb_sublevel_unique_names <- scoot_object[layer_present] %>%
        lapply(methods::slot, name = type) %>%
        lapply("[[", gb) %>%
        lapply(names) %>%
        unlist() %>%
        unique()
      for (i in gb_sublevel_unique_names) {
        df <- lapply(X = layer_present,
                     function(x) {
                       if (!is.null(scoot_object[[x]]@composition[[gb]][[i]])) {
                         scoot_object[[x]]@composition[[gb]][[i]] %>%
                           mutate(scoot_sample = x)
                       }
                     })
        df <- data.table::rbindlist(df, use.names=TRUE, fill=TRUE)
        comp_prop[[gb]][[i]] <- df
      }
    }


    # Aggregated_profile
    message("Merging aggregated profiles of " , gb, "...")

    type <- "pseudobulk"

    # Aggregated_profile Pseudobulk
    celltypes <- lapply(names(scoot_object),
                        function(x) {
                          colnames(scoot_object[[x]]@aggregated_profile[[type]][[gb]])
                        }) %>%
      unlist() %>% unique()

    for (ct in celltypes) {
      ct_present <- lapply(names(scoot_object),
                           function(x) {
                             ct %in% colnames(scoot_object[[x]]@aggregated_profile[[type]][[gb]]) }) %>%
        unlist()
      layer_ct_present <- names(scoot_object)[(names(scoot_object) %in% layer_present) & ct_present]
      df <- lapply(X = layer_ct_present,
                   function(x) {
                     scoot_object[[x]]@aggregated_profile[[type]][[gb]][, ct]
                   })
      df <- do.call(cbind, df)
      df <- Matrix::Matrix(df, sparse = TRUE)
      colnames(df) <- layer_ct_present
      avg_expr[[gb]][[ct]] <- df
    }

    # join each matrix per celltype into a single matrix changing the colnames to accommodate sample source
    if (tolower(pseudobulk_matrix) == "unique") {
      avg_expr[[gb]] <- lapply(names(avg_expr[[gb]]),
                               function(x) {
                                 mat <- avg_expr[[gb]][[x]]
                                 colnames(mat) <- paste(x, colnames(mat), sep = "__")
                                 mat <- mat %>% as.data.frame() %>%
                                   tibble::rownames_to_column("gene")
                               }) %>%
        reduce(full_join, by = "gene") %>%
        # convert NA to 0
        mutate_if (is.numeric, ~ifelse(is.na(.), 0, .)) %>%
        tibble::column_to_rownames("gene") %>%
        as.matrix() %>%
        Matrix::Matrix(., sparse = TRUE)
    }


    # Aggregated_profile Signatures

    type <- "signatures"

    df <- lapply(X = layer_present,
                 function(x) {
                   if (!is.null(scoot_object[[x]]@aggregated_profile[[type]][[gb]])) {
                     scoot_object[[x]]@aggregated_profile[[type]][[gb]] %>%
                       mutate(scoot_sample = x)
                   }
                 })
    df <- data.table::rbindlist(df,
                                use.names = TRUE,
                                fill = TRUE)
    aggr_sig[[gb]] <- df
  }

  scoot <- methods::new("scoot",
                        metadata = metadata,
                        composition = comp_prop,
                        aggregated_profile = list("pseudobulk" = avg_expr,
                                                  "signatures" = aggr_sig)
  )
  return(scoot)
}


#' Get metrics of cell types pseudobulk clustering
#'
#' @param scoot_object A scoot class object obtained with \link{scoot_object} or pseudobulk raw count matrix (\code{scootObject@aggregated_profile$pseudobulk$layer_1})
#' @param cluster_by Vector indicating the variable for clustering, default is celltype (for the annotation) and sample
#' @param min_samples Minimum number of samples to calculate scores (do not go lower than 5)
#' @param cluster_by_drop_na Whether to keep (FALSE) or drop (TRUE) NAs present in cluster_by column.
#' @param batching Vector indicating the variable for batching to allow calculating scores per batch, to account for batch effect
#' @param scores scores to compute clustering of samples based on cell type prediction.
#' @param knn_k Number of k-nearest neighbors used for the calculation of modularity score and leiden clustering. Default "auto" uses the square root of the total number of samples.
#' @param n_clust Number of clusters for unsupervised clustering. Use "auto" to automatically determine it. Alternatively, put an integer to force a specific number of clusters.
#' @param dist_method Method to compute distance between celltypes, default is euclidean.
#' @param NbClust_method Hierarchical clustering method for fastNbClust (improved version of NbClust::NbClust).
#' @param ntests Number of shuffling events to calculate p-value for scores
#' @param seed Set seed for random shuffling events to calculate p-value for scores
#' @param pval.combine Method for combining p-values if calculated using batching. Default is "zmethod" weighted (Stouffer's) Z-method, weighting by sample size per batch. Alternatively, Fisher's method can be used with "fisher".
#' @param pca_comps_labs_invisible Parameter to pass to fviz_pca defining which labels to plot for composition PCA, see documentation of fviz_pca for details.
#' @param pca_pb_labs_invisible Parameter to pass to fviz_pca defining which labels to plot for pseudobulk PCA, see documentation of fviz_pca for details.
#' @param pca_sig_labs_invisible Parameter to pass to fviz_pca defining which labels to plot for signatures PCA, see documentation of fviz_pca for details.
#' @param nvar_genes Number of variable genes to assess samples. Default is 500.
#' @param black_list List of genes to discard from clustering, if "default" object "default_black_list" object is used. Alternative black listed genes can be provided as a vector or list.
#' @param pca_n_hvg Number of most highly variable genes to plot on the PCA (for pseudobulks).
#' @param ncores The number of cores to use, by default, all available cores - 2.
#' @param bparam A \code{BiocParallel::bpparam()} object that tells how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param progressbar Whether to show a progressbar or not

#' @importFrom tidyr separate pivot_wider
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallelly availableCores
#' @importFrom data.table rbindlist
#' @importFrom dplyr mutate mutate_if filter %>% coalesce mutate_all full_join row_number
#' @importFrom tibble rownames_to_column column_to_rownames remove_rownames
#' @importFrom ggplot2 aes geom_point guides theme geom_col labs guide_legend annotate theme_bw ggtitle geom_ribbon element_text set_last_plot
#' @importFrom MatrixGenerics rowVars rowMins
#' @importFrom BiocGenerics counts
#' @importFrom stringr str_to_title
#' @importFrom stats prcomp na.omit formula pnorm t.test
#' @importFrom metap sumlog sumz


#' @return Metrics of cell types pseudobulk clustering
#' @export get_cluster_score
#'


get_cluster_score <- function(scoot_object = NULL,
                              cluster_by = NULL,
                              min_samples = 5,
                              cluster_by_drop_na = TRUE,
                              batching = NULL,
                              scores = c("silhouette_isolated",
                                         # "silhouette",
                                         "modularity"),
                              dist_method = "euclidean",
                              clust_method = c("pam_clust", "h_clust"),
                              knn_k = "auto",
                              max_nc = 5,
                              n_clust = "auto",
                              hdbscan_min_pts = 3,
                              NbClust_method = "ward.D2",

                              # For scores p-value calculation
                              ntests = 100, # number of shuffling events
                              seed = 22, # seed for random shuffling
                              pval_combine_method = "weighted_zmethod",

                              # For PCA
                              pca_comps_labs_invisible = c("quali"),
                              pca_pb_labs_invisible = c("quali"),
                              pca_sig_labs_invisible = c("quali"),

                              # Pseudobulk params
                              nvar_genes = 500,
                              black_list = NULL,
                              pca_n_hvg = 20,

                              ncores = round(parallelly::availableCores() / 2), # Otherwise can lead to memory issues
                              bparam = NULL,
                              progressbar = TRUE) {

  # Drop last ggplot2 from memory to prevent memory garbage collection from main environment
  ggplot2::set_last_plot(NULL)

  if (is.null(scoot_object) ||
      !inherits(scoot_object, "scoot")) {
    stop("Please provide a scoot class object or a count matrix.")
  }
  if (is.null(cluster_by)) {
    stop("Please provide a metadata column name to cluster by.")
  }

  if (!any(cluster_by %in% names(scoot_object@metadata))) {
    stop("cluster_by variables ", paste(cluster_by, collapse = ", "), " not found in metadata")
  } else if (!all(cluster_by %in% names(scoot_object@metadata))) {
    cluster_by <- cluster_by[cluster_by %in% names(scoot_object@metadata)]
    message("Only found ", paste(cluster_by, collapse = ", ") , " as grouping variable for scoot Object.")
  }

  if (!is.numeric(min_samples) ||
      !length(min_samples) == 1) {
    stop("Please provide a single numeric for min_samples.")
  }
  min_samples <- round(min_samples)
  if (min_samples < 5) {
    min_samples <- 5
    message("min_samples samples should not be less than 5 and was automatically set to 5")
  }

  if (!n_clust == "auto" &&
      (!is.numeric(n_clust) || !length(n_clust) == 1 || n_clust < 2)) {
    stop("n_clust must be a single numeric >= 2 or 'auto'.")
  }

  # Need to replace special characters
  colnames(scoot_object@metadata) <- make.names(colnames(scoot_object@metadata))
  cluster_by <- make.names(cluster_by)

  for (i in cluster_by) {
    if (length(unique(scoot_object@metadata[[i]])) == 1) {
      stop("All values are the same in cluster_by ", i, ". Please provide a metadata column with at least two different groups.")
    }
  }

  # Convert NA to factor level or drop
  for (i in cluster_by) {
    if (cluster_by_drop_na) {
      scoot_object@metadata <- scoot_object@metadata[!is.na(scoot_object@metadata[[i]]), ]
    } else {
      scoot_object@metadata[[i]] <- as.character(scoot_object@metadata[[i]])
      scoot_object@metadata[[i]][is.na(scoot_object@metadata[[i]])] <- "NA"
    }
    # cluster_by column must be factor
    scoot_object@metadata[[i]] <- as.factor(scoot_object@metadata[[i]])
  }

  # set parallelization parameters
  param <- set_parallel_params(ncores = ncores,
                               bparam = bparam,
                               progressbar = progressbar)

  # empty list to fill in the loop
  results <- list()


  # Process data ###############################################

  ## Pre-process data ###############################################
  data <- get_cluster_score_pre_proc(scoot_object,
                                     cluster_by_drop_na,
                                     min_samples,
                                     dist_method,
                                     nvar_genes = nvar_genes,
                                     cluster_by = 1,
                                     black_list)


  ## Unsupervised ###############################################
  # message("\nUnsupervised analysis\n")
  #
  # clust_args <- list(max_nc_adj = max_nc_adj,
  #                    fviz_nbclust_method = "silhouette",
  #                    nclust = "auto")
  #
  # message("\nProcessing cell type composition\n")
  #
  # comp_layers <- names(data[["composition"]])
  #
  # for (layer in comp_layers) {
  #   if (layer == comp_layers[1]) {
  #     df <- data[["composition"]][["feature_matrix"]]
  #
  #     cluster_labels <- do.call(clust_method, clust_args)
  #   }
  # }


  ## Supervised ###############################################
  message("\nSupervised analysis\n")

  get_scores_args <- list(scores = scores,
                          dist_method = dist_method,
                          knn_k = knn_k,
                          max_nc = max_nc,
                          n_clust = n_clust,
                          NbClust_method = NbClust_method,
                          ntests = ntests,
                          seed = seed)

  for (cluster_col in cluster_by) {
    message("Processing ", cluster_col)

    ## Process celltype composition ###############################################
    message("\nProcessing cell type composition\n")

    type <- "composition"

    comp_layers <- names(scoot_object@composition)

    for (layer in comp_layers) {

      get_scores_args[c("title", "invisible")] <- list(title = paste(cluster_col,
                                                                     stringr::str_to_title(type),
                                                                     layer),
                                                       invisible = pca_comps_labs_invisible)

      if (inherits(scoot_object@composition[[layer]], "data.frame")) {
        mat <- scoot_object@composition[[layer]][, c("celltype", "clr", "scoot_sample"), with = FALSE]

        if (cluster_by_drop_na) {
          mat <- mat %>% filter(scoot_sample %in% scoot_object@metadata[["scoot_sample"]])
        }
        mat <- mat %>%
          tidyr::pivot_wider(names_from = scoot_sample,
                             values_from = clr) %>%
          stats::na.omit() %>%
          tibble::column_to_rownames(var = "celltype") %>%
          t() %>%
          scale(center = TRUE,
                scale = FALSE)

        if (nrow(mat) >= min_samples) {

          if (is.null(batching))  {
            cluster_labels <- scoot_object@metadata %>%
              dplyr::filter(scoot_sample %in% row.names(mat)) %>%
              .[[cluster_col]]

            get_scores_args[c("feat_mat", "cluster_labels")] <- list(feat_mat = mat,
                                                                     cluster_labels = cluster_labels)

            results[[cluster_col]][[type]][[layer]] <- do.call(get_scores, get_scores_args)
          }

          if (!is.null(batching))  {
            for (b_var in batching) {
              if (b_var != cluster_col) {
                b_var_res_summary <- list()
                for (b in unique(scoot_object@metadata[[b_var]])) {
                  meta <- scoot_object@metadata %>%
                    dplyr::filter(get(b_var) == b) %>%
                    dplyr::filter(scoot_sample %in% row.names(mat))

                  cluster_labels <- meta[[cluster_col]]

                  m <- mat[ , row.names(mat) %in% as.character(meta[["scoot_sample"]])] %>%
                    scale(center = TRUE,
                          scale = FALSE)

                  get_scores_args[c("feat_mat", "cluster_labels")] <- list(feat_mat = mat,
                                                                           cluster_labels = cluster_labels)

                  results[[cluster_col]][[type]][[layer]][[b_var]][[b]] <- do.call(get_scores, get_scores_args)

                  for (score in scores) {
                    b_var_res_summary[[score]][["summary"]] <- c(
                      b_var_res_summary[[score]][["summary"]],
                      results[[cluster_col]][[type]][[layer]][[b_var]][[b]][["scores"]][[score]][["summary"]])
                    b_var_res_summary[[score]][["n"]] <- c(
                      b_var_res_summary[[score]][["n"]],
                      results[[cluster_col]][[type]][[layer]][[b_var]][[b]][["scores"]][[score]][["n"]])
                    b_var_res_summary[[score]][["p_value"]] <- c(
                      b_var_res_summary[[score]][["p_value"]],
                      results[[cluster_col]][[type]][[layer]][[b_var]][[b]][["scores"]][[score]][["p_value"]])
                  }
                }

                for (score in scores) {
                  b_var_res_summary[[score]][["summary"]] <- stats::na.omit(b_var_res_summary[[score]][["summary"]])
                  b_var_res_summary[[score]][["summary"]] <-
                    mean(b_var_res_summary[[score]][["summary"]])

                  # Requires the number of samples per batch, so run before summing n
                  p_values_not_na <- stats::na.omit(b_var_res_summary[[score]][["p_value"]])
                  p_values_not_na_len <- length(p_values_not_na)
                  if (p_values_not_na_len > 1) {
                    b_var_res_summary[[score]] <- combine_pvals(b_var_res_summary[[score]],
                                                                pval_combine_method)
                  } else if (p_values_not_na_len == 1) {
                    b_var_res_summary[[score]][["p_value"]] <- p_values_not_na
                  } else {
                    b_var_res_summary[[score]][["p_value"]] <- NULL
                    b_var_res_summary[[score]][["summary"]] <- NULL
                  }

                  b_var_res_summary[[score]][["n"]] <-
                    sum(b_var_res_summary[[score]][["n"]])
                }
                results[[cluster_col]][[type]][[layer]][[b_var]][["all"]][["scores"]] <- b_var_res_summary
              }
            }
          }
        } else {
          results[[cluster_col]][[type]][[layer]] <- NULL
        }

      } else if (is.list(scoot_object@composition[[layer]])) {
        ggplot2::set_last_plot(NULL)
        results[[cluster_col]][[type]][[layer]] <-
          BiocParallel::bplapply(
            X = names(scoot_object@composition[[layer]]),
            BPPARAM = param,
            function(i) {
              mat <- scoot_object@composition[[layer]][[i]][, c("celltype", "clr", "scoot_sample"), with = F]

              get_scores_args[["title"]] <- paste(cluster_col,
                                                  stringr::str_to_title(type),
                                                  layer,
                                                  i)

              if (cluster_by_drop_na) {
                mat <- mat %>% filter(scoot_sample %in% scoot_object@metadata[["scoot_sample"]])
              }
              mat <- mat %>%
                tidyr::pivot_wider(names_from = scoot_sample,
                                   values_from = clr) %>%
                tibble::column_to_rownames(var = "celltype") %>%
                t() %>%
                scale(center = TRUE,
                      scale = FALSE)

              if (nrow(mat) >= min_samples) {

                if (is.null(batching))  {
                  cluster_labels <- scoot_object@metadata %>%
                    dplyr::filter(scoot_sample %in% row.names(mat)) %>%
                    .[[cluster_col]]

                  if (nrow(mat) > 1) {
                    get_scores_args[c("feat_mat", "cluster_labels")] <-
                      list(feat_mat = mat, cluster_labels = cluster_labels)
                    res <- do.call(get_scores, get_scores_args)
                    return(res)
                  } else {
                    return(NULL)
                  }
                }

                if (!is.null(batching)) {
                  res <- list()
                  for (b_var in batching) {
                    if (b_var != cluster_col) {
                      b_var_res_summary <- list()
                      for (b in unique(scoot_object@metadata[[b_var]])) {
                        meta <- scoot_object@metadata %>%
                          dplyr::filter(get(b_var) == b) %>%
                          dplyr::filter(scoot_sample %in% row.names(mat))

                        cluster_labels <- meta[[cluster_col]]

                        m <- mat[ , row.names(mat) %in% as.character(meta[["scoot_sample"]])] %>%
                          scale(center = TRUE,
                                scale = FALSE)

                        if (nrow(m) > 1) {
                          get_scores_args[c("feat_mat", "cluster_labels")] <-
                            list(feat_mat = m, cluster_labels = cluster_labels)
                          res[[b_var]][[b]] <-
                            do.call(get_scores, get_scores_args)
                        }

                        for (score in scores) {
                          b_var_res_summary[[score]][["summary"]] <- c(
                            b_var_res_summary[[score]][["summary"]],
                            res[[b_var]][[b]][["scores"]][[score]][["summary"]])
                          b_var_res_summary[[score]][["n"]] <- c(
                            b_var_res_summary[[score]][["n"]],
                            res[[b_var]][[b]][["scores"]][[score]][["n"]])
                          b_var_res_summary[[score]][["p_value"]] <- c(
                            b_var_res_summary[[score]][["p_value"]],
                            res[[b_var]][[b]][["scores"]][[score]][["p_value"]])
                        }
                      }

                      for (score in scores) {
                        b_var_res_summary[[score]][["summary"]] <- stats::na.omit(b_var_res_summary[[score]][["summary"]])
                        b_var_res_summary[[score]][["summary"]] <-
                          mean(b_var_res_summary[[score]][["summary"]])

                        # Requires the number of samples per batch, so run before summing n
                        p_values_not_na <- stats::na.omit(b_var_res_summary[[score]][["p_value"]])
                        p_values_not_na_len <- length(p_values_not_na)
                        if (p_values_not_na_len > 1) {
                          b_var_res_summary[[score]] <- combine_pvals(b_var_res_summary[[score]],
                                                                      pval_combine_method)
                        } else if (p_values_not_na_len == 1) {
                          b_var_res_summary[[score]][["p_value"]] <- p_values_not_na
                        } else {
                          b_var_res_summary[[score]][["p_value"]] <- NULL
                          b_var_res_summary[[score]][["summary"]] <- NULL
                        }

                        b_var_res_summary[[score]][["n"]] <-
                          sum(b_var_res_summary[[score]][["n"]])
                      }
                      res[[b_var]][["all"]][["scores"]] <- b_var_res_summary
                    }
                  }
                  if(length(res) == 0) {
                    return(NULL)
                  } else {
                    return(res)
                  }
                }
              } else {
                return(NULL)
              }
              ggplot2::set_last_plot(NULL)
            }
          )
        names(results[[cluster_col]][[type]][[layer]]) <-
          names(scoot_object@composition[[layer]])
      }
    }


    ## Process pseudobulk ###############################################
    message("\nProcessing Pseudobulks\n")

    type <- "pseudobulk"

    pb_layers <- names(scoot_object@aggregated_profile[[type]])

    # Get black list
    if (is.null(black_list)) {
      black_list <- default_black_list
    }
    black_list <- unlist(black_list)

    for (layer in pb_layers) {
      ggplot2::set_last_plot(NULL)
      results[[cluster_col]][[type]][[layer]] <- BiocParallel::bplapply(
        X = names(scoot_object@aggregated_profile[[type]][[layer]]),
        BPPARAM = param,
        function(i) {
          mat <- as.matrix(scoot_object@aggregated_profile[[type]][[layer]][[i]])

          get_scores_args[c("title", "invisible", "select_var")] <-
            list(title = paste(cluster_col,
                               stringr::str_to_title(type),
                               layer,
                               i),
                 invisible = pca_pb_labs_invisible,
                 select_var = list(cos2 = pca_n_hvg))

          if (ncol(mat) >= min_samples) {

            # Remove black listed genes from the matrix
            mat <- mat[!row.names(mat) %in% black_list,]

            meta <- scoot_object@metadata %>%
              dplyr::filter(scoot_sample %in% colnames(mat))
            cluster_labels <- meta[[cluster_col]]
            if (cluster_by_drop_na) {
              mat <- mat[, colnames(mat) %in% as.character(meta[["scoot_sample"]])]
            }

            if (is.null(batching))  {
              if (length(unique(cluster_labels)) > 1) {
                mat <- DESeq2.normalize(matrix = mat,
                                        metadata = meta,
                                        cluster_by = cluster_col,
                                        nvar_genes = nvar_genes,
                                        black_list = black_list)

                get_scores_args[c("feat_mat", "cluster_labels")] <-
                  list(feat_mat = t(mat), cluster_labels = cluster_labels)
                res <- do.call(get_scores, get_scores_args)
                return(res)
              } else {
                return(NULL)
              }
            }

            if (!is.null(batching)) {
              res <- list()
              for (b_var in batching) {
                if (b_var != cluster_col) {
                  b_var_res_summary <- list()
                  for (b in unique(scoot_object@metadata[[b_var]])) {
                    met <- scoot_object@metadata %>%
                      dplyr::filter(get(b_var) == b) %>%
                      dplyr::filter(scoot_sample %in% colnames(mat))

                    cluster_labels <- met[[cluster_col]]

                    m <- mat[ , colnames(mat) %in% met[["scoot_sample"]]]

                    if (length(unique(cluster_labels)) > 1) {
                      m <- DESeq2.normalize(matrix = m,
                                            metadata = met,
                                            cluster_by = cluster_col,
                                            nvar_genes = nvar_genes,
                                            black_list = black_list)

                      if (nrow(m) > 1) {
                        get_scores_args[c("feat_mat", "cluster_labels")] <-
                          list(feat_mat = t(m), cluster_labels = cluster_labels)
                        res[[b_var]][[b]] <- do.call(get_scores, get_scores_args)
                      }
                    }
                    for (score in scores) {
                      b_var_res_summary[[score]][["summary"]] <- c(
                        b_var_res_summary[[score]][["summary"]],
                        res[[b_var]][[b]][["scores"]][[score]][["summary"]])
                      b_var_res_summary[[score]][["n"]] <- c(
                        b_var_res_summary[[score]][["n"]],
                        res[[b_var]][[b]][["scores"]][[score]][["n"]])
                      b_var_res_summary[[score]][["p_value"]] <- c(
                        b_var_res_summary[[score]][["p_value"]],
                        res[[b_var]][[b]][["scores"]][[score]][["p_value"]])
                    }
                  }

                  for (score in scores) {
                    b_var_res_summary[[score]][["summary"]] <- stats::na.omit(b_var_res_summary[[score]][["summary"]])
                    b_var_res_summary[[score]][["summary"]] <-
                      mean(b_var_res_summary[[score]][["summary"]])

                    # Requires the number of samples per batch, so run before summing n
                    p_values_not_na <- stats::na.omit(b_var_res_summary[[score]][["p_value"]])
                    p_values_not_na_len <- length(p_values_not_na)
                    if (p_values_not_na_len > 1) {
                      b_var_res_summary[[score]] <- combine_pvals(b_var_res_summary[[score]],
                                                                  pval_combine_method)
                    } else if (p_values_not_na_len == 1) {
                      b_var_res_summary[[score]][["p_value"]] <- p_values_not_na
                    } else {
                      b_var_res_summary[[score]][["p_value"]] <- NULL
                      b_var_res_summary[[score]][["summary"]] <- NULL
                    }

                    b_var_res_summary[[score]][["n"]] <-
                      sum(b_var_res_summary[[score]][["n"]])
                  }
                  res[[b_var]][["all"]][["scores"]] <- b_var_res_summary
                }
              }
              if(length(res) == 0) {
                return(NULL)
              } else {
                return(res)
              }
            }
          } else {
            return(NULL)
          }
          ggplot2::set_last_plot(NULL)
        }
      )
      names(results[[cluster_col]][[type]][[layer]]) <-
        names(scoot_object@aggregated_profile[[type]][[layer]])
    }
    get_scores_args[["select_var"]] <- NULL


    ## Process signatures ###############################################
    message("\nProcessing Signatures\n")

    type <- "signatures"

    comp_layers <- names(scoot_object@aggregated_profile[[type]])

    if (!is.null(comp_layers)) {
      for (layer in comp_layers) {
        cols <- colnames(scoot_object@aggregated_profile[[type]][[layer]])
        signatures <- cols[!cols %in% c("celltype", "scoot_sample")]

        results[[cluster_col]][[type]][[layer]] <- BiocParallel::bplapply(
          X = signatures,
          BPPARAM = param,
          function(i) {
            mat <- scoot_object@aggregated_profile[[type]][[layer]][, c("celltype", i, "scoot_sample"), with = FALSE]

            get_scores_args[c("title", "invisible")] <-
              list(title = paste(cluster_col,
                                 stringr::str_to_title(type),
                                 layer,
                                 i),
                   invisible = pca_sig_labs_invisible)

            if (cluster_by_drop_na) {
              mat <- mat %>% filter(scoot_sample %in% scoot_object@metadata[["scoot_sample"]])
            }
            mat <- mat %>%
              tidyr::pivot_wider(names_from = scoot_sample, values_from = i) %>%
              tibble::column_to_rownames(var = "celltype") %>%
              replace(is.na(.), 0) %>%
              t() %>%
              scale(center = TRUE,
                    scale = TRUE)

            if (nrow(mat) >= min_samples) {
              if (is.null(batching))  {

                # Drop columns containing NAN, caused by scaling zero variance columns
                mat <- mat[ , colSums(is.nan(mat)) == 0]

                cluster_labels <- scoot_object@metadata %>%
                  filter(scoot_sample %in% row.names(mat)) %>%
                  .[[cluster_col]]

                get_scores_args[c("feat_mat", "cluster_labels")] <-
                  list(feat_mat = mat, cluster_labels = cluster_labels)
                res <- do.call(get_scores, get_scores_args)
                return(res)
              }

              if (!is.null(batching)) {
                res <- list()
                for (b_var in batching) {
                  if (b_var != cluster_col) {
                    b_var_res_summary <- list()
                    for (b in unique(scoot_object@metadata[[b_var]])) {
                      meta <- scoot_object@metadata %>%
                        dplyr::filter(get(b_var) == b) %>%
                        dplyr::filter(scoot_sample %in% row.names(mat))

                      cluster_labels <- meta[[cluster_col]]

                      m <- mat[ , row.names(mat) %in% as.character(meta[["scoot_sample"]])] %>%
                        scale(center = TRUE,
                              scale = TRUE)

                      # Drop columns containing NAN, caused by scaling zero variance columns
                      m <- m[ , colSums(is.nan(m)) == 0]

                      if (ncol(m) >= min_samples) {
                        get_scores_args[c("feat_mat", "cluster_labels")] <-
                          list(feat_mat = m, cluster_labels = cluster_labels)
                        res[[b_var]][[b]] <- do.call(get_scores, get_scores_args)

                        for (score in scores) {
                          b_var_res_summary[[score]][["summary"]] <- c(
                            b_var_res_summary[[score]][["summary"]],
                            res[[b_var]][[b]][["scores"]][[score]][["summary"]])
                          b_var_res_summary[[score]][["n"]] <- c(
                            b_var_res_summary[[score]][["n"]],
                            res[[b_var]][[b]][["scores"]][[score]][["n"]])
                          b_var_res_summary[[score]][["p_value"]] <- c(
                            b_var_res_summary[[score]][["p_value"]],
                            res[[b_var]][[b]][["scores"]][[score]][["p_value"]])
                        }
                      }

                      for (score in scores) {
                        b_var_res_summary[[score]][["summary"]] <- stats::na.omit(b_var_res_summary[[score]][["summary"]])
                        b_var_res_summary[[score]][["summary"]] <-
                          mean(b_var_res_summary[[score]][["summary"]])

                        # Requires the number of samples per batch, so run before summing n
                        p_values_not_na <- stats::na.omit(b_var_res_summary[[score]][["p_value"]])
                        p_values_not_na_len <- length(p_values_not_na)
                        if (p_values_not_na_len > 1) {
                          b_var_res_summary[[score]] <- combine_pvals(b_var_res_summary[[score]],
                                                                      pval_combine_method)
                        } else if (p_values_not_na_len == 1) {
                          b_var_res_summary[[score]][["p_value"]] <- p_values_not_na
                        } else {
                          b_var_res_summary[[score]][["p_value"]] <- NULL
                          b_var_res_summary[[score]][["summary"]] <- NULL
                        }

                        b_var_res_summary[[score]][["n"]] <-
                          sum(b_var_res_summary[[score]][["n"]])
                      }
                      res[[b_var]][["all"]][["scores"]] <- b_var_res_summary
                    }
                  }
                }
                if(length(res) == 0) {
                  return(NULL)
                } else {
                  return(res)
                }
              }
            } else {
              return(NULL)
            }
            ggplot2::set_last_plot(NULL)
          }
        )
        names(results[[cluster_col]][[type]][[layer]]) <- signatures
      }
    }
    ggplot2::set_last_plot(NULL)
    invisible(gc())
  }

  results[["data"]] <- data

  # Save user set parameters for summarize_cluster_scores
  results[["params"]][["cluster_by"]] <- cluster_by
  results[["params"]][["batching"]] <- batching
  results[["params"]][["scores"]] <- scores

  return(results)
}




#' Summarize scores and plot heat map
#'
#' @param data Output from get_cluster_scores
#' @param topN Integer indicating number of topN highest scoring (most discriminating) features ranked for each score
#' @param create_plots Boolean indicating whether to create and show plots or not
#' @param p_adjustment Whether to adjust p-value columns or not
#' @param p_adjust_method Method for adjusting p-values (see stats::p.adjust for methods)
#' @param p_value_cutoff p-value (mean of all p-value columns) cutoff to filter out non-significant results
#'
#' @importFrom dplyr mutate filter
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom stats p.adjust na.omit quantile
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom ggplotify as.ggplot
#'
#' @return Average silhouette widths per grouping variable. Optionally, a heatmap plot for visualization
#' @export summarize_cluster_scores
#'


summarize_cluster_scores <- function(data = NULL,
                                     topN = NULL,
                                     create_plots = TRUE,
                                     p_adjustment = TRUE,
                                     p_adjust_method = "fdr",
                                     p_value_cutoff = 0.05) {

  show_variables <- c(".summary", ".n", ".p_value")

  if (is.null(data)) {
    message("Please provide scores object (output from get_cluster_scores)")
  }

  cluster_by <- data[["params"]][["cluster_by"]]
  batching <- data[["params"]][["batching"]]
  scores <- data[["params"]][["scores"]]

  data[["params"]] <- NULL # Remove params

  data_conts <- unlist(data)
  df_pars <- list()

  for (par in show_variables) {
    data_conts_temp <- data_conts[endsWith(names(data_conts), par)]

    if (!is.null(batching)) {
      data_conts_temp <- data_conts_temp[grepl(".all.scores.",
                                               names(data_conts_temp))]
      names(data_conts_temp) <- gsub("all.", "", names(data_conts_temp))
    }

    if (length(data_conts_temp) == 0) {
      stop("Error happened at ", show_variables, ". No values left after filtering.")
    } else {
      target_string <- paste0(scores, collapse = "|")
      target_string <- paste0("(", target_string, "+)")
      data_conts_names_split <- strsplit(gsub(target_string,"~\\1",
                                              names(data_conts_temp)), "~")

      df <- data.frame(t(do.call(cbind, data_conts_names_split))) %>%
        dplyr::mutate(X3 = sapply(data_conts_temp, "[[", 1)) %>%
        tidyr::pivot_wider(names_from = "X1",
                           values_from = "X3") %>%
        tibble::column_to_rownames(var = "X2") %>%
        t() %>%
        as.data.frame()

      row.names(df) <- gsub(".scores.",
                            "",
                            row.names(df))

      df_pars[[par]] <- df
    }
  }

  # Check if all columns with sample numbers are equal
  df <- as.matrix(df_pars[[".n"]])
  all_n_cols_equal <- all(apply(df, 2, identical, df[,1]))
  if (all_n_cols_equal) {
    df_pars[[".n"]] <- df_pars[[".n"]][, 1, drop = FALSE]
    colnames(df_pars[[".n"]]) <- "n"
  }

  df <- merge(df_pars[[show_variables[1]]],
              df_pars[[show_variables[2]]],
              by = "row.names",
              all = TRUE)
  row.names(df) <- df$Row.names
  df$Row.names <- NULL

  df <- merge(df,
              df_pars[[show_variables[3]]],
              by = "row.names",
              all = TRUE)
  row.names(df) <- df$Row.names
  df$Row.names <- NULL

  df_cluster_by_list <- list()
  for (c in cluster_by) {
    df_cluster_by_list[[c]] <- df %>%
      dplyr::filter(row.names(.) %>%
                      startsWith(c))
    row.names(df_cluster_by_list[[c]]) <- gsub(paste0(c, "."),
                                               "",
                                               row.names(df_cluster_by_list[[c]]))
    # Remove batching variable name
    if (!is.null(batching)) {
      batch_names <- batching[!batching %in% c]
      batch_names <- batch_names %>%
        paste0(".", .) %>%
        paste(collapse = "|")
      row.names(df_cluster_by_list[[c]]) <- gsub(batch_names,
                                                 "",
                                                 row.names(df_cluster_by_list[[c]]))
    }

    # Adjust p-values and filter
    if (p_adjustment) {
      p_val_cols <- which(grepl(".p_value",
                                colnames(df_cluster_by_list[[c]])))
      for (i in p_val_cols) {
        df_cluster_by_list[[c]][, i] <-
          stats::p.adjust(df_cluster_by_list[[c]][, i],
                          method = p_adjust_method)
      }

      df_cluster_by_list[[c]] <-
        df_cluster_by_list[[c]][rowMeans(df_cluster_by_list[[c]][, p_val_cols]) <= p_value_cutoff, ]

      df_cluster_by_list[[c]] <- stats::na.omit(df_cluster_by_list[[c]])

      if (nrow(df_cluster_by_list[[c]]) == 0) {
        df_cluster_by_list[[c]] <- NULL
        message(paste("For ", c, " no separation was found after p-value cutoff. You can try to set it higher."))

        next
      }
    }

    # Get topN scored results
    if (!is.null(topN) &&
        is.numeric(topN) &&
        length(topN) == 1) {
      summary_score_cols <- which(grepl(".summary",
                                        colnames(df_cluster_by_list[[c]])))
      topN_rownames <- c()
      for (i in summary_score_cols) {
        topN_rownames_i <- row.names(df_cluster_by_list[[c]])[
          order(df_cluster_by_list[[c]][, i], decreasing = TRUE)][1:topN]
        topN_rownames <- c(topN_rownames,
                           topN_rownames_i)
      }
      df_cluster_by_list[[c]] <- df_cluster_by_list[[c]][unique(topN_rownames), ]
    }
  }

  # Remove NULL elements
  df_cluster_by_list <- Filter(Negate(is.null), df_cluster_by_list)

  # Check if all df_cluster_by_list were NULL
  if (length(df_cluster_by_list) == 0) {

    message("No significant separation found for cluster_by provided")

  } else {

    if (create_plots) {
      df_list <- list()
      plot_list <- list()

      for (n in names(df_cluster_by_list)) {
        df <- stats::na.omit(df_cluster_by_list[[n]])

        n_breaks <- min(100, dim(df)[1] * dim(df)[2])
        quantiles <- seq(0, 1, 1/n_breaks)

        # Check that variance is not zero
        nonzero_var_cols <- unlist(lapply(df, function(x) !length(unique(x))==1))
        breaks <- df[, nonzero_var_cols] %>%
          scale() %>%
          stats::quantile(., quantiles) %>%
          unique()

        color_breaks <- length(breaks)

        color <- grDevices::colorRampPalette(
          rev(RColorBrewer::brewer.pal(n = 7,
                                       name = "RdYlBu")))(color_breaks)

        if (nrow(df) > 1) {
          scale <- "column"
        } else {
          scale <- "none"
        }

        df_pval_invert <- df
        df_pval <- df_pval_invert[, grepl(".p_value", names(df))]
        df_pval <- log10(1 / (df_pval + 1e-16))
        df_pval_invert[, grepl(".p_value", names(df))] <- df_pval

        plot_list[[n]] <- pheatmap::pheatmap(df_pval_invert,
                                             main = n,
                                             angle_col = 315,
                                             scale = scale,
                                             display_numbers = round(df, 3),
                                             number_color = "black",
                                             color = color,
                                             breaks = breaks,
                                             legend_breaks = 0,
                                             legend_labels = "",
                                             cluster_cols = FALSE,
                                             cluster_rows = FALSE,
                                             silent = TRUE)[[4]]
      }
      df_cluster_by_list[["plots"]][["plot_list"]] <- plot_list

      g <- gridExtra::grid.arrange(
        gridExtra::arrangeGrob(grobs = plot_list,
                               ncol = length(plot_list))
      )
      df_cluster_by_list[["plots"]][["summary_plot"]] <- ggplotify::as.ggplot(g)
    }

    return(df_cluster_by_list)
  }
}




#' Render plots summarizing celltype proportions and distribution in samples
#'
#'
#' @param obj_list List of Seurat objects
#' @param annot_col Metadata column(s) containing the cell type annotations
#' @param bottom_mar Adjustable bottom margin for long sample names

#' @importFrom stats setNames

#' @return Get percentage of not annotated cells per sample and plot it.
#' @export nas_per_sample
#'

nas_per_sample <- function(obj_list = NULL,
                           annot_col = c("scGate_multi"),
                           return_plot = TRUE,
                           bottom_mar = 10.2) {
  if (is.null(obj_list) &
      !is.list(obj_list) &
      !all(lapply(obj_list, inherits, "Seurat"))) {
    stop("Please provide a list of seurat objects")
  }

  na_list <- list()

  for (col in annot_col) {
    na_perc_per_sample <- c()
    for (i in names(obj_list)) {
      if (col %in% names(obj_list[[i]]@meta.data)) {
        percs <- prop.table(table(obj_list[[i]]@meta.data[[col]], useNA = "ifany"))*100
        nas_perc <- unname(percs[is.na(names(percs))])
        na_perc_per_sample <- c(na_perc_per_sample, stats::setNames(nas_perc, i))
      } else {
        stop(paste(col, " not found in obj_list item ", i))
      }
    }
    par(mar = c(bottom_mar, 4.1, 4.1, 2.1))
    if (return_plot) {
      barplot(na_perc_per_sample,
              main = col,
              ylab = "Percentage of NA values",
              las=2)
    }
    na_list[[col]] <- na_perc_per_sample
  }
  return(na_list)
}




#' Render bar plots summarizing celltype proportions and distribution in samples and groups
#'
#'
#' @param scoot_object A scoot class object (typically after applying merge_scoot_objects onto a list of scootObjects)
#' @param layer Default "layer_1" if you have one cell type annotation layer in your scoot_object. Alternatively "layer_2" etc. if you have multiple layers of annotation depths.
#' @param return_plot_to_var Optionally, you can save the ggplots to a variable if you would like to further modify and adapt the plots on your own.
#' @param facet_by This allows you to pass a metadata column name present in your scoot_object$metadata to show your samples in facets with ggplot facet_grid, for example by "condition".

#' @importFrom ggplot2 ggplot aes geom_bar theme element_text ggtitle facet_grid
#' @importFrom stats reformulate
#' @importFrom patchwork wrap_plots

#' @return Plotting function to show the cell type composition from scoot object across different samples.
#' @export composition_barplot
#'

composition_barplot <- function(scoot_object = NULL,
                                sample_col = NULL,
                                layer = "layer_1",
                                return_plot_to_var = FALSE,
                                facet_by = NULL) {

  # Need to replace special characters
  colnames(scoot_object@metadata) <- make.names(colnames(scoot_object@metadata))

  sample_col <- "scoot_sample"

  if (is.null(scoot_object)) {
    stop("Please provide input scoot_object")
  }
  if (!length(layer) == 1 || !is.character(layer)) {
    stop("Please provide one character string for layer parameter")
  }
  if (!is.null(facet_by)) {
    if (!is.character(facet_by)) {
      stop("Please provide a character string or a vector of character strings for the facet_by parameter")
    }
    facet_by <- make.names(facet_by)
    facet_by_in_colnames <- facet_by %in% names(scoot_object@metadata)
    if (!all(facet_by_in_colnames)) {
      facet_by_not_in_colnames <- facet_by[!facet_by_in_colnames]
      stop(paste("facet_by ", facet_by_not_in_colnames, " not found in scoot_object@metadata column names"))
    }
  }


  comps <- scoot_object@composition[[layer]]
  meta <- scoot_object@metadata

  if (!is.null(facet_by)) {
    facet_by_reformulate <- reformulate(facet_by)
  }

  ylab <- "Relative abundance (%)"

  if (is.data.frame(comps)) {
    comp <- merge(comps, meta[, c(sample_col, facet_by), drop=FALSE], by = sample_col)

    p <- ggplot(comp, aes(x = scoot_sample, y = freq, fill = celltype)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlab("") + ylab(ylab)

    if (!is.null(facet_by)) {

      p <- p + facet_grid(facet_by_reformulate,
                          space  = "free",
                          scales = "free")
    }

    if (return_plot_to_var) {
      return(p)
    } else {
      print(p)
    }
  }

  else {
    p_list <- list()
    for (ct in names(comps)) {
      comp <- scoot_object@composition[[layer]][[ct]]
      comp <- merge(comp, meta[, c(sample_col, facet_by), drop=FALSE], by = sample_col)

      p_list[["plot_list"]][[ct]] <- ggplot(comp, aes(x = scoot_sample, y = freq, fill = celltype)) +
        geom_bar(stat = "identity") +
        xlab("") + ylab(ylab) +
        theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        ggtitle(ct)

      if (!is.null(facet_by)) {
        p_list[["plot_list"]][[ct]] <- p_list[["plot_list"]][[ct]] +
          facet_grid(facet_by_reformulate,
                     space  = "free",
                     scales = "free")
      }

    }
    p_list[["arranged_plots"]] <- patchwork::wrap_plots(p_list[["plot_list"]])

    if (return_plot_to_var) {
      return(p_list)
    } else {
      print(p_list[["arranged_plots"]])
    }
  }
}


#' Render box plots summarizing celltype proportions and distribution in samples and groups
#'
#'
#' @param scoot_object A scoot class object (typically after applying merge_scoot_objects onto a list of scootObjects)
#' @param plot_var Column in the scoot_object$composition: either "freq" for cell type relative abundance in percent or "clr" (for centered log-ratio transformed). Default: "clr" as it is better suited for statistical analysis and is better able to also show low abundant cell types.
#' @param layer Default "layer_1" if you have one cell type annotation layer in your scoot_object. Alternatively "layer_2" etc. if you have multiple layers of annotation depths.
#' @param return_plot_to_var Optionally, you can save the ggplots to a variable if you would like to further modify and adapt the plots on your own.
#' @param group_by This allows you to pass a metadata column name present in your scoot_object$metadata to show your samples in groups, for example by "condition".
#' @param facet_by This allows you to pass a metadata column name present in your scoot_object$metadata to show your samples in facets with ggplot facet_grid, for example by "condition".
#' @param pval_method Specify method how to calculate the p-value. Wilcoxon test is recommended as compositional data might not fullfill the assumption of a gaussian distribution. For alternatives, see documentation of ggpubr::stat_pwc.
#' @param p_adjust_method Method for adjusting p-values (see ggpubr::stat_pwc for available methods)
#' @param palette Choose a palette of your liking. For available palettes, see ggsci package. Default: "lancet"
#' @param legend_position Where to put the legend. Possible options: "top", "right", "bottom", "left"

#' @importFrom ggplot2 ggplot aes geom_boxplot theme element_text ggtitle facet_grid position_jitterdodge
#' @importFrom ggpubr stat_pwc ggboxplot
#' @importFrom stats reformulate
#' @importFrom patchwork wrap_plots

#' @return Plotting function to show the cell type composition from scoot object across different samples.
#' @export composition_boxplot
#'

composition_boxplot <- function(scoot_object = NULL,
                                plot_var = "clr",
                                layer = "layer_1",
                                return_plot_to_var = FALSE,
                                group_by = NULL,
                                facet_by = NULL,
                                pval_method = "wilcox.test",
                                p_adjust_method = "BH",
                                palette = "lancet",
                                legend_position = "right") {

  # Need to replace special characters
  colnames(scoot_object@metadata) <- make.names(colnames(scoot_object@metadata))

  sample_col <- "scoot_sample"

  if (is.null(scoot_object)) {
    stop("Please provide input scoot_object")
  }
  if (!length(plot_var) == 1 ||
      !is.character(plot_var) ||
      !plot_var %in% c("freq", "clr")) {
    stop("Please provide one character string for plot_var, either 'freq' or 'clr'")
  }
  if (!length(layer) == 1 || !is.character(layer)) {
    stop("Please provide one character string for layer parameter")
  }
  if (!is.null(group_by)) {
    if (!length(group_by) == 1 || !is.character(group_by)) {
      stop("Please provide one character string for the group_by parameter")
    }
    group_by <- make.names(group_by)
    group_by_gg <- sym(group_by)
    nr_of_boxplots <- scoot_object@metadata[[group_by]] %>%
      unique() %>%
      length() %>%
      "*"(1.5) %>%
      round()
    scoot_object@metadata[group_by] <- lapply(scoot_object@metadata[group_by], as.factor)
  }
  if (!is.null(facet_by)) {
    if (!is.character(facet_by)) {
      stop("Please provide a character string or a vector of character strings for the facet_by parameter")
    }
    facet_by <- make.names(facet_by)

    facet_by_in_colnames <- facet_by %in% names(scoot_object@metadata)
    if (!all(facet_by_in_colnames)) {
      facet_by_not_in_colnames <- facet_by[!facet_by_in_colnames]
      stop(paste("facet_by ", facet_by_not_in_colnames, " not found in scoot_object@metadata column names"))
    }
  }

  comps <- scoot_object@composition[[layer]]
  meta <- scoot_object@metadata

  plot_var_gg <- sym(plot_var)

  if (plot_var == "clr") {
    ylab <- "Relative abundance (clr)"
  }
  if (plot_var == "freq") {
    ylab <- "Relative abundance (%)"
  }
  if (plot_var == "freq") {
    ylab <- "Cell counts"
  }

  if (is.data.frame(comps)) {
    comp <- merge(comps, meta[, c(sample_col, group_by, facet_by), drop=FALSE], by = sample_col)

    # Need to check if group_by is NULL
    # Due to a presumed bug, if group_by is passed as variable to ggboxplot, even if it is assigned NULL, it throws an error
    if (is.null(group_by)) {
      p <- ggboxplot(comp,
                     x = "celltype",
                     y = plot_var,
                     xlab = "",
                     ylab = ylab,
                     outlier.shape = NA,
                     palette = palette,
                     facet.by = facet_by,
                     legend = legend_position) +
        geom_jitter(width = 0.2, size = 1)
    } else {
      p <- ggboxplot(comp,
                     x = "celltype",
                     y = plot_var,
                     xlab = "",
                     ylab = ylab,
                     color = group_by,
                     outlier.shape = NA,
                     palette = palette,
                     facet.by = facet_by,
                     legend = legend_position) +
        geom_jitter(mapping = aes(color = !!group_by_gg), position=position_jitterdodge(jitter.width = 1/nr_of_boxplots), size = 1) +
        stat_pwc(aes(group = !!group_by_gg),
                 label = "p.signif",
                 method = pval_method,
                 p.adjust.method = p_adjust_method,
                 p.adjust.by = "panel",
                 tip.length = 0,
                 hide.ns = TRUE)
    }

    p <- p +
      theme(axis.text.x = element_text(angle = 45, hjust=1)) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

    if (return_plot_to_var) {
      return(p)
    } else {
      print(p)
    }
  }

  else {
    p_list <- list()
    for (ct in names(comps)) {
      comp <- scoot_object@composition[[layer]][[ct]]
      comp <- merge(comp, meta[, c(sample_col, group_by, facet_by), drop=FALSE], by = sample_col)

      # Need to check if group_by is NULL
      # Due to a presumed bug, if group_by is passed as variable to ggboxplot, even if it is assigned NULL, it throws an error
      if (is.null(group_by)) {
        p_list[["plot_list"]][[ct]] <- ggboxplot(comp,
                                                 x = "celltype",
                                                 y = plot_var,
                                                 xlab = "",
                                                 ylab = ylab,
                                                 outlier.shape = NA,
                                                 palette = palette,
                                                 facet.by = facet_by,
                                                 legend = legend_position) +
          geom_jitter(width = 0.2, size = 1)
      } else {
        p_list[["plot_list"]][[ct]] <- ggboxplot(comp,
                                                 x = "celltype",
                                                 y = plot_var,
                                                 xlab = "",
                                                 ylab = ylab,
                                                 color = group_by,
                                                 outlier.shape = NA,
                                                 palette = palette,
                                                 facet.by = facet_by,
                                                 legend = legend_position) +
          geom_jitter(mapping = aes(color = !!group_by_gg), position=position_jitterdodge(jitter.width = 1/nr_of_boxplots), size = 1) +
          stat_pwc(aes(group = !!group_by_gg),
                   label = "p.signif",
                   method = pval_method,
                   p.adjust.method = p_adjust_method,
                   p.adjust.by = "panel",
                   tip.length = 0,
                   hide.ns = TRUE)
      }

      p_list[["plot_list"]][[ct]] <- p_list[["plot_list"]][[ct]] +
        theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

    }
    p_list[["arranged_plots"]] <- patchwork::wrap_plots(p_list[["plot_list"]])

    if (return_plot_to_var) {
      return(p_list)
    } else {
      print(p_list[["arranged_plots"]])
    }
  }
}
