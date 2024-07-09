# Function get.HiTObject
#' Get a HiT object summarizing cell type classification and aggregated profiles.
#'
#'
#'
#' @param object A Seurat object or a list of Seurat objects
#' @param group.by List or vector with one or multiple Seurat object metadata columns with cell type predictions to group by (e.g. layer 1 cell type classification)
#' @param split.by A Seurat object metadata column to split by compositional data (e.g. sample names) \link{get.celltype.composition}
#' @param min.cells.composition Parameter for internal \link{get.celltype.composition}. Minimum number of cells annotated in a group.by parameter to render the cell type composition. If the number of cells annotated for a certain group.by parameter is less, compositional data will not be rendered. Default value is 10.
#' @param min.cells.aggregated Parameter for internal \link{get.aggregated.profile} and \link{get.aggregated.signature}. Minimum number of cells per sample and cell type to aggregate data for pseudobulk or aggregated signatures. Aggregated data for cell types with less than that value will not be returned in aggregated data. Default value is 10.
#' @param name.additional.signatures Names of additional signatures as found in object metadata to take into account.
#' @param useNA logical whether to return aggregated profile for NA (undefined) cell types, default is FALSE.
#' @param clr_zero_impute_perc Parameter for internal \link{get.celltype.composition}.
#' @param layers_links Parameter for internal \link{get.celltype.composition}
#' @param assay Parameter for internal \link{get.aggregated.profile}.
#' @param layer1_link Column of metadata linking layer1 prediction (e.g. scGate ~ CellOntology_ID) in order to perform subsetting for second layer classification.
#'
#' @importFrom methods setClass new
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallelly availableCores
#'
#' @import SeuratObject
#'
#' @return HiT object summarizing cell type classification and aggregated profiles.
#' @export get.HiTObject
#'


get.HiTObject <- function(object,
                          group.by = list("layer1" = c("scGate_multi"),
                                          "layer2" = c("functional.cluster")
                          ),
                          split.by = NULL,
                          min.cells.composition = 10,
                          min.cells.aggregated = 10,
                          name.additional.signatures = NULL,
                          useNA = FALSE,
                          clr_zero_impute_perc = 1,
                          layers_links = c("scGate_multi" = "functional.cluster"),
                          assay = "RNA",
                          layer1_link = "CellOntology_ID",

                          ncores = parallelly::availableCores() - 2,
                          bparam = NULL,
                          progressbar = TRUE) {

  args <- list(group.by,
               split.by,
               min.cells.composition,
               min.cells.aggregated,
               name.additional.signatures,
               useNA,
               clr_zero_impute_perc,
               layers_links,
               assay,
               layer1_link)

  # if object is a single Seurat object, turn into a list
  if (!is.list(object)) {
    hit.list <- get.HiTObject.helper(object, group.by,
                                     split.by,
                                     min.cells.composition,
                                     min.cells.aggregated,
                                     name.additional.signatures,
                                     useNA,
                                     clr_zero_impute_perc,
                                     layers_links,
                                     assay,
                                     layer1_link)
  } else {
    # set parallelization parameters
    param <- set_parallel_params(ncores = ncores,
                                 bparam = bparam,
                                 progressbar = progressbar)

    hit.list <- BiocParallel::bplapply(X = object,
                                       BPPARAM = param,
                                       function(x) {
                                         do.call(get.HiTObject.helper,
                                                 c(x, args)
                                         )
                                       })
  }

  return(hit.list)
}

get.HiTObject.helper <- function(object,
                                 group.by = list("layer1" = c("scGate_multi"),
                                                 "layer2" = c("functional.cluster")
                                 ),
                                 split.by = NULL,
                                 min.cells.composition = 10,
                                 min.cells.aggregated = 10,
                                 name.additional.signatures = NULL,
                                 useNA = FALSE,
                                 clr_zero_impute_perc = 1,
                                 layers_links = c("scGate_multi" = "functional.cluster"),
                                 assay = "RNA",
                                 layer1_link = "CellOntology_ID") {

  if (is.null(object)) {
    stop("Please provide a Seurat object")
  } else if (!inherits(object, "Seurat")) {
    stop("Not Seurat object included, cannot be processed.\n")
  }


  if (is.null(group.by)) {
    stop("Please provide at least one grouping variable")
  }

  if (!is.list(group.by)) {
    group.by <- as.list(group.by)
  }

  if (!any(group.by %in% names(object@meta.data))) {
    stop("Group.by variables ", paste(group.by, collapse = ", "), " not found in metadata")
  } else if (!all(group.by %in% names(object@meta.data))) {
    group.by <- group.by[group.by %in% names(object@meta.data)]
    message("Only found ", paste(group.by, collapse = ", ") , " as grouping variable for HiT Object.")
  }

  for (v in seq_along(group.by)) {
    if (is.null(names(group.by)[[v]]) || is.na(names(group.by)[[v]])) {
      names(group.by)[[v]] <- paste0("layer", v)
    }
  }


  # Make list of list for layers for predictions slot
  pred.list <- list()
  for (a in names(group.by)) {
    pred.list[[a]] <- list()
    pred.list[[a]][[group.by[[a]]]] <- object@meta.data[,group.by[[a]], drop = FALSE]
  }



  # Run extended if object got scGate info
  # Extract values from misc slot from object
  if (any(group.by == "scGate_multi")) {
    layer.scgate <- names(group.by[group.by == "scGate_multi"])
    if (is.null(name.additional.signatures)) {
      name.additional.signatures <- object@misc$layer1_param$additional.signatures
    }
    scgate.models <- object@misc$layer1_param$scGate_models


    if (!is.null(name.additional.signatures)) {
      sig <- grep(paste(name.additional.signatures, collapse = "|"),
                  names(object@meta.data), value = TRUE)
      sig.df <- object@meta.data[, sig]
      pred.list[[layer.scgate]][["additional_signatures"]] <- sig.df
    }

    if (!is.null(scgate.models)) {
      scgate <- grep(paste(paste0(scgate.models, "$"), collapse = "|"),
                     names(object@meta.data), value = TRUE)
      scgate.df <- object@meta.data[, scgate]
      pred.list[[layer.scgate]][["scGate_is.pure"]] <- scgate.df
    }

    if (!is.null(name.additional.signatures) && !is.null(scgate.models)) {
      ucell <- names(object@meta.data)[!names(object@meta.data) %in% c(sig, scgate)] %>%
        grep("_UCell$", ., value = TRUE)

      ucell.df <- object@meta.data[, ucell]

      pred.list[[layer.scgate]][["UCell_scores"]] <- ucell.df
    }
  }

  # Run extended if object got ProjecTILs info
  # Extract values from misc slot from object
  if (any(group.by == "functional.cluster")) {
    layer.pt <- names(group.by[group.by == "functional.cluster"])
    pred.list[[layer.pt]][["functional.cluster"]] <- object@meta.data[,"functional.cluster", drop = FALSE]
    pred.list[[layer.pt]][["functional.cluster.conf"]] <- object@meta.data[,"functional.cluster.conf", drop = FALSE]
  }


  # Compute proportions
  # message("\nComputing cell type composition...\n")

  comp.prop <- get.celltype.composition(object,
                                        group.by.composition = group.by,
                                        min.cells.composition = min.cells.composition,
                                        split.by = split.by,
                                        useNA = useNA,
                                        clr_zero_impute_perc = clr_zero_impute_perc,
                                        layer1_link = layer1_link,
                                        layers_links = layers_links)
  # Compute avg expression
  # message("\nComputing aggregated profile...\n")

  avg.expr <- get.aggregated.profile(object,
                                     group.by.aggregated = group.by,
                                     min.cells.aggregated = min.cells.aggregated,
                                     assay = assay,
                                     useNA = useNA)

  aggr.signature <- get.aggregated.signature(object,
                                             group.by.aggregated = group.by,
                                             min.cells.aggregated = min.cells.aggregated,
                                             name.additional.signatures = name.additional.signatures,
                                             useNA = useNA)

  hit <- methods::new("HiT",
                      metadata = object@meta.data,
                      predictions = pred.list,
                      aggregated_profile = list("pseudobulk" = avg.expr,
                                                "signatures" = aggr.signature),
                      composition = comp.prop
  )
  return(hit)
}


#' Calculate cell type composition or frequencies
#'
#' @param object A seurat object
#' @param group.by.composition The Seurat object metadata column(s) containing celltype annotations (provide as character vector, containing the metadata column name(s))
#' @param split.by A Seurat object metadata column to split by (e.g. sample names)
#' @param min.cells.composition Set a minimum threshold for number of cells to calculate relative abundance (e.g. less than 10 cells -> no relative abundance will be calculated)
#' @param useNA Whether to include not annotated cells or not (labelled as "NA" in the group.by.composition). Can be defined separately for each group.by.composition (provide single boolean or vector of booleans)
#' @param layers_links Named vector indicating the relation of multiple layer classification, default is \code{c("scGate_multi", "functional.cluster")}. Name of the vector element ought to be layer1 and element layer2. This parameter is used to compute relative compositional data of layer2 within layer1 classes.
#' @param clr_zero_impute_perc To calculate the clr-transformed relative abundance ("clr_freq"), zero values are not allowed and need to be imputed (e.g. by adding a pseudo cell count). Instead of adding a pseudo cell count of flat +1, here a pseudo cell count of +1% of the total cell count will be added to all cell types, to better take into consideration the relative abundance ratios (e.g. adding +1 cell to a total cell count of 10 cells would have a different, i.e. much larger effect, than adding +1 to 1000 cells).
#' @param layer1_link Column of metadata linking layer1 prediction (e.g. scGate ~ CellOntology_ID) in order to perform subsetting for second layer classification.

#' @importFrom Hotelling clr
#' @importFrom dplyr group_by summarize filter ungroup mutate select left_join n coalesce bind_rows across all_of
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom rrapply rrapply
#'
#' @return Cell type compositions as a list of data.frames containing cell counts, relative abundance (freq) and clr-transformed freq (freq_clr), respectively.
#' @export get.celltype.composition
#'
#' @examples
#' devtools::install_github('satijalab/seurat-data')
#' library(SeuratData)
#' options(timeout = max(300, getOption("timeout")))
#' InstallData("panc8")
#' data("panc8")
#' panc8 = UpdateSeuratObject(object = panc8)
#' # Calculate overall composition
#' celltype.compositions.overall <- get.celltype.composition(object = panc8, group.by.composition = "celltype")
#'
#' # Calculate sample-wise composition

#' celltype.compositions.sample_wise <- get.celltype.composition(object = panc8, group.by.composition = "celltype", split.by = "orig.ident")
#'
get.celltype.composition <- function(object = NULL,
                                     group.by.composition = NULL,
                                     split.by = NULL,
                                     min.cells.composition = 10,
                                     useNA = FALSE,
                                     layers_links = c("scGate_multi" = "functional.cluster"),
                                     clr_zero_impute_perc = 1,
                                     layer1_link = "CellOntology_ID") {

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

  if (is.null(group.by.composition)) {
    stop("Please specificy a group.by.composition variable")
  }

  # Assess wheter split.by variable is in metadata
  if (!is.null(split.by) &&
      !split.by %in% names(meta.data)) {
    stop("Split.by variable not found in meta.data!\n")
  }

  if (length(useNA) != 1 &&
      length(useNA) != length(group.by.composition)) {
    stop("useNA variable must be of length 1 or the same length as group.by.composition (group.by)")
  }

  # convert group.by.composition to list if not
  if (!is.list(group.by.composition)) {
    group.by.composition <- as.list(group.by.composition)
  }


  # Rename group.by.composition if not indicated
  for (v in seq_along(group.by.composition)) {
    if (is.null(names(group.by.composition)[[v]]) ||
        is.na(names(group.by.composition)[[v]])) {
      names(group.by.composition)[[v]] <- group.by.composition[[v]]
    }

  }



  # Keep only grouping variables in metadata
  if (!any(group.by.composition %in% names(meta.data))) {
    stop("Group.by variables ", paste(group.by.composition, collapse = ", "), " not found in metadata")
  } else if (!all(group.by.composition %in% names(meta.data))) {
    group.present <- group.by.composition %in% names(meta.data)
    # accommodate useNA vector
    if (length(useNA) == length(group.by.composition)) {
      useNA <- useNA[group.present]

    }
    # keep only grouping variables found in metadata
    group.by.composition <- group.by.composition[group.present]
    message("Only found ", paste(group.by.composition, collapse = ", ") , " as grouping variables.")
  }

  # Factorize cell type to keep all in compositional data
  meta.data <- meta.data %>%
    dplyr::mutate_at(unlist(group.by.composition), as.factor)

  # evaluate which layers are included in group.by.composition and layers_links
  for (l in seq_along(layers_links)) {
    if (!all(c(names(layers_links[l]), layers_links[l]) %in%
             group.by.composition)) {
      message("Not computing relative proportion values of layer2 within layer1 for ",
              paste0(names(layers_links[l]), " -> ", layers_links[l]))
    }
  }


  # Rep useNa parameters according to group.by.composition variable
  if (length(useNA) == 1) {
    useNA <- rep(useNA, length(group.by.composition))
  }

  # list to store compositions
  celltype.compositions <- list()


  for (i in seq_along(group.by.composition)) {

    suppressMessages(
      {
        ctable <- compositional_data(data = meta.data,
                                     split.by = split.by,
                                     group.by.1 = group.by.composition[[i]],
                                     useNA = useNA[i],
                                     clr_zero_impute_perc = clr_zero_impute_perc)

        # stop if very few cells
        if (sum(ctable$cell_counts) < min.cells.composition) {
          # Return empty data.frame
          ctable <- ctable[0,]
        }
        # get proportion relative to layer1 types
        if (group.by.composition[[i]] %in% layers_links) {
          # stop if very few cells
          if (sum(ctable$cell_counts) < min.cells.composition) {
            # Return empty list
            ctable <- list()
          } else {
            lay <- layers_links[which(layers_links == group.by.composition[[i]])]
            if (lay %in% names(object@misc$layer2_param)) {
              meta.split <- split(meta.data,
                                  meta.data[[names(lay)]])
              # switch list names to cell ontology ID
              cellonto_dic <- lapply(meta.split, function(x) {
                nam <- unique(x[[names(lay)]])
                val <- unique(x[[layer1_link]])
                names(val) <- nam
                return(val)
              }) %>%
                unname() %>%
                unlist()

              levs <- object@misc$layer2_param[[lay]]$levels2_per_levels1
              names(levs) <- names(cellonto_dic[match(names(levs), cellonto_dic)])

              # If a celltype was not detected, drop it
              levs <- levs[!is.na(names(levs))]

              # Filter only layer1 cells with representation in layer2
              meta.split <- meta.split[names(levs)]

              # add factor to group by layer1
              for (f in names(meta.split)) {
                meta.split[[f]]$functional.cluster <- factor(meta.split[[f]]$functional.cluster,
                                                             levels = levs[[f]])
              }

              ctable.split <- lapply(meta.split,
                                     compositional_data,
                                     split.by = split.by,
                                     group.by.1 = group.by.composition[[i]],
                                     useNA = useNA[i])

              # If not enough cells, return empty dataframe
              ctable.split <- lapply(ctable.split,
                                     function (x) {
                                       if (sum(x$cell_counts) < min.cells.composition) {
                                         x[0,]
                                       } else {x}
                                     })

              ctable <- c(ctable.split)
            }
          }
        }
      })

    ## Append
    celltype.compositions[[names(group.by.composition)[i]]] <- ctable
  }

  return(celltype.compositions)

}



#' Compute aggregated gene expression
#'
#' Function to compute aggregated expression (pseudobulk, i.e. sum counts per ident), and average expression by indicated cell type or grouping variable.
#'
#'
#' @param object A seurat object or a list of seurat objects
#' @param group.by.aggregated The Seurat object metadata column(s) containing celltype annotations
#' @param min.cells.aggregated Minimum number of cells per sample and cell type to aggregate data for pseudobulk or aggregated signatures. Aggregated data for cell types with less than that value will not be returned in aggregated data. Default value is 10.
#' @param assay Assay to retrieve information. By default "RNA".
#' @param useNA logical whether to return aggregated profile for NA (undefined) cell types, default is FALSE.
#' @param ... Extra parameters for internal Seurat functions: AverageExpression, AggregateExpression, FindVariableFeatures

#' @importFrom Seurat AverageExpression AggregateExpression FindVariableFeatures
#' @importFrom dplyr filter
#' @return Average and aggregated expression as a list of matrices for all genes and indicated gene lists filtering.
#' @export get.aggregated.profile


get.aggregated.profile <- function(object,
                                   group.by.aggregated = NULL,
                                   min.cells.aggregated = 10,
                                   assay = "RNA",
                                   useNA = FALSE,
                                   ...) {

  if (is.null(object)) {
    stop("Please provide a Seurat object")
    if (!inherits(object, "Seurat")) {
      stop("Please provide a Seurat object")
    }
  }

  if (is.null(group.by.aggregated)) {
    stop("Please specificy a group.by.aggregated variable")
  }

  # convert group.by.aggregated to list if not
  if (!is.list(group.by.aggregated)) {
    group.by.aggregated <- as.list(group.by.aggregated)
  }

  # Rename group.by.aggregated if not indicated
  for (v in seq_along(group.by.aggregated)) {
    if (is.null(names(group.by.aggregated)[[v]]) || is.na(names(group.by.aggregated)[[v]])) {
      names(group.by.aggregated)[[v]] <- group.by.aggregated[[v]]
    }
  }

  if (!any(group.by.aggregated %in% names(object@meta.data))) {
    stop("Group.by variables ", paste(group.by.aggregated, collapse = ", "), " not found in metadata")
  } else if (!all(group.by.aggregated %in% names(object@meta.data))) {
    group.by.aggregated <- group.by.aggregated[group.by.aggregated %in% names(object@meta.data)]
    message("Only found ", paste(group.by.aggregated, collapse = ", ") , " as grouping variables.")
  }


  avg.exp <- list()

  # loop over different grouping
  for (i in names(group.by.aggregated)) {

    # compute pseudobulk
    suppressWarnings({

      # Calculate total (sample) pseudobulks
      if (i == names(group.by.aggregated)[1]) {

        # Calculate pseudobulk for ALL cells in sample
        avg.exp[[i]] <- object@assays[["RNA"]]["counts"]
        row_names <- row.names(avg.exp[[i]])
        avg.exp[[i]] <- Matrix::Matrix(rowSums(avg.exp[[i]]))
        row.names(avg.exp[[i]]) <- row_names
        colnames(avg.exp[[i]]) <- "all"

        # Calculate pseudobulk for only annotated cells in sample

        # Subset to remove not annotated cells
        object <- object[, which(!is.na(object[[group.by.aggregated[[i]]]]))]

        # Calculate pseudobulk of all annotated cells
        mat <- object@assays[["RNA"]]["counts"]
        row_names <- row.names(mat)
        mat <- Matrix::Matrix(rowSums(mat))
        row.names(mat) <- row_names
        colnames(mat) <- "all.annotated_only"

        avg.exp[[i]] <- cbind(avg.exp[[i]], mat)
      }

      # remove from aggregated data cell with less than min.cells.aggregated
      cnts <- compositional_data(object@meta.data,
                                 group.by.1 = group.by.aggregated[[i]],
                                 only.counts = TRUE,
                                 useNA = useNA)

      if (nrow(cnts) == 0) {
        next
      }

      # Remove cell types with less than min.cells.aggregated
      keep <- cnts[cnts[["cell_counts"]] > min.cells.aggregated, 1] %>%
        unlist()
      object@meta.data$fltr <- object@meta.data[[group.by.aggregated[[i]]]]
      object <- object[, as.character(object$fltr) %in% keep]

      # Handle not annotated cells being labelled as NA
      # object@meta.data[[group.by.aggregated[[i]]]] <- as.character(object@meta.data[[group.by.aggregated[[i]]]])
      # object@meta.data[[group.by.aggregated[[i]]]][is.na(object@meta.data[[group.by.aggregated[[i]]]])] <- "NA"
      # object@meta.data[[group.by.aggregated[[i]]]] <- as.factor(object@meta.data[[group.by.aggregated[[i]]]])

      if (length(unique(object@meta.data[[group.by.aggregated[[i]]]])) >= 2) {
        mat <-
          Seurat::AggregateExpression(object,
                                      group.by = group.by.aggregated[[i]],
                                      assays = assay,
                                      verbose = FALSE,
                                      ...)[[assay]]
        avg.exp[[i]] <- cbind(avg.exp[[i]], mat)
      } else {
        # Handle case if there is only one cell type
        col_name <- as.character(unique(object@meta.data[[group.by.aggregated[[i]]]]))
        colnames(avg.exp[[i]]) <- col_name
      }
    })

    if (useNA == FALSE) {
      # Drop NA column
      avg.exp[[i]] <- avg.exp[[i]][, !colnames(avg.exp[[i]]) %in% "NA", drop = FALSE]
    }
  }
  return(avg.exp)
}


#' Compute aggregated additional signatures by cell type
#'
#' Function to compute aggregated signatures of predicted cell types.
#'
#'
#' @param object A seurat object or metadata data frame.
#' @param group.by.aggregated The Seurat object metadata column(s) containing celltype annotations (idents).
#' @param min.cells.aggregated Minimum number of cells per sample and cell type to aggregate data for pseudobulk or aggregated signatures. Aggregated data for cell types with less than that value will not be returned in aggregated data. Default value is 10.
#' @param name.additional.signatures Names of additional signatures to compute the aggregation per cell type.
#' @param fun Function to aggregate the signature, e.g. mean or sum.
#' @param useNA logical whether to return aggregated signatures for NA (undefined) cell types, default is FALSE.

#' @importFrom dplyr group_by summarize_at filter
#' @return Aggregated signature score for each indicated cell type grouping Results is NULL of not additional signatures are indicated or present in metadata.
#' @export get.aggregated.signature


get.aggregated.signature <- function(object,
                                     group.by.aggregated = NULL,
                                     min.cells.aggregated = 10,
                                     name.additional.signatures = NULL,
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


  if (is.null(group.by.aggregated)) {
    stop("Please specificy a group.by.aggregated variable")
  }

  # convert group.by.aggregated to list if not
  if (!is.list(group.by.aggregated)) {
    group.by.aggregated <- as.list(group.by.aggregated)
  }

  # Rename group.by.aggregated if not indicated
  for (v in seq_along(group.by.aggregated)) {
    if (is.null(names(group.by.aggregated)[[v]]) ||
        is.na(names(group.by.aggregated)[[v]])) {
      names(group.by.aggregated)[[v]] <- group.by.aggregated[[v]]
    }
  }

  if (!any(group.by.aggregated %in% names(meta.data))) {
    stop("Group.by variables ", paste(group.by.aggregated, collapse = ", "), " not found in metadata")
  } else if (!all(group.by.aggregated %in% names(meta.data))) {
    group.by.aggregated <- group.by.aggregated[group.by.aggregated %in% names(meta.data)]
    message("Only found ", paste(group.by.aggregated, collapse = ", ") , " as grouping variable.")
  }

  if (is.null(name.additional.signatures)) {
    name.additional.signatures <- object@misc$layer1_param$additional.signatures
  }

  if (is.null(name.additional.signatures)) {
    message("No additional signatures indicated. Returning NULL")
    aggr.sig <- NULL
  } else {

    if (!any(grepl(paste(name.additional.signatures, collapse = "|"),
                   names(meta.data)))) {
      stop("No additional signatures found in this object metadata")
    }

    add.sig.cols <- grep(paste(name.additional.signatures, collapse = "|"),
                         names(meta.data), value = TRUE)

    if (length(add.sig.cols) > length(name.additional.signatures)) {
      for (i in name.additional.signatures) {
        if (sum(grep(i, names(meta.data))) > 1) {
          meta_cols <- names(meta.data)[grep(i, names(meta.data))]
          message(paste("Signatue", i, "was found in multiple metadata columns:", meta_cols))
        }
      }
      stop("The name of at least one signature provided was found in multiple metadata columns.
           Please give them a more unique name, e.g. by appending '_signature' to the name")
    }

    aggr.sig <- list()

    for (e in names(group.by.aggregated)) {
      # remove from aggregated data cell with less than min.cells.aggregated
      cnts <- compositional_data(meta.data,
                                 group.by.1 = group.by.aggregated[[e]],
                                 only.counts = TRUE,
                                 useNA = useNA)

      keep <- cnts[cnts[["cell_counts"]] > min.cells.aggregated, 1] %>%
        unlist()


      aggr.sig[[e]] <- meta.data %>%
        dplyr::filter(.data[[group.by.aggregated[[e]]]] %in% keep) %>%
        dplyr::group_by(.data[[group.by.aggregated[[e]]]]) %>%
        dplyr::summarize_at(add.sig.cols, fun, na.rm = TRUE)

      # filter out NA if useNA=F
      if (!useNA) {
        aggr.sig[[e]] <- aggr.sig[[e]] %>%
          dplyr::filter(!is.na(.data[[group.by.aggregated[[e]]]]))
      }

      colnames(aggr.sig[[e]])[1] <- "celltype"
    }
  }
  return(aggr.sig)
}



#' Merge HiTObjects
#'
#' @param hit.object List of HiTObjects
#' @param group.by If only merging for certain layers of annotation is intended, layers names can be indicated here as vector. Otherwise all layers present in all HiT object will be merged.
#' @param metadata.vars Variables to keep as metadata. (Default: NULL, keeping unique metadata columns per sample, dropping single-cell metadata)
#' @param pseudobulk.matrix Paramater to determine whther obtain the pseudobulk matrix as a single matrix (\code{"unique"}), or as one matrix for each cell type in the layer (\code{"list"})
#' @param ncores The number of cores to use, by default, all available cores - 2.
#' @param bparam A \code{BiocParallel::bpparam()} object that tells how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param progressbar Whether to show a progressbar or not
#' @param verbose Whether to show optional messages or not

#' @importFrom dplyr mutate mutate_if filter %>% mutate_all full_join
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallelly availableCores
#' @importFrom data.table rbindlist
#' @importFrom methods slot

#' @return Merged HiTObject
#' @export merge.HiTObjects
#'

merge.HiTObjects <- function(hit.object = NULL,
                             group.by = NULL,
                             metadata.vars = NULL,
                             pseudobulk.matrix = "list",
                             ncores = parallelly::availableCores() - 2,
                             bparam = NULL,
                             progressbar = FALSE,
                             verbose = FALSE) {

  if (is.null(hit.object) ||
      !is.list(hit.object)) {
    stop("Please provide a list of HiT object containing more than one sample")
    if (suppressWarnings(!all(lapply(hit.object, function(x) {inherits(x, "HiT")})))) {
      stop("Not all components of the list are HiT objects.")
    }
  }

  # TODO: delete? -------------------------------
  # # give name to list of hit objects
  for (v in seq_along(hit.object)) {
    if (is.null(names(hit.object)[v]) ||
        is.na(names(hit.object)[v])) {
      names(hit.object)[v] <- paste0("Sample", v)
    }
  }

  # fetch which layers are included in all HiT objects
  layers_in_hit <- lapply(hit.object, function(x) {
    a <- names(x@composition)
    b <- names(x@aggregated_profile$pseudobulk)
    c <- names(x@aggregated_profile$signatures)
    u <- unique(c(a,b,c))
    return(u)
  })

  # if group.by is not NULL retrieve the layers present in HiT object
  if (is.null(group.by)) {
    group.by <- unique(unlist(layers_in_hit))
  }


  if (suppressWarnings(!all(lapply(layers_in_hit, function(x) {any(group.by %in% x)})))) {
    stop("None of the supplied HiT object contain at least one of these layers: ", paste(group.by, collapse = ", "),
         ".\nPlease make sure to indicate the name of the layer.")
  } else {
    present <- sapply(group.by,
                      function(char) {
                        unlist(lapply(layers_in_hit,
                                      function(df) {char %in% df}
                        ))
                      }
    )

    present.sum <- colSums(present)

    for (l in names(group.by)) {
      message("** ", l, " present in ", present.sum[[l]], " / ", length(hit.object), " HiT objects.")
    }
  }

  if (!is.null(metadata.vars)) {
    if (verbose) {message("\n#### Metadata ####")}
    if (suppressWarnings(!all(lapply(hit.object, function(x) {any(metadata.vars %in% names(x@metadata))})))) {
      message("Not all supplied HiT object contain ", paste(metadata.vars, collapse = ", "),
              "metadata elements in their metadata")
    }
    if (verbose) {
      in.md <- sapply(metadata.vars,
                      function(char) {
                        unlist(lapply(hit.object,
                                      function(df) {char %in% colnames(df@metadata)}
                        ))
                      }
      )
      in.md.sum <- colSums(in.md)

      for (l in metadata.vars) {
        message("** ", l, " present in ",
                in.md.sum[[l]], " / ", length(hit.object),
                " HiT objects.")
      }
    }
  }

  # set parallelization parameters
  param <- set_parallel_params(ncores = ncores,
                               bparam = bparam,
                               progressbar = progressbar)

  # Join data from all HiT objects

  # Metadata
  if (is.null(metadata.vars)) {
    # Per HiTObject from input, get metadata column names which have all the same value.
    # E.g. each HiTObject is from a specific sample, so each HiTObject has metadata column "Sample" which contains all the same value (e.g. "sample1" for sample 1, "sample2" for sample 2, etc.)
    # Other metadata columns, e.g. scGate_multi can be dropped, as the merged HiTObject is a summary per sample, so single-cell metadata columns don't make sense.
    all.names <- list()
    for (x in names(hit.object)) {
      all.names[[x]] <- names(hit.object[[x]]@metadata)
    }
    all.names_in.all <- Reduce(intersect, all.names)

    umc <- list()
    for (x in names(hit.object)) {
      umc[[x]] <- apply(hit.object[[x]]@metadata[all.names_in.all], 2, function(y) length(unique(y))) == 1
    }
    umc <- umc %>% Reduce("&", .)
    metadata.vars <- names(umc)[umc == TRUE]
  }
  metadata <- lapply(names(hit.object),
                     function (x) {hit.object[[x]]@metadata[1, metadata.vars] %>%
                         mutate(hitme.sample = x)
                     }
  )
  metadata <- data.table::rbindlist(metadata, use.names=TRUE, fill=TRUE)

  comp.prop <- list()
  avg.expr <- list()
  aggr.signature <- list()

  for (gb in group.by) {

    layer_present <- row.names(present)[present[,gb]]

    # Composition
    message("Merging compositions of " , gb, "...")

    type <- "composition"

    # assess that all HiT objects in the list contain the composition data for that layer
    is_df_check <- hit.object[layer_present] %>%
      lapply(methods::slot, name = type) %>%
      lapply("[[", gb) %>%
      lapply(is.data.frame) %>%
      unlist()
    is_list_check <- hit.object[layer_present] %>%
      lapply(methods::slot, name = type) %>%
      lapply("[[", gb) %>%
      lapply(function(x) {is.list(x) &!inherits(x, "data.frame")}) %>%
      unlist()

    if (all(is_df_check)) {
      df <- BiocParallel::bplapply(X = layer_present,
                                   BPPARAM = param,
                                   function(x) {
                                     hit.object[[x]]@composition[[gb]] %>%
                                       mutate(hitme.sample = x)
                                   })
      df <- data.table::rbindlist(df, use.names=TRUE, fill=TRUE)
      comp.prop[[gb]] <- df

    } else if (all(is_list_check)) {
      gb_sublevel_unique_names <- hit.object[layer_present] %>%
        lapply(methods::slot, name = type) %>%
        lapply("[[", gb) %>%
        lapply(names) %>%
        unlist() %>%
        unique()
      for (i in gb_sublevel_unique_names) {
        df <- BiocParallel::bplapply(X = layer_present,
                                     BPPARAM = param,
                                     function(x) {
                                       if (!is.null(hit.object[[x]]@composition[[gb]][[i]])) {
                                         hit.object[[x]]@composition[[gb]][[i]] %>%
                                           mutate(hitme.sample = x)
                                       }
                                     })
        df <- data.table::rbindlist(df, use.names=TRUE, fill=TRUE)
        comp.prop[[gb]][[i]] <- df
      }
    }


    # Aggregated_profile
    message("Merging aggregated profiles of " , gb, "...")

    type <- "pseudobulk"

    # Aggregated_profile Pseudobulk
    celltypes <- lapply(names(hit.object),
                        function(x) {
                          colnames(hit.object[[x]]@aggregated_profile[[type]][[gb]])
                        }) %>%
      unlist() %>% unique()

    for (ct in celltypes) {
      ct_present <- lapply(names(hit.object),
                           function(x) {
                             ct %in% colnames(hit.object[[x]]@aggregated_profile[[type]][[gb]]) }) %>%
        unlist()
      layer_ct_present <- names(hit.object)[(names(hit.object) %in% layer_present) & ct_present]
      df <- BiocParallel::bplapply(X = layer_ct_present,
                                   BPPARAM = param,
                                   function(x) {
                                     hit.object[[x]]@aggregated_profile[[type]][[gb]][, ct]
                                   })
      df <- do.call(cbind, df)
      df <- Matrix::Matrix(df, sparse = TRUE)
      colnames(df) <- layer_ct_present
      avg.expr[[gb]][[ct]] <- df
    }

    # join each matrix per celltype into a single matrix changing the colnames to accommodate sample source
    if (tolower(pseudobulk.matrix) == "unique") {
      avg.expr[[gb]] <- lapply(names(avg.expr[[gb]]),
                               function(x) {
                                 mat <- avg.expr[[gb]][[x]]
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

    df <- BiocParallel::bplapply(X = layer_present,
                                 BPPARAM = param,
                                 function(x) {
                                   if (!is.null(hit.object[[x]]@aggregated_profile[[type]][[gb]])) {
                                     hit.object[[x]]@aggregated_profile[[type]][[gb]] %>%
                                       mutate(hitme.sample = x)
                                   }
                                 })
    df <- data.table::rbindlist(df,
                                use.names = TRUE,
                                fill = TRUE)
    aggr.signature[[gb]] <- df
  }

  hit <- methods::new("HiT",
                      metadata = metadata,
                      composition = comp.prop,
                      aggregated_profile = list("pseudobulk" = avg.expr,
                                                "signatures" = aggr.signature)
  )
  return(hit)
}


#' Get metrics of cell types pseudobulk clustering
#'
#' @param hit.object A Hit class object obtained with \link{get.HiTObject} or pseudobulk raw count matrix (\code{HitObject@aggregated_profile$pseudobulk$layer1})
#' @param cluster.by Vector indicating the variable for clustering, default is celltype (for the annotation) and sample
#' @param cluster.by.drop.na Whether to keep (FALSE) or drop (TRUE) NAs present in cluster.by column.
#' @param batching Vector indicating the variable for batching to allow calculating scores per batch, to account for batch effect
#' @param scores Scores to compute clustering of samples based on cell type prediction.
#' @param modularity.k Number of k-nearest neighbours to use to build graph for the calculation of the modularity score
#' @param dist.method Method to compute distance between celltypes, default is euclidean.
#' @param hclust.method Hierarchical clustering method for hclust. Options are: "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC). (See hclust package for details)
#' @param ntests Number of shuffling events to calculate p-value for scores
#' @param seed Set seed for random shuffling events to calculate p-value for scores
#' @param pval.combine Method for combining p-values if calculated using batching. Default is "zmethod" weighted (Stouffer's) Z-method, weighting by sample size per batch. Alternatively, Fisher's method can be used with "fisher".
#' @param pca_comps_labs_invisible Parameter to pass to fviz_pca defining which labels to plot for composition PCA, see documentation of fviz_pca for details.
#' @param pca_pb_labs_invisible Parameter to pass to fviz_pca defining which labels to plot for pseudobulk PCA, see documentation of fviz_pca for details.
#' @param pca_sig_labs_invisible Parameter to pass to fviz_pca defining which labels to plot for signatures PCA, see documentation of fviz_pca for details.
#' @param ndim Number of dimensions to be use for PCA clustering metrics. Default is 10.
#' @param nVarGenes Number of variable genes to assess samples. Default is 500.
#' @param black.list List of genes to discard from clustering, if "default" object "default_black_list" object is used. Alternative black listed genes can be provided as a vector or list.
#' @param ncores The number of cores to use, by default, all available cores - 2.
#' @param bparam A \code{BiocParallel::bpparam()} object that tells how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param progressbar Whether to show a progressbar or not

#' @importFrom tidyr separate pivot_wider
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallelly availableCores
#' @importFrom data.table rbindlist
#' @importFrom dplyr mutate mutate_if filter %>% coalesce mutate_all full_join row_number
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom tibble rownames_to_column column_to_rownames remove_rownames
#' @importFrom caret nearZeroVar
#' @importFrom ggplot2 aes geom_point guides theme geom_col labs guide_legend annotate theme_bw ggtitle geom_ribbon element_text
#' @importFrom DESeq2 DESeqDataSetFromMatrix vst estimateSizeFactors
#' @importFrom MatrixGenerics rowVars rowMins
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics counts
#' @importFrom ggdendro ggdendrogram
#' @importFrom stringr str_to_title
#' @importFrom stats prcomp na.omit formula pnorm t.test
#' @importFrom factoextra fviz_pca
#' @importFrom scran buildKNNGraph
#' @importFrom igraph modularity set_vertex_attr layout_nicely V
#' @importFrom patchwork wrap_plots plot_layout plot_annotation
#' @importFrom ggraph ggraph geom_edge_link geom_node_point
#' @importFrom metap sumlog sumz

#' @return Metrics of cell types pseudobulk clustering
#' @export get.cluster.score
#'


get.cluster.score <- function(hit.object = NULL,
                              cluster.by = NULL,
                              cluster.by.drop.na = TRUE,
                              batching = NULL,
                              scores = c("Silhouette_isolated", "Silhouette", "Modularity"),
                              modularity.k = 3,
                              dist.method = "euclidean",
                              hclust.method = "complete",

                              # For scores p-value calculation
                              ntests = 100, # number of shuffling events
                              seed = 22, # seed for random shuffling
                              pval.combine.method = "weighted_zmethod",

                              # For PCA
                              pca_comps_labs_invisible = c("quali"),
                              pca_pb_labs_invisible = c("var", "quali"),
                              pca_sig_labs_invisible = c("quali"),

                              # Pseudobulk params
                              ndim = 10,
                              nVarGenes = 500,
                              black.list = NULL,

                              ncores = round(parallelly::availableCores() / 3), # to reduce memory load
                              bparam = NULL,
                              progressbar = TRUE) {

  if (is.null(hit.object) ||
      !inherits(hit.object, "HiT")) {
    stop("Please provide a Hit class object or a count matrix.")
  }
  if (is.null(cluster.by)) {
    stop("Please provide a metadata column name to cluster by.")
  }

  if (!any(cluster.by %in% names(hit.object@metadata))) {
    stop("Group.by variables ", paste(cluster.by, collapse = ", "), " not found in metadata")
  } else if (!all(cluster.by %in% names(hit.object@metadata))) {
    cluster.by <- cluster.by[cluster.by %in% names(hit.object@metadata)]
    message("Only found ", paste(cluster.by, collapse = ", ") , " as grouping variable for HiT Object.")
  }

  # Need to replace special characters
  colnames(hit.object@metadata) <- make.names(colnames(hit.object@metadata))
  cluster.by <- make.names(cluster.by)

  for (i in cluster.by) {
    if (length(unique(hit.object@metadata[[i]])) == 1) {
      stop("All values are the same in cluster.by ", i, ". Please provide a metadata column with at least two different groups.")
    }
  }

  # Convert NA to factor level or drop
  for (i in cluster.by) {
    if (cluster.by.drop.na) {
      hit.object@metadata <- hit.object@metadata[!is.na(hit.object@metadata[[i]]), ]
    } else {
      hit.object@metadata[[i]] <- as.character(hit.object@metadata[[i]])
      hit.object@metadata[[i]][is.na(hit.object@metadata[[i]])] <- "NA"
    }
    # cluster.by column must be factor
    hit.object@metadata[[i]] <- as.factor(hit.object@metadata[[i]])
  }

  scores <- str_to_title(tolower(scores))

  # set parallelization parameters
  param <- set_parallel_params(ncores = ncores,
                               bparam = bparam,
                               progressbar = progressbar)

  # empty list to fill in the loop
  results <- list()


  # Process data ###############################################
  for (cluster_col in cluster.by) {
    message("Processing ", cluster_col)

    ## Process celltype composition ###############################################
    message("\nProcessing cell type composition\n")

    type <- "composition"

    comp_layers <- names(hit.object@composition)

    for (layer in comp_layers) {
      if (inherits(hit.object@composition[[layer]], "data.frame")) {
        mat <- hit.object@composition[[layer]][, c("celltype", "clr", "hitme.sample"), with = FALSE]
        if (cluster.by.drop.na) {
          mat <- mat %>% filter(hitme.sample %in% hit.object@metadata[["hitme.sample"]])
        }
        mat <- mat %>%
          tidyr::pivot_wider(names_from = hitme.sample,
                             values_from = clr) %>%
          stats::na.omit() %>%
          tibble::column_to_rownames(var = "celltype") %>%
          scale(center = TRUE,
                scale = FALSE)

        if (is.null(batching))  {
          cluster_labels <- hit.object@metadata %>%
            dplyr::filter(hitme.sample %in% colnames(mat)) %>%
            .[[cluster_col]]

          results[[cluster_col]][[type]][[layer]] <-
            get.scores(matrix = mat,
                       cluster_labels = cluster_labels,
                       scores = scores,
                       modularity.k = modularity.k,
                       dist.method = dist.method,
                       ntests = ntests,
                       seed = seed,
                       title = paste(cluster_col,
                                     stringr::str_to_title(type),
                                     layer),
                       invisible = pca_comps_labs_invisible)
        }

        if (!is.null(batching))  {
          for (b_var in batching) {
            if (b_var != cluster_col) {
              b_var_res_summary <- list()
              for (b in unique(hit.object@metadata[[b_var]])) {
                meta <- hit.object@metadata %>%
                  dplyr::filter(get(b_var) == b) %>%
                  dplyr::filter(hitme.sample %in% colnames(mat))

                cluster_labels <- meta[[cluster_col]]

                m <- mat[ , colnames(mat) %in% as.character(meta[["hitme.sample"]])] %>%
                  scale(center = TRUE,
                        scale = FALSE)

                results[[cluster_col]][[type]][[layer]][[b_var]][[b]] <-
                  get.scores(matrix = m,
                             cluster_labels = cluster_labels,
                             scores = scores,
                             modularity.k = modularity.k,
                             dist.method = dist.method,
                             ntests = ntests,
                             seed = seed,
                             title = paste(cluster_col,
                                           stringr::str_to_title(type),
                                           layer),
                             invisible = pca_comps_labs_invisible)
                for (score in scores) {
                  b_var_res_summary[[score]][["summary"]] <- c(
                    b_var_res_summary[[score]][["summary"]],
                    results[[cluster_col]][[type]][[layer]][[b_var]][[b]][["Scores"]][[score]][["summary"]])
                  b_var_res_summary[[score]][["n"]] <- c(
                    b_var_res_summary[[score]][["n"]],
                    results[[cluster_col]][[type]][[layer]][[b_var]][[b]][["Scores"]][[score]][["n"]])
                  b_var_res_summary[[score]][["p_value"]] <- c(
                    b_var_res_summary[[score]][["p_value"]],
                    results[[cluster_col]][[type]][[layer]][[b_var]][[b]][["Scores"]][[score]][["p_value"]])
                }
              }

              for (score in scores) {
                b_var_res_summary[[score]][["summary"]] <-
                  mean(b_var_res_summary[[score]][["summary"]])

                # Requires the number of samples per batch, so run before summing n
                p_values_not_na <- stats::na.omit(b_var_res_summary[[score]][["p_value"]])
                p_values_not_na_len <- length(p_values_not_na)
                if (p_values_not_na_len > 1){
                  b_var_res_summary[[score]] <- combine_pvals(b_var_res_summary[[score]],
                                                              pval.combine.method)
                } else if (p_values_not_na_len == 1) {
                  b_var_res_summary[[score]][["p_value"]] <- p_values_not_na
                } else {
                  b_var_res_summary[[score]][["p_value"]] <- NULL
                  b_var_res_summary[[score]][["summary"]] <- NULL
                }

                b_var_res_summary[[score]][["n"]] <-
                  sum(b_var_res_summary[[score]][["n"]])
              }
              results[[cluster_col]][[type]][[layer]][[b_var]][["all"]][["Scores"]] <- b_var_res_summary
            }
          }
        }

      } else if (is.list(hit.object@composition[[layer]])) {
        results[[cluster_col]][[type]][[layer]] <-
          BiocParallel::bplapply(
            X = names(hit.object@composition[[layer]]),
            BPPARAM = param,
            function(i){
              mat <- hit.object@composition[[layer]][[i]][, c("celltype", "clr", "hitme.sample"), with = F]
              if (cluster.by.drop.na) {
                mat <- mat %>% filter(hitme.sample %in% hit.object@metadata[["hitme.sample"]])
              }
              mat <- mat %>%
                tidyr::pivot_wider(names_from = hitme.sample,
                                   values_from = clr) %>%
                tibble::column_to_rownames(var = "celltype")

              mat <- scale(mat, center = TRUE,
                           scale = FALSE)

              if (is.null(batching))  {
                cluster_labels <- hit.object@metadata %>%
                  dplyr::filter(hitme.sample %in% colnames(mat)) %>%
                  .[[cluster_col]]

                if (nrow(mat) > 1) {
                  res <- get.scores(matrix = mat,
                                    cluster_labels = cluster_labels,
                                    scores = scores,
                                    modularity.k = modularity.k,
                                    dist.method = dist.method,
                                    ntests = ntests,
                                    seed = seed,
                                    title = paste(cluster_col,
                                                  stringr::str_to_title(type),
                                                  layer,
                                                  i),
                                    invisible = pca_comps_labs_invisible)
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
                    for (b in unique(hit.object@metadata[[b_var]])) {
                      meta <- hit.object@metadata %>%
                        dplyr::filter(get(b_var) == b) %>%
                        dplyr::filter(hitme.sample %in% colnames(mat))

                      cluster_labels <- meta[[cluster_col]]

                      m <- mat[ , colnames(mat) %in% as.character(meta[["hitme.sample"]])] %>%
                        scale(center = TRUE, scale = FALSE)

                      if (nrow(m) > 1) {
                        res[[b_var]][[b]] <-
                          get.scores(matrix = m,
                                     cluster_labels = cluster_labels,
                                     scores = scores,
                                     modularity.k = modularity.k,
                                     dist.method = dist.method,
                                     ntests = ntests,
                                     seed = seed,
                                     title = paste(cluster_col,
                                                   stringr::str_to_title(type),
                                                   layer,
                                                   i),
                                     invisible = pca_comps_labs_invisible)
                      }

                      for (score in scores) {
                        b_var_res_summary[[score]][["summary"]] <- c(
                          b_var_res_summary[[score]][["summary"]],
                          res[[b_var]][[b]][["Scores"]][[score]][["summary"]])
                        b_var_res_summary[[score]][["n"]] <- c(
                          b_var_res_summary[[score]][["n"]],
                          res[[b_var]][[b]][["Scores"]][[score]][["n"]])
                        b_var_res_summary[[score]][["p_value"]] <- c(
                          b_var_res_summary[[score]][["p_value"]],
                          res[[b_var]][[b]][["Scores"]][[score]][["p_value"]])
                      }
                    }

                    for (score in scores) {
                      b_var_res_summary[[score]][["summary"]] <-
                        mean(b_var_res_summary[[score]][["summary"]])

                      # Requires the number of samples per batch, so run before summing n
                      p_values_not_na <- stats::na.omit(b_var_res_summary[[score]][["p_value"]])
                      p_values_not_na_len <- length(p_values_not_na)
                      if (p_values_not_na_len > 1) {
                        b_var_res_summary[[score]] <- combine_pvals(b_var_res_summary[[score]],
                                                                    pval.combine.method)
                      } else if (p_values_not_na_len == 1) {
                        b_var_res_summary[[score]][["p_value"]] <- p_values_not_na
                      } else {
                        b_var_res_summary[[score]][["p_value"]] <- NULL
                        b_var_res_summary[[score]][["summary"]] <- NULL
                      }

                      b_var_res_summary[[score]][["n"]] <-
                        sum(b_var_res_summary[[score]][["n"]])
                    }
                    res[[b_var]][["all"]][["Scores"]] <- b_var_res_summary
                  }
                }
                if(length(res) == 0) {
                  return(NULL)
                } else {
                  return(res)
                }
              }
            }
          )
        names(results[[cluster_col]][[type]][[layer]]) <-
          names(hit.object@composition[[layer]])
      }
    }


    ## Process pseudobulk ###############################################
    message("\nProcessing Pseudobulks\n")

    type <- "pseudobulk"

    pb_layers <- names(hit.object@aggregated_profile[[type]])

    for (layer in pb_layers) {
      results[[cluster_col]][[type]][[layer]] <- BiocParallel::bplapply(
        X = names(hit.object@aggregated_profile[[type]][[layer]]),
        BPPARAM = param,
        function(i){
          mat <- hit.object@aggregated_profile[[type]][[layer]][[i]]
          meta <- hit.object@metadata %>%
            dplyr::filter(hitme.sample %in% colnames(mat))
          cluster_labels <- meta[[cluster_col]]
          if (cluster.by.drop.na) {
            mat <- mat[, colnames(mat) %in% as.character(meta[["hitme.sample"]])]
          }

          if (is.null(batching))  {
            if (length(unique(cluster_labels)) > 1) {
              mat <- preproc_pseudobulk(matrix = mat,
                                        metadata = meta,
                                        cluster.by = cluster_col,
                                        nVarGenes = nVarGenes,
                                        black.list = black.list)

              res <- get.scores(matrix = mat,
                                cluster_labels = cluster_labels,
                                scores = scores,
                                modularity.k = modularity.k,
                                dist.method = dist.method,
                                ntests = ntests,
                                seed = seed,
                                title = paste(cluster_col,
                                              stringr::str_to_title(type),
                                              layer,
                                              i),
                                invisible = pca_pb_labs_invisible)
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
                for (b in unique(hit.object@metadata[[b_var]])) {
                  met <- hit.object@metadata %>%
                    dplyr::filter(get(b_var) == b) %>%
                    dplyr::filter(hitme.sample %in% colnames(mat))

                  cluster_labels <- met[[cluster_col]]

                  m <- mat[ , colnames(mat) %in% met[["hitme.sample"]]]

                  if (length(unique(cluster_labels)) > 1) {
                    m <- preproc_pseudobulk(matrix = m,
                                            metadata = met,
                                            cluster.by = cluster_col,
                                            nVarGenes = nVarGenes,
                                            black.list = black.list)

                    if (nrow(m) > 1) {
                      res[[b_var]][[b]] <-
                        get.scores(matrix = m,
                                   cluster_labels = cluster_labels,
                                   scores = scores,
                                   modularity.k = modularity.k,
                                   dist.method = dist.method,
                                   ntests = ntests,
                                   seed = seed,
                                   title = paste(cluster_col,
                                                 stringr::str_to_title(type),
                                                 layer,
                                                 i),
                                   invisible = pca_pb_labs_invisible)
                    }
                  }
                  for (score in scores) {
                    b_var_res_summary[[score]][["summary"]] <- c(
                      b_var_res_summary[[score]][["summary"]],
                      res[[b_var]][[b]][["Scores"]][[score]][["summary"]])
                    b_var_res_summary[[score]][["n"]] <- c(
                      b_var_res_summary[[score]][["n"]],
                      res[[b_var]][[b]][["Scores"]][[score]][["n"]])
                    b_var_res_summary[[score]][["p_value"]] <- c(
                      b_var_res_summary[[score]][["p_value"]],
                      res[[b_var]][[b]][["Scores"]][[score]][["p_value"]])
                  }
                }

                for (score in scores) {
                  b_var_res_summary[[score]][["summary"]] <-
                    mean(b_var_res_summary[[score]][["summary"]])

                  # Requires the number of samples per batch, so run before summing n
                  p_values_not_na <- stats::na.omit(b_var_res_summary[[score]][["p_value"]])
                  p_values_not_na_len <- length(p_values_not_na)
                  if (p_values_not_na_len > 1) {
                    b_var_res_summary[[score]] <- combine_pvals(b_var_res_summary[[score]],
                                                                pval.combine.method)
                  } else if (p_values_not_na_len == 1) {
                    b_var_res_summary[[score]][["p_value"]] <- p_values_not_na
                  } else {
                    b_var_res_summary[[score]][["p_value"]] <- NULL
                    b_var_res_summary[[score]][["summary"]] <- NULL
                  }

                  b_var_res_summary[[score]][["n"]] <-
                    sum(b_var_res_summary[[score]][["n"]])
                }
                res[[b_var]][["all"]][["Scores"]] <- b_var_res_summary
              }
            }
            if(length(res) == 0) {
              return(NULL)
            } else {
              return(res)
            }
          }
        }
      )
      names(results[[cluster_col]][[type]][[layer]]) <-
        names(hit.object@aggregated_profile[[type]][[layer]])
    }


    ## Process signatures ###############################################
    message("\nProcessing Signatures\n")

    type <- "signatures"

    comp_layers <- names(hit.object@aggregated_profile[[type]])

    if (!is.null(comp_layers)) {
      for (layer in comp_layers) {
        cols <- colnames(hit.object@aggregated_profile[[type]][[layer]])
        signatures <- cols[! cols %in% c("celltype", "hitme.sample")]

        results[[cluster_col]][[type]][[layer]] <- BiocParallel::bplapply(
          X = signatures,
          BPPARAM = param,
          function(i){
            mat <- hit.object@aggregated_profile[[type]][[layer]][, c("celltype", i, "hitme.sample"), with = FALSE]
            if (cluster.by.drop.na) {
              mat <- mat %>% filter(hitme.sample %in% hit.object@metadata[["hitme.sample"]])
            }
            mat <- mat %>%
              tidyr::pivot_wider(names_from = hitme.sample, values_from = i) %>%
              tibble::column_to_rownames(var = "celltype") %>%
              replace(is.na(.), 0)

            if (is.null(batching))  {
              mat <- scale(mat,
                           center = TRUE,
                           scale = TRUE)

              # Drop columns containing NAN, caused by scaling zero variance columns
              mat <- mat[ , colSums(is.nan(mat)) == 0]

              cluster_labels <- hit.object@metadata %>%
                filter(hitme.sample %in% colnames(mat)) %>%
                .[[cluster_col]]

              res <- get.scores(matrix = mat,
                                cluster_labels = cluster_labels,
                                scores = scores,
                                modularity.k = modularity.k,
                                dist.method = dist.method,
                                ntests = ntests,
                                seed = seed,
                                title = paste(cluster_col,
                                              stringr::str_to_title(type),
                                              layer,
                                              i),
                                invisible = pca_sig_labs_invisible)
              return(res)
            }

            if (!is.null(batching)) {
              res <- list()
              for (b_var in batching) {
                if (b_var != cluster_col) {
                  b_var_res_summary <- list()
                  for (b in unique(hit.object@metadata[[b_var]])) {
                    meta <- hit.object@metadata %>%
                      dplyr::filter(get(b_var) == b) %>%
                      dplyr::filter(hitme.sample %in% colnames(mat))

                    cluster_labels <- meta[[cluster_col]]

                    m <- mat[ , colnames(mat) %in% as.character(meta[["hitme.sample"]])] %>%
                      scale(center = TRUE,
                            scale = TRUE)

                    # Drop columns containing NAN, caused by scaling zero variance columns
                    m <- m[ , colSums(is.nan(m)) == 0]

                    res[[b_var]][[b]] <-
                      get.scores(matrix = m,
                                 cluster_labels = cluster_labels,
                                 scores = scores,
                                 modularity.k = modularity.k,
                                 dist.method = dist.method,
                                 ntests = ntests,
                                 seed = seed,
                                 title = paste(cluster_col,
                                               stringr::str_to_title(type),
                                               layer,
                                               i),
                                 invisible = pca_sig_labs_invisible)
                    for (score in scores) {
                      b_var_res_summary[[score]][["summary"]] <- c(
                        b_var_res_summary[[score]][["summary"]],
                        res[[b_var]][[b]][["Scores"]][[score]][["summary"]])
                      b_var_res_summary[[score]][["n"]] <- c(
                        b_var_res_summary[[score]][["n"]],
                        res[[b_var]][[b]][["Scores"]][[score]][["n"]])
                      b_var_res_summary[[score]][["p_value"]] <- c(
                        b_var_res_summary[[score]][["p_value"]],
                        res[[b_var]][[b]][["Scores"]][[score]][["p_value"]])
                    }
                  }

                  for (score in scores) {
                    b_var_res_summary[[score]][["summary"]] <-
                      mean(b_var_res_summary[[score]][["summary"]])

                    # Requires the number of samples per batch, so run before summing n
                    p_values_not_na <- stats::na.omit(b_var_res_summary[[score]][["p_value"]])
                    p_values_not_na_len <- length(p_values_not_na)
                    if (p_values_not_na_len > 1) {
                      b_var_res_summary[[score]] <- combine_pvals(b_var_res_summary[[score]],
                                                                  pval.combine.method)
                    } else if (p_values_not_na_len == 1) {
                      b_var_res_summary[[score]][["p_value"]] <- p_values_not_na
                    } else {
                      b_var_res_summary[[score]][["p_value"]] <- NULL
                      b_var_res_summary[[score]][["summary"]] <- NULL
                    }

                    b_var_res_summary[[score]][["n"]] <-
                      sum(b_var_res_summary[[score]][["n"]])
                  }
                  res[[b_var]][["all"]][["Scores"]] <- b_var_res_summary
                }
              }
              if(length(res) == 0) {
                return(NULL)
              } else {
                return(res)
              }
            }
          }
        )
        names(results[[cluster_col]][[type]][[layer]]) <- signatures
      }
    }
  }

  # Save user set parameters for summarize.cluster.scores
  results[["params"]][["cluster.by"]] <- cluster.by
  results[["params"]][["batching"]] <- batching
  results[["params"]][["scores"]] <- scores

  return(results)
}




#' Summarize scores and plot heat map
#'
#' @param data Output from get.cluster.scores
#' @param topN Integer indicating number of topN highest scoring (most discriminating) features ranked for each score
#' @param create_plots Boolean indicating whether to create and show plots or not
#' @param p.adjustment Whether to adjust p-value columns or not
#' @param p.adjust.method Method for adjusting p-values (see stats::p.adjust for methods)
#' @param p.value_cutoff p-value (mean of all p-value columns) cutoff to filter out non-significant results
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
#' @export summarize.cluster.scores
#'


summarize.cluster.scores <- function(data = NULL,
                                     topN = NULL,
                                     create_plots = TRUE,
                                     p.adjustment = TRUE,
                                     p.adjust.method = "fdr",
                                     p.value_cutoff = 0.05) {

  show.variables <- c(".summary", ".n", ".p_value")

  if (is.null(data)) {
    message("Please provide scores object (output from get.cluster.scores)")
  }

  cluster.by <- data[["params"]][["cluster.by"]]
  batching <- data[["params"]][["batching"]]
  scores <- data[["params"]][["scores"]]

  data[["params"]] <- NULL # Remove params

  data_conts <- unlist(data)
  df_pars <- list()

  for (par in show.variables) {
    data_conts_temp <- data_conts[endsWith(names(data_conts), par)]

    if (!is.null(batching)) {
      data_conts_temp <- data_conts_temp[grepl(".all.Scores.",
                                               names(data_conts_temp))]
      names(data_conts_temp) <- gsub("all.", "", names(data_conts_temp))
    }

    if (length(data_conts_temp) == 0) {
      stop("Error happened at ", show.variables, ". No values left after filtering.")
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

      row.names(df) <- gsub(".Scores.",
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

  df <- merge(df_pars[[show.variables[1]]],
              df_pars[[show.variables[2]]],
              by = "row.names",
              all = TRUE)
  row.names(df) <- df$Row.names
  df$Row.names <- NULL

  df <- merge(df,
              df_pars[[show.variables[3]]],
              by = "row.names",
              all = TRUE)
  row.names(df) <- df$Row.names
  df$Row.names <- NULL

  df_cluster.by_list <- list()
  for (c in cluster.by) {
    df_cluster.by_list[[c]] <- df %>%
      dplyr::filter(row.names(.) %>%
               startsWith(c))
    row.names(df_cluster.by_list[[c]]) <- gsub(paste0(c, "."),
                                               "",
                                               row.names(df_cluster.by_list[[c]]))
    # Remove batching variable name
    if (!is.null(batching)) {
      batch_names <- batching[!batching %in% c]
      batch_names <- batch_names %>%
        paste0(".", .) %>%
        paste(collapse = "|")
      row.names(df_cluster.by_list[[c]]) <- gsub(batch_names,
                                                 "",
                                                 row.names(df_cluster.by_list[[c]]))
    }

    # Adjust p-values and filter
    if (p.adjustment) {
      p_val_cols <- which(grepl(".p_value",
                                colnames(df_cluster.by_list[[c]])))
      for (i in p_val_cols){
        df_cluster.by_list[[c]][, i] <-
          stats::p.adjust(df_cluster.by_list[[c]][, i],
                          method = p.adjust.method)
      }

      df_cluster.by_list[[c]] <-
        df_cluster.by_list[[c]][rowMeans(df_cluster.by_list[[c]][, p_val_cols]) <= p.value_cutoff, ]

      df_cluster.by_list[[c]] <- stats::na.omit(df_cluster.by_list[[c]])

      if (nrow(df_cluster.by_list[[c]]) == 0) {
        df_cluster.by_list[[c]] <- NULL
        message(paste("For ", c, " no separation was found after p-value cutoff. You can try to set it higher."))

        next
      }
    }

    # Get topN scored results
    if (!is.null(topN) &&
        is.numeric(topN) &&
        length(topN) == 1) {
      summary_score_cols <- which(grepl(".summary",
                                        colnames(df_cluster.by_list[[c]])))
      topN_rownames <- c()
      for (i in summary_score_cols) {
        topN_rownames_i <- row.names(df_cluster.by_list[[c]])[
          order(df_cluster.by_list[[c]][, i], decreasing = TRUE)][1:topN]
        topN_rownames <- c(topN_rownames,
                           topN_rownames_i)
      }
      df_cluster.by_list[[c]] <- df_cluster.by_list[[c]][unique(topN_rownames), ]
    }
  }

  # Remove NULL elements
  df_cluster.by_list <- Filter(Negate(is.null), df_cluster.by_list)

  # Check if all df_cluster.by_list were NULL
  if (length(df_cluster.by_list) == 0) {

    message("No significant separation found for cluster.by provided")

  } else {

    if (create_plots) {
      df_list <- list()
      plot_list <- list()

      for (n in names(df_cluster.by_list)) {
        df <- stats::na.omit(df_cluster.by_list[[n]])

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
        df_pval <- 1 / (df_pval + 1e-16)
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
      df_cluster.by_list[["plots"]][["plot_list"]] <- plot_list

      g <- gridExtra::grid.arrange(
        gridExtra::arrangeGrob(grobs = plot_list,
                               ncol=length(plot_list))
      )
      df_cluster.by_list[["plots"]][["summary_plot"]] <- ggplotify::as.ggplot(g)
    }

    return(df_cluster.by_list)
  }
}




#' Render plots summarizing celltype proportions and distribution in samples
#'
#'
#' @param obj.list List of Seurat objects
#' @param annot.col Metadata column(s) containing the cell type annotations
#' @param bottom.mar Adjustable bottom margin for long sample names

#' @importFrom stats setNames

#' @return Get percentage of not annotated cells per sample and plot it.
#' @export nas.per.sample
#'

nas.per.sample <- function (obj.list = NULL,
                            annot.col = c("scGate_multi"),
                            return.plot = TRUE,
                            bottom.mar = 10.2) {
  if (is.null(obj.list) &
      !is.list(obj.list) &
      !all(lapply(obj.list, inherits, "Seurat"))) {
    stop("Please provide a list of seurat objects")
  }

  na_list <- list()

  for (col in annot.col) {
    na_perc_per_sample <- c()
    for (i in names(obj.list)) {
      if (col %in% names(obj.list[[i]]@meta.data)) {
        percs <- prop.table(table(obj.list[[i]]@meta.data[[col]], useNA = "ifany"))*100
        nas_perc <- unname(percs[is.na(names(percs))])
        na_perc_per_sample <- c(na_perc_per_sample, stats::setNames(nas_perc, i))
      } else {
        stop(paste(col, " not found in obj.list item ", i))
      }
    }
    par(mar = c(bottom.mar, 4.1, 4.1, 2.1))
    if (return.plot) {
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
#' @param hit.object A Hit class object (typically after applying merge.HiTObjects onto a list of HiTObjects)
#' @param layer Default "layer1" if you have one cell type annotation layer in your hit.object. Alternatively "layer2" etc. if you have multiple layers of annotation depths.
#' @param return.plot.to.var Optionally, you can save the ggplots to a variable if you would like to further modify and adapt the plots on your own.
#' @param facet.by This allows you to pass a metadata column name present in your hit.object$metadata to show your samples in facets with ggplot facet_grid, for example by "condition".

#' @importFrom ggplot2 ggplot aes geom_bar theme element_text ggtitle facet_grid
#' @importFrom stats reformulate
#' @importFrom patchwork wrap_plots

#' @return Plotting function to show the cell type composition from HiTME object across different samples.
#' @export composition.barplot
#'

composition.barplot <- function (hit.object = NULL,
                                 sample.col = NULL,
                                 layer = "layer1",
                                 return.plot.to.var = FALSE,
                                 facet.by = NULL) {

  # Need to replace special characters
  colnames(hit.object@metadata) <- make.names(colnames(hit.object@metadata))

  sample.col <- "hitme.sample"

  if (is.null(hit.object)) {
    stop("Please provide input hit.object")
  }
  if (!length(layer) == 1 || !is.character(layer)) {
    stop("Please provide one character string for layer parameter")
  }
  if (!is.null(facet.by)) {
    if (!is.character(facet.by)) {
      stop("Please provide a character string or a vector of character strings for the facet.by parameter")
    }
    facet.by <- make.names(facet.by)
    facet.by.in.colnames <- facet.by %in% names(hit.object@metadata)
    if (!all(facet.by.in.colnames)) {
      facet.by.not.in.colnames <- facet.by[!facet.by.in.colnames]
      stop(paste("facet.by ", facet.by.not.in.colnames, " not found in hit.object@metadata column names"))
    }
  }


  comps <- hit.object@composition[[layer]]
  meta <- hit.object@metadata

  if (!is.null(facet.by)) {
    facet.by_reformulate <- reformulate(facet.by)
  }

  if (is.data.frame(comps)) {
    comp <- merge(comps, meta[, c(sample.col, facet.by), drop=FALSE], by = sample.col)

    p <- ggplot(comp, aes(x = hitme.sample, y = freq, fill = celltype)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust=1))

    if (!is.null(facet.by)) {

      p <- p + facet_grid(facet.by_reformulate,
                          space  = "free",
                          scales = "free")
    }

    if (return.plot.to.var) {
      return(p)
    } else {
      print(p)
    }
  }

  else {
    p_list <- list()
    for (ct in names(comps)) {
      comp <- hit.object@composition[[layer]][[ct]]
      comp <- merge(comp, meta[, c(sample.col, facet.by), drop=FALSE], by = sample.col)

      p_list[["plot_list"]][[ct]] <- ggplot(comp, aes(x = hitme.sample, y = freq, fill = celltype)) +
        geom_bar(stat = "identity") +
        ggtitle(ct) +
        theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        ggtitle(ct)

      if (!is.null(facet.by)) {
        p_list[["plot_list"]][[ct]] <- p_list[["plot_list"]][[ct]] +
          facet_grid(facet.by_reformulate,
                     space  = "free",
                     scales = "free")
      }

    }
    p_list[["arranged_plots"]] <- patchwork::wrap_plots(p_list[["plot_list"]])

    if (return.plot.to.var) {
      return(p_list)
    } else {
      print(p_list[["arranged_plots"]])
    }
  }
}


#' Render box plots summarizing celltype proportions and distribution in samples and groups
#'
#'
#' @param hit.object A Hit class object (typically after applying merge.HiTObjects onto a list of HiTObjects)
#' @param plot.var Column in the hit.object$composition: either "freq" for cell type relative abundance in percent or "clr" (for centered log-ratio transformed). Default: "clr" as it is better suited for statistical analysis and is better able to also show low abundant cell types.
#' @param layer Default "layer1" if you have one cell type annotation layer in your hit.object. Alternatively "layer2" etc. if you have multiple layers of annotation depths.
#' @param return.plot.to.var Optionally, you can save the ggplots to a variable if you would like to further modify and adapt the plots on your own.
#' @param group.by This allows you to pass a metadata column name present in your hit.object$metadata to show your samples in groups, for example by "condition".
#' @param facet.by This allows you to pass a metadata column name present in your hit.object$metadata to show your samples in facets with ggplot facet_grid, for example by "condition".
#' @param pval.method Specify method how to calculate the p-value. Wilcoxon test is recommended as compositional data might not fullfill the assumption of a gaussian distribution. For alternatives, see documentation of ggpubr::stat_pwc.
#' @param p.adjust.method Method for adjusting p-values (see ggpubr::stat_pwc for available methods)
#' @param palette Choose a palette of your liking. For available palettes, see ggsci package. Default: "lancet"
#' @param legend.position Where to put the legend. Possible options: "top", "right", "bottom", "left"

#' @importFrom ggplot2 ggplot aes geom_boxplot theme element_text ggtitle facet_grid position_jitterdodge
#' @importFrom ggpubr stat_pwc ggboxplot
#' @importFrom stats reformulate
#' @importFrom patchwork wrap_plots

#' @return Plotting function to show the cell type composition from HiTME object across different samples.
#' @export composition.boxplot
#'

composition.boxplot <- function (hit.object = NULL,
                                 plot.var = "clr",
                                 layer = "layer1",
                                 return.plot.to.var = FALSE,
                                 group.by = NULL,
                                 facet.by = NULL,
                                 pval.method = "wilcox.test",
                                 p.adjust.method = "BH",
                                 palette = "lancet",
                                 legend.position = "right") {

  # Need to replace special characters
  colnames(hit.object@metadata) <- make.names(colnames(hit.object@metadata))

  sample.col <- "hitme.sample"

  if (is.null(hit.object)) {
    stop("Please provide input hit.object")
  }
  if (!length(plot.var) == 1 ||
      !is.character(plot.var) ||
      !plot.var %in% c("freq", "clr")) {
    stop("Please provide one character string for plot.var, either 'freq' or 'clr'")
  }
  if (!length(layer) == 1 || !is.character(layer)) {
    stop("Please provide one character string for layer parameter")
  }
  if (!is.null(group.by)) {
    if (!length(group.by) == 1 || !is.character(group.by)) {
      stop("Please provide one character string for the group.by parameter")
    }
    group.by <- make.names(group.by)
    group.by.gg <- sym(group.by)
    nr_of_boxplots <- hit.object@metadata[[group.by]] %>%
      unique() %>%
      length() %>%
      "*"(1.5) %>%
      round()
    hit.object@metadata[group.by] <- lapply(hit.object@metadata[group.by], as.factor)
  }
  if (!is.null(facet.by)) {
    if (!is.character(facet.by)) {
      stop("Please provide a character string or a vector of character strings for the facet.by parameter")
    }
    facet.by <- make.names(facet.by)

    facet.by.in.colnames <- facet.by %in% names(hit.object@metadata)
    if (!all(facet.by.in.colnames)) {
      facet.by.not.in.colnames <- facet.by[!facet.by.in.colnames]
      stop(paste("facet.by ", facet.by.not.in.colnames, " not found in hit.object@metadata column names"))
    }
  }

  comps <- hit.object@composition[[layer]]
  meta <- hit.object@metadata

  plot.var.gg <- sym(plot.var)

  if (is.data.frame(comps)) {
    comp <- merge(comps, meta[, c(sample.col, group.by, facet.by), drop=FALSE], by = sample.col)

    # Need to check if group.by is NULL
    # Due to a presumed bug, if group.by is passed as variable to ggboxplot, even if it is assigned NULL, it throws an error
    if (is.null(group.by)) {
      p <- ggboxplot(comp,
                     x = "celltype",
                     y = plot.var,
                     outlier.shape = NA,
                     palette = palette,
                     facet.by = facet.by,
                     legend = legend.position) +
        geom_jitter(width = 0.2, size = 1)
    } else {
      p <- ggboxplot(comp,
                     x = "celltype",
                     y = plot.var,
                     color = group.by,
                     outlier.shape = NA,
                     palette = palette,
                     facet.by = facet.by,
                     legend = legend.position) +
        geom_jitter(mapping = aes(color = !!group.by.gg), position=position_jitterdodge(jitter.width = 1/nr_of_boxplots), size = 1) +
        stat_pwc(aes(group = !!group.by.gg),
                 label = "p.signif",
                 method = pval.method,
                 p.adjust.method = p.adjust.method,
                 p.adjust.by = "panel",
                 tip.length = 0,
                 hide.ns = TRUE)
    }

    p <- p +
      theme(axis.text.x = element_text(angle = 45, hjust=1)) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

    if (return.plot.to.var) {
      return(p)
    } else {
      print(p)
    }
  }

  else {
    p_list <- list()
    for (ct in names(comps)) {
      comp <- hit.object@composition[[layer]][[ct]]
      comp <- merge(comp, meta[, c(sample.col, group.by, facet.by), drop=FALSE], by = sample.col)

      # Need to check if group.by is NULL
      # Due to a presumed bug, if group.by is passed as variable to ggboxplot, even if it is assigned NULL, it throws an error
      if (is.null(group.by)) {
        p_list[["plot_list"]][[ct]] <- ggboxplot(comp,
                                                 x = "celltype",
                                                 y = plot.var,
                                                 outlier.shape = NA,
                                                 palette = palette,
                                                 facet.by = facet.by,
                                                 legend = legend.position) +
          geom_jitter(width = 0.2, size = 1)
      } else {
        p_list[["plot_list"]][[ct]] <- ggboxplot(comp,
                                                 x = "celltype",
                                                 y = plot.var,
                                                 color = group.by,
                                                 outlier.shape = NA,
                                                 palette = palette,
                                                 facet.by = facet.by,
                                                 legend = legend.position) +
          geom_jitter(mapping = aes(color = !!group.by.gg), position=position_jitterdodge(jitter.width = 1/nr_of_boxplots), size = 1) +
          stat_pwc(aes(group = !!group.by.gg),
                   label = "p.signif",
                   method = pval.method,
                   p.adjust.method = p.adjust.method,
                   p.adjust.by = "panel",
                   tip.length = 0,
                   hide.ns = TRUE)
      }

      p_list[["plot_list"]][[ct]] <- p_list[["plot_list"]][[ct]] +
        theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

    }
    p_list[["arranged_plots"]] <- patchwork::wrap_plots(p_list[["plot_list"]])

    if (return.plot.to.var) {
      return(p_list)
    } else {
      print(p_list[["arranged_plots"]])
    }
  }
}
