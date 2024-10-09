# scoot class definition ##############################################################################

## Build S4 objects to store data
methods::setClass(Class = "scoot",
                  slots = list(
                    data = "list",
                    metadata = "data.frame",
                    version = "package_version"
                  )
)





set_parallel_params <- function(ncores,
                                bparam,
                                progressbar)
{
  if (is.null(ncores)) {
    ncores <- 1
  }

  if (ncores > parallelly::availableCores()) {
    ncores <- parallelly::availableCores()
    message("Using more cores available in this computer, reducing number of cores to ", ncores)
  }

  # set parallelization parameters
  if (is.null(bparam)) {
    if (ncores > 1) {
      param <- BiocParallel::MulticoreParam(workers =  ncores,
                                            progressbar = progressbar)
    } else {
      param <- BiocParallel::SerialParam(progressbar = progressbar)
    }
  } else {
    param <- bparam
  }
  return(param)
}





compositional_data <- function(data,
                               split_by = NULL,
                               ann_layer_col_1 = NULL,
                               useNA = FALSE,
                               clr_zero_impute_perc = 1,
                               only_counts = FALSE) {

  if (all(is.na(data[[ann_layer_col_1]]))) {
    if (!only_counts) {
      ctable <- data.frame("celltype" = character(),
                           "cell_counts" = integer(),
                           "freq" = numeric(),
                           "clr" = numeric())
    } else {
      ctable <- data.frame("celltype" = character(),
                           "cell_counts" = integer())
    }
    return(ctable)
  } else {
    # set grouping variables
    gr_vars <- c(split_by, ann_layer_col_1)
    gr_vars2 <- c(split_by)

    ctable <- data %>%
      # drop = F keeps all levels of the factor
      dplyr::group_by(dplyr::across(dplyr::all_of(gr_vars)), .drop = F) %>%
      dplyr::summarize(cell_counts = dplyr::n()) %>%
      dplyr::ungroup()

    colnames(ctable)[1] <- "celltype"

    if (!only_counts) {
      ctable <- ctable %>%
        dplyr::filter(if (!useNA) !is.na(.data[["celltype"]])
                      else rep(TRUE, n())) %>%
        dplyr::group_by(across(all_of(gr_vars2))) %>%
        dplyr::mutate(freq = cell_counts/sum(cell_counts) * 100,
                      !!"celltype" := coalesce(.data[["celltype"]], "NA")) %>%
        as.data.frame() %>%
        na.omit()

      if (nrow(ctable) > 0) {
        # compute clr
        clr_df <- ctable %>%
          dplyr::select(-cell_counts) %>%
          # add pseudocount
          dplyr::mutate(freq = freq + clr_zero_impute_perc) %>%
          tidyr::pivot_wider(names_from = "celltype",
                             values_from = "freq")

        # accommodate df for clr transformation
        ## Remove character columns
        num_cols_bool_idx <- sapply(clr_df, is.numeric)
        num_cols <- names(clr_df)[num_cols_bool_idx]
        chr_cols <- names(clr_df)[!num_cols_bool_idx]
        clr_df_ref <- clr_df %>% dplyr::select(all_of(num_cols))

        clr <- Hotelling::clr(clr_df_ref)

        # add extra cols (if any)
        clr <- cbind(clr_df[,chr_cols],clr)  %>%
          tidyr::pivot_longer(-chr_cols,
                              names_to = "celltype",
                              values_to = "clr")
        # join clr df to main dataframe
        ctable <- dplyr::left_join(ctable,
                                   clr,
                                   by = c(chr_cols, "celltype"))
      }
    }

    return(ctable)
  }
}


# get_cluster_score helper functions ##############################################################################

filter_feat_dist <- function(df,
                             min_samples,
                             dist_method) {
  if (nrow(df) >= min_samples) {
    d <- list()
    d[["feature_matrix"]] <- df
    dist_mat <- stats::dist(df, method = dist_method)
    d[["distance_matrix"]] <- as.matrix(dist_mat)
    return(d)
  } else {
    return(NULL)
  }
}


pre_proc_comp <- function(df_comp,
                          min_samples,
                          dist_method) {
  df <- df_comp %>%
    tidyr::pivot_wider(names_from = scoot_sample,
                       values_from = clr) %>%
    stats::na.omit() %>%
    tibble::column_to_rownames(var = "celltype") %>%
    t() %>%
    scale(center = TRUE,
          scale = FALSE)

  feat_dist_list <- filter_feat_dist(df, min_samples, dist_method)
  return(feat_dist_list)
}

get_cluster_score_pre_proc <- function(scoot_object,
                                       cluster_by_drop_na,
                                       min_samples,
                                       dist_method,
                                       nvar_genes,
                                       cluster_by,
                                       black_list) {

  data <- list()

  ## Process celltype composition ###############################################
  type <- "composition"

  comp_layers <- names(scoot_object@data[[type]])

  for (layer in comp_layers) {

    if (inherits(scoot_object@data[[type]][[layer]], "data.frame")) {
      df <- scoot_object@data[[type]][[layer]][, c("celltype", "clr", "scoot_sample"), with = FALSE]
      if (cluster_by_drop_na) {
        df <- df %>% filter(scoot_sample %in% scoot_object@metadata[["scoot_sample"]])
      }
      data[[type]][[layer]] <- pre_proc_comp(as.data.frame(df), min_samples, dist_method)

    } else if (is.list(scoot_object@data[[type]][[layer]])) {
      data[[type]][[layer]] <- sapply(
        names(scoot_object@data[[type]][[layer]]),
        function(i) {
          df <- scoot_object@data[[type]][[layer]][[i]][, c("celltype", "clr", "scoot_sample"), with = F]
          if (cluster_by_drop_na) {
            df <- df %>% filter(scoot_sample %in% scoot_object@metadata[["scoot_sample"]])
          }
          pre_proc_comp(as.data.frame(df), min_samples, dist_method)
        }, simplify = FALSE, USE.NAMES = TRUE
      )
    }
  }


  ## Process pseudobulk ###############################################
  type <- "pseudobulk"

  pb_layers <- names(scoot_object@data[[type]])

  # Get black list
  if (is.null(black_list)) {
    black_list <- default_black_list
  }
  black_list <- unlist(black_list)

  for (layer in pb_layers) {
    data[[type]][[layer]] <- sapply(
      names(scoot_object@data[[type]][[layer]]),
      function(i) {
        mat <- as.matrix(scoot_object@data[[type]][[layer]][[i]])
        if (cluster_by_drop_na) {
          mat <- mat[, colnames(mat) %in% as.character(scoot_object@metadata[["scoot_sample"]])]
        }
        if (length(mat) == 0 || !is.matrix(mat)) {
          return(NULL)
        }
        # Remove black listed genes from the matrix
        mat <- mat[!row.names(mat) %in% black_list,]

        meta <- scoot_object@metadata %>%
          dplyr::filter(scoot_sample %in% colnames(mat))

        mat <- DESeq2.normalize(matrix = mat,
                                metadata = meta,
                                cluster_by = cluster_by,
                                nvar_genes = nvar_genes,
                                black_list = black_list)

        mat <- t(mat)
        d <- filter_feat_dist(mat, min_samples, dist_method)
        return(d)
      }, simplify = FALSE, USE.NAMES = TRUE
    )
  }


  ## Process signatures ###############################################
  type <- "signatures"

  sig_layers <- names(scoot_object@data[[type]])

  if (!is.null(sig_layers)) {
    for (layer in sig_layers) {
      cols <- colnames(scoot_object@data[[type]][[layer]])
      signatures <- cols[!cols %in% c("celltype", "scoot_sample")]

      data[[type]][[layer]] <- lapply(
        signatures,
        function(i) {
          df <- scoot_object@data[[type]][[layer]][, c("celltype", i, "scoot_sample"), with = FALSE]
          if (cluster_by_drop_na) {
            df <- df %>% filter(scoot_sample %in% scoot_object@metadata[["scoot_sample"]])
          }
          df <- df %>%
            tidyr::pivot_wider(names_from = scoot_sample,
                               values_from = i) %>%
            tibble::column_to_rownames(var = "celltype") %>%
            replace(is.na(.), 0) %>%    # Replace NAs with zero
            t() %>%
            scale(center = TRUE,
                  scale = TRUE) %>%
            as.data.frame() %>%
            select_if(~ !any(is.na(.)))
          filter_feat_dist(df, min_samples, dist_method)
        }
      )
      names(data[[type]][[layer]]) <- signatures
    }
  }

  return(data)
}



#' @importFrom DESeq2 DESeqDataSetFromMatrix vst estimateSizeFactors
#' @importFrom SummarizedExperiment assay

DESeq2.normalize <- function(matrix,
                             metadata,
                             cluster_by,
                             nvar_genes = 500,
                             black_list = NULL) {

  suppressMessages({
    suppressWarnings({

      # Normalize pseudobulk data using DESeq2
      # do formula for design with the cluster_by elements in order
      dformula <-  stats::formula(paste("~", cluster_by))
      matrix <- DESeq2::DESeqDataSetFromMatrix(countData = matrix,
                                               colData = metadata,
                                               design = dformula)
      matrix <- DESeq2::estimateSizeFactors(matrix)

      nsub <- min(1000, sum(rowMeans(BiocGenerics::counts(matrix, normalized=TRUE)) > 10 ))

      # transform counts using vst
      matrix <- DESeq2::vst(matrix, blind = T, nsub = nsub)
      matrix <- SummarizedExperiment::assay(matrix)

      # get top variable genes
      rv <- MatrixGenerics::rowVars(matrix)
      select <- order(rv, decreasing=TRUE)[seq_len(min(nvar_genes, length(rv)))]
      select <- row.names(matrix)[select]

      matrix <- matrix[select[select %in% row.names(matrix)],]

    })
  })

  return(matrix)
}



#' @importFrom cluster pam
#' @importFrom ggpubr get_legend as_ggplot
#' @importFrom factoextra fviz_pca fviz_nbclust
#' @importFrom scran buildKNNGraph
#' @importFrom igraph modularity set_vertex_attr layout_nicely V strength gorder
#' @importFrom patchwork wrap_plots plot_layout plot_annotation wrap_elements plot_spacer
#' @importFrom ggraph ggraph geom_edge_link geom_node_point
#' @importFrom mclust Mclust mclustBIC
#' @importFrom pheatmap pheatmap
#' @importFrom ggplotify as.ggplot
#' @importFrom ComplexHeatmap Heatmap rowAnnotation

get_scores_unsup <- function(data,
                             metadata,
                             scores,
                             knn_k,
                             max_nc,
                             n_clust,
                             NbClust_method,
                             ntests,
                             seed,
                             title = "", # Title for summary plot
                             select.var, # Number of pca loadings to show
                             show_clustering = TRUE,
                             heatmap_metadata_cols
) {

  feat_mat <- data[["feature_matrix"]]
  dist_mat <- data[["distance_matrix"]]

  nsamples <- dim(feat_mat)[1]

  # Need at least 4 samples, otherwise clustering doe not make sense
  if (nsamples >= 6) {
    results <- list()
    clus_scores_list <- list()
    pca_plot_list <- list()

    # Plot PCA ###############################################
    results[["plots"]][["pca_loadings"]] <- plot_pca(feat_mat,
                                                     label = "var",
                                                     labelsize = 3,
                                                     select.var = select.var) +
      ggplot2::ggtitle("PCA showing loadings")

    results[["plots"]][["pca_samples"]] <- plot_pca(feat_mat,
                                                    invisible = "var",
                                                    select.var = select.var) +
      ggplot2::ggtitle("PCA showing sample names")


    # Unsupervised clustering ###############################################

    # Set maximum number of clusters (cannot be more than samples - 1)
    if (nsamples <= max_nc) {
      max_nc <- length(nsamples) - 1
    } else {
      max_nc <- max_nc
    }

    # Use rule of thumb to determine k nearest neighbours to consider
    if (knn_k == "auto") {
      knn_k <- round(sqrt(nsamples))
      if (knn_k < 3) {
        knn_k <- 3
      }
    }

    # Args for calc_score
    args_list <- list(feat_mat = feat_mat,
                      dist_mat = dist_mat,
                      scores = scores,
                      ntests = ntests,
                      seed = seed,
                      knn_k = knn_k)


    ## Hierarchical clustering ###############################################

    # Methods commonly throwing errors or not creating specific Best.nc are commented out
    indeces <- c("kl", "ch", "hartigan",
                 # "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin",
                 "cindex", "db", "silhouette", "duda", "pseudot2",
                 # "beale",
                 "ratkowsky", "ball", "ptbiserial", "gap",
                 # "frey",
                 "mcclain", "gamma", "gplus", "tau", "dunn",
                 # "hubert",
                 "sdindex",
                 # "dindex",
                 "sdbw")

    # Determine the optimal number of clusters
    clust_results <- sapply(indeces, function(x) {
      tryCatch({
        fastNbClust(data = feat_mat,
                    min.nc = 1,
                    max.nc = max_nc,
                    method = NbClust_method,
                    index = x)
      }, error = function(e) {
        return(NULL)
      })
    })

    clust_results[sapply(clust_results, is.null)] <- NULL

    nclust_table <- clust_results %>%
      sapply(function(x) x[["Best.nc"]][["Number_clusters"]]) %>%
      unlist() %>%
      table()

    nclust <- nclust_table %>%
      which.max() %>%
      nclust_table[.] %>%
      names() %>%
      as.numeric()

    # From the methods that agree on the number of clusters
    clust_results_filtered <- clust_results[sapply(clust_results, function(x) x[["Best.nc"]][["Number_clusters"]] == nclust)]

    # Get clustering for each method
    cluster_labels_df <- clust_results_filtered %>%
      sapply(function(x) x[["Best.partition"]])

    if (nclust > 1) {
      # For each method (column) get consensus cluster assignment
      cluster_labels_unsup <- apply(cluster_labels_df, 1, function(x) names(which.max(table(x)))) %>%
        as.factor()

      args_list[["cluster_labels"]] <- cluster_labels_unsup
      clus_scores <- do.call(calc_scores, args_list)
      clus_scores_summary <- lapply(clus_scores[["scores"]], function (x) {x[["summary"]]})
      clus_scores_summary <- unlist(clus_scores_summary[!is.null(clus_scores_summary)])
      score_title <- paste0(paste(scores, "=", round(clus_scores_summary, 3)), collapse = "\n")

      clus_scores_list[["hclust"]] <- clus_scores
      results[["clustering"]][["hclust_labels"]] <- cluster_labels_unsup
      pca_plot_list[["hclust"]] <- plot_pca(feat_mat,
                                            color_cluster_by = cluster_labels_unsup,
                                            invisible = "var",
                                            label = "none",
                                            addEllipses = TRUE) +
        ggplot2::ggtitle("Hierarchical", score_title) +
        ggplot2::theme(legend.position="none")
    }


    ## PAM clustering ###############################################

    if (n_clust == "auto") {
      # Find number of clusters with PAM and silhouette method

      p <-
        factoextra::fviz_nbclust(x = feat_mat,
                                 FUNcluster = cluster::pam,
                                 method = "silhouette",
                                 k.max = max_nc,
                                 print.summary = TRUE) +
        theme_minimal() +
        ggtitle("Clusters suggested by\nK-Medoids (Partitioning Around Medoids) and silhouette method")
      # results[["plots"]][["number_of_clusters"]] <- p
      nclust <- as.numeric(p$data$clusters[which.max(p$data$y)])
    } else {
      if (length(cluster_labels) < n_clust) {
        nclust <- length(cluster_labels)
      } else {
        nclust <- n_clust
      }
    }

    cluster_labels_unsup <- cluster::pam(feat_mat, nclust, cluster.only=TRUE, nstart = 30)

    if (nclust > 1) {
      args_list[["cluster_labels"]] <- cluster_labels_unsup
      clus_scores <- do.call(calc_scores, args_list)
      clus_scores_summary <- lapply(clus_scores[["scores"]], function (x) {x[["summary"]]})
      clus_scores_summary <- unlist(clus_scores_summary[!is.null(clus_scores_summary)])
      score_title <- paste0(paste(scores, "=", round(clus_scores_summary, 3)), collapse = "\n")

      clus_scores_list[["pam"]] <- clus_scores
      results[["clustering"]][["pam_labels"]] <- cluster_labels_unsup
      pca_plot_list[["pam"]] <- plot_pca(feat_mat,
                                         color_cluster_by = cluster_labels_unsup,
                                         invisible = "var",
                                         label = "none",
                                         addEllipses = TRUE) +
        ggplot2::ggtitle("PAM", score_title) +
        ggplot2::theme(legend.position="none")
    }


    ## GMM clustering ###############################################

    # use nclust from PAM
    cluster_labels_unsup <- mclust::Mclust(feat_mat,
                                           G = nclust,
                                           verbose = FALSE)$classification

    if (nclust > 1) {
      args_list[["cluster_labels"]] <- cluster_labels_unsup
      clus_scores <- do.call(calc_scores, args_list)
      clus_scores_summary <- lapply(clus_scores[["scores"]], function (x) {x[["summary"]]})
      clus_scores_summary <- unlist(clus_scores_summary[!is.null(clus_scores_summary)])
      score_title <- paste0(paste(scores, "=", round(clus_scores_summary, 3)), collapse = "\n")

      clus_scores_list[["gmm"]] <- clus_scores
      results[["clustering"]][["gmm_labels"]] <- cluster_labels_unsup
      pca_plot_list[["gmm"]] <- plot_pca(feat_mat,
                                         color_cluster_by = cluster_labels_unsup,
                                         invisible = "var",
                                         label = "none",
                                         addEllipses = TRUE) +
        ggplot2::ggtitle("GMM", score_title) +
        ggplot2::theme(legend.position="none")
    }


    ## Leiden clustering ###############################################

    adj_graph <- igraph::as.undirected(cccd::nng(feat_mat, k = knn_k))
    r <- quantile(igraph::strength(adj_graph))[2] / (igraph::gorder(adj_graph) - 1) / 4
    cluster_labels_unsup <- igraph::cluster_leiden(adj_graph, resolution_parameter = r)$membership

    nclust <- length(unique(cluster_labels_unsup))

    # Try other resolutions if only one or more than 3 clusters were found
    if (nclust == 1) {
      r <- quantile(igraph::strength(adj_graph))[2] / (igraph::gorder(adj_graph) - 1) / 2
      cluster_labels_unsup <- igraph::cluster_leiden(adj_graph, resolution_parameter = r)$membership
    }
    if (nclust > 3) {
      r <- quantile(igraph::strength(adj_graph))[2] / (igraph::gorder(adj_graph) - 1) / 8
      cluster_labels_unsup <- igraph::cluster_leiden(adj_graph, resolution_parameter = r)$membership

      if (length(unique(cluster_labels_unsup))) {
        r <- quantile(igraph::strength(adj_graph))[2] / (igraph::gorder(adj_graph) - 1) / 4
        cluster_labels_unsup <- igraph::cluster_leiden(adj_graph, resolution_parameter = r)$membership
      }
    }

    names(cluster_labels_unsup) <- row.names(feat_mat)

    nclust <- length(unique(cluster_labels_unsup))

    if (nclust > 1) {
      args_list[["cluster_labels"]] <- cluster_labels_unsup
      clus_scores <- do.call(calc_scores, args_list)
      clus_scores_summary <- lapply(clus_scores[["scores"]], function (x) {x[["summary"]]})
      clus_scores_summary <- unlist(clus_scores_summary[!is.null(clus_scores_summary)])
      score_title <- paste0(paste(scores, "=", round(clus_scores_summary, 3)), collapse = "\n")

      clus_scores_list[["leiden"]] <- clus_scores
      results[["clustering"]][["leiden_labels"]] <- cluster_labels_unsup
      pca_plot_list[["leiden"]] <- plot_pca(feat_mat,
                                            color_cluster_by = cluster_labels_unsup,
                                            invisible = "var",
                                            label = "none",
                                            addEllipses = TRUE) +
        ggplot2::ggtitle("Leiden", score_title) +
        ggplot2::theme(legend.position="none")
    }

    ## Combine scores


    ## Combine clustering plots ###############################################

    results[["plots"]][["combined"]] <-
      patchwork::wrap_elements(
        patchwork::wrap_plots(pca_plot_list,
                              ncol = 2) +
          patchwork::plot_layout(widths = 1) +
          patchwork::plot_annotation("Unsupervised clustering",
                                     theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                                                               hjust = 0,
                                                                                               size = 20)))
      )


    # Plot hclust heatmap ###############################################

    # # Subset metadata to samples present in feat_mat
    # row.names(metadata) <- metadata[["scoot_sample"]]
    # metadata <- metadata[row.names(feat_mat), ]
    #
    # if (show_clustering) {
    #   # Add clustering labels
    #   for (c in names(results[["clustering"]])) {
    #     metadata[, c] <- results[["clustering"]][[c]]
    #   }
    #   heatmap_cols <- c(heatmap_cols, names(results[["clustering"]]))
    #
    #   metadata <- metadata[, heatmap_cols]
    #
    #   # Remove metadata columns which have different values for each sample
    #   dropcols <- apply(metadata, 2, function(y) length(unique(y))) == nrow(metadata)
    #   dropcols_names <- names(dropcols)[dropcols == TRUE]
    #   keepcols_names <- !(colnames(metadata) %in% dropcols_names)
    #   metadata <- metadata[, keepcols_names]
    #
    #
    #   # Reorder metadata for annotation_rows
    #   annotations <- metadata[row.names(corMat), ]
    #   ha <- ComplexHeatmap::rowAnnotation(df = annotations)
    # } else {
    #   ha <- NULL
    # }

    # ha <- NULL
    #
    # # Plot correlation heatmap
    # results[["plots"]][["hclust_hmap"]] <- ComplexHeatmap::Heatmap(scale(t(feat_mat)),
    #                                                                top_annotation = ha
    # ) %>%
    #   ggplotify::as.ggplot() +
    #   ggplot2::ggtitle("Heatmap hierarchically clustered")



#     corMat <- cor(t(feat_mat), method = "pearson")
#
#     if (show_clustering) {
#       # Reorder metadata for annotation_rows
#       annotations <- metadata[row.names(corMat), ]
#       ha <- ComplexHeatmap::rowAnnotation(df = annotations)
#     }
#
#     # Plot correlation heatmap
#     results[["plots"]][["corr_hmap"]] <- ComplexHeatmap::Heatmap(corMat,
#                                                                  left_annotation = ha
#     ) %>%
#       ggplotify::as.ggplot() +
#       ggplot2::ggtitle("Correlation heatmap")
#
#     # Plot correlation heatmap sorted by dist_matrix
#     results[["plots"]][["corr_hmap_distsort"]] <- ComplexHeatmap::Heatmap(corMat,
#                                                                           clustering_distance_columns = as.dist(dist_mat),
#                                                                           clustering_distance_rows = as.dist(dist_mat),
#                                                                           left_annotation = ha
#     ) %>%
#       ggplotify::as.ggplot() +
#       ggplot2::ggtitle("Correlation heatmap sorted by distance matrix")


    # Combine plots ###############################################

    results[["plots"]][["summary_plot"]] <-
      patchwork::wrap_plots(results[["plots"]][["pca_loadings"]],
                            results[["plots"]][["pca_samples"]],
                            # results[["plots"]][["hclust_hmap"]],
                            results[["plots"]][["combined"]],
                            ncol = 2) +
      patchwork::plot_layout(widths = 1) +
      patchwork::plot_annotation(title,
                                 theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                                                           hjust = 0,
                                                                                           size = 20)))
    ggplot2::set_last_plot(NULL)


    # Combine scores ###############################################

    for (s in scores) {
      avg_summary <- clus_scores_list %>%
        lapply(function (x) {x[["scores"]][[s]][["summary"]]}) %>%
        unlist() %>%
        mean(na.rm = TRUE)
      avg_n <- clus_scores_list %>%
        lapply(function (x) {x[["scores"]][[s]][["n"]]}) %>%
        unlist() %>%
        mean(na.rm = TRUE)
      avg_p_value <- clus_scores_list %>%
        lapply(function (x) {x[["scores"]][[s]][["p_value"]]}) %>%
        unlist() %>%
        mean(na.rm = TRUE)

      results[["scores"]][[s]] <- list("summary" = avg_summary,
                                       "n" = avg_n,
                                       "p_value" = avg_p_value)
    }


    return(results)

  } else {

    return(NULL)
  }
}



#' @importFrom cluster pam
#' @importFrom ggpubr get_legend as_ggplot
#' @importFrom factoextra fviz_pca fviz_nbclust
#' @importFrom scran buildKNNGraph
#' @importFrom igraph modularity set_vertex_attr layout_nicely V strength gorder
#' @importFrom patchwork wrap_plots plot_layout plot_annotation wrap_elements plot_spacer
#' @importFrom ggraph ggraph geom_edge_link geom_node_point

get_scores_sup <- function(data,
                           cluster_labels,
                           scores,
                           knn_k,
                           ntests,
                           seed,
                           title = "", # Title for summary plot
                           select.var,
                           heatmap_annotations_title) {

  feat_mat <- data[["feature_matrix"]]
  dist_mat <- data[["distance_matrix"]]

  if (knn_k == "auto") {
    knn_k <- round(sqrt(nrow(feat_mat)))
    if (knn_k < 3) {
      knn_k <- 3
    }
  }

  # Check if there are at least 2 clusters,
  # that there are at least 6 samples and
  # that there are more samples than clusters (not each sample is one separate cluster)
  if (length(unique(cluster_labels)) > 1 &
      length(cluster_labels) >= 6 &
      nrow(feat_mat) > length(unique(cluster_labels))) {

    # Supervised: Calculate scores + plots ###############################################

    # Args for calc_score
    args_list <- list(feat_mat = feat_mat,
                      dist_mat = dist_mat,
                      cluster_labels = cluster_labels,
                      scores = scores,
                      ntests = ntests,
                      seed = seed,
                      knn_k = knn_k,
                      show_plot = TRUE)

    results <- do.call(calc_scores, args_list)


    # Plot PCA ###############################################

    results[["plots"]][["pca_loadings"]] <- plot_pca(feat_mat,
                                                     color_cluster_by = cluster_labels,
                                                     label = "var",
                                                     labelsize = 3,
                                                     select.var = select.var) +
      ggplot2::ggtitle("PCA showing loadings")

    results[["plots"]][["pca_samples"]] <- plot_pca(feat_mat,
                                                    color_cluster_by = cluster_labels,
                                                    invisible = "var",
                                                    select.var = select.var) +
      ggplot2::ggtitle("PCA showing sample names")


    # # Plot hclust heatmap ###############################################
    #
    # ha <- data.frame(column1 = cluster_labels)
    # row.names(ha) <- names(cluster_labels)
    # colnames(ha) <- heatmap_annotations_title
    #
    # # Plot correlation heatmap
    # results[["plots"]][["hclust_hmap"]] <- ComplexHeatmap::Heatmap(scale(t(feat_mat)),
    #                                                                top_annotation = ha
    # ) %>%
    #   ggplotify::as.ggplot() +
    #   ggplot2::ggtitle("Heatmap hierarchically clustered")


    # Combine plots ###############################################

    results[["plots"]][["summary_plot"]] <- patchwork::wrap_plots(results[["plots"]][["pca_loadings"]],
                                                                  results[["plots"]][["pca_samples"]],
                                                                  # results[["plots"]][["hclust_hmap"]],
                                                                  results[["plots"]][["combined"]],
                                                                  ncol = 2) +
      patchwork::plot_layout(widths = 1) +
      patchwork::plot_annotation(title,
                                 theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                                                           hjust = 0,
                                                                                           size = 20)))
    ggplot2::set_last_plot(NULL)

    return(results)

  } else {

    return(NULL)
  }
}





## get_scores helpers ##############################################################################

### Calculate silhouette score of group vs all others (instead of nearest other group)

calc_sil_iso <- function(dist_mat,
                         labels,
                         return_mean = TRUE,
                         ...) {

  sils <- lapply(
    unique(labels),
    function(a) {

      clus1 <- ifelse(labels == a, 1, 2)

      sil <- cluster::silhouette(clus1, dist_mat)

      # silhouette score for each sample
      # change names back to character
      sil_clus1 <- as.data.frame(sil) %>%
        dplyr::filter(cluster == 1) %>%
        dplyr::mutate(cluster = a) %>%
        dplyr::arrange(dplyr::desc(sil_width))


      size <- nrow(sil_clus1)

      sil.sumA <- data.frame(group = a,
                             size = size,
                             avg_sil_width = mean(sil_clus1$sil_width))

      return(list("samples" = sil_clus1, # silhouette score for each sample
                  "avg_per_group" = sil.sumA))
    }
  )

  # join results
  sils <- lapply(sils, function(x) {x[["samples"]]}) %>%
    data.table::rbindlist() %>%
    dplyr::rename(group = cluster) %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(dplyr::desc(sil_width), .by_group = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rowid = dplyr::row_number())

  if (return_mean) {
    return(mean(sils[["sil_width"]]))
  } else {
    return(sils)
  }
}


### Calculate (classical) silhouette score

calc_sil <- function(dist_mat,
                     labels,
                     return_mean = TRUE,
                     ...) {

  sils <- cluster::silhouette(x = as.numeric(factor(labels)),
                              dist = dist_mat) %>%
    as.data.frame() %>%
    dplyr::rename(group = cluster) %>%
    dplyr::mutate(group = labels) %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(dplyr::desc(sil_width), .by_group = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rowid = dplyr::row_number())

  if (return_mean) {
    return(mean(sils[["sil_width"]]))
  } else {
    return(sils)
  }
}



### Calculate modularity score

calc_modularity <- function(feat_mat,
                            labels,
                            k = 3,
                            ...) {

  results <- list()

  g <- scran::buildKNNGraph(feat_mat,
                            transposed = TRUE,
                            k = k)

  # Calculate modularity score
  results[["score"]] <- igraph::modularity(g,
                                           membership = as.numeric(factor(labels)))
  results[["graph"]] <- g

  return(results)
}


### Permutation test to calculate p-value
# Permutation testing allows for calculating the p-value, e.g. for the silhouette score
# to calculate how likely it is to observe a specific labelling
# For permutation testing labels are randomly shuffled to calculate the random distribution of labels as the null hypothesis

perm_test <- function(fun,
                      args_list,
                      obs,
                      ntests,
                      seed) {

  if (!(length(ntests) == 1 &&
        is.numeric(ntests) &&
        ntests >= 0)) {
    stop("Please provide a number for ntests permutation testing")
  }

  if (!(length(ntests) == 1 &&
        is.numeric(seed))) {
    stop("Please provide a number for setting seed")
  }

  labels <- args_list[["labels"]]

  # perform the shuffling
  if (ntests > 0) {

    random_values <- c()

    for (u in 1:ntests) {

      # seeding for reproducibility
      set.seed(seed + u)
      random_labels <- sample(labels) %>%
        as.factor() %>%
        as.numeric()
      args_list[["labels"]] <- random_labels

      rand_val <- do.call(fun, args_list) %>%
        .[[1]]

      random_values <- c(random_values, rand_val)
    }

    # compute p-value
    p_val <- p_val_zscore(obs = obs,
                          random_values = random_values)

    return(p_val)

  } else if (ntests == 0) {
    return(NULL)
  }
}



### Compute p-value based on z score

p_val_zscore <- function(obs,
                         random_values) {

  mean <- mean(random_values)
  sd <- sd(random_values)
  z_score <- (obs - mean) / sd
  p_val <- 1 - stats::pnorm(z_score)
  return(p_val)
}



### Calculate all scores

calc_scores <- function(feat_mat,
                        dist_mat,
                        cluster_labels,
                        scores,
                        ntests,
                        seed,
                        knn_k,
                        show_plot = FALSE) {

  results <- list()

  suppressMessages({
    suppressWarnings({

      for (s in scores) {

        ## Silhouette_isolated (new) ###############################################
        if (s == "silhouette_isolated") {
          sils <- calc_sil_iso(dist_mat = dist_mat,
                               labels = cluster_labels,
                               return_mean = FALSE)

          avg_per_group <- sils %>%
            dplyr::group_by(group) %>%
            dplyr::summarize(avg_sil_width = mean(sil_width))

          avg_sil <- mean(sils[["sil_width"]])

          ci_intervals <- stats::t.test(sils[["sil_width"]])$conf.int

          args_list <- list(dist_mat = dist_mat,
                            labels = cluster_labels)

          p_val <- perm_test(fun = calc_sil_iso,
                             args_list = args_list,
                             obs = avg_sil,
                             ntests = ntests,
                             seed = seed)


          results[["scores"]][[s]] <- list("samples" = sils,
                                           "avg_per_group" = avg_per_group,
                                           "summary" = avg_sil,
                                           "conf_int" = ci_intervals,
                                           "n" = nrow(sils),
                                           "p_value" = p_val)

          if (show_plot) {
            results[["plots"]][[s]] <- plot_silhouette(sil_scores = results[["scores"]][[s]],
                                                       title = "Silhouette (isolated) plot")
          }
        }

        ## Silhouette (original) ###############################################
        if (s == "silhouette") {
          sils <- calc_sil(dist_mat = dist_mat,
                           labels = cluster_labels,
                           return_mean = FALSE)

          avg_per_group <- sils %>%
            dplyr::group_by(group) %>%
            dplyr::summarize(avg_sil_width = mean(sil_width))

          avg_sil <- mean(sils[["sil_width"]])

          ci_intervals <- t.test(sils[["sil_width"]])$conf.int

          args_list <- list(dist_mat = dist_mat,
                            labels = cluster_labels)

          p_val <- perm_test(fun = calc_sil,
                             args_list = args_list,
                             obs = avg_sil,
                             ntests = ntests,
                             seed = seed)

          results[["scores"]][[s]] <- list("samples" = sils,
                                           "avg_per_group" = avg_per_group,
                                           "summary" = avg_sil,
                                           "conf_int" = ci_intervals,
                                           "n" = nrow(sils),
                                           "p_value" = p_val)

          if (show_plot) {
            results[["plots"]][[s]] <- plot_silhouette(sil_scores = results[["scores"]][[s]],
                                                       title = "Silhouette plot")
          }
        }

        ## Modularity ###############################################
        if (s == "modularity") {

          if (length(cluster_labels) > (knn_k)) {
            mod_res <- calc_modularity(feat_mat = feat_mat,
                                       labels = cluster_labels,
                                       k = knn_k)
            g <- mod_res[["graph"]]

            # Plotting the graph
            g <- igraph::set_vertex_attr(g, "name", value = cluster_labels)

            args_list <- list(feat_mat = feat_mat,
                              labels = cluster_labels,
                              k = knn_k)

            p_val <- perm_test(fun = calc_modularity,
                               args_list = args_list,
                               obs = mod_res[["score"]],
                               ntests = ntests,
                               seed = seed)

            results[["scores"]][[s]] <- list("igraph" = g,
                                             "summary" = mod_res[["score"]],
                                             "n" = length(g),
                                             "p_value" = p_val)

            if (show_plot) {
              # Using ggraph layout = "kk" instead of default "stress",
              # as "kk" handles disconnected communities much better (showing them separately, instead of fur balling them like "stress")
              results[["plots"]][[s]] <- ggraph::ggraph(g, layout = 'kk') +
                ggraph::geom_edge_link(color = "grey", edge_width = 0.2) +
                ggraph::geom_node_point(ggplot2::aes(fill = as.factor(names(igraph::V(g)))),
                                        shape = 21,
                                        color = "black",
                                        size = 3) +
                ggplot2::ggtitle(paste("kNN plot with k = ", knn_k,
                                       "\nModularity score = ", round(mod_res[["score"]], 3),
                                       ifelse(!is.null(p_val),
                                              paste("\np-value:",
                                                    format.pval(p_val, digits = 3)), ""))) +
                ggplot2::labs(fill = "Groups") +
                ggplot2::theme(panel.background = element_rect(fill = "white"))
            }
          }
        }
      }
    })
  })

  if (show_plot) {
    ## Combine scores plots ###############################################

    # Extract the legend and convert to ggplot
    leg <- ggpubr::get_legend(results[["plots"]][[1]])
    results[["plots"]] <- lapply(results[["plots"]],
                                 function (x) x + ggplot2::theme(legend.position="none"))
    results[["plots"]][["leg_plot"]] <- ggpubr::as_ggplot(leg)

    results[["plots"]][["combined"]] <-
      patchwork::wrap_elements(
        patchwork::wrap_plots(results[["plots"]],
                              ncol = 2) +
          patchwork::plot_layout(widths = 1) +
          patchwork::plot_annotation("Supervised scores",
                                     theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                                                               hjust = 0,
                                                                                               size = 20)))
      )
  }

  return(results)
}


### get_scores plot helpers

plot_pca <- function(matrix,
                     color_cluster_by = "none",
                     invisible = "none",
                     label = "all",
                     labelsize = 2,
                     pointsize = 1.5,
                     select.var = NULL,
                     alpha.var = 0.3,
                     addEllipses = FALSE,
                     repel=FALSE) {

  # Remove constant columns with variance = 0
  constant_columns <- apply(matrix, 2, function(col) var(col) == 0)
  constant_columns_indices <- which(constant_columns)
  if (length(constant_columns_indices) > 0) {
    matrix <- matrix[, -constant_columns_indices, drop=FALSE]
  }

  # Otherwise, stats::prcomp would fail with error

  # If there are still columns left
  if (ncol(matrix) > 0) {
    res.pca <- stats::prcomp(matrix,
                             center = TRUE,
                             scale. = FALSE)
    suppressWarnings(
      suppressMessages(
        p <- factoextra::fviz_pca(res.pca,
                                  habillage = color_cluster_by,
                                  invisible = invisible,
                                  label = label,
                                  labelsize = labelsize,
                                  pointsize = pointsize,
                                  select.var = select.var,
                                  alpha.var = alpha.var,
                                  addEllipses = addEllipses,
                                  repel=repel,
                                  mean.point = FALSE) +
          # do not rescale x and y-axes, so scale does not get distorted and represents actual distance better
          coord_equal() +
          # Remove shapes added by fviz_pca
          scale_shape_manual(values = c(rep(19, length(unique(color_cluster_by)))))
      )
    )

    return(p)
  }
}


plot_silhouette <- function(sil_scores,
                            title = "title") {

  xend <- length(sil_scores[["samples"]][["sil_width"]])
  ci <- sil_scores[["conf_int"]]
  m <- sil_scores[["summary"]]

  p <- ggplot2::ggplot(sil_scores[["samples"]],
                       ggplot2::aes(x = rowid,
                                    y = sil_width,
                                    fill = group)) +
    ggplot2::geom_bar(position = "dodge", stat = "identity") +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::ggtitle(paste0(title,
                            "\nAvg sil width = ", round(m, 3),
                            "\n95% CI: [", round(ci[1], 3), ", ", round(ci[2], 3), "]",
                            ifelse(!is.null(sil_scores[["p_value"]]),
                                   paste("\np-value: ",
                                         format.pval(sil_scores[["p_value"]],
                                                     digits = 3)), ""))) +
    ggplot2::labs(fill = "Groups") +
    ggplot2::annotate("ribbon", x = 1:xend, ymin = ci[1], ymax = ci[2], alpha = 0.15) +
    ggplot2::annotate("segment", x = 1, xend = xend, y = ci[1], lty = 1, alpha = 0.15) +
    ggplot2::annotate("segment", x = 1, xend = xend, y = ci[2], lty = 1, alpha = 0.15) +
    ggplot2::annotate("segment", x = 1, xend = xend, y = m, lty = 2) +
    ggplot2::theme_bw()

  return(p)
}



### Combine multiple p-values from get.cluster.scores if batching

combine_pvals <- function(list_pvals_n,
                          pval_combine_method) {

  use_fallback_method <- FALSE

  # To prevent error from metap, replace zeros with some small value
  list_pvals_n[["p_value"]][list_pvals_n[["p_value"]] == 0] <- 1e-30

  if (!"p_value" %in% names(list_pvals_n)) {
    return(list_pvals_n)
  } else
    if (pval_combine_method == "weighted_zmethod") {
      if (!length(list_pvals_n[["p_value"]]) == length(list_pvals_n[["n"]])) {
        use_fallback_method <- TRUE
      } else {
        list_pvals_n[["p_value"]] <- metap::sumz(p = list_pvals_n[["p_value"]],
                                                 weights = list_pvals_n[["n"]])[["p"]][1, 1]
      }
    } else
      if (pval_combine_method == "fisher" |
          use_fallback_method) {
        list_pvals_n[["p_value"]] <- metap::sumlog(p = list_pvals_n[["p_value"]])[["p"]]
      } else {
        stop("pval_combine_method not recognized.
         Please see documentation for possible methods.")
      }

  list_pvals_n[["n"]] <- sum(list_pvals_n[["n"]])

  return(list_pvals_n)
}




# feature_corr_plot helper functions ##############################################################################

#' @importFrom stats cor

# Function to calculate Pearson correlation between distance matrices without dropping or imputing
correlation_between_dms <- function(dm1, dm2) {
  # Extract row and column names
  rows_dm1 <- row.names(dm1)
  cols_dm1 <- colnames(dm1)
  rows_dm2 <- row.names(dm2)
  cols_dm2 <- colnames(dm2)

  # Find common samples
  common_samples <- intersect(rows_dm1, rows_dm2)

  if (length(common_samples) >= 10) {
    # Ensure common samples are in the same order in both matrices
    dm1_common <- dm1[common_samples, common_samples]
    dm2_common <- dm2[common_samples, common_samples]

    corrs <- c()

    # Convert rows to ranks
    dm1_common_ranks <- t(apply(data.frame(dm1_common), 1, rank))
    dm2_common_ranks <- t(apply(data.frame(dm2_common), 1, rank))

    # Flatten matrices into vectors
    dm1_flat <- as.vector(dm1_common_ranks)
    dm2_flat <- as.vector(dm2_common_ranks)

    # Calculate Pearson correlation
    corrs[i] <- stats::cor(dm1_flat, dm2_flat)

    return(corrs)
  }
}


get_dist_mats <- function(scores) {
  data <- scores[["data"]]

  dm_list <- list()

  for (m in names(data)) {
    for (l in names(data[[m]])) {
      if ("distance_matrix" %in% names(data[[m]][[l]])) {
        dm <- data[[m]][[l]][["distance_matrix"]]
        if (!is.null(dm)) {
          dm_list[[paste0(c(m,l), collapse = "_")]] <- dm
        }
      }
      else {
        for (s in names(data[[m]][[l]])) {
          dm <- data[[m]][[l]][[s]][["distance_matrix"]]
          if (!is.null(dm)) {
            dm_list[[paste0(c(m,l,s), collapse = "_")]] <- dm
          }
        }
      }
    }
  }

  return(dm_list)
}




# Plot L1, L2, and L3 summary PCA ##############################################################################

scores_composite_pca <- function(scores,
                                 scoot_summary
) {

  cluster_by <- names(scores)[!names(scores) %in% c("unsupervised", "params", "data")]

  results <- list()

  for (c in cluster_by) {

    # L1 ###############################################

    mat <- scores[["data"]][["composition"]][["layer_1"]][["distance_matrix"]]
    results[[c]][["matrices"]][["L1"]] <- mat

    clust_labels <- scoot_summary@metadata[match(row.names(mat), scoot_summary@metadata[["scoot_sample"]]), ][[c]]
    sil_iso <- calc_sil_iso(dist_mat = as.dist(mat),
                            labels = clust_labels,
                            return_mean = TRUE) %>% round(3)

    results[[c]][["plots"]][["L1"]] <- plot_pca(mat,
                                                label = "none",
                                                invisible = c("var", "quali"),
                                                color_cluster_by = clust_labels) +
      ggtitle(paste("L1 - PCA - low-res cell type composition\nSilhouette (iso) score:", sil_iso))


    # L2 ###############################################

    mats <- list()
    for (i in names(scores[["data"]][["composition"]][["layer_2"]])) {
      mats[[i]] <- scores[["data"]][["composition"]][["layer_2"]][[i]][["distance_matrix"]]
    }

    mat <- get_avg_matrix(mats)

    results[[c]][["matrices"]][["L2"]] <- mat

    clust_labels <- scoot_summary@metadata[match(row.names(mat), scoot_summary@metadata[["scoot_sample"]]), ][[c]]
    sil_iso <- calc_sil_iso(dist_mat = as.dist(mat),
                            labels = clust_labels,
                            return_mean = TRUE) %>% round(3)

    results[[c]][["plots"]][["L2"]] <- plot_pca(mat,
                                                label = "none",
                                                invisible = c("var", "quali"),
                                                color_cluster_by = clust_labels) +
      ggtitle(paste("L2 - PCA - hi-res cell type composition\nSilhouette (iso) score:", sil_iso))


    # L3 ###############################################

    mats <- list()
    for (i in names(scores[["data"]][["signatures"]][["layer_1"]])) {
      mats[[paste0("layer_1_", i)]] <- scores[["data"]][["signatures"]][["layer_1"]][[i]][["distance_matrix"]]
    }
    for (i in names(scores[["data"]][["signatures"]][["layer_2"]])) {
      mats[[paste0("layer_2_", i)]] <- scores[["data"]][["signatures"]][["layer_1"]][[i]][["distance_matrix"]]
    }

    mat <- get_avg_matrix(mats)

    results[[c]][["matrices"]][["L3"]] <- mat

    clust_labels <- scoot_summary@metadata[match(row.names(mat), scoot_summary@metadata[["scoot_sample"]]), ][[c]]
    sil_iso <- calc_sil_iso(dist_mat = as.dist(mat),
                            labels = clust_labels,
                            return_mean = TRUE) %>% round(3)

    results[[c]][["plots"]][["L3"]] <- plot_pca(mat,
                                                label = "none",
                                                invisible = c("var", "quali"),
                                                color_cluster_by = clust_labels) +
      ggtitle(paste("L3 - PCA - signature expression per cell type\nSilhouette (iso) score:", sil_iso))


    # L1_L2_L3_combined ###############################################

    mat <- get_avg_matrix(results[[c]][["matrices"]])

    results[[c]][["matrices"]][["L1_L2_L3_combined"]] <- mat

    clust_labels <- scoot_summary@metadata[match(row.names(mat), scoot_summary@metadata[["scoot_sample"]]), ][[c]]
    sil_iso <- calc_sil_iso(dist_mat = as.dist(mat),
                            labels = clust_labels,
                            return_mean = TRUE) %>% round(3)

    results[[c]][["plots"]][["L1_L2_L3_combined"]] <- plot_pca(mat,
                                                               label = "none",
                                                               invisible = c("var", "quali"),
                                                               color_cluster_by = clust_labels) +
      ggtitle(paste("L1_L2_L3_combined - PCA\nSilhouette (iso) score:", sil_iso))


    # B1 ###############################################

    mat <- scores[["data"]][["pseudobulk"]][["layer_1"]][["all"]][["distance_matrix"]]
    results[[c]][["matrices"]][["B1"]] <- mat

    clust_labels <- scoot_summary@metadata[match(row.names(mat), scoot_summary@metadata[["scoot_sample"]]), ][[c]]
    sil_iso <- calc_sil_iso(dist_mat = as.dist(mat),
                            labels = clust_labels,
                            return_mean = TRUE) %>% round(3)

    results[[c]][["plots"]][["B1"]] <- plot_pca(mat,
                                                label = "none",
                                                invisible = c("var", "quali"),
                                                color_cluster_by = clust_labels) +
      ggtitle(paste("B1 - pseudobulk PCA\nSilhouette (iso) score:", sil_iso))


    # Combine plots ###############################################

    results[[c]][["plots"]][["summary_plot"]] <-
      patchwork::wrap_plots(results[[c]][["plots"]][["B1"]],
                            results[[c]][["plots"]][["L1_L2_L3_combined"]],
                            patchwork::plot_spacer(),
                            results[[c]][["plots"]][["L1"]],
                            results[[c]][["plots"]][["L2"]],
                            results[[c]][["plots"]][["L3"]],
                            nrow = 2) +
      patchwork::plot_layout(widths = 1) +
      patchwork::plot_annotation(c,
                                 theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                                                           hjust = 0,
                                                                                           size = 20)))
    print(results[[c]][["plots"]][["summary_plot"]])
  }

  return(results)
}


## scores_composite_pca helper  ###############################################

get_avg_matrix <- function(mats) {

  # Identify all unique row and column names
  all_row_names <- mats %>% lapply(row.names) %>% unlist() %>% unique()
  all_col_names <- mats %>% lapply(colnames) %>% unlist() %>% unique()

  # Extend each matrix to have all row and column names
  extended_matrices <- lapply(mats,
                              extend_matrix,
                              all_rows = all_row_names,
                              all_cols = all_col_names)

  # Sum all extended matrices
  sum_matrix <- Reduce("+", extended_matrices)

  # Normalize sum_matrix by number of sample appearance
  # E.g. one sample might have values for all matrices (e.g. have values for all 4 subtypes in layer 2)
  # but another sample might be missing one subtype, so divide this one only by 3 (subtypes in layer 2)
  avg_matrix <- sum_matrix
  for (i in colnames(sum_matrix)) {
    len <- 1
    for (j in 1:length(extended_matrices)) {
      if (!all(extended_matrices[[j]][i, ] == 0)) {
        len <- len + 1
      }
    }
    avg_matrix[i, ] <- avg_matrix[i, ] / len
    # avg_matrix[, i] <- avg_matrix[, i] / len
  }

  return(avg_matrix)
}


# Because not all sub cell types are present in all samples, we need to to create a common dist matrix
# Define a function to extend a matrix to include all row and column names
extend_matrix <- function(mat,
                          all_rows,
                          all_cols) {

  # Create a matrix of zeros with the complete set of row and column names
  extended_mat <- matrix(0,
                         nrow = length(all_rows),
                         ncol = length(all_cols),
                         dimnames = list(all_rows,
                                         all_cols))

  # Assign values from the original matrix to the extended matrix
  common_rows <- intersect(rownames(mat), all_rows)
  common_cols <- intersect(colnames(mat), all_cols)
  extended_mat[common_rows, common_cols] <- mat[common_rows, common_cols]

  return(extended_mat)
}


# Convenience functions  ###############################################

#' Save object list to disk, in parallel
#'
#' @param obj.list A list of Seurat objects
#' @param dir File location
#' @param ncores Number of CPU cores to use for parallelization
#' @param progressbar Whether to show a progressbar or not
#'
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallelly availableCores
#' @return A list of Seurat objects saved to disk as separate .rds files
#' @export save_objs
#'
#' @examples
#' save_objs(obj.list, "./output/samples")

save_objs <- function(obj.list,
                      dir,
                      ncores = parallelly::availableCores() - 2,
                      progressbar = TRUE) {

  BiocParallel::bplapply(
    X = names(obj.list),
    BPPARAM =  BiocParallel::MulticoreParam(workers = ncores,
                                            progressbar = progressbar),
    function(x) {
      file_name <- file.path(dir, sprintf("%s.rds", x))
      saveRDS(obj.list[[x]], file_name)
    })
}



#' Read all .rds files in a directory and return list of Seurat objects, in parallel
#'
#' @param dir File location
#' @param file.list A list of files (with full pathname)
#' @param ncores Number of CPU cores to use for parallelization
#' @param progressbar Whether to show a progressbar or not
#'
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom stringr str_remove_all
#' @importFrom parallelly availableCores
#' @return A list of Seurat objects read into R
#' @export read_objs
#'
#' @examples
#' obj.list <- read_objs("./output/samples")

read_objs <- function(dir = NULL,
                      file.list = NULL,
                      ncores = parallelly::availableCores() - 2,
                      progressbar = TRUE) {

  if (!is.null(dir) & is.null(file.list)) {
    file_names <- list.files(dir)
    file_paths <- file.path(dir, file_names)
  } else if (is.null(dir) & !is.null(file.list)) {
    file_paths <- file.list[endsWith(files, '.rds')]
    file_names <- gsub("^.*/", "", file_paths)
  }
  obj.list <- BiocParallel::bplapply(
    X = file_paths,
    BPPARAM =  BiocParallel::MulticoreParam(workers = ncores,
                                            progressbar = progressbar),
    function(x) {
      readRDS(file.path(x))
    })
  names(obj.list) <- stringr::str_remove_all(file_names, '.rds')
  return(obj.list)
}
