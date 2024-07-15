# scoot class definition ##############################################################################

## Build S4 objects to store data
methods::setClass(Class = "scoot",
                  slots = list(
                    metadata = "data.frame",
                    composition = "list",
                    aggregated_profile = "list",
                    version = "list"
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
      param <- BiocParallel::SerialParam()
    }
  } else {
    param <- bparam
  }
  return(param)
}





compositional_data <- function(data,
                               split_by = NULL,
                               group_by_1 = NULL,
                               useNA = FALSE,
                               clr_zero_impute_perc = 1,
                               only_counts = FALSE) {

  if (all(is.na(data[[group_by_1]]))) {
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
    gr_vars <- c(split_by, group_by_1)
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


# get.cluster.score helper functions ##############################################################################

### Pre-process pseudobulk count data

preproc_pseudobulk <- function(matrix,
                               metadata,
                               cluster_by,
                               nvar_genes = 500,
                               black_list = NULL) {

  suppressMessages({
    suppressWarnings({
      matrix <- DESeq2.normalize(matrix = matrix,
                                 metadata = metadata,
                                 cluster_by = cluster_by,
                                 nvar_genes = nvar_genes,
                                 black_list = black_list)
    })
  })

  return(matrix)
}



### Just a helper for preproc_pseudobulk, additional pre-processing steps are needed

DESeq2.normalize <- function(matrix,
                             metadata,
                             cluster_by,
                             nvar_genes = 500,
                             black_list = NULL) {

  # Get black list
  if (is.null(black_list)) {
    utils::data("default_black_list")
  }
  black_list <- unlist(black_list)

  # Normalize pseudobulk data using DESeq2
  # do formula for design with the cluster_by elements in order
  dformula <-  stats::formula(paste("~", cluster_by))
  data <- DESeq2::DESeqDataSetFromMatrix(countData = matrix,
                                         colData = metadata,
                                         design = dformula)
  data <- DESeq2::estimateSizeFactors(data)

  nsub <- min(1000, sum(rowMeans(BiocGenerics::counts(data, normalized=TRUE)) > 10 ))

  # transform counts using vst
  data <- DESeq2::vst(data, blind = T, nsub = nsub)
  data <- SummarizedExperiment::assay(data)

  # Remove black listed genes from the matrix
  data <- data[!row.names(data) %in% black_list,]

  # get top variable genes
  rv <- MatrixGenerics::rowVars(data)
  select <- order(rv, decreasing=TRUE)[seq_len(min(nvar_genes, length(rv)))]
  select <- row.names(data)[select]

  data <- data[select[select %in% row.names(data)],]

  return(data)
}




get_scores <- function(matrix,
                       cluster_labels,
                       scores,
                       modularity_k,
                       dist_method = "euclidean",
                       ntests = 100, # number of shuffling events
                       seed = 22, # seed for random shuffling
                       title = "", # Title for summary
                       # For PCA
                       invisible = c("var", "quali")) {

  matrix <- t(matrix)

  results <- list()

  # Check if there are at least 2 clusters,
  # that there are more than 4 samples and
  # that there are more samples than clusters (not each sample is one separate cluster)
  if (length(unique(cluster_labels)) > 1 &
      length(cluster_labels) > 4 &
      nrow(matrix) > length(unique(cluster_labels))) {

    results[["feature_matrix"]] <- matrix
    results[["distance_matrix"]] <- stats::dist(matrix, method = dist_method)
    mat <- as.matrix(results[["distance_matrix"]])

    # Plot PCA ###############################################
    results[["plots"]][["pca"]] <- plot_pca(matrix,
                                            color_cluster_by = cluster_labels,
                                            label = "var",
                                            invisible = invisible) +
      ggplot2::ggtitle("PCA")


    # Clustering ###############################################

    ## PAM clustering ###############################################

    # Find number of clusters with PAM and silhouette method
    if (length(cluster_labels) <= 10) {
      k_max <- length(cluster_labels) - 1
    } else {
      k_max <- 10
    }

    p <-
      factoextra::fviz_nbclust(x = mat,
                               FUNcluster = cluster::pam,
                               method = "silhouette",
                               k.max = k_max,
                               print.summary = TRUE) +
      theme_minimal() +
      ggtitle("Clusters suggested by\nK-Medoids (Partitioning Around Medoids) and silhouette method")
    # results[["plots"]][["number_of_clusters"]] <- p
    nclust <- as.numeric(p$data$clusters[which.max(p$data$y)])

    clusters <- cluster::pam(matrix, nclust, cluster.only=TRUE, nstart = 30)

    # Plot PCA
    if (nclust == 1) {
      results[["plots"]][["pca_pam_clustered"]] <- plot_pca(matrix,
                                                            label = "var",
                                                            invisible = invisible)
    } else {
      results[["plots"]][["pca_pam_clustered"]] <- plot_pca(matrix,
                                                            label = "var",
                                                            invisible = invisible,
                                                            color_cluster_by = clusters,
                                                            add_ellipses = TRUE) +
        ggplot2::ggtitle("PCA - PAM clustered")
    }


    ## Leiden clustering ###############################################

    adj_graph <- igraph::as.undirected(cccd::nng(matrix, k = 3))
    r <- quantile(strength(adj_graph))[2] / (gorder(adj_graph) - 1) / 4
    leiden_clusters <- igraph::cluster_leiden(adj_graph, resolution_parameter = r)$membership

    # Plot PCA
    if (length(unique(leiden_clusters)) == 1) {
      results[["plots"]][["pca_leiden_clustered"]] <- plot_pca(matrix,
                                                               label = "var",
                                                               invisible = invisible)
    } else {
      results[["plots"]][["pca_leiden_clustered"]] <- plot_pca(matrix,
                                                               label = "var",
                                                               invisible = invisible,
                                                               color_cluster_by = leiden_clusters,
                                                               add_ellipses = TRUE) +
        ggplot2::ggtitle("PCA - Leiden clustered")
    }

    # Calculate scores + plots ###############################################

    suppressMessages({
      suppressWarnings({

        for (s in scores) {

          ## Silhouette_isolated (new) ###############################################
          if (s == "silhouette_isolated") {
            sils <- calc_sil_onelabel(labels = cluster_labels,
                                      dist = results[["distance_matrix"]],
                                      return_mean_for_permtest = FALSE)

            avg_per_group <- sils %>%
              dplyr::group_by(group) %>%
              dplyr::summarize(avg_sil_width = mean(sil_width))

            avg_sil <- mean(sils[["sil_width"]])

            ci_intervals <- stats::t.test(sils[["sil_width"]])$conf.int

            p_val <- perm_test(fun = calc_sil_onelabel,
                               data = results[["distance_matrix"]],
                               labels = cluster_labels,
                               obs = avg_sil,
                               ntests = ntests,
                               seed = seed)


            results[["scores"]][[s]] <- list("samples" = sils,
                                             "avg_per_group" = avg_per_group,
                                             "summary" = avg_sil,
                                             "conf_int" = ci_intervals,
                                             "n" = nrow(sils),
                                             "p_value" = p_val)

            p <- plot_silhouette(sil_scores = results[["scores"]][[s]],
                                 title = "Silhouette (isolated) plot")

            results[["plots"]][[s]] <- p
          }

          ## Silhouette (original) ###############################################
          if (s == "silhouette") {
            sils <- calc_sil(labels = cluster_labels,
                             dist = results[["distance_matrix"]],
                             return_mean_for_permtest = FALSE)

            avg_per_group <- sils %>%
              dplyr::group_by(group) %>%
              dplyr::summarize(avg_sil_width = mean(sil_width))

            avg_sil <- mean(sils[["sil_width"]])

            ci_intervals <- t.test(sils[["sil_width"]])$conf.int

            p_val <- perm_test(fun = calc_sil,
                               data = results[["distance_matrix"]],
                               labels = cluster_labels,
                               obs = avg_sil,
                               ntests = ntests,
                               seed = seed)

            results[["scores"]][[s]] <- list("samples" = sils,
                                             "avg_per_group" = avg_per_group,
                                             "summary" = avg_sil,
                                             "conf_int" = ci_intervals,
                                             "n" = nrow(sils),
                                             "p_value" = p_val)

            p <- plot_silhouette(sil_scores = results[["scores"]][[s]],
                                 title = "Silhouette plot")

            results[["plots"]][[s]] <- p
          }

          ## Modularity ###############################################
          if (s == "modularity") {

            if (length(cluster_labels) >= (modularity_k + 1)) {
              g <- scran::buildKNNGraph(matrix,
                                        transposed = TRUE,
                                        k = modularity_k)

              # Calculate modularity score
              modularity_score <- igraph::modularity(g, membership = as.numeric(factor(cluster_labels)))

              # Plotting the graph
              g <- igraph::set_vertex_attr(g, "name", value = cluster_labels)

              p_val <- perm_test(fun = calc_modularity,
                                 data = matrix,
                                 labels = cluster_labels,
                                 obs = modularity_score,
                                 ntests = ntests,
                                 seed = seed)

              results[["scores"]][[s]] <- list("igraph" = g,
                                               "summary" = modularity_score,
                                               "n" = length(g),
                                               "p_value" = p_val)

              # Using ggraph layout = "kk" instead of default "stress",
              # as "kk" handles disconnected communities much better (showing them separately, instead of fur balling them like "stress")
              p <- ggraph::ggraph(g, layout = 'kk') +
                ggraph::geom_edge_link(color = "grey", edge_width = 0.2) +
                ggraph::geom_node_point(ggplot2::aes(fill = as.factor(names(igraph::V(g)))),
                                        shape = 21,
                                        color = "black",
                                        size = 5) +
                ggplot2::ggtitle(paste("KNN plot with k = ", modularity_k,
                                       "\nModularity score = ", round(modularity_score, 3),
                                       ifelse(!is.null(p_val),
                                              paste("\np-value:",
                                                    format.pval(p_val, digits = 3)), ""))) +
                ggplot2::labs(fill = "Groups") +
                ggplot2::theme(panel.background = element_rect(fill = "white"))

              results[["plots"]][[s]] <- p
            }
          }
        }
      })
    })

    # Plot dendrogram  (TODO NEEDS REWORK) ###############################################
    # pc2 <- pc$x[,1:ndim] %>%
    #   as.data.frame() %>%
    #   mutate(celltype = metadata[[df.score[x, 1]]])
    #
    # # Calculate the row-wise average grouping by row names
    # pc2 <- aggregate(. ~ celltype, data = pc2, FUN = mean) %>%
    #   tibble::column_to_rownames("celltype") %>%
    #   as.matrix()
    #
    # dist_group <- stats::dist(pc2,
    #                           method = df.score[x, 2])
    # hclust <- stats::hclust(dist_group,
    #                         method = hclust_method)
    # dendo <- ggdendro::ggdendrogram(as.dendrogram(hclust)) +
    #   ggtitle(paste0("Hierarchical clustering dendrogram - ",
    #                  hclust_method))


    # Combine plots ###############################################
    results[["plots"]][["summary_plot"]] <- patchwork::wrap_plots(results[["plots"]],
                                                                  ncol = 2) +
      patchwork::plot_layout(widths = 1) +
      patchwork::plot_annotation(title,
                                 theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                                                           hjust = 0,
                                                                                           size = 20)))
    return(results)

  } else {

    return(NULL)
  }
}





## get_scores helpers ##############################################################################

### Calculate silhouette score of group vs all others (instead of nearest other group)

calc_sil_onelabel <- function(labels,
                              dist,
                              return_mean_for_permtest = TRUE) {

  sils <- lapply(
    unique(labels),
    function(a) {

      clus1 <- ifelse(labels == a, 1, 2)

      sil <- cluster::silhouette(clus1, dist)

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

  if (return_mean_for_permtest) {
    return(mean(sils[["sil_width"]]))
  } else {
    return(sils)
  }
}


### Calculate (classical) silhouette score

calc_sil <- function(labels,
                     dist,
                     return_mean_for_permtest = TRUE) {

  sils <- cluster::silhouette(x = as.numeric(factor(labels)),
                              dist = dist) %>%
    as.data.frame() %>%
    dplyr::rename(group = cluster) %>%
    dplyr::mutate(group = labels) %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(dplyr::desc(sil_width), .by_group = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rowid = dplyr::row_number())

  if (return_mean_for_permtest) {
    return(mean(sils[["sil_width"]]))
  } else {
    return(sils)
  }
}



### Calculate modularity score

calc_modularity <- function(labels,
                            matrix,
                            k = 3) {

  g <- scran::buildKNNGraph(matrix,
                            transposed = TRUE,
                            k = k)

  # Calculate modularity score
  modularity_score <- igraph::modularity(g, membership = as.numeric(factor(labels)))

  return(modularity_score)
}


### Permutation test to calculate p-value
# Permutation testing allows for calculating the p-value, e.g. for the silhouette score
# to calculate how likely it is to observe a specific labelling
# For permutation testing labels are randomly shuffled to calculate the random distribution of labels as the null hypothesis

perm_test <- function(fun,
                      data,
                      labels,
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

  # perform the shuffling
  if (ntests > 0) {

    random_values <- c()

    for (u in 1:ntests) {

      # seeding for reproducibility
      set.seed(seed + u)
      random_labels <- sample(labels) %>%
        as.factor() %>%
        as.numeric()

      rand_val <- fun(random_labels, data) %>%
        mean()

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




### get_scores plot helpers

plot_pca <- function(matrix,
                     scale = FALSE,
                     label = "all",
                     invisible = c("var", "quali"),
                     geom_var = c("arrow", "text"),
                     col_var = "steelblue",
                     alpha_var = 0.3,
                     repel = TRUE,
                     color_cluster_by = "none",
                     add_ellipses = FALSE) {

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
                             scale. = scale)
    suppressWarnings(
      suppressMessages(
        p <- factoextra::fviz_pca(res.pca,
                                  habillage = color_cluster_by,
                                  addEllipses = add_ellipses,
                                  label = label,
                                  pointsize = 3,
                                  invisible = invisible,
                                  geom.var = geom_var,
                                  col.var = col_var,
                                  alpha.var = alpha_var,
                                  repel = repel) +
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
                            "\nAverage silhouette width = ", round(m, 3),
                            "   95% CI: [", round(ci[1], 3), ", ", round(ci[2], 3), "]",
                            ifelse(!is.null(sil_scores[["p_value"]]),
                                   paste("\np-value: ",
                                         format.pval(sil_scores[["p_value"]],
                                                     digits = 3)), ""))) +
    ggplot2::labs(fill = "Groups") +
    ggplot2::geom_ribbon(aes(x = 1:xend,
                             ymin = ci[1],
                             ymax = ci[2]),
                         alpha = 0.15,
                         inherit.aes = FALSE) +
    ggplot2::annotate("segment",
                      x = 1,
                      xend = xend,
                      y = ci[1],
                      lty = 1,
                      alpha = 0.15) +
    ggplot2::annotate("segment",
                      x = 1,
                      xend = xend,
                      y = ci[2],
                      lty = 1,
                      alpha = 0.15) +
    ggplot2::annotate("segment",
                      x = 1,
                      xend = xend,
                      y = m,
                      lty = 2) +
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
