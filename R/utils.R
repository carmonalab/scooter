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


# get_cluster_score helper functions ##############################################################################

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



### Just a helper for preproc_pseudobulk, if additional pre-processing steps are needed

#' @importFrom DESeq2 DESeqDataSetFromMatrix vst estimateSizeFactors
#' @importFrom SummarizedExperiment assay

DESeq2.normalize <- function(matrix,
                             metadata,
                             cluster_by,
                             nvar_genes = 500,
                             black_list = NULL) {

  # Get black list
  if (is.null(black_list)) {
    data(default_black_list, envir=environment())
  }
  black_list <- unlist(black_list)

  # Remove black listed genes from the matrix
  matrix <- matrix[!row.names(matrix) %in% black_list,]

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

  # get top variable genes
  rv <- MatrixGenerics::rowVars(data)
  select <- order(rv, decreasing=TRUE)[seq_len(min(nvar_genes, length(rv)))]
  select <- row.names(data)[select]

  data <- data[select[select %in% row.names(data)],]

  return(data)
}


#' @importFrom cluster pam
#' @importFrom NbClust NbClust
#' @importFrom ggpubr get_legend as_ggplot
#' @importFrom factoextra fviz_pca fviz_nbclust
#' @importFrom scran buildKNNGraph
#' @importFrom igraph modularity set_vertex_attr layout_nicely V strength gorder
#' @importFrom patchwork wrap_plots plot_layout plot_annotation wrap_elements plot_spacer
#' @importFrom ggraph ggraph geom_edge_link geom_node_point

get_scores <- function(matrix,
                       cluster_labels,
                       scores,
                       modularity_k,
                       max_nc,
                       pam_nclusters,
                       hdbscan_min_pts,
                       dist_method = "euclidean",
                       ntests = 100, # number of shuffling events
                       seed = 22, # seed for random shuffling
                       title = "", # Title for summary
                       # For PCA
                       invisible = c("var", "quali"),
                       select_var = NULL) {

  matrix <- t(matrix)

  results <- list()

  # Check if there are at least 2 clusters,
  # that there are more than 4 samples and
  # that there are more samples than clusters (not each sample is one separate cluster)
  if (length(unique(cluster_labels)) > 1 &
      length(cluster_labels) > 4 &
      nrow(matrix) > length(unique(cluster_labels))) {

    results[["feature_matrix"]] <- matrix
    results[["distance_matrix"]] <- as.matrix(stats::dist(matrix, method = dist_method))
    # Do not store dist object, as this can lead to the following error when trying to view it in RStudio (depending on which packages are loaded):
    # Error in .Primitive("[")(x, 1:6, , drop = FALSE) : incorrect number of dimensions

    # Plot PCA ###############################################
    results[["plots"]][["pca_loadings"]] <- plot_pca(matrix,
                                                     color_cluster_by = cluster_labels,
                                                     label = "var",
                                                     invisible = invisible,
                                                     select_var = select_var) +
      ggplot2::ggtitle("PCA showing loadings")

    results[["plots"]][["pca_samples"]] <- plot_pca(matrix,
                                                    color_cluster_by = cluster_labels) +
      ggplot2::ggtitle("PCA showing sample names")

    ## Plot PCA scree plot ###############################################
    #
    # # Perform PCA
    # pca_result <- stats::prcomp(matrix,
    #                             center = TRUE,
    #                             scale. = FALSE)
    #
    # # Variance explained by each principal component
    # variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
    # pc_numbers <- 1:length(variance_explained)
    # scree_data <- data.frame(PC = pc_numbers, VarianceExplained = variance_explained)
    #
    # # Plot scree plot
    # results[["plots"]][["pca_scree"]] <- ggplot(scree_data,
    #                                             aes(x = PC,
    #                                                 y = VarianceExplained)) +
    #   geom_point(color = "blue", size = 3) +
    #   geom_line(color = "blue", size = 1) +
    #   labs(x = "Principal Component",
    #        y = "Proportion of Variance Explained",
    #        title = "Scree Plot") +
    #   theme_minimal()


    # Supervised: Calculate scores + plots ###############################################

    scores_plot_list <- list()

    suppressMessages({
      suppressWarnings({

        for (s in scores) {

          ## Silhouette_isolated (new) ###############################################
          if (s == "silhouette_isolated") {
            sils <- calc_sil_onelabel(labels = cluster_labels,
                                      dist = as.dist(results[["distance_matrix"]]),
                                      return_mean_for_permtest = FALSE)

            avg_per_group <- sils %>%
              dplyr::group_by(group) %>%
              dplyr::summarize(avg_sil_width = mean(sil_width))

            avg_sil <- mean(sils[["sil_width"]])

            ci_intervals <- stats::t.test(sils[["sil_width"]])$conf.int

            p_val <- perm_test(fun = calc_sil_onelabel,
                               data = as.dist(results[["distance_matrix"]]),
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

            scores_plot_list[[s]] <- p
          }

          ## Silhouette (original) ###############################################
          if (s == "silhouette") {
            sils <- calc_sil(labels = cluster_labels,
                             dist = as.dist(results[["distance_matrix"]]),
                             return_mean_for_permtest = FALSE)

            avg_per_group <- sils %>%
              dplyr::group_by(group) %>%
              dplyr::summarize(avg_sil_width = mean(sil_width))

            avg_sil <- mean(sils[["sil_width"]])

            ci_intervals <- t.test(sils[["sil_width"]])$conf.int

            p_val <- perm_test(fun = calc_sil,
                               data = as.dist(results[["distance_matrix"]]),
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

            scores_plot_list[[s]] <- p
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
                                        size = 3) +
                ggplot2::ggtitle(paste("kNN plot with k = ", modularity_k,
                                       "\nModularity score = ", round(modularity_score, 3),
                                       ifelse(!is.null(p_val),
                                              paste("\np-value:",
                                                    format.pval(p_val, digits = 3)), ""))) +
                ggplot2::labs(fill = "Groups") +
                ggplot2::theme(panel.background = element_rect(fill = "white"))

              scores_plot_list[[s]] <- p
            }
          }
        }
      })
    })

    ## Combine scores plots ###############################################

    # Extract the legend and convert to ggplot
    leg <- ggpubr::get_legend(scores_plot_list[[1]])
    scores_plot_list <- lapply(scores_plot_list, function (x) x + ggplot2::theme(legend.position="none"))
    scores_plot_list[["leg_plot"]] <- ggpubr::as_ggplot(leg)

    results[["plots"]][["scores_combined"]] <-
      patchwork::wrap_elements(
        patchwork::wrap_plots(scores_plot_list,
                              ncol = 2) +
          patchwork::plot_layout(widths = 1) +
          patchwork::plot_annotation("Supervised scores",
                                     theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                                                               hjust = 0,
                                                                                               size = 20)))
      )


    # Unsupervised clustering ###############################################

    cluster_plot_list <- list()

    # Set maximum number of clusters (cannot be more than samples - 1)
    if (length(cluster_labels) <= max_nc) {
      max_nc_adj <- length(cluster_labels) - 1
    } else {
      max_nc_adj <- max_nc
    }


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
        NbClust::NbClust(data = matrix,
                         min.nc = 1,
                         max.nc = max_nc_adj,
                         method = "ward.D2",
                         index = x)
      }, error=function(e) NULL)
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


    # Plot PCA
    if (nclust == 1) {
      cluster_plot_list[["pca_hclust_clustered"]] <- plot_pca(matrix,
                                                              label = "var",
                                                              pointsize = 1.5) +
        ggplot2::ggtitle("Hierarchical\nMethod: ward.D2")
    } else {
      # For each method (column) get consensus cluster assignment
      hclust_clusters <- apply(cluster_labels_df, 1, function(x) names(which.max(table(x)))) %>%
        as.factor()

      cluster_plot_list[["pca_hclust_clustered"]] <- plot_pca(matrix,
                                                              label = "var",
                                                              pointsize = 1.5,
                                                              color_cluster_by = hclust_clusters,
                                                              add_ellipses = TRUE) +
        ggplot2::ggtitle("Hierarchical\nMethod: ward.D2") +
        ggplot2::theme(legend.position="none")
    }


    ## PAM clustering ###############################################

    if (pam_nclusters == "auto") {
      # Find number of clusters with PAM and silhouette method

      p <-
        factoextra::fviz_nbclust(x = matrix,
                                 FUNcluster = cluster::pam,
                                 method = "silhouette",
                                 k.max = max_nc_adj,
                                 print.summary = TRUE) +
        theme_minimal() +
        ggtitle("Clusters suggested by\nK-Medoids (Partitioning Around Medoids) and silhouette method")
      # results[["plots"]][["number_of_clusters"]] <- p
      nclust <- as.numeric(p$data$clusters[which.max(p$data$y)])
    } else {
      if (length(cluster_labels) < pam_nclusters) {
        nclust <- length(cluster_labels)
      } else {
        nclust <- pam_nclusters
      }
    }

    pam_clusters <- cluster::pam(matrix, nclust, cluster.only=TRUE, nstart = 30)

    # Plot PCA
    if (nclust == 1) {
      cluster_plot_list[["pca_pam_clustered"]] <- plot_pca(matrix,
                                                           label = "var",
                                                           pointsize = 1.5) +
        ggplot2::ggtitle(paste0("PAM\nk: ", nclust))
    } else {
      cluster_plot_list[["pca_pam_clustered"]] <- plot_pca(matrix,
                                                           label = "var",
                                                           pointsize = 1.5,
                                                           color_cluster_by = pam_clusters,
                                                           add_ellipses = TRUE) +
        ggplot2::ggtitle(paste0("PAM\nk: ", nclust)) +
        ggplot2::theme(legend.position="none")
    }


    ## Leiden clustering ###############################################

    adj_graph <- igraph::as.undirected(cccd::nng(matrix, k = 3))
    r <- quantile(igraph::strength(adj_graph))[2] / (igraph::gorder(adj_graph) - 1) / 4
    leiden_clusters <- igraph::cluster_leiden(adj_graph, resolution_parameter = r)$membership

    # Plot PCA
    if (length(unique(leiden_clusters)) == 1) {
      cluster_plot_list[["pca_leiden_clustered"]] <- plot_pca(matrix,
                                                              label = "var",
                                                              pointsize = 1.5) +
        ggplot2::ggtitle(paste0("Leiden\nResolution: ", round(r, 4)))
    } else {
      cluster_plot_list[["pca_leiden_clustered"]] <- plot_pca(matrix,
                                                              label = "var",
                                                              pointsize = 1.5,
                                                              color_cluster_by = leiden_clusters,
                                                              add_ellipses = TRUE) +
        ggplot2::ggtitle(paste0("Leiden\nResolution: ", round(r, 4))) +
        ggplot2::theme(legend.position="none")
    }


    ## HDBSCAN clustering ###############################################

    cl <- dbscan::hdbscan(matrix, minPts = hdbscan_min_pts)
    hdbscan_clusters <- cl$cluster+1

    # Plot PCA
    if (length(unique(hdbscan_clusters)) == 1) {
      cluster_plot_list[["pca_hdbscan_clustered"]] <- plot_pca(matrix,
                                                               label = "var",
                                                               pointsize = 1.5) +
        ggplot2::ggtitle(paste0("HDBSCAN\nminPts: ", hdbscan_min_pts))
    } else {
      cluster_plot_list[["pca_hdbscan_clustered"]] <- plot_pca(matrix,
                                                               label = "var",
                                                               pointsize = 1.5,
                                                               color_cluster_by = hdbscan_clusters,
                                                               add_ellipses = TRUE) +
        ggplot2::ggtitle(paste0("HDBSCAN\nminPts: ", hdbscan_min_pts)) +
        ggplot2::theme(legend.position="none")
    }


    ## Combine clustering plots ###############################################

    results[["plots"]][["clustering_combined"]] <-
      patchwork::wrap_elements(
        patchwork::wrap_plots(cluster_plot_list,
                              ncol = 2) +
          patchwork::plot_layout(widths = 1) +
          patchwork::plot_annotation("Unsupervised clustering",
                                     theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                                                               hjust = 0,
                                                                                               size = 20)))
      )


    # Combine plots ###############################################

    results[["plots"]][["summary_plot"]] <- patchwork::wrap_plots(results[["plots"]],
                                                                  ncol = 2) +
      patchwork::plot_layout(widths = 1) +
      patchwork::plot_annotation(title,
                                 theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                                                           hjust = 0,
                                                                                           size = 20)))

    results[["plots"]][["scores"]] <- scores_plot_list
    results[["plots"]][["clustering"]] <- cluster_plot_list

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
                     label = "all",
                     pointsize = 3,
                     invisible = c("var", "quali"),
                     select_var = NULL,
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
                             scale. = FALSE)
    suppressWarnings(
      suppressMessages(
        p <- factoextra::fviz_pca(res.pca,
                                  habillage = color_cluster_by,
                                  addEllipses = add_ellipses,
                                  label = label,
                                  pointsize = pointsize,
                                  invisible = invisible,
                                  select.var = select_var,
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
                            "\nAvg sil width = ", round(m, 3),
                            "\n95% CI: [", round(ci[1], 3), ", ", round(ci[2], 3), "]",
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


# Temporary functions ##############################################################################

scores_composite_pca <- function(scores,
                                 scoot_summary
) {

  cluster_by <- names(scores)[!names(scores) %in% "params"]

  results <- list()

  for (c in cluster_by) {

    # L1 ###############################################

    mat <- scores[[c]][["composition"]][["layer_1"]][["distance_matrix"]]
    results[[c]][["matrices"]][["L1"]] <- mat

    res.pca <- prcomp(mat)

    clust_labels <- scoot_summary@metadata[match(row.names(mat), scoot_summary@metadata[["scoot_sample"]]), ][[c]]
    sil <- cluster::silhouette(as.numeric(clust_labels), mat)
    mean_sil <- sil[, 3] %>% mean() %>% round(3)

    results[[c]][["plots"]][["L1"]] <- factoextra::fviz_pca(res.pca,
                                                            habillage = clust_labels,
                                                            label = "var",
                                                            pointsize = 3,
                                                            invisible = c("var", "quali"),
                                                            geom = "point") +
      ggtitle(paste("L1 - PCA based on low-res cell type composition\n Silhouette score:", mean_sil)) +
      coord_equal()


    # L2 ###############################################

    mats <- list()
    for (i in names(scores[[c]][["composition"]][["layer_2"]])) {
      mats[[i]] <- scores[[c]][["composition"]][["layer_2"]][[i]][["distance_matrix"]]
    }

    avg_matrix <- get_avg_matrix(mats)

    results[[c]][["matrices"]][["L2"]] <- avg_matrix

    res.pca <- prcomp(avg_matrix)

    clust_labels <- scoot_summary@metadata[match(row.names(avg_matrix), scoot_summary@metadata[["scoot_sample"]]), ][[c]]
    sil <- cluster::silhouette(as.numeric(clust_labels), avg_matrix)
    mean_sil <- sil[, 3] %>% mean() %>% round(3)

    results[[c]][["plots"]][["L2"]] <- factoextra::fviz_pca(res.pca,
                                                            habillage = clust_labels,
                                                            pointsize = 3,
                                                            invisible = c("var", "quali"),
                                                            geom = "point") +
      ggtitle(paste("L2 - PCA based on hi-res cell type composition\nSilhouette score:", mean_sil)) +
      coord_equal()


    # L3 ###############################################

    mats <- list()
    for (i in names(scores[[c]][["signatures"]][["layer_1"]])) {
      mats[[paste0("layer_1_", i)]] <- scores[[c]][["signatures"]][["layer_1"]][[i]][["distance_matrix"]]
    }
    for (i in names(scores[[c]][["signatures"]][["layer_2"]])) {
      mats[[paste0("layer_2_", i)]] <- scores[[c]][["signatures"]][["layer_1"]][[i]][["distance_matrix"]]
    }

    avg_matrix <- get_avg_matrix(mats)

    results[[c]][["matrices"]][["L3"]] <- avg_matrix

    res.pca <- prcomp(avg_matrix)

    clust_labels <- scoot_summary@metadata[match(row.names(avg_matrix), scoot_summary@metadata[["scoot_sample"]]), ][[c]]
    sil <- cluster::silhouette(as.numeric(clust_labels), avg_matrix)
    mean_sil <- sil[, 3] %>% mean() %>% round(3)

    results[[c]][["plots"]][["L3"]] <- factoextra::fviz_pca(res.pca,
                                                            habillage = clust_labels,
                                                            pointsize = 3,
                                                            invisible = c("var", "quali"),
                                                            geom = "point") +
      ggtitle(paste("L3 - PCA based on signature expression per cell type\nSilhouette score:", mean_sil)) +
      coord_equal()


    # L1_L2_L3_combined ###############################################

    avg_matrix <- get_avg_matrix(results[[c]][["matrices"]])

    results[[c]][["matrices"]][["L1_L2_L3_combined"]] <- avg_matrix

    res.pca <- prcomp(avg_matrix)

    clust_labels <- scoot_summary@metadata[match(row.names(avg_matrix), scoot_summary@metadata[["scoot_sample"]]), ][[c]]
    sil <- cluster::silhouette(as.numeric(clust_labels), avg_matrix)
    mean_sil <- sil[, 3] %>% mean() %>% round(3)

    results[[c]][["plots"]][["L1_L2_L3_combined"]] <- factoextra::fviz_pca(res.pca,
                                                                           habillage = clust_labels,
                                                                           pointsize = 3,
                                                                           invisible = c("var", "quali"),
                                                                           geom = "point") +
      ggtitle(paste("L1_L2_L3_combined - PCA \nSilhouette score:", mean_sil)) +
      coord_equal()


    # B1 ###############################################

    mat <- scores[[c]][["pseudobulk"]][["layer_1"]][["all"]][["distance_matrix"]]
    results[[c]][["matrices"]][["B1"]] <- mat

    res.pca <- prcomp(mat)

    clust_labels <- scoot_summary@metadata[match(row.names(mat), scoot_summary@metadata[["scoot_sample"]]), ][[c]]
    sil <- cluster::silhouette(as.numeric(clust_labels), mat)
    mean_sil <- sil[, 3] %>% mean() %>% round(3)

    results[[c]][["plots"]][["B1"]] <- factoextra::fviz_pca(res.pca,
                                                            habillage = clust_labels,
                                                            label = "var",
                                                            pointsize = 3,
                                                            invisible = c("var", "quali"),
                                                            geom = "point") +
      ggtitle(paste("B1 - pseudobulk PCA\n Silhouette score:", mean_sil)) +
      coord_equal()


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
