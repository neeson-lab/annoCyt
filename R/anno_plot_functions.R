#' Create diagram showing flow of lineage table
#'
#' Generates a basic visualisation of your lineage table in the form of a tree diagram. This will give you an idea of the annotation logic.
#'
#' @param lineagedf  Dataframe storing user-defined cell type lineage. Must contain columns "label" (name of the cell type), "parent" (the parent cell type) and "value" (the markers that define a cell type).
#'
#' @return ggraph tree diagram.
#' @export
#'
#' @examples
#' annoDiagram(eg_lineageTable)
annoDiagram = function(lineagedf) {

  lineagedf = .updateLineageTable(lineagedf)

  lineagedf = dplyr::filter(lineagedf, parent != "special")

  glt =  dplyr::select(lineagedf, label, parent) |> dplyr::relocate(label, .after = parent)
  colnames(glt) = c("from", "to")
  mygraph = igraph::graph_from_data_frame(glt)

  diagram = ggraph::ggraph(mygraph, layout = 'tree', circular = FALSE) +
    ggraph::geom_edge_link(ggplot2::aes(start_cap = ggraph::label_rect(node1.name), end_cap = ggraph::label_rect(node2.name)), arrow = grid::arrow(length = grid::unit(4, 'mm'))) +
    ggraph::geom_node_label(ggplot2::aes( label=name)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")

  return(diagram)
}


#' Plot expression per cluster with annotation labels
#'
#' Heatmap of the median scaled expression per cluster labelled with the cluster ID and annotation assignment. Uses CATALYST::plotExprHeatmap() to generate heatmap.
#'
#' @param sceobj A SingleCellExperiment object
#' @param lineage_df DataFrame of cell type hierarchy.
#' @param cluster_code A string of the clustering code used for annotation (usually metaX or som100); must be available in names(cluster_codes(sceobj))
#' @param anno_name String of the name of annotation to be used in SCE; will be present in names(cluster_codes(sceobj))
#' @param output_matrix Logical; if TRUE, outputs the matrix used to generate heatmap instead of the heatmap itself.
#' @param marker_order Optional string vector of order in which to plot markers in the heatmap. Must contain all markers present in expression data of SCE.
#' @param ... Optional additions to ComplexHeatmap heatmap
#'
#' @return Heatmap of median scaled expression per cluster; if output_matrix == TRUE, outputs matrix used to generate heatmap instead
#' @export
#'
#' @examples
#' plotClusterHeatmap(eg_annotatedSCE, eg_lineageTable,  anno_name = "ct2 ",output_matrix = FALSE)
plotClusterHeatmap = function(sceobj, lineage_df, cluster_code = "som100", anno_name = "celltype1",
                              output_matrix = F, marker_order = NULL,  ...){

  S4Vectors::metadata(sceobj)$cluster_codes$numberedAnno = as.factor(paste(
    S4Vectors::metadata(sceobj)$cluster_codes[[cluster_code]],
    S4Vectors::metadata(sceobj)$cluster_codes[[anno_name]],
    sep = " "
  ))

  x0 = CATALYST::plotExprHeatmap(sceobj,  by="cluster_id", k = "numberedAnno", row_anno=F, perc = T)
  x = x0@matrix
  # order labels by DFS
  labels = unique(S4Vectors::metadata(sceobj)$cluster_codes[[anno_name]])

  lineage_df = .updateLineageTable(lineage_df)
  lineage_df = dplyr::filter(lineage_df, parent != "special")

  glt =  dplyr::select(lineage_df, label, parent) |>
    dplyr::relocate(label, .after = parent)
  colnames(glt) = c("from", "to")
  mygraph = igraph::graph_from_data_frame(glt)
  dfs_order <- unlist(lapply("root", .dfs_paths, graph = mygraph))
  dfs_order = gsub("root/", "", dfs_order)
  # Get leaf nodes (no outgoing edges)
  leaf_names <- igraph::V(mygraph)$name[igraph::degree(mygraph, mode = "out") == 0]
  # Identify which dfs_order paths are non-leaves
  is_leaf_path <- sapply(dfs_order, function(path) {
    utils::tail(strsplit(path, "/")[[1]], 1) %in% leaf_names
  })

  # Append identifier to non-leaf paths in dfs_order
  dfs_order_display <- ifelse(is_leaf_path, dfs_order, paste0(dfs_order, "/unknown"))
  rowlabels = rownames(x)[order(.match_partial(rownames(x), dfs_order_display))]

  roworder = order(factor(rownames(x), levels = rowlabels))

  if (is.vector(marker_order)){
    if ((missing <- setdiff(marker_order, colnames(x))) |> length() > 0) {
      stop("Missing values: ", paste(missing, collapse = ", "))
    }
    cbels = colnames(x)
    colorder = order(factor(cbels, levels = marker_order))
    y = x[roworder, colorder, drop=FALSE]
  } else {
    y = x[roworder, , drop=FALSE]
  }


  yhm = ComplexHeatmap::Heatmap(
    matrix = y,
    col = circlize::colorRamp2(seq(min(y), max(y), l = n <- 100), grDevices::colorRampPalette( rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(n)),
    rect_gp = grid::gpar(col = "white"),
    cluster_rows = F ,
    cluster_columns = F,
    show_row_dend = F,
    show_column_dend = F,
    ...
  )

  if (output_matrix){
    return(x@matrix)
  } else {
  return(yhm)
    }
}

#' Plot values of a metric per cluster
#'
#' Generates a basic dot plot of metric values per cluster. Will generate a plot for every column in cluster_values_df.
#' Will specify the the annotation labels (as well as cluster ID) in cluster_annotation_df if provided.
#'
#' Designed to assist with annotation rule refinement.
#' This function is often useful when determining expression or singscore thresholds for the rules e.greaterthan/e.lessthan/s.greaterthan/s.lessthan .
#'
#' @param cluster_values_df A dataframe containing columns with values per cluster (usually either expression data or singscore values)
#' @param cluster_annotation_df Optional; a dataframe containing annotation label per cluster; output of annotateClusters().
#' @param as.image Logical; if TRUE, will export each plot as a .png into output_dir
#' @param output_dir String containing directory path of where to store images; only required when as.image is TRUE
#'
#' @return List of plots per metric in cluster_values_df. If as.image is TRUE, check output_dir for plots.
#' @export
#'
#' @examples
#' plotClusterValues(cluster_values_df = eg_expr, 
#'     cluster_annotation_df = eg_annoTbl, as.image = FALSE)
plotClusterValues = function(cluster_values_df,
                             cluster_annotation_df = NULL,
                             as.image = FALSE,
                             output_dir = NULL) {

  if (!"cluster_id" %in% colnames(cluster_values_df)) {
    stop("cluster_values_df must contain 'cluster_id'")
  }


  total_cols <- colnames(cluster_values_df)[!(colnames(cluster_values_df) == "cluster_id")]

  tmp <- dplyr::mutate(cluster_values_df, cluster_id = as.numeric(cluster_id))
  tmp <- dplyr::arrange(tmp, cluster_id)
  tmp <- dplyr::mutate(tmp, cluster_id = as.character(cluster_id))

  if (!is.null(cluster_annotation_df)){
    if (!"cluster_id" %in% colnames(cluster_annotation_df)) {
      stop("cluster_annotation_df must contain 'cluster_id'")
    }
    ann <- dplyr::mutate(cluster_annotation_df, cluster_id = as.character(cluster_id))
    tmp <- dplyr::left_join(tmp, ann, by = "cluster_id")
    tmp <- dplyr::mutate(tmp, numlabel = paste(cluster_id, annotation))
  } else {
    tmp <- dplyr::mutate(tmp, numlabel = as.character(cluster_id))
  }

  scores.tmp <- dplyr::select(tmp, numlabel, tidyr::all_of(total_cols))

  plots = lapply(total_cols, function(col) {
    ggplot2::ggplot(scores.tmp, ggplot2::aes(
      x = forcats::fct_reorder(as.factor(numlabel), -.data[[col]]),
      y = .data[[col]]
    )) +
      ggplot2::geom_point() +
      ggplot2::labs(x = "cluster", y = col, title = col) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "bottom",
        strip.background = ggplot2::element_rect(fill = "white"),
        strip.text.x = ggplot2::element_text(size = 30),
        axis.text.x = ggplot2::element_text(size = 10, angle=90, hjust=0.95, vjust=0.2),
        axis.text.y = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_blank()
      )
  })

  names(plots) = total_cols

  if (as.image) {
    if (is.null(output_dir)) {
      output_dir = file.path(getwd(), "plots")
      message("No output directory specified; saving plots to ", output_dir)
    } else {
      message("Saving plots to: ", output_dir)
    }

    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    lapply(names(plots), function(nm) {
      ggplot2::ggsave(
        filename = file.path(output_dir, paste0(nm, ".png")),
        plot = plots[[nm]],
        width = 20,
        height = 10
      )
    })
  } else {
    return(plots)

  }
}


#' Plot singscore values of sibling cell types
#'
#' A plot to compare the singscore values of sibling cell types across clusters.
#' Designed to assist with annotation rule refinement. This plot is useful for plotting cell types specified in the max() rule and to finetune thresholds.
#'
#'
#'
#' @param ranking_scores_df A dataframe of scores for each cluster and cell type; often the output of generateScores()
#' @param compare.labels A string vector specifying the cell types to compare. Must be present in the label column of your lineage table
#' @param cluster_annotation_df A dataframe containing annotation label per cluster; output of annotateClusters().
#' @param all.clusters Optional; vector containing clusters to show in plot.
#'
#' @return ggplot2 plot
#' @export
#'
#' @examples
#' plotSiblings(eg_singscores, compare.labels = c("Tcell", "NK"),
#'     cluster_annotation_df = eg_annoTbl)
plotSiblings = function(ranking_scores_df,
                        compare.labels,
                        cluster_annotation_df,
                        all.clusters = NULL) {

  score_cols <- paste0(compare.labels, "_TotalScore")
  missing_cols <- setdiff(score_cols, colnames(ranking_scores_df))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The following compare.labels are missing in ranking_scores_df: ",
        paste(compare.labels[score_cols %in% missing_cols], collapse = ", ")
      )
    )
  }

  scores <- dplyr::arrange(ranking_scores_df, as.numeric(cluster_id))
  ann <- dplyr::mutate(cluster_annotation_df,
                       cluster_id = as.character(cluster_id))
  scores <- dplyr::left_join(scores, ann, by = "cluster_id")

  # choose clusters
  if (is.null(all.clusters)) {
    keep_clusters <- ann$cluster_id[
      vapply(
        strsplit(ann$annotation, "/"),
        function(x) any(x %in% compare.labels),
        logical(1)
      )
    ]
  } else {
    keep_clusters <- ann$cluster_id
  }

  ss <- scores[, c("cluster_id", "annotation", score_cols), drop = FALSE]
  ss <- ss[ss$cluster_id %in% keep_clusters, , drop = FALSE]
  ss$ds <- paste(ss$cluster_id, ss$annotation)
  ss <- tidyr::pivot_longer(
    ss,
    cols = dplyr::all_of(score_cols),
    names_to = "diffstate",
    values_to = "score"
  )

  ss <- dplyr::arrange(ss, annotation)
  ss$ds <- factor(ss$ds, levels = unique(ss$ds))


  plot1 = ggplot2::ggplot(ss) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
    ggplot2::geom_point(ggplot2::aes(x = ds, y = score,  group = diffstate, colour = diffstate),  size = 3) +
    ggplot2::geom_line((ggplot2::aes(x = ds, y = score,  group = diffstate, colour = diffstate)),  lwd=1, alpha = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "white"),
      strip.text.x = ggplot2::element_text(size = 30),
      axis.text.x = ggplot2::element_text(size = 10, angle=90, hjust=0.95, vjust=0.2),
      axis.text.y = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_blank()
    ) + ggplot2::labs(title = paste(compare.labels, collapse = " v "), y = "singscore value", x = "annotation")

  return(plot1)
}




#' Check validity of lineage table
#'
#' Checks the following in the lineage table:
#' \itemize{
#'  \item{All markers are present in exprdf (expression dataframe)}
#'  \item{Correct syntax}
#'  \item{All parent labels are present in "label" column}
#'  \item{Only one root parent}
#' }

#'
#' @param lineagedf Dataframe of cell type hierarchy; user specified
#' @param exprdf Dataframe of expression values per cluster; often the output of scaleClusteExpr()
#'
#' @return If invalid syntax in lineagedf, will specify what to fix.
#' @export
#'
#' @examples
#' checkLineageTable(eg_lineageTable, eg_expr)
checkLineageTable = function(lineagedf, exprdf) {

  valid_markers = colnames(exprdf)
  valid_markers = valid_markers[valid_markers != "cluster_id"]

  ## ---- 1. Marker checks ----
  for (i in seq_along(lineagedf$value)) {

    markers = strsplit(lineagedf$value[i], ",")[[1]]
    markers = trimws(markers)

    for (m in markers) {

      # Invalid neg usage
      if (grepl("neg", m) && !grepl("neg$", m)) {
        stop(
          sprintf(
            "Invalid 'neg' usage in marker '%s' (row %d)",
            m, i
          ),
          call. = FALSE
        )
      }

      # Remove trailing neg
      marker_clean = sub("neg$", "", m)

      if (!marker_clean %in% valid_markers) {
        stop(
          sprintf(
            "Marker '%s' not found in valid_markers (row %d)",
            marker_clean, i
          ),
          call. = FALSE
        )
      }
    }
  }

  ## ---- 2. Parent must exist in label (except root/special) ----
  parents_to_check = setdiff(unique(lineagedf$parent), c("root", "special"))

  missing_labels = setdiff(parents_to_check, lineagedf$label)

  if (length(missing_labels) > 0) {
    stop(
      sprintf(
        "Parent value(s) not found in label column: %s",
        paste(missing_labels, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  ## ---- 3. Exactly one root in parent ----
  n_root = sum(lineagedf$parent == "root", na.rm = TRUE)

  if (n_root != 1) {
    stop(
      sprintf(
        "Expected exactly one 'root' in parent column, found %d",
        n_root
      ),
      call. = FALSE
    )
  }

  ## ---- 4. No forward slashes in label ----
  invalid_idx = which(grepl("/", lineagedf$label))

  if (length(invalid_idx) > 0) {
    stop(
      sprintf(
        "'/' character not allowed in label; please replace with another symbol. Problem at row(s): %s",
        paste(invalid_idx, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  return(message("Valid syntax in lineage table. Good job!"))
}

#' Check validity of rules table
#'
#' @details Checks the following in the rules table:
#' \itemize{
#'  \item{all markers are present in exprdf (expression dataframe)}
#'  \item{Labels are present in lineage table (lineagedf)}
#'  \item{Name of each rule is allowed}
#'  \item{Threshold are numerical and only present for rules requiring it}
#' }
#'
#' @param rulesdf Dataframe of rules; created by user
#' @param lineagedf Dataframe of lineage table; created by user
#' @param exprdf Dataframe of expression values per cluster; likely the output of scaleClusterExpr()
#'
#' @return If invalid syntax in rulesdf, will specify what to fix.
#' @export
#'
#' @examples
#' checkRulesTable(eg_rulesTable, eg_lineageTable, eg_expr)
checkRulesTable = function(rulesdf, lineagedf, exprdf) {

  valid_markers = colnames(exprdf)
  valid_markers = valid_markers[valid_markers != "cluster_id"]


  ## ---- 1. label must exist in lineagedf$label ----
  missing_labels = setdiff(unique(rulesdf$label), lineagedf$label)

  if (length(missing_labels) > 0) {
    stop(
      sprintf(
        "rulesdf$label value(s) not found in lineagedf$label: %s",
        paste(missing_labels, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  ## ---- helpers ----
  is_es_rule = function(rule) {
    grepl("^(?:[A-Za-z0-9]+)?(e|s)\\.(greaterthan|lessthan)$", rule)
  }

  base_es_rule = function(rule) {
    grepl("^(e|s)\\.(greaterthan|lessthan)$", rule)
  }

  marker_es_rule = function(rule) {
    is_es_rule(rule) && !base_es_rule(rule)
  }

  ## ---- 2. rule validation ----
  for (i in seq_len(nrow(rulesdf))) {

    rule = rulesdf$rule[i]

    ## A. e./s. rules
    if (is_es_rule(rule)) {

      if (marker_es_rule(rule)) {
        marker = sub("(e|s)\\.(greaterthan|lessthan)$", "", rule)

        if (!marker %in% valid_markers) {
          stop(
            sprintf(
              "Invalid marker '%s' in rule '%s' (row %d)",
              marker, rule, i
            ),
            call. = FALSE
          )
        }
      }

      ## B. max(x,y,...) — must have ≥ 2 values
    } else if (grepl("^max\\([^()]+\\)$", rule)) {

      args = strsplit(
        gsub("^max\\(|\\)$", "", rule),
        ","
      )[[1]]

      args = trimws(args)

      if (length(args) < 2) {
        stop(
          sprintf(
            "max() must contain at least two values (row %d)",
            i
          ),
          call. = FALSE
        )
      }

      missing = setdiff(args, lineagedf$label)

      if (length(missing) > 0) {
        stop(
          sprintf(
            "Unknown label(s) in rule '%s': %s (row %d)",
            rule,
            paste(missing, collapse = ", "),
            i
          ),
          call. = FALSE
        )
      }

      ## C. not.x / isct.x
    } else if (grepl("^(not|isct)\\.", rule)) {

      arg = sub("^(not|isct)\\.", "", rule)

      if (!arg %in% lineagedf$label) {
        stop(
          sprintf(
            "Unknown label '%s' in rule '%s' (row %d)",
            arg, rule, i
          ),
          call. = FALSE
        )
      }

    } else {
      stop(
        sprintf(
          "Invalid rule syntax '%s' (row %d)",
          rule, i
        ),
        call. = FALSE
      )
    }

    ## ---- 3. threshold validation ----
    if (is_es_rule(rule)) {
      if (!is.numeric(rulesdf$threshold[i]) || is.na(rulesdf$threshold[i])) {
        stop(
          sprintf(
            "Rule '%s' requires numeric threshold (row %d)",
            rule, i
          ),
          call. = FALSE
        )
      }
    } else {
      if (!is.na(rulesdf$threshold[i])) {
        stop(
          sprintf(
            "Rule '%s' must NA in threshold column (row %d)",
            rule, i
          ),
          call. = FALSE
        )
      }
    }
  }

  return(message("Valid syntax in rules table. Well done!"))
}


.dfs_paths <- function(graph, node, current_path = "") {
  full_path <- if (current_path == "") node else paste(current_path, node, sep = "/")

  children <- igraph::neighbors(graph, node, mode = "out")
  child_names <- igraph::V(graph)$name[children]

  if (length(child_names) == 0) {
    return(full_path)
  } else {
    return(c(full_path, unlist(lapply(child_names, .dfs_paths, graph = graph, current_path = full_path))))
  }
}

.match_partial <- function(labels, order_vec) {
  sapply(labels, function(lbl) {
    lbl_clean <- sub("^\\d+\\.?\\d* ", "", lbl)
    which(order_vec == lbl_clean)[1]
  })
}
