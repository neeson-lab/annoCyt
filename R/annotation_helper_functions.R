
.getSingscoreValues = function(rankdata, gs.name, gs.up, gs.down = NULL, centerScore = TRUE, knownDirection = TRUE, dispersionFun = mad) {
  set.up = GSEABase::GeneSet(gs.up, geneIdType=GSEABase::SymbolIdentifier(), setName=gs.name)
  res = NULL
  if (!is.null(gs.down)) {
    set.down =  GSEABase::GeneSet(gs.down, geneIdType=GSEABase::SymbolIdentifier(), setName=gs.name)
    res = singscore::simpleScore(rankData=rankdata, upSet = set.up, downSet = set.down, knownDirection = knownDirection, centerScore = centerScore, dispersionFun = dispersionFun)
    res = res |> tibble::rownames_to_column("cluster_id") |> dplyr::select("cluster_id", "TotalScore", "TotalDispersion")
    cnames = colnames(res)
    gs.name2 = make.names(gs.name)
    cnames[2] = paste(gs.name2, cnames[2], sep="_")
    cnames[3] = paste(gs.name2, cnames[3], sep="_")
    colnames(res) = cnames

  } else {
    res = singscore::simpleScore(rankData=rankdata, upSet = set.up, knownDirection = knownDirection, centerScore = centerScore, dispersionFun = dispersionFun)
    res = res |> tibble::rownames_to_column("cluster_id") |> dplyr::select("cluster_id", "TotalScore", "TotalDispersion")
    cnames = colnames(res)
    gs.name2 = make.names(gs.name)
    cnames[2] = paste(gs.name2, cnames[2], sep="_")
    cnames[3] = paste(gs.name2, cnames[3], sep="_")
    colnames(res) = cnames

  }

  return(res[,1:2])
}


#' Generate 'singscore' values for each cluster and label
#'
#' Derives a score for every cluster and label based on marker expression. This score represents how likely a cluster is to be a particular cell type based on the markers specified in the lineage table for each cell type.
#' Scores will be generated for any entry in the lineage table with all positive markers OR a mix of positive and negative markers, but NOT for entries with only negative markers.
#'
#' @param lineagedf Dataframe storing user-defined cell type lineage. Must contain columns "label" (name of the cell type), "parent" (the parent cell type) and "value" (the markers that define a cell type).
#' @param exprdf  Dataframe storing expression values per marker per cluster
#'
#' @return A dataframe of scores for each cluster and cell type.
#' @export
#'
#' @examples
#' generateScores(eg_lineageTable, eg_expr)
generateScores = function(lineagedf, exprdf) {

  lineagedf = .updateLineageTable(lineagedf)
  lineagedf = lineagedf |> dplyr::ungroup() |> dplyr::distinct(label, value, .keep_all = TRUE)

  exprdf = tibble::column_to_rownames(exprdf, var = "cluster_id")
  rankdata = singscore::rankGenes(exprdf)

  ranking.scores = data.frame()

  for (i in 1:nrow(lineagedf)) {

    scoremarkers = NA
    negmarkers = NA

    scorename = lineagedf[i, "label"][[1]]
    scoremarkers = lineagedf[i, "value"][[1]][[1]]

    if (any(grepl("neg", scoremarkers))) {
      containsneg = grepl("neg", scoremarkers)

      posmarkers = scoremarkers[!containsneg]
      negmarkers = scoremarkers[containsneg]
      negmarkers = gsub("neg", "", negmarkers)

      if (length(scoremarkers) == length(negmarkers)){
        next
      }

    } else {
      posmarkers = scoremarkers
      negmarkers = NULL
    }

    # if a cell type has all negative markers, don't generate a singscore
    if (is.null(negmarkers)) {
      ascore = .getSingscoreValues(t(rankdata), gs.name = scorename, gs.up = posmarkers)
    } else {
      ascore = .getSingscoreValues(
        t(rankdata),
        gs.name = scorename,
        gs.up = posmarkers,
        gs.down = negmarkers
      )
    }

    if(nrow(ranking.scores) == 0){
      ranking.scores = ascore
    } else {
      ranking.scores = dplyr::full_join(ranking.scores, ascore, by = "cluster_id")
    }
  }

  return(ranking.scores)
}


.updateLineageTable = function(lineagedf) {

  newlt = data.frame(label = lineagedf$label, parent = lineagedf$parent)

  ms = strsplit(lineagedf$value, ",")

  newlt = lineagedf |> dplyr::mutate(value = stringr::str_split(lineagedf$value, ","))

  return(newlt)
}


#' Get cell counts per cluster per sample
#'
#' @param sceobj SingleCellExperiment object containing expression data and metadata
#' @param cluster_code Get cluster counts per this level of clustering. Must be available in names(cluster_codes(sce))
#' @param anno_name Name of annotation specified in annoToSCE(). Must be available in names(cluster_codes(sce))
#' @param thr_clusters_df Optional dataframe from output of compareToThreshold()
#'
#' @return Dataframe of cell counts per cluster and sample
#' @export
#'
#' @examples
#' getCountsPerCluster(eg_annotatedSCE, anno_name = "ct2")
#' getCountsPerCluster(eg_annotatedSCE, anno_name = "ct2", thr_clusters_df = eg_indepAnno)
getCountsPerCluster = function(sceobj, cluster_code = "som100", anno_name = "annoLabel", thr_clusters_df = NULL){

  if (!inherits(sceobj, "SingleCellExperiment")) {
    stop("sceobj must be a SingleCellExperiment object.")
  }

  if (!(cluster_code %in% names(CATALYST::cluster_codes(sceobj)))){
    stop("cluster_code not valid; must be present in names(cluster_codes(x))")
  }

  anno.per.cell = as.data.frame(SummarizedExperiment::colData(sceobj))

  celltype.meta = data.frame(cluster_id = S4Vectors::metadata(sceobj)$cluster_codes[[cluster_code]],
                             annotation = S4Vectors::metadata(sceobj)$cluster_codes[[anno_name]])

  if (cluster_code != "som100"){
  data.clustercodes = S4Vectors::metadata(sceobj)$cluster_codes |> dplyr::select("som100", !!cluster_code)
  anno.per.cell = dplyr::left_join(anno.per.cell, data.clustercodes, by = c( "cluster_id" = "som100"))
  anno.per.cell = dplyr::left_join(anno.per.cell, celltype.meta, by = stats::setNames("cluster_id", cluster_code))
  anno.per.cell = anno.per.cell |> dplyr::select(-cluster_id) |> dplyr::rename(cluster_id = !!cluster_code)
  } else {
    anno.per.cell = dplyr::left_join(anno.per.cell, celltype.meta, by = "cluster_id")

  }


  cluster.counts.df = anno.per.cell |>
    dplyr::group_by(sample_id, annotation, cluster_id) |> dplyr::summarise(n_cluster=dplyr::n())

  ncols = max(stringr::str_count(cluster.counts.df$annotation,"/")) + 1

  cluster.counts.df = tidyr::separate_wider_delim(
    cluster.counts.df,
    cols = annotation,
    delim = "/",
    names = paste0("label", 1:ncols),
    too_few = "align_start",
    cols_remove = FALSE
  )

  if (is.data.frame(thr_clusters_df)){
    with.markers = TRUE
    cluster.counts.df =  dplyr::left_join(cluster.counts.df, thr_clusters_df, by = "cluster_id" )
  }

  cluster.counts.df = dplyr::relocate(cluster.counts.df, sample_id, cluster_id, annotation)

  return(cluster.counts.df)
}



#' Get counts per label of annotation hierarchy per sample
#'
#' Combines the count of clusters with the same annotation label to generate a counts dataframe per annotation label per sample.
#'
#' @param cluster.counts.df Dataframe from the output of getCountsPerCluster()
#' @param add.proportions Logical; if TRUE, will add the proportion of each label level using the previous label level as the denominator (p_labelX = n_labelX / n_label(X-1))
#'
#' @return Dataframe of counts per label per sample
#' @export
#'
#' @examples
#' getCountsPerLabel(eg_clusterCounts)
#' 
#' @importFrom rlang :=
getCountsPerLabel = function(cluster.counts.df, add.proportions = TRUE) {
  # Identify dynamically which columns start with "label"
  label_cols = grep("^label", names(cluster.counts.df), value = TRUE)

  # Start with the base dataframe
  label.counts = cluster.counts.df

  # STEP 1: Create n_labelX columns dynamically
  for (i in seq_along(label_cols)) {
    level_cols = c("sample_id", label_cols[1:i])
    new_col_name = paste0("n_", label_cols[i])

    label.counts = label.counts |>
      dplyr::group_by(dplyr::across(dplyr::all_of(level_cols))) |>
      dplyr::mutate(!!new_col_name := dplyr::case_when(is.na(.data[[label_cols[i]]]) ~ NA_real_,
                                                       TRUE ~ sum(n_cluster))
      ) |>
      dplyr::ungroup()
  }

  # STEP 2: Create total count per sample for p_label1 denominator
  label.counts = label.counts |>
    dplyr::group_by(sample_id) |>
    dplyr::mutate(total_sample = sum(n_cluster)) |>
    dplyr::ungroup()

  # STEP 3: Create proportion columns dynamically (INCLUDING p_label1)
  if (add.proportions){
  for (i in seq_along(label_cols)) {
    current_n = paste0("n_", label_cols[i])
    new_prop = paste0("p_", label_cols[i])

    if (i == 1) {
      # p_label1 = n_label1 / total_sample
      label.counts = label.counts |>
        dplyr::mutate(!!new_prop := dplyr::case_when(
          is.na(.data[[current_n]]) ~ NA_real_,
          TRUE ~ .data[[current_n]] / total_sample
        ) )
    } else {
      # p_labelX = n_labelX / n_label(X-1)
      prev_n = paste0("n_", label_cols[i-1])
      label.counts = label.counts |>
        dplyr::mutate(
          !!new_prop := dplyr::case_when(
            is.na(.data[[current_n]]) ~ NA_real_,
            is.na(.data[[prev_n]]) ~ NA_real_,
            TRUE ~ .data[[current_n]] / .data[[prev_n]]
          )
        )
    }
  }
    }


  # Final dataframe
  tokeep = colnames(label.counts)
  tokeep = tokeep[grepl("sample_id|annotation|label", tokeep)]

  label.counts = label.counts |>
    dplyr::select(dplyr::all_of(tokeep)) |>
    dplyr::ungroup() |>
    dplyr::group_by(sample_id) |>
    dplyr::distinct(annotation, .keep_all = TRUE)

  return(label.counts)
}

.getAllLabels = function(lineagedf){
  # from lineage table, get all possible full labels
  lineagedf = dplyr::filter(lineagedf, parent != "special")

  glt = lineagedf |> dplyr::select(label, parent) |> dplyr::relocate(label, .after = parent)
  colnames(glt) = c("from", "to")
  # mygraph = graph_from_data_frame(glt)

  # Create directed graph
  g = igraph::graph_from_data_frame(glt, directed = TRUE)

  # Identify root nodes (no incoming edges)
  roots = setdiff(igraph::V(g)$name, glt$to)

  # Identify leaf nodes (no outgoing edges)
  leaves = setdiff(igraph::V(g)$name, glt$from)

  # Function to get all simple paths from roots to leaves
  paths = unlist(
    lapply(roots, function(r)
      igraph::all_simple_paths(g, from = r, to = leaves)
    ),
    recursive = FALSE
  )

  # Convert paths to readable label strings
  path_strings = sapply(paths, function(p) paste(igraph::V(g)[p]$name, collapse = "/"))


  path_strings = gsub("root/", "", path_strings)

  return(path_strings)
}



#' Add missing cell types to cluster counts dataframe
#'
#' Some cell types may be present in the lineage dataframe that were not detected during the annotation step. Alternatively, some cell types may be not be detected in specific samples.
#' Users may be find it useful to know which cell types are not present at all or in specific samples.
#' This function adds the missing cell types per sample to the cluster counts dataframe with a count of 0.
#'
#' @param cluster_counts Dataframe of cluster counts per sample; output of getCountsPerCluster()
#' @param cluster_annotation_df Dataframe of annotation assignment per cluster; output of annotateClusters()
#' @param lineagedf Dataframe of lineage table of cell type hierarchy.
#'
#' @return Dataframe of cluster counts per sample with missing cell types added.
#' @export
#'
#' @examples
#' addMissing(eg_clusterCounts, eg_annoTbl, eg_lineageTable)
addMissing = function(cluster_counts, cluster_annotation_df, lineagedf) {

  if (!all(c("annotation", "cluster_id") %in% colnames(cluster_annotation_df))) {
    stop(
      "Columns missing from cluster_annotation_df. Please ensure there is a column called `annotation` and `cluster_id`."
    )
  }

  # all possible labels
  potential.labels <- .getAllLabels(lineagedf)

  # ensure all labels exist in assignment
  missing.labels <- setdiff(potential.labels, cluster_annotation_df$annotation)

  if (length(missing.labels) > 0) {
    nclusters <- max(as.numeric(cluster_annotation_df$cluster_id))

    additional.assignment <- tidyr::tibble(
      cluster_id = as.character(seq(nclusters + 1, length.out = length(missing.labels))),
      annotation = missing.labels
    )

    cluster_annotation_df <- dplyr::bind_rows(cluster_annotation_df, additional.assignment)
  }

  ncols <- max(stringr::str_count(cluster_annotation_df$annotation, "/")) + 1

  cluster_annotation_df <- cluster_annotation_df |>
    dplyr::select(cluster_id, annotation) |>
    tidyr::separate_wider_delim(
      cols = annotation,
      delim = "/",
      names = paste0("label", seq_len(ncols)),
      too_few = "align_start",
      cols_remove = FALSE
    ) |>
    dplyr::relocate(annotation, .after = cluster_id)

  full_grid <- tidyr::expand_grid(
    sample_id = unique(cluster_counts$sample_id),
    cluster_id = unique(cluster_annotation_df$cluster_id)
  )

  cluster.counts.all <- full_grid |>
    dplyr::left_join(cluster_counts, by = c("sample_id", "cluster_id")) |>
    dplyr::mutate(n_cluster = dplyr::coalesce(n_cluster, 0)) |>
    dplyr::select(!tidyr::matches("label|annotation")) |>
    dplyr::left_join(cluster_annotation_df,  by = "cluster_id")

  return(cluster.counts.all)
}


#' Compare metric to threshold per cluster
#' 
#' 
#' Input a matrix and threshold values determine if column value passes the provided threshold.
#' Threshold is user-defined.
#' Generalised to take any matrix containing cluster-level values (i.e, not just expression data).
#'
#' @param scores_matrix Dataframe or matrix containing values of some metric per cluster
#' @param threshold_matrix Named vector or matrix containing threshold values per metric. Must be present in column names of scores_matrix
#' @param clusters Optional vector of clusters to run comparison on. Must be present in row names of scores_matrix
#'
#' @return Dataframe of presence/absence of metric per cluster.
#' @export
#'
#' @examples
#' compareToThreshold(eg_expr, c( "KI67" = 0.5))
compareToThreshold = function(scores_matrix, threshold_matrix, clusters = NULL) {
  scores_matrix = as.data.frame(scores_matrix)



  # Subset clusters
  if (!is.null(clusters)) {
    missing_clusters = setdiff(clusters, rownames(scores_matrix))
    if (length(missing_clusters) > 0) {
      stop(paste("clusters not found in expression data:", paste(missing_clusters, collapse = ", ")))
    }
    scores_matrix = scores_matrix[clusters, , drop = FALSE]
  }

  # Handle threshold as vector
  if (is.vector(threshold_matrix)) {
    matched_markers = intersect(names(threshold_matrix), colnames(scores_matrix))

    if (length(matched_markers) == 0) {
      stop("No matching markers between scores_matrix and threshold_matrix.")
    } else if (length(matched_markers) != length(threshold_matrix)) {
        unmatched_markers = setdiff(names(threshold_matrix), colnames(scores_matrix))
        message(paste("Warning:", unmatched_markers ,
                      "not found in scores_matrix; not included in output.\n"))
    }

    scores_matrix = scores_matrix[, matched_markers, drop = FALSE]

    threshold_matrix = as.data.frame(matrix(
      rep(threshold_matrix[matched_markers], each = nrow(scores_matrix)),
      nrow = nrow(scores_matrix),
      ncol = length(matched_markers)
    ))

    rownames(threshold_matrix) = rownames(scores_matrix)
    colnames(threshold_matrix) = matched_markers

  } else {
    threshold_matrix = as.data.frame(threshold_matrix)

    matched_clusters = intersect(rownames(threshold_matrix), rownames(scores_matrix))
    matched_markers = intersect(colnames(threshold_matrix), colnames(scores_matrix))

    scores_matrix = scores_matrix[matched_clusters, matched_markers, drop = FALSE]
    threshold_matrix = threshold_matrix[matched_clusters, matched_markers, drop = FALSE]

    if (!all(dim(scores_matrix) == dim(threshold_matrix))) {
      stop("Expression and threshold data frames must match in dimensions after filtering.")
    }
  }

  # Build result matrix
  result = matrix(nrow = nrow(scores_matrix), ncol = ncol(scores_matrix))
  rownames(result) = rownames(scores_matrix)
  colnames(result) = colnames(scores_matrix)

  # Loop over markers (columns now)
  for (j in seq_len(ncol(scores_matrix))) {
    marker_name = colnames(scores_matrix)[j]
    result[, j] = ifelse(scores_matrix[, j] > threshold_matrix[, j],
                         paste0(marker_name, "+"),
                         paste0(marker_name, "-"))
  }

  result_df = as.data.frame(result) |>
    tibble::rownames_to_column(var = "cluster_id")

  return(result_df)
}



#' Add cluster annotation assignment to SCE object
#'
#'
#' @param sceobj A SingleCellExperiment object
#' @param cluster_annotation_df Dataframe of annotation labels per cluster; output of annotateClusters(). Must contain the columns "cluster_id" and "annotation"
#' @param cluster_code String of the clustering code used for annotation (usually metaX or som100); must be available in names(cluster_codes(sceobj))
#' @param anno_name String of the name of annotation to be used in SCE; will be present in names(cluster_codes(sceobj))
#' @param overwrite Logical; whether to update annotation under anno_name if anno_name already present
#'
#' @return An updated SingleCellExperiment object with annotation added to cluster_codes() and colData
#' @export
#'
#' @examples
#' \dontrun{
#' sce = CATALYST::cluster(eg_sce_raw, features = "type", 
#'     xdim=10, ydim=10, maxK=80, verbose=FALSE, seed=1234)
#' sce = annoToSCE(sce, eg_annoTbl, anno_name = "ct2")
#' }
#' 
#' 
#' @importFrom rlang :=
annoToSCE = function(sceobj, cluster_annotation_df, cluster_code = "som100", anno_name = "annoLabel", overwrite = FALSE){

  if (!inherits(sceobj, "SingleCellExperiment")) {
    stop("sceobj must be a SingleCellExperiment object.")
  }

   if (!(cluster_code %in% names(CATALYST::cluster_codes(sceobj)))){
    stop("cluster_code not valid; must be present in names(cluster_codes(x))")
   }

  # add to cluster codes
  data.clustercodes = S4Vectors::metadata(sceobj)$cluster_codes
  cluster_annotation_df_min = cluster_annotation_df[ , c("cluster_id", "annotation")]

  merging.tbl = cluster_annotation_df_min |>
    dplyr::filter(!is.na(cluster_id)) |>
    dplyr::mutate(!!cluster_code:=as.character(cluster_id), !!anno_name:=as.character(annotation)) |>
    dplyr::select(!!cluster_code, !!anno_name)

  sceobj = CATALYST::mergeClusters(sceobj, k = cluster_code, table = merging.tbl, id = anno_name, overwrite = overwrite)


# add to coldata
  data.clustercodes =  S4Vectors::metadata(sceobj)$cluster_codes

    match_clusters = data.clustercodes[ , c("som100", anno_name)]
    colnames(match_clusters)[2] = "annotation"
    SummarizedExperiment::colData(sceobj)[[anno_name]] = match_clusters$annotation[match(SummarizedExperiment::colData(sceobj)$cluster_id, match_clusters$som100)]

  return(sceobj)
}
