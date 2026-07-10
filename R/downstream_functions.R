#' Get dimension reduction information from SCE object
#'
#'
#'
#' @param sceobj A SingleCellExperiment object
#' @param cluster_annotation_df Dataframe of annotation labels per cluster; output of annotateClusters(). Must contain the columns "cluster_id" and "annotation"
#' @param dr_method String specifying dimension reduction method
#' @param anno_name String of the name of annotation as in SCE; will be present in names(cluster_codes(sceobj))
#'
#' @return Dataframe of dimension reduction coordinates
#' @export
#'
#' @examples
#' \dontrun{
#' sce = CATALYST::runDR(sce, dr = "UMAP", cells = 5000)
#'  umapdata = getDRdata(sce, eg_annoTbl, dr_method = "UMAP", anno_name = "ct2")
#' }
#'
getDRdata = function(sceobj, cluster_annotation_df, dr_method = "UMAP", anno_name) {
  xy =SingleCellExperiment::reducedDim(sceobj, dr_method)
  colnames(xy) = c("x", "y")
  drdata =
    data.frame(SummarizedExperiment::colData(sceobj), xy, check.names = FALSE) |>
    dplyr::left_join(cluster_annotation_df, by = c("cluster_id",stats::setNames("annotation", anno_name))) |>
    dplyr::filter(!is.na(x) & !is.na(y)) |>
    dplyr::relocate(c("x", "y"), .before = "sample_id")
  return(drdata)
}



#' Downsample and combine batch SCEs and perform dimension reduction
#'
#' When cytometry analyses were run per batch with the same panel, use this function to combine batch-level expression data prior to dimension reduction visualisation.
#'
#' @details
#' Process is follows:
#' \itemize{
#'  \item{Expression data will first be scaled per batch.}
#'  \item{Downsample per sample}
#'  \item{Expression data combined across batches}
#'  \item{Dimension reduction performed}
#' }
#'
#' @param input.list List of SCE objects. Expects one SCE per batch.
#' @param scale.margin Value to scale values by. Default value is 2.
#' @param n_cells_per_sample Number of cells to keep per sample in each batch. Defaults to 5000 cells per sample.
#' @param dr_method String specifying dimension reduction method
#' @param cluster_code Get cluster counts per this level of clustering. Must be the same in all SCE objects
#' @param ... Parameters passed on to scater's runX() where X is a dimension reduction method
#'
#' @return A dataframe containing dimension reduction coordinates, as well as colData (sample metadata) and scaled expression values per cell.
#' @export
#'
#' @examples
#' \dontrun{ 
#' combineBatchDR(list(sce.batch1, sce.batch2))}

combineBatchDR = function(input.list,
                          scale.margin = 2,
                           n_cells_per_sample = 5000,
                           dr_method = "UMAP",
                           cluster_code = "som100",
                           ...) {

  # validate SCEs
  sce_list = lapply(seq_along(input.list), function(i) {
    x = input.list[[i]]
    if (is.character(x))
      x = readRDS(x)
    if (!inherits(x, "SingleCellExperiment"))
      stop(i, " is not an SCE object")
    x
  })

  # scale expression and attach metadata
  data_list = lapply(sce_list, function(sce) {
    col_df = as.data.frame(SummarizedExperiment::colData(sce))
    col_df$cluster_id = as.character(col_df$cluster_id)

    #cluster_codes = as.data.frame(S4Vectors::metadata(sce)$cluster_codes)
    #col_df = dplyr::inner_join(col_df, cluster_codes, by = c("cluster_id" = cluster_code))

    expr_df = as.data.frame(t(SummarizedExperiment::assay(sce, "exprs")))
    expr_scaled = as.data.frame(CATALYST:::.scale_exprs(expr_df, margin = scale.margin))

    cbind(expr_scaled, col_df)
  })

  # downsample per sample
  ds_list <- lapply(data_list, function(df) {
    df <- df |>
      dplyr::group_by(sample_id) |>
      dplyr::slice_sample(n = n_cells_per_sample) |>
      dplyr::ungroup() |>
      dplyr::mutate(row_id = paste(run_id, cell.id, sep = "_"))

  })

  # combine
  total_df = dplyr::bind_rows(ds_list)
  total_df = tibble::column_to_rownames(total_df, "row_id")
  total_df = total_df[, colSums(is.na(total_df)) == 0]

  # select marker columns
  markers = rownames(sce_list[[1]])
  expr_mat = as.matrix(total_df[, markers])

  # find the scater run* function
  dr_method_upper = toupper(dr_method)

  # get all exported scater functions starting with run
  run_methods = getNamespaceExports("scater")
  run_methods = run_methods[grepl("^run", run_methods)]

  run_methods_upper = toupper(sub("^run", "", run_methods))

  if (!(dr_method_upper %in% run_methods_upper)) {
    stop("Unknown dimension reduction method in scater: ", dr_method)
  }

  # find corresponding function
  dr_fun_name = run_methods[which(run_methods_upper == dr_method_upper)]
  dr_fun = get(dr_fun_name, envir = asNamespace("scater"))

  # run DR with additional arguments
  temp_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = t(expr_mat))  # scater expects genes x cells
  )

  dr_result = dr_fun(temp_sce, exprs_values = "counts" , ...)

  # combine
  dr_df <- as.data.frame(
    SingleCellExperiment::reducedDim(dr_result, dr_method_upper)
  )
  n_dim = ncol(dr_df)

  if (dr_method_upper != "PCA"){
  colnames(dr_df)[1:min(2, n_dim)] = paste0("Dim", 1:min(2, n_dim))}

  return(cbind(dr_df, total_df))
}
