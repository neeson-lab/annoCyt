.cleanMarkers <- function(channelmarker.tbl, format) {

  # Locate keyword positions
  channel_pos <- regexpr("channel", format)
  marker_pos  <- regexpr("marker",  format)

  if (channel_pos == -1 || marker_pos == -1) {
    stop("Format must contain both 'channel' and 'marker', e.g. 'channel_marker' or 'marker-channel'.")
  }

  # Determine order
  channel_first <- channel_pos < marker_pos

  # Extract separator
  if (channel_first) {
    sep <- substr(format, channel_pos + nchar("channel"), marker_pos - 1)
  } else {
    sep <- substr(format, marker_pos + nchar("marker"), channel_pos - 1)
  }

  if (nchar(sep) == 0) {
    stop("Format must contain a separator between 'channel' and 'marker', e.g. 'channel_marker'.")
  }

  # Extract prefix and suffix
  first_keyword_pos <- min(channel_pos, marker_pos)
  prefix <- if (first_keyword_pos > 1) substr(format, 1, first_keyword_pos - 1) else ""

  last_keyword_end <- if (channel_first) {
    marker_pos + nchar("marker") - 1
  } else {
    channel_pos + nchar("channel") - 1
  }
  suffix <- if (last_keyword_end < nchar(format)) substr(format, last_keyword_end + 1, nchar(format)) else ""

  # Strip prefix and suffix
  if (nchar(prefix) > 0) channelmarker.tbl$markers <- sub(paste0("^", .regexEscape(prefix)), "", channelmarker.tbl$markers)
  if (nchar(suffix) > 0) channelmarker.tbl$markers <- sub(paste0(.regexEscape(suffix), "$"), "", channelmarker.tbl$markers)

  # Filter to entries containing the separator
  channelmarker.tbl <- channelmarker.tbl[grepl(sep, channelmarker.tbl$markers, fixed = TRUE), ]

  # Split into channel and marker
  # parts <- strsplit(channelmarker.tbl$markers, regexEscape(sep), perl = TRUE)

  channelmarker.tbl$markers <- if (channel_first) {
    sub(paste0(".*", sep), "", channelmarker.tbl$markers)
  } else {
    sub(paste0(sep, ".*"), "", channelmarker.tbl$markers)
  }

  return(channelmarker.tbl)
}

# Helper to escape special regex characters in user-supplied strings
.regexEscape <- function(x) {
  gsub("([.\\\\|*?+()\\[\\]{}^$#])", "\\\\\\1", x, perl = TRUE)
}

.customRead.ncdFlowSet = function(filename, isWriteSlice=FALSE, transformation=FALSE, common_channels=NULL, channel2marker_LOOKUP=NULL,  clean.marker_FUN=NULL, clean.marker_mapping=NULL,...) {

  if (!is.character(filename)) {
    stop("Expecting only a single .fcs file path")
  }
  nfs = ncdfFlow::read.ncdfFlowSet(filename, isWriteSlice = isWriteSlice, transformation=transformation,...)

  #apply common channels filtering if given
  if (!is.null(common_channels)) {
    nfs = nfs[,common_channels]
  }

  channelmarker.tbl = channel2marker_LOOKUP |> dplyr::mutate(clean.markers=markers)


  if (!is.null(clean.marker_FUN)) {
    ##-- rename markers naming into the cleaned version.
    #(a)create channel-marker-clean.marker table
    channelmarker.tbl = channelmarker.tbl |>
      dplyr::mutate(clean.markers=clean.marker_FUN(markers))
  }


  #(b) apply marker mapping if given
  if (!is.null(clean.marker_mapping)) {
    #convert cleanmarkers mapping from list to dataframe.
    mm = data.frame(clean.marker_mapping, check.names = FALSE) |> t() |> data.frame(check.names = FALSE)
    colnames(mm) = c("alias")
    mm = mm |> tibble::rownames_to_column("clean.markers") |> dplyr::mutate_if(is.factor, as.character)

    channelmarker.tbl = channelmarker.tbl |>
      dplyr::left_join(mm, by="clean.markers") |>
      dplyr::mutate(final.clean.name=dplyr::coalesce(alias, clean.markers)) |>
      dplyr::mutate(clean.markers = final.clean.name)
  }
  #(c) apply clean.marker names
  channel.marker = channelmarker.tbl |> dplyr::select(channel, clean.markers)
  channel.marker.data = channel.marker$clean.markers
  names(channel.marker.data) = channel.marker$channel
  flowCore::markernames(nfs) = channel.marker.data
  return(nfs)
}


#' Creates dataframe mapping channel to marker
#'
#' @param metadata Dataframe containing sample-level metadata information
#' @param exclude.FSC.SSC Logical specifying whether to exclude FSC and SSC
#' @param rm.markers Character vector specifying which markers to remove, if any
#' @param rm.channels Character vector specifying which channels to remove, if any
#' @param format String specifying the format of the channel and marker as it was exported in the .fcs files. This will separate them accordingly
#'
#' @return Dataframe specifying which channel each marker is mapped to
#' @export
#'
#' @examples
#' eg_sampleMetadata = data.frame(
#' file_name=system.file("extdata",
#' "eg_rawdata.fcs",
#' package = "annoCyt"),
#' sample_name="control",
#' condition="control",
#' patient_id="control",
#' run_id="run1",
#' sample_id="run1_control")
#' 
#' createChannelMarkerTable(eg_sampleMetadata, format = "channel_marker", 
#'     rm.markers = c("Cisplatin", "DNA1", "DNA2", "FVS"))
createChannelMarkerTable = function(metadata, exclude.FSC.SSC = T, rm.markers = NULL, rm.channels = NULL, format = NULL) {
  distinctRun.info = metadata[1,c("file_name", "run_id")]
  nfs = suppressMessages(ncdfFlow::read.ncdfFlowSet(distinctRun.info$file_name[1], isWriteSlice = FALSE, transformation=FALSE))

  channelmarker.tbl = nfs  |>
    flowCore::markernames() |>
    as.list() |>
    utils::stack() |>
    dplyr::rename(markers = values, channel = ind) |>
    dplyr::relocate(channel, .before = markers)

    if (!is.null(format)) {
      channelmarker.tbl = .cleanMarkers(channelmarker.tbl, format)
    }

  if (!is.null(rm.markers)) {
    channelmarker.tbl = dplyr::filter(channelmarker.tbl, !(markers %in% rm.markers))
  }

  if (!is.null(rm.channels)) {
    channelmarker.tbl = dplyr::filter(channelmarker.tbl, !(channel %in% rm.channels))
  }

  return(channelmarker.tbl)
}


#' Create flowset of .fcs files
#'
#' @param sample.annot Dataframe of sample-level metadata information
#' @param channelmarker.tbl Dataframe specifying which channel each marker is mapped to
#' @param seed Set seed for reproducibility
#'
#' @return A flowSet object.
#' @export
#'
#' @examples
#' eg_sampleMetadata = data.frame(
#' file_name=system.file("extdata",
#' "eg_rawdata.fcs",
#' package = "annoCyt"),
#' sample_name="control",
#' condition="control",
#' patient_id="control",
#' run_id="run1",
#' sample_id="run1_control")
#' createFlowset(eg_sampleMetadata, eg_channelMarkerTbl)
createFlowset = function(sample.annot, channelmarker.tbl, seed  =42){
  set.seed(seed)

  # run info
  allRun.info = sample.annot |> dplyr::select(file_name, run_id)
  distinctRun.info = sample.annot |>
    dplyr::group_by(run_id) |>
    dplyr::filter(dplyr::row_number() == 1) |>
    dplyr::ungroup() |>
    dplyr::select(file_name, run_id)

  ## markers mapping, per run, only only if applicable.
  runs.markers.mapping=NULL

  #1) Check the first file of each distinct runs, and get common channels
  common.channels = NULL
  for(fcspath in distinctRun.info$file_name) {
    afcs = .customRead.ncdFlowSet(fcspath, isWriteSlice = FALSE, channel2marker_LOOKUP=channelmarker.tbl, clean.marker_FUN=NULL) |> colnames()
    if (is.null(common.channels)) {
      common.channels = afcs
    } else {
      common.channels = intersect(common.channels, afcs)
    }
  }

  #1b) validation check: check for markers discrepancy
  # channelmarker.tbl=NULL
  # for(ridx in 1:nrow(distinctRun.info)) {
  #   rowdata = distinctRun.info[ridx,]
  #   fcs.runId = rowdata$run_id
  #   fcs.metadata.common = customRead.ncdFlowSet(rowdata$file_name, isWriteSlice = FALSE, channel2marker_LOOKUP=def.channel.marker, clean.marker_FUN=NULL, common_channels = common.channels, clean.marker_mapping=NULL)
  #
  #   fcs.channelmarker.tbl = getChannelMarkerTable(fcs.metadata.common) |>
  #                           mutate(!!fcs.runId := markers) |>
  #
  #                             dplyr::select(-markers)
  #   if(is.null(channelmarker.tbl)) {
  #     channelmarker.tbl = fcs.channelmarker.tbl
  #   }
  #   else {
  #     #should not happen as using common channels
  #     if (sum(channelmarker.tbl$channel != fcs.channelmarker.tbl$channel)>0) {
  #       stop("Error: channel order must be the same among .fcs files")
  #     }
  #     channelmarker.tbl = channelmarker.tbl |> left_join(fcs.channelmarker.tbl, by="channels")
  #   }
  # }
  #
  # ## Stops if there is markers discrepancies
  # count.run.diff.data =  channelmarker.tbl |>
  #                           select(-channels) |>
  #                           apply(MARGIN = 1, FUN=n_distinct)
  # channelmarker.tbl.filtered = channelmarker.tbl |>
  #                                 mutate(ndistinct=count.run.diff.data-1) |>
  #                                 filter(ndistinct > 0) |>
  #                                 dplyr::select(-ndistinct)
  # if (nrow(channelmarker.tbl.filtered) > 0) {
  #   print(channelmarker.tbl.filtered)
  #   stop("Unresolved markers discrepancy across runs.")
  # }
  #
  # #stops if the expected canonical panel marker is missing
  # unidentified.markers = lineage_channels  |> setdiff(channelmarker.tbl[,2] )
  # if (length(unidentified.markers) > 0) {
  #   stop(si("Unidentified markers: `paste(unidentified.markers, collapse=';')`"))
  # }



  #2) Check the first file of each distinct runs, and get common channels
  #allframes= vector("list", length(allRun.fcspath))
  #i=1

  allNcdflowset=NULL
  for(ridx in 1:nrow(allRun.info)) {
    rowdata = allRun.info[ridx,]
    fcs.runId = rowdata$run_id
    fobj.subset = .customRead.ncdFlowSet(rowdata$file_name, isWriteSlice = TRUE, transformation=FALSE,
                                        channel2marker_LOOKUP=channelmarker.tbl,
                                        clean.marker_FUN=NULL,
                                        common_channels = common.channels,
                                        clean.marker_mapping=runs.markers.mapping[[fcs.runId]],
                                        truncate_max_range = FALSE)

    #fobj = flowCore::read.FCS(fcspath)
    #fobj = ncdfFlow::read.ncdfFlowSet(fcspath)
    #fobj.subset = fobj[,common.channels]

    if (is.null(allNcdflowset)) {
      allNcdflowset = fobj.subset
    }
    else {
      allNcdflowset = methods::rbind2(allNcdflowset, fobj.subset)
    }

    #allframes[[i]] = fobj.subset
    #i = i+1
  }
  #fs = flowSet(allframes)
  global.fset = allNcdflowset

  return(global.fset)
}


# ### get emrging1 assignment per cell
# getMerging1Assg = function(sceObj) {
#   rd = colData(sceObj) |>
#     as.data.frame() |>
#     dplyr::rename(cluster_som_id=cluster_id) |>
#     dplyr::mutate(cluster_som_id = as.character(cluster_som_id)) |>
#     dplyr::select(cluster_som_id)
#   clustering.somID = metadata(sceObj)$cluster_codes[,"som100"] |> as.character()
#   clustering.merging1 = metadata(sceObj)$cluster_codes[,"merging1"] |> as.character()
#   clustering.tbl = data.frame(cluster_som_id=clustering.somID, merging1=clustering.merging1, stringsAsFactors = FALSE)
#
#   res = rd |> left_join(clustering.tbl, by="cluster_som_id")
#   return(res$merging1)
# }
#
# getMetaClusterAssg = function(sceObj,k=26) {
#   rd = colData(sceObj) |>
#     as.data.frame() |>
#     dplyr::rename(cluster_som_id=cluster_id) |>
#     dplyr::mutate(cluster_som_id = as.character(cluster_som_id)) |>
#     dplyr::select(cluster_som_id)
#   clustering.somID = metadata(sceObj)$cluster_codes[,"som100"] |> as.character()
#   clustering.metaCluster = metadata(sceObj)$cluster_codes[,si("meta`k`")] |> as.character()
#   clustering.tbl = data.frame(cluster_som_id=clustering.somID, metaCluster=clustering.metaCluster, stringsAsFactors = FALSE)
#
#   res = rd |> left_join(clustering.tbl, by="cluster_som_id")
#   return(res$metaCluster)
# }


#' Create SingleCellExperiment object for your cytometry data
#'
#' @param channel.marker.table Dataframe mapping the marker with the cytometry channel
#' @param sample.annot Dataframe of sample-level metadata information
#' @param global.fset flowSet object
#' @param cofactor.val Cofactor value for arcsinh-transformation
#' @param ... Passed to CATALYST::prepData()
#'
#' @return A SingleCellExperiment object containing expression data and metadata
#' @export
#'
#' @examples
#' \dontrun{
#' sce = createSCE(eg_channelMarkerTbl, eg_sampleMetadata, eg_fset, cofactor.val = 5)
#' }

createSCE = function(channel.marker.table,
                    sample.annot,
                    global.fset,
                    cofactor.val = NULL,
                  ...
){

  global.catalyst.panel = channel.marker.table |>
    dplyr::rename(fcs_colname = channel, antigen = markers) |>
    dplyr::mutate(marker_class="type")

  #use external annot dataset, more flexible to change..
  global.catalyst.md =  dplyr::mutate(sample.annot, fn=basename(file_name))

  md.factors = colnames(sample.annot)
  md.factors = md.factors[!md.factors %in% c("file_name", "sample_id")]

  sceobj = CATALYST::prepData(global.fset, global.catalyst.panel, global.catalyst.md, cofactor=cofactor.val, FACS=T,
                          panel_cols = list(channel = "fcs_colname", antigen = "antigen"),
                          md_cols = list(file = "fn", id = "sample_id",
                                         factors = md.factors), ...)

    #assign cellId from sample_id and sequential number
    md = SummarizedExperiment::colData(sceobj) |>
      as.data.frame() |>
      dplyr::group_by(sample_id) |>
      dplyr::mutate(cellid= dplyr::row_number()) |>
      dplyr::mutate(cellid=paste(sample_id, cellid, sep="_"))

    SummarizedExperiment::colData(sceobj)$cell.id = md$cellid

  return(sceobj)
}


sampleMedianExprHeatmap = function(sceobj, markers=NULL, ...) {
  expr = sceobj@assays@data$exprs |> t() |> as.data.frame()
  rdata = SingleCellExperiment::colData(sceobj) |> as.data.frame()

  effmarkers =  markers
  if(is.null(markers)) {
    effmarkers =  colnames(expr)
  }

  expr.med.bySmplId = expr[,effmarkers, drop=FALSE] |>
    dplyr::bind_cols(dplyr::select(rdata, "sample_id")) |>
    dplyr::group_by(sample_id) |>
    dplyr::summarize_all(median) |>
    dplyr::ungroup() |>
    tibble::column_to_rownames("sample_id")

  palette = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  hm_cols = grDevices::colorRampPalette(palette)(100)
  hm1 = ComplexHeatmap::Heatmap(matrix = t(expr.med.bySmplId),
                                 col = hm_cols,
                                 name =  "Normalised Expression",
                                 column_names_gp = grid::gpar(fontsize = 8),
                                 row_names_gp = grid::gpar(fontsize = 8),
                                 rect_gp = grid::gpar(col = "white"), na_col = "lightgrey",
                                 cluster_rows = TRUE, cluster_columns = TRUE,
                                 clustering_distance_columns = "pearson",
                                 show_row_dend = FALSE, show_column_dend = FALSE,
                                 column_names_side = "top", row_names_side = "left", ...)
  return(hm1)
}

# ## Note: hacked in cluster() of daFrame, taking expr matrix (rather than flowset), and return FlowSOM som id directly (no concensusclustering run)
# cluster.expr.somIds = function(exprMatrix, cols_to_use = NULL, xdim = 10, ydim = 10,  seed = 1) {
#
#   sel.cols = ifelse(is.null(cols_to_use), colnames(exprMatrix), cols_to_use)
#   fsom = ReadInput(flowFrame(as.matrix(exprMatrix)))
#   set.seed(seed)
#   som = BuildSOM(fsom, colsToUse = cols_to_use, silent = TRUE,
#                   xdim = xdim, ydim = ydim)
#
#   clusterId = factor(som$map$mapping[, 1])
#   return(clusterId)
# }



#' Get unscaled median expression of each marker per cluster
#'
#' @param sceobj A SingleCellExperiment object
#' @param cluster_code A cluster code available names(cluster_codes(sce))
#'
#' @return A dataframe of median expression per marker per cluster
#' @export
#'
#' @examples
#' \dontrun{
#' getRawExpr(sce)
#' }
#' 
getRawExpr = function(sceobj, cluster_code = "som100"){

  data.coldata = as.data.frame(SingleCellExperiment::colData(sceobj))
  data.coldata$cluster_id = as.character(data.coldata$cluster_id)
  data.coldata = data.coldata |>
    dplyr::inner_join(
      as.data.frame(S4Vectors::metadata(sceobj)$cluster_codes)[,c("som100", cluster_code)],
      by = c("cluster_id" = "som100")
    )

  if (cluster_code != "som100"){
    data.coldata = dplyr::select(data.coldata, -cluster_id) |>
      dplyr::rename(cluster_id = !!cluster_code)
  }

  data.expr = as.data.frame(t(sceobj@assays@data$exprs))

  data.expr.median = data.expr |>
    dplyr::bind_cols(data.coldata) |>
    dplyr::group_by(eff.metacluster) |>
    dplyr::summarize_if(is.numeric,.funs=median) |>
    reshape2::melt(id.vars="eff.metacluster", variable.name="marker", value.name="median.expr")

  return(data.expr.median)
}


#' Get scaled median expression of each marker per cluster
#'
#' @param sceobj A SingleCellExperiment object
#' @param cluster_code A cluster code available in names(cluster_codes(sce))
#' @param scale.margin Value to scale values by
#'
#' @return A dataframe of scaled median expression per marker per cluster
#' @export
#'
#' @examples
#' \dontrun{
#' df.expr = scaleClusterExpr(sce) # run post clustering
#' }

scaleClusterExpr = function(sceobj, cluster_code = "som100", scale.margin=2) {

  if (!inherits(sceobj, "SingleCellExperiment")) {
    stop("sceobj must be a SingleCellExperiment object.")
  }

  if (!(cluster_code %in% names(CATALYST::cluster_codes(sceobj)))){
    stop("cluster_code not valid; must be present in names(cluster_codes(x))")
  }


  data.coldata = as.data.frame(SingleCellExperiment::colData(sceobj))
  data.coldata$cluster_id = as.character(data.coldata$cluster_id)
  data.coldata = data.coldata |>
    dplyr::inner_join(
      as.data.frame(S4Vectors::metadata(sceobj)$cluster_codes)[,c("som100", cluster_code)],
      by = c("cluster_id" = "som100")
    )

  if (cluster_code != "som100"){
    data.coldata = dplyr::select(data.coldata, -cluster_id) |>
    dplyr::rename(cluster_id = !!cluster_code)
    }

  data.expr = as.data.frame(t(sceobj@assays@data$exprs))

  data.expr.median.scaled = data.expr |> CATALYST:::.scale_exprs(margin = scale.margin) |> as.data.frame() |>
    dplyr::bind_cols(data.coldata) |>
    dplyr::group_by(cluster_id) |>
    dplyr::summarize_if(is.numeric, .funs =  median) |>
    dplyr::arrange(as.numeric(cluster_id))

  return(data.expr.median.scaled)

}
