#' Example channel-marker table dataframe
#' 
#' Dataframe mapping channel to marker from the cytometry experiment. 
#' Empty channels and irrelevant markers/channels have been filtered out with createChannelMarkerTable(). 
#'
#' @format A data frame with 42 rows and 2 variables:
#' \describe{
#'   \item{channel}{Name of channel}
#'   \item{markers}{Name of marker}
#' }
"eg_channelMarkerTbl"

#' Example of raw (unfiltered) channel-marker table dataframe
#' 
#' Dataframe mapping channel to marker from the cytometry experiment. 
#' Comes directly from ncdfFlow::read/ncdfFlowSet(). 
#' Has not been filtered to remove irrelevant/empty channels. 
#' Depending on how user set up experiment prior to computational processing, both the channel and marker may appear in the marker channel.
#'
#' @format A data frame with 85 rows and 2 variables:
#' \describe{
#'   \item{channel}{Name of channel}
#'   \item{markers}{Name of marker}
#' }
"eg_channelMarkerTblUnfilt"


#' Example flowSet
#' 
#' Example flowSet containing example fcs files
#'
#' @format A flowSet with two fcs files.
"eg_fset"


#' Example raw SingleCellExperiment object
#' 
#' Example SingleCellExperiment object generated from example flowSet.
#'
#' @format A SingleCellExperimentObject
"eg_sce_raw"

#' Example SingleCellExperiment object with annotation labels
#' 
#' Example SingleCellExperiment object with clustering and annotation labels.
#'
#' @format A SingleCellExperimentObject
"eg_annotatedSCE"


#' Example lineage table
#' 
#' Example lineage table specifying hierarchy of cell types
#'
#' @format A dataframe
#' \describe{
#'   \item{label}{Name of cell type}
#'   \item{parent}{Name of parent of cell type}
#'   \item{value}{Markers composing a cell type separated by a comma. Markers negative for a cell type are noted as marker name + "neg"}
#' }
"eg_lineageTable"

#' Example rules table
#' 
#' Example rules table specifying rules for how to call each cell type
#'
#' @format A dataframe
#' \describe{
#'   \item{label}{Name of cell type}
#'   \item{rule}{Name of rule used to call the cell type}
#'   \item{threshold}{Value denoting positive/negative for a metric; only applicable for certain rules}
#' }
"eg_rulesTable"

#' Example singscore values
#' 
#' Example singscore values for each cluster
#'
#' @format A dataframe
#' \describe{Generated for each row in lineage table (unless only negative markers in a cell type)}
"eg_singscores"


#' Example expression values
#' 
#' Example scaled median expression values for each cluster
#'
#' @format A dataframe
#' \describe{Generated for each cluster}
"eg_expr"

#' Example output of annotateClusters()
#' 
#' Example annotation dataframe detailing annotation label per cluster.
#'
#' @format A dataframe
#' \describe{Generated for each cluster based on lineage table and rules}
"eg_annoTbl"

#' Example annotated SCE object
#' 
#' Example SCE object post-clustering and post-annotation steps.
#'
#' @format SingleCellExperiment object
"eg_annotatedSCE"


#' Example dataframe of independent annotation
#' 
#' Dataframe of detailing if clusters passed the marker threshold.
#'
#' @format A dataframe
"eg_indepAnno"


#' Example dataframe of counts per cluster and sample
#'
#' @format A dataframe
"eg_clusterCounts"


#' Example dataframe of counts per label and sample
#' 
#' Counts concatenated across clusters with same annotation label.
#'
#' @format A dataframe
"eg_labelCounts"