#' @include internal.R
#'
NULL

#' Seurat V3 Pre Process
#'
#' Pre-process Seurat V3 object for DoubletDecon
#'
#' @param object Seurat object which has been pre-processed with clustering and dimensional reduction
#' @param assay Name of assay to pull the count data from; default is 'RNA'.
#' @param markers Output of FindAllMarkers. If NULL, will recalculate it with default parameters.
#' @param top.n Number of markers for the top_n function. Default is 50.
#' @param group.by Name of one metadata columns to group cells by; If NULL, will use current identity.
#' @param output.dir If set, save the output files to .txt format. Defauly is NULL.
#'
#' @return a list with three elements:
#' \code{newExpressionFile} - Seurat expression file in ICGS format (ICGS genes)
#' \code{newFullExpressionFile} - Seurat expression file in ICGS format (all genes)
#' \code{newGroupsFile} - Groups file ICGS format
#'
#' @references DePasquale, E. A., Schnell, D. J., Van Camp, P. J., Valiente-Alandí, Í., Blaxall, B. C., Grimes, H. L., ... & Salomonis, N. (2019). DoubletDecon: Deconvoluting Doublets from Single-Cell RNA-Sequencing Data. Cell reports, 29(6), 1718-1727.
#'
#' @import Seurat
#' @import dplyr
#'
#' @keywords Seurat
#'
#' @export
#'
PrepDoubletDecon <- function(
  object,
  assay = NULL,
  markers = NULL,
  top.n = 50,
  group.by = NULL,
  output.dir = NULL
) {
  assay <- assay %||% DefaultAssay(object = object)

  #extract expression file
  counts <- as.data.frame(GetAssayData(object = object, assay = assay, slot = 'counts'))

  #extract clusters
  if (!is.null(group.by)) {
    Idents(object = object) <- object[[group.by]]
  }
  clusters <- as.data.frame(Idents(object = object))

  #extract top 50 markers
  if (is.null(markers)) {
    markers <- FindAllMarkers(object = object, assay = assay, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  }
  cluster.check <-  levels(clusters[, 1]) %in% levels(markers[, 'cluster'])
  if (any(cluster.check == FALSE)) {
    stop('Grouping variable is not matching FindAllMarkers output')
  }
  markers <- markers %>% group_by(cluster) %>% top_n(n = top.n, wt = avg_logFC)

  #Find and replace "-"
  colnames(counts) <- gsub("-",".", colnames(x = counts))
  clusters[, 1] <- gsub("-", ".", clusters[, 1])

  #Start cluster numbers at 1
  if (class(x = markers$cluster) == 'factor'){
    if (min(x = as.numeric(x = as.character(x = markers$cluster))) == 0){
      markers$cluster <- as.numeric(x = as.character(x = markers$cluster)) + 1
      clusters[, 1] <- as.numeric(x = clusters[, 1]) + 1
    } else if (min(x = as.numeric(x = as.character(x = markers$cluster))) == 1){
      markers$cluster <- as.numeric(x = as.character(x = markers$cluster))
      clusters[, 1] <- as.numeric(x = clusters[,1])
    } else {
      stop('Unexpected cluster numbering scheme. Cluster numbers are expected to be continuous numbers starting from either 0 or 1. Please check conversion for correctness following this function.')
    }
  }

  #Reorder clusters file based on cluster number (this will be the new order for the cells)
  clusters2 <- clusters[order(clusters[, 1]), , drop = FALSE]

  #Reorder columns
  counts <- counts[, row.names(x = clusters2)]

  #Reorder genes file based on gene cluster number (this will be the new order for the genes)
  markers2 <- markers[order(markers$cluster), ]

  #Subset expression file for genes in the genes file
  full.counts <- counts
  counts <- counts[row.names(x = counts) %in% as.character(x = markers$gene),]

  #Reorder genes
  geneOrder <- intersect(x = markers2$gene, y = as.character(x = row.names(x = counts)))
  counts <- counts[match(geneOrder, row.names(x = counts)), ]

  #Add columnn_clusters-flat
  full.counts <- rbind(clusters2[, 1], full.counts)
  counts <- rbind(clusters2[,1], counts)
  row.names(x = full.counts)[1] <- 'column_clusters-flat'
  row.names(x = counts)[1] <- 'column_clusters-flat'

  #Add row_clusters-flat
  markers3 <- markers2[match(geneOrder, markers2$gene), ]
  rowToAdd <- c(NA, markers3$cluster)
  rowToAdd2 <- rep(NA, nrow(x = full.counts))
  counts <- cbind(rowToAdd, counts)
  full.counts <- cbind(rowToAdd2, full.counts)
  colnames(counts)[1] <- 'row_clusters-flat'
  colnames(full.counts)[1] <- 'row_clusters-flat'

  #Make groups file
  groups <- cbind(as.numeric(x = counts[1, 2:ncol(x = counts)]), as.numeric(x = counts[1, 2:ncol(x = counts)]))
  row.names(x = groups) <- as.character(x = colnames(x = counts)[2:ncol(x = counts)])

  if (!is.null(output.dir)) {
    write.table(counts, file = file.path(output.dir, 'ICGS_expression.txt'), sep = '\t')
    write.table(full.counts, file = file.path(output.dir, 'ICGS_fullExpression.txt'), sep = '\t')
    write.table(groups, file = file.path(output.dir, 'ICGS_groups.txt'), sep = '\t', col.names = F)
  }

  return(list(newExpressionFile = counts,
              newFullExpressionFile = full.counts,
              newGroupsFile = groups))
}
