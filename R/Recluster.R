#' Recluster
#'
#' This function reclusters doublets and non doublets seperately based on deconvolution analysis and returns a new expression file, groups file, and DeconFreq table for downstream analyses.
#'
#' @param isADoublet isADoublet data.frame from IsDoublet.
#' @param data Processed data from CleanUpInput (or RemoveCellCycle).
#' @param groups Processed groups file from CleanUpInput.
#' @param log_file_name used for saving run notes to log file
#'
#' @return newData2 - processed expression and groups file, reordered, with correct new cluster numbers.
#' @return decon - DeconCalledFreq table with all doublets 100 percent doublet and all non doublets at 0 percent doublet frequency.
#'
#' @keywords recluster
#'
#' @export
#'
Recluster <- function(
  isADoublet,
  data,
  groups,
  log_file_name
) {
  #Get list of doublet samples
  doubletCells <- row.names(subset(isADoublet, isADoublet == TRUE))
  notDoubletCells <- row.names(subset(isADoublet, isADoublet == FALSE))

  #Make new expression tables for doublets and non-doublets individually
  doubletCellsData <- data[2:nrow(x = data), colnames(x = data) %in% doubletCells]
  notDoubletCellsData <- data[2:nrow(x = data), colnames(x = data) %in% notDoubletCells]
  doubletCellsData <- t(as.matrix(doubletCellsData))
  notDoubletCellsData <- t(as.matrix(notDoubletCellsData))

  #Use "clustering" by deconvolution predicted groups.
  nondoublets <- groups[notDoubletCells, ] #non doublets groups
  doublets <- NULL #doublets groups
  isADoublet2 <- subset(isADoublet, isADoublet == TRUE)
  isADoublet2 <- isADoublet2$Cell_Types
  partialGroups <- cbind(isADoublet2, isADoublet2)
  row.names(x = partialGroups) <- row.names(x = subset(isADoublet, isADoublet == TRUE))
  colnames(x = partialGroups) <- colnames(x = groups)
  doublets <- partialGroups[order(partialGroups[, 1]), ]
  isADoublet3 <- doublets[, 2]
  doublets[, 1] <- as.integer(x = as.factor(x = doublets[, 1]))
  doublets[, 2] <- as.integer(x = as.factor(x = doublets[, 2]))
  doublets[, 1] <- as.integer(x = doublets[, 1]) + length(x = unique(nondoublets[, 1]))
  doublets[, 2] <- isADoublet3

  #merge these files together and create new groups classification
  newGroups <- rbind(nondoublets, doublets)

  #Make Decon frequency table to return with 0% and 100%
  uniqueClusters <- as.character(x = unique(newGroups[, 2]))
  DeconCalledFreq <- as.data.frame(matrix(nrow = length(x = uniqueClusters), ncol = 1), row.names = uniqueClusters)
  DeconCalledFreq[1:length(x = unique(nondoublets[, 1])), 1] <- 0
  DeconCalledFreq[(length(x = unique(nondoublets[, 1])) + 1):nrow(x = DeconCalledFreq), 1] <- 100

  #create new reordered expression file for return
  newData <- data[, match(row.names(x = newGroups), colnames(x = data)) ]
  if (colnames(x = data)[1] %in% 'row_clusters.flat' || colnames(x = data)[1] %in% 'row_clusters-flat') {
    newData2 <- CleanUpInput(newData, newGroups, data[2:nrow(x = data), 1], log_file_name = log_file_name)
  } else {
    newData2 <- CleanUpInput(newData, newGroups, rowClusters = NULL, log_file_name = log_file_name)
  }
  return(list(newData2 = newData2, decon = DeconCalledFreq))
}
