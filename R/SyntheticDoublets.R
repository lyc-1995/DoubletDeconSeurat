#' Synthetic Doublets
#'
#' This function creates synthetic doublets by averaging gene expression profiles from each combination of clusters to generate deconvolution profiles for each type of doublet.
#'
#' @param data Processed data from Clean_Up_Input (or Remove_Cell_Cycle).
#' @param groups Processed groups file from Clean_Up_Input.
#' @param groupsMedoids New groups file based on blacklisted clusters from Blacklist_Groups
#' @param newMedoids New combined medoids from Blacklist_groups.
#' @param num_doubs The user defined number of doublets to make for each pair of clusters
#' @param only50 Use only synthetic doublets created with 50\%/50\% mix of parent cells, as opposed to the extended option of 30\%/70\% and 70\%/30\%, default is TRUE.
#' @param location Directory where output should be stored
#' @param seed Set a random seed, default is 100.
#'
#' @return averagesAverages = average deconvolution profiles for each combination of cell types.
#' @return doubletCellsInput2 = deconvolution profiles for synthetic doublet cells (for quality control).
#'
#' @keywords synthetic
#'
#' @export
#'
SyntheticDoublets <- function(
  data,
  groups,
  groupsMedoids,
  newMedoids,
  num_doubs,
  only50,
  location,
  seed = 100
){
  set.seed(seed = seed)

  #Override the original groups for making synthetics with groups based on the blacklisted clusters
  groups <- groupsMedoids

  #Number of doublets created per cluster pair
  ndub <- num_doubs

  #Make data frame to hold new doublets (30 per combination of clusters)
  pairs <- combn(unique(groups[, 2]), 2)
  doubletCellsInput <- as.data.frame(matrix(ncol = 2, nrow = ((length(x = pairs)/2)*ndub)))
  doubletCellsInput2 <- as.data.frame(matrix(ncol = 4, nrow = ((length(x = pairs)/2)*ndub))) #for info

  #For each combination of medoids, select 30 pairs
  for (pair in 0:((length(x = pairs)/2) - 1)) {
    for (synth in 1:ndub) {
      doubletCellsInput[synth + (ndub*pair), 1] <- sample(row.names(x = subset(groups, groups[, 2] == pairs[1, (pair + 1)])), 1, replace = FALSE)
      doubletCellsInput2[synth + (num_doubs*pair), 1] <- doubletCellsInput[synth + (num_doubs*pair), 1]
      doubletCellsInput[synth + (ndub*pair), 2] <- sample(row.names(x = subset(groups, groups[, 2] == pairs[2, (pair + 1)])), 1, replace = FALSE)
      doubletCellsInput2[synth + (num_doubs*pair), 2] <- doubletCellsInput[synth + (num_doubs*pair), 2]
      doubletCellsInput2[synth + (num_doubs*pair), 3] <- pairs[1, (pair + 1)]
      doubletCellsInput2[synth + (num_doubs*pair), 4] <- pairs[2, (pair + 1)]
    }
  }

  doubletAverages <- as.data.frame(matrix(ncol = nrow(x = doubletCellsInput), nrow = nrow(x = data)))
  row.names(x = doubletAverages) <- row.names(x = data)
  doubletAverages[1, ] <- rep((length(x = unique(groups[, 1])) + 1), ncol(x = doubletAverages))

  doubletAverages <- as.data.frame(matrix(ncol = nrow(x = doubletCellsInput)*3, nrow = nrow(x = data)))
  row.names(x = doubletAverages) <- row.names(x = data)

  doubletAverages[1, ] <- rep((length(unique(groups[, 1])) + 1), ncol(x = doubletAverages)*3)

  for (doublet in 1:(ncol(x = doubletAverages)/3)) {
    cell1 <- as.character(x = doubletCellsInput[doublet, 1])
    cell2 <- as.character(x = doubletCellsInput[doublet, 2])
    expression1 <- data[2:nrow(x = data), which(colnames(x = data) == cell1)]
    expression2 <- data[2:nrow(x = data), which(colnames(x = data) == cell2)]
    temp <- cbind(expression1, expression2)
    newExpression <- apply(temp, 1, weighted.mean, c(0.5, 0.5))
    newExpression_a <- apply(temp, 1, weighted.mean, c(0.70, 0.30))
    newExpression_b <- apply(temp, 1, weighted.mean, c(0.30, 0.70))
    doubletAverages[2:nrow(x = doubletAverages), doublet] <- newExpression
    doubletAverages[2:nrow(x = doubletAverages), doublet + (ncol(x = doubletAverages)/3)] <- newExpression_a
    doubletAverages[2:nrow(x = doubletAverages), doublet + ((ncol(x = doubletAverages)/3)*2)] <- newExpression_b
    colnames(x = doubletAverages)[doublet] <- paste0(cell1, '-', cell2, '-even')
    colnames(x = doubletAverages)[doublet + (ncol(x = doubletAverages)/3)] <- paste0(cell1, '-', cell2, '-one')
    colnames(x = doubletAverages)[doublet + ((ncol(x = doubletAverages)/3)*2)] <- paste0(cell1, '-', cell2, '-two')
  }

  #50%/50% references only or 30%/70% and 70%/30% included
  if (only50 == TRUE) {
    mult <- 1
  } else {
    mult <- 3
  }

  results <- DeconRNASeq(doubletAverages[2:nrow(x = doubletAverages), ], newMedoids)
  resultsreadable <- round(results$out.all*100,2)
  write.table(resultsreadable, file = file.path(location, 'resultsreadable_synths.txt'), sep = '\t')
  row.names(x = resultsreadable) <- colnames(x = doubletAverages)

  averagesAverages <- as.data.frame(matrix(ncol = ncol(x = resultsreadable), nrow = (length(x = pairs)/2)*mult))
  colnames(x = averagesAverages) <- colnames(x = resultsreadable)
  i <- 1
  for (clust in 1:nrow(x = averagesAverages)) {
    averagesAverages[clust, ] <- apply(resultsreadable[i:(i + (num_doubs - 1)), ], 2, mean)
    if (clust %in% 1:(length(pairs)/2)) {
      row.names(x = averagesAverages)[clust] <- paste0(pairs[1, clust], '-', pairs[2, clust], '-even')
    } else if (clust %in% ((length(x = pairs)/2) + 1):((length(x = pairs)/2)*2)) {
      row.names(x = averagesAverages)[clust] <- paste0(pairs[1, clust - (length(x = pairs)/2)], '-', pairs[2, clust - (length(x = pairs)/2)], '-one')
    } else {
      row.names(x = averagesAverages)[clust] <- paste0(pairs[1, clust - length(x = pairs)], '-', pairs[2, clust - length(x = pairs)], '-two')
    }
    i <- i + num_doubs
  }
  return(list(averagesAverages = averagesAverages, doubletCellsInput2 = doubletCellsInput2))
}

