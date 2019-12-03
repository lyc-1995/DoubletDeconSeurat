#' Is A Doublet
#'
#' This function uses deconvolution analysis (DeconRNASeq) to evaluate each cell for equal contribution from blacklisted clusters.
#'
#' @param data Processed data from CleanUpInput (or RemoveCellCycle).
#' @param newMedoids New combined medoids from BlacklistGroups.
#' @param groups Processed groups file from CleanUpInput.
#' @param synthProfiles Average profiles of synthetic doublets from SyntheticDoublets.
#' @param log_file_name used for saving run notes to log file
#'
#' @return isADoublet - data.frame with each cell as a row and whether it is called a doublet by deconvolution analysis.
#' @return resultsreadable - data.frame with results of deconvolution analysis (cell by cluster) in percentages.
#'
#' @keywords doublet deconvolution decon
#'
#' @export
#'
IsDoublet <- function(
  data,
  newMedoids,
  groups,
  synthProfiles,
  log_file_name
) {
  #create data frame to store doublets table
  isADoublet <- data.frame(matrix(ncol = 4, nrow = (ncol(x = data) - 1)))
  rownames(x = isADoublet) <- colnames(x = data)[2:ncol(x = data)]
  rownames(x = newMedoids) <- rownames(x = data)[2:nrow(x = data)]

  #run DeconRNASeq with new medoids and data
  results <- DeconRNASeq(data[2:nrow(x = data), 2:ncol(x = data)], newMedoids)
  resultsreadable <- round(results$out.all*100, 2)
  rownames(x = resultsreadable) <- rownames(x = isADoublet) #make an easily readable results table

  #get average profiles for cell clusters
  averagesReal <- as.data.frame(matrix(ncol = ncol(x = resultsreadable), nrow = length(x = unique(groups[, 2]))))
  colnames(x = averagesReal) <- colnames(x = resultsreadable)
  for (clust in 1:length(x = unique(groups[, 2]))) {
    cells <- row.names(x = subset(groups, groups[, 1] == clust))
    subsetResults <- resultsreadable[row.names(x = resultsreadable) %in% cells,]
    averagesReal[clust, ] <- apply(subsetResults, 2, mean)
  }

  #create a table with average profiles of cell clusters and synthetic combinations
  allProfiles <- rbind(averagesReal, synthProfiles)

  #this section determines the profile with the highest correlation to the given cell and determines if it is one of the doublet profiles
  for (cell in 1:nrow(x = isADoublet)) {
    if (ncol(x = resultsreadable) == 2) { #If there are only 2 groups, correlation won't work, so I use minimum euclidean distance instead
      a <- rbind(allProfiles, resultsreadable[cell, ])
      b <- as.matrix(dist(a))
      c <- b[nrow(x = b), 1:(ncol(x = b) - 1)]
      chosenCorrelation <- c[c %in% min(x = c)]
      isADoublet[cell, 1] <- 100 - chosenCorrelation #100-euclidean distance
      isADoublet[cell, 2] <- names(chosenCorrelation)
      if (names(chosenCorrelation) %in% unique(groups[, 2])) { #it is an original cluster
        isADoublet[cell, 3] <- FALSE
      } else {
        isADoublet[cell, 3] <- TRUE
      }
    } else {
      #correlations=apply(allProfiles, 1, cor, resultsreadable[cell,])
      correlations <- apply(allProfiles, 1, cor, resultsreadable[cell, ])
      sortCorrelations <- sort(correlations, decreasing = TRUE)[1:2]
      maxCorrelation1 <- which(correlations == sortCorrelations[1])
      maxCorrelation2 <- which(correlations == sortCorrelations[2])
      chosenCorrelation <- maxCorrelation1
      isADoublet[cell, 1] <- correlations[chosenCorrelation]
      correlatedCluster <- row.names(x = allProfiles)[chosenCorrelation]
      isADoublet[cell, 2] <- correlatedCluster
      if (chosenCorrelation > length(x = unique(groups[, 2]))) {
        isADoublet[cell, 3] <- TRUE
      } else {
        isADoublet[cell, 3] <- FALSE
      }
    }
  }
  isADoublet[, 4] <- groups[, 2]
  colnames(x = isADoublet) <- c('Distance','Cell_Types', 'isADoublet', 'Group_Cluster')

  message(paste0(length(which(isADoublet$isADoublet == TRUE)), '/', nrow(x = isADoublet),  ' possible doublets removed'))
  cat(paste0(length(which(isADoublet$isADoublet == TRUE)), '/', nrow(x = isADoublet),  ' possible doublets removed'), file = log_file_name, append = TRUE, sep = '\n')
  return(list(isADoublet = isADoublet, resultsreadable = resultsreadable))
}
