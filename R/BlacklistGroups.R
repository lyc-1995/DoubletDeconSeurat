#' Blacklist Groups
#'
#' This function is used for the calculation of blacklisted medoids and groups. It calculates medoids, performs medoid correlations and creates the binary correlation table between medoids, otherwise known as the blacklist. It makes a blacklist heatmap. It makes new medoids based on blacklisted combined clusters and a new groups file based on the blacklisted combined clusters.
#'
#' @param data  Processed data from CleanUpInput (or RemoveCellCycle)
#' @param groups Processed groups file from CleanUpInput
#' @param rhop x in mean+x*SD to determine upper cutoff for correlation in the blacklist. Default is 1.
#' @param centroid_flag use centroids if data is sparse
#' @param log_file_name used for saving run notes to log file
#'
#' @return newMedoids - new medoids data.frame for the new combined blacklisted clusters.
#' @return newGroups - new groups file containing cluster assignment based on new combined blacklisted clusters.
#'
#' @keywords blacklist correlation
#'
#' @export
#'
BlacklistGroups <- function(
  data,
  groups,
  rhop,
  centroid_flag,
  log_file_name
) {
  #Step 1: calculate medoids
  medoids <- data.frame(rep(NA, nrow(x = data) - 1))
  clusters <- length(x = unique(groups[, 1]))
  for (cluster in 1:clusters) {
    if (centroid_flag == TRUE) {
      medoids <- cbind(medoids, apply(data[2:nrow(x = data), which(data[1, ] == cluster), drop = FALSE], 1, mean))
    } else {
      medoids <- cbind(medoids, apply(data[2:nrow(x = data), which(data[1, ] == cluster), drop = FALSE], 1, median))
    }
  }
  medoids <- medoids[, -1]
  colnames(x = medoids) <- unique(groups[, 2])

  #Step 2: medoid correlations
  cormedoids <- cor(medoids, method = "pearson")

  #Step 3: create blacklist (binary correlation table between medoids)
  blacklist <- data.frame(matrix(ncol = ncol(x = medoids), nrow = ncol(x = medoids)))
  cutoff <- mean(cormedoids) + (rhop)*sd(cormedoids) #based on mean + 1SD * user provided multiplier
  for (rrow in 1:nrow(cormedoids)) {
    for (ccol in 1:ncol(cormedoids)) {
      if (cormedoids[rrow, ccol] > cutoff) {
        blacklist[rrow, ccol] <- 1
      } else {
        blacklist[rrow, ccol] <- 0
      }
    }
  }
  blacklist_original_order <- blacklist #keep original order with no names so I can pull the correct clusters out when trying to make medoids
  row.names(x = blacklist) <- row.names(x = cormedoids)
  colnames(x = blacklist) <- row.names(x = cormedoids)

  #Step 4: make blacklist heatmap
  BLheatmap <- heatmap.2(as.matrix(blacklist_original_order),
                         Colv = TRUE, # clustering of columns
                         Rowv = TRUE, # clustering of rows
                         xlab = 'Cell Types', #x axis title
                         ylab =  'Cell Types', #y axis title
                         trace = 'none',
                         main = 'Cluster Merge') #main title
  blacklist_original_order <- blacklist_original_order[BLheatmap$rowInd, BLheatmap$colInd]
  blacklist <- blacklist[BLheatmap$rowInd, BLheatmap$colInd]

  #Step 5: make new medoids
  blacklistCluster <- try(mcl(blacklist, addLoops = FALSE)$Cluster)
  if (class(blacklistCluster) == "try-error") {
    message("Unable to perform mcl function for blacklist clustering, please try a different rhop.")
    stop()
  }
  i <- -1 #if the cluster is 0 (meaning no combining) change the name of the cluster to make it unique
  for (cluster in 1:length(x = blacklistCluster)) {
    if (blacklistCluster[cluster] == 0) {
      blacklistCluster[cluster] <- i
      i <- i-1
    }
  }
  uniquelist <- unique(blacklistCluster) #list of unique clusters
  nunique <- length(x = uniquelist) #number of unique clusters
  newMedoids <- data.frame(matrix(ncol = nunique, nrow = (nrow(x = data) - 1)))
  for (cluster in 1:nunique) { #for each new cluster, assign medoid to new data.frame
    temp <- which(blacklistCluster == uniquelist[cluster])
    if (centroid_flag == TRUE) {
      newMedoids[, cluster] <- apply(data[2:nrow(x = data), (data[1, ] %in% rownames(x = blacklist_original_order)[temp]), drop = FALSE], 1, mean) #this is where the original order is critical
    } else {
      newMedoids[, cluster] <- apply(data[2:nrow(x = data), (data[1, ] %in% rownames(x = blacklist_original_order)[temp]), drop = FALSE], 1, median) #this is where the original order is critical
    }
    colnames(x = newMedoids)[cluster] <- paste(rownames(x = blacklist)[temp], collapse = '-') #need to give the columns meaningful names (combination names of combined clusters)
  }
  row.names(x = newMedoids) <- row.names(x = data)[2:nrow(x = data)]

  #Step 6: make new combined groups file
  newGroups <- data.frame(matrix(ncol = ncol(x = groups) + 1, nrow = 1))
  for (cluster in 1:nunique) {
    temp <- which(blacklistCluster == uniquelist[cluster])
    temp1.5 <- as.integer(row.names(x = blacklist_original_order)[temp])
    temp2 <- groups[groups[, 1] %in% temp1.5, , drop = FALSE]
    temp3 <- cbind(temp2, rep(paste(rownames(x = blacklist)[temp], collapse = '-'), nrow(x = temp2)))
    colnames(x = temp3) <- colnames(x = newGroups)
    newGroups <- rbind(newGroups, temp3)
  }
  newGroups <- newGroups[-1, ]
  newGroups[, 2] <- as.integer(x = as.factor(x = newGroups[, 3]))
  newGroups <- newGroups[, -1]

  message(paste0('New blacklisted clusters: ', paste(unique(newGroups[, 2]), sep = "' '", collapse = ", ")))
  cat(paste0('New blacklisted clusters: ', paste(unique(newGroups[, 2]), sep = "' '", collapse = ", ")), file = log_file_name, append = TRUE, sep = "\n")
  return(list(newMedoids = newMedoids, newGroups = newGroups))
}
