#' Main DoubletDecon v1.0.1
#'
#' This is the main function. This function identifies clusters of doublets with a combination of deconvolution analysis and unique gene expression and individual doublet cells with deconvolution analysis.
#'
#' @param rawDataFile Name of file containing ICGS or Seurat expression data (gene by cell)
#' @param groupsFile Name of file containing group assignments (3 column: cell, group(numeric), group(numeric or character))
#' @param filename Unique filename to be incorporated into the names of outputs from the functions.
#' @param location Directory where output should be stored
#' @param fullDataFile Name of file containing full expression data (gene by cell). Default is NULL.
#' @param removeCC Remove cell cycle gene cluster by KEGG enrichment. Default is FALSE.
#' @param species Species as scientific species name, KEGG ID, three letter	species abbreviation, or NCBI ID. Default is "mmu".
#' @param rhop x in mean+x*SD to determine upper cutoff for correlation in the blacklist. Default is 1.
#' @param write Write output files as .txt files. Default is TRUE.
#' @param PMF Use step 2 (unique gene expression) in doublet determination criteria. Default is TRUE.
#' @param useFull Use full gene list for PMF analysis. Requires fullDataFile. Default is FALSE.
#' @param heatmap Boolean value for whether to generate heatmaps. Default is TRUE. Can be slow to datasets larger than ~3000 cells.
#' @param centroids Use centroids as references in deconvolution instead of the default medoids.
#' @param num_doubs The user defined number of doublets to make for each pair of clusters. Default is 100.
#' @param only50 use only synthetic doublets created with 50\%/50\% mix of parent cells, as opposed to the extended option of 30\%/70\% and 70\%/30\%, default is FALSE.
#' @param min_uniq minimum number of unique genes required for a cluster to be rescued
#' @param seed Set a random seed. If set NULL, will not use random seed.
#'
#' @return data_processed = new expression file (cleaned).
#' @return groups_processed = new groups file (cleaned).
#' @return PMF_results = pseudo marker finder t-test results (gene by cluster).
#' @return DRS_doublet_table = each cell and whether it is called a doublet by deconvolution analysis.
#' @return DRS_results = results of deconvolution analysis (cell by cluster) in percentages.
#' @return Decon_called_freq = percentage of doublets called in each cluster by deconvolution analysis.
#' @return Final_doublets_groups = new groups file containing only doublets.
#' @return Final_nondoublets_groups = new groups file containing only non doublets.
#'
#' @keywords doublet decon main
#'
#' @export
#'
DoubletDecon <- function(
  rawDataFile,
  groupsFile,
  filename,
  location,
  fullDataFile = NULL,
  removeCC = FALSE,
  species = "mmu",
  rhop = 1,
  write = TRUE,
  PMF = TRUE,
  useFull = FALSE,
  heatmap = TRUE,
  centroids = FALSE,
  num_doubs = 100,
  only50 = FALSE,
  min_uniq = 4,
  seed = NULL
){

  #load required packages
  message("Loading packages")
  suppressMessages(require(DeconRNASeq))
  suppressMessages(require(gplots))
  suppressMessages(require(plyr))
  suppressMessages(require(MCL))
  suppressMessages(require(clusterProfiler))
  suppressMessages(require(mygene))
  suppressMessages(require(tidyr))
  suppressMessages(require(R.utils)) #for new PMF
  suppressMessages(require(dplyr)) #for new PMF
  suppressMessages(require(foreach)) #for new PMF
  suppressMessages(require(doParallel)) #for new PMF
  suppressMessages(require(stringr)) #for new PMF

  #Check variables
  if(is.character(rawDataFile) != TRUE & is.data.frame(rawDataFile) != TRUE) {stop('ERROR: rawDataFile must be a character string!')}
  if(is.character(groupsFile) != TRUE & is.data.frame(groupsFile) != TRUE & is.matrix(groupsFile) != TRUE) {stop('ERROR: groupsFile must be a character string!')}
  if(is.character(filename) != TRUE) {stop('ERROR: filename must be a character string!')}
  if(is.character(location) != TRUE) {stop('ERROR: location must be a character string!')}
  if(is.character(fullDataFile) != TRUE & is.null(fullDataFile) != TRUE & is.data.frame(fullDataFile) != TRUE) {stop('ERROR: fullDataFile must be a character string or NULL!')}
  if(is.logical(removeCC) != TRUE) {stop('ERROR: removeCC must be TRUE or FALSE!')}
  if(is.character(species) != TRUE) {stop('ERROR: species must be a character string!')}
  if(is.numeric(rhop) != TRUE) {stop('ERROR: rhop must be numeric!')}
  if(is.logical(write) != TRUE) {stop('ERROR: write must be TRUE or FALSE!')}
  if(is.logical(PMF) != TRUE) {stop('ERROR: PMF must be TRUE or FALSE!')}
  if(is.logical(useFull) != TRUE) {stop('ERROR: useFull must be TRUE or FALSE!')}
  if(is.logical(heatmap) != TRUE) {stop('ERROR: heatmap must be TRUE or FALSE!')}
  if(is.logical(centroids) != TRUE) {stop('ERROR: centroids must be TRUE or FALSE!')}
  if(is.numeric(num_doubs) != TRUE) {stop('ERROR: numdoubs must be numeric!')}
  if(is.logical(only50) != TRUE) {stop('ERROR: only50 must be TRUE or FALSE!')}
  if(is.numeric(min_uniq) != TRUE) {stop('ERROR: min_uniq must be numeric!')}

  #Read in data
  message('Reading data')

  ICGS2_flag <- FALSE #set for checking if the input file is in ICGS2 format

  if (class(rawDataFile) == 'character') {
    #NEW: test for ICGS2
    rawDataHeader <- read.table(rawDataFile, sep = '\t', header = FALSE, row.names = 1, nrows = 1, stringsAsFactors = FALSE)
    if (length(x = grep(':', rawDataHeader[2])) == 1) {
      ICGS2_flag <- TRUE
      ICGS2 <- ICGS2toICGS1(rawDataFile = rawDataFile, groupsFile = groupsFile)
      rawData <- ICGS2$rawData
    } else {
      rawData <- read.table(rawDataFile, sep = '\t', header = T, row.names = 1)
    }
  } else {
    warning("if using ICGS2 file input, please import 'rawDataFile' and 'groupsFile' as path/location instead of an R object.")
    rawData <- rawDataFile
  }

  if (class(groupsFile) == 'character') {
    if (ICGS2_flag == TRUE) {
      groups <- ICGS2$groups
    } else {
      groups <- read.table(groupsFile, sep = '\t', header = FALSE, row.names = 1)
    }
  } else {
    groups <- groupsFile
  }

  #Clean up data and groups file
  message('Processing raw data')
  data <- CleanUpInput(rawData = rawData, groups = groups)
  og_processed_data <- data$processed
  groups <- data$groups

  #Centroids or medoids?
  if (centroids == TRUE) {
    centroid_flag <- TRUE
  } else {
    centroid_flag <- FALSE
  }

  #Original data heatmap
  if (heatmap == TRUE) {
    message('Creating original data heatmap')
    breaks <- seq(0, #start point of color key
                  as.numeric(x = quantile(data.matrix(data$processed[2:nrow(x = data$processed), 2:ncol(x = data$processed)]), 0.99)),  #end point of color key
                  by = 0.05) #length of sub-division
    mycol <- colorpanel(n = length(x = breaks) - 1, low = "black", high = "yellow") #heatmap colors
    suppressWarnings(DDheatmap(data.matrix(data$processed[2:nrow(x = data$processed), 2:ncol(x = data$processed)]), #the data matrix
                               Colv = FALSE, # No clustering of columns
                               Rowv = FALSE, #no clustering of rows
                               dendrogram = "none", #do not generate dendrogram
                               col = mycol, #colors used in heatmap
                               ColSideColors = as.color(Renumber(data$processed[1,2:ncol(x = data$processed)]), alpha = 1, seed = 4), #column color bar
                               RowSideColors = as.color(Renumber(data$processed[2:nrow(x = data$processed), 1]), alpha = 1, seed = 2), # row color bar
                               breaks = breaks, #color key details
                               trace = "none", #no trace on map
                               na.rm = TRUE, #ignore missing values
                               margins = c(5, 5), # size and layout of heatmap window
                               labRow = NA, #turn off gene labels
                               labCol = NA, #turn off cell labels
                               xlab = "Samples", #x axis title
                               ylab =  "Genes", # y axis title
                               main = paste0("Original data: ", filename))) #main title
  }

  #Remove cell cycle gene cluster (optional)
  if (removeCC == TRUE) {
    message('Removing cell cycle clusters')
    data <- RemoveCellCycle(data = data$processed, species = species)
  } else {
    data <- data$processed
  }
  if (write == TRUE) {
    write.table(data, file = file.path(location, paste0('data_processed_', filename, '.txt')), sep = '\t')
    write.table(groups, file = file.path(location, paste0('groups_processed_', filename, '.txt')), sep = '\t')
  }

  #Calculate medoids, medoid correlations, blacklist to create new combine medoids
  message('Combining similar clusters')
  BL <- BlacklistGroups(data = data, groups = groups, rhop = rhop, centroid_flag = centroid_flag)
  newMedoids <- BL$newMedoids
  groupsMedoids <- BL$newGroups

  #Create synthetic doublets to get average synthetic profiles
  message('Creating synthetic doublet profiles')
  if (.Platform$OS.type == "unix") {
    sink("/dev/null") #hides DeconRNASeq output
    synthProfilesx <- SyntheticDoublets(data = data, groups = groups, groupsMedoids = groupsMedoids, newMedoids = newMedoids, num_doubs = num_doubs, only50 = only50, location = location, seed = seed)
    sink()
  } else {
    synthProfilesx <- SyntheticDoublets(data = data, groups = groups, groupsMedoids = groupsMedoids, newMedoids = newMedoids, num_doubs = num_doubs, only50 = only50, location = location, seed = seed)
  }
  synthProfiles <- synthProfilesx$averagesAverages
  doubletCellsInput2 <- synthProfilesx$doubletCellsInput2
  if (write == TRUE) {
    write.table(doubletCellsInput2, file = file.path(location, paste0('Synth_doublet_info_', filename, '.txt')), sep = '\t')
  }

  #Calculate doublets using DeconRNASeq
  message('Step 1: Removing possible doublets')
  if (.Platform$OS.type == "unix") {
    sink("/dev/null") #hides DeconRNASeq output
    doubletTable <- IsDoublet(data = data, newMedoids = newMedoids, groups = groups, synthProfiles = synthProfiles)
    sink()
  } else {
    doubletTable <- IsDoublet(data = data, newMedoids = newMedoids, groups = groups, synthProfiles = synthProfiles)
  }
  if (write == TRUE) {
    write.table(doubletTable$isADoublet, file = file.path(location, paste0('DRS_doublet_table_', filename, '.txt')), sep = '\t')
    write.table(doubletTable$resultsreadable, file = file.path(location, paste0('DRS_results_', filename, '.txt')), sep = "\t")
  }

  #Recluster doublets and non-doublets
  message('Step 2: Re-clustering possible doublets')
  reclusteredData <- Recluster(isADoublet = doubletTable$isADoublet, data = data, groups = groups)
  data <- reclusteredData$newData2$processed
  groups <- reclusteredData$newData2$groups
  write.table(data, file = file.path(location, paste0('data_processed_reclust_', filename, '.txt')), sep = '\t', col.names = NA, quote = FALSE)
  write.table(groups, file = file.path(location, paste0('groups_processed_reclust_', filename, '.txt')), sep = '\t')

  #Run Pseudo Marker Finder to identify clusters with no unique gene expression
  if (PMF == FALSE) {
    message('SKIPPING Step 3: Rescuing cells with unique gene expression')
    PMFresults <- NULL
  } else {
    message('Step 3: Rescuing cells with unique gene expression')
    if (useFull == TRUE) {
      PMFresults <- PseudoMarkerFinder(groups = as.data.frame(groups), redu_data2 = file.path(location, paste0('data_processed_reclust_', filename, '.txt')), full_data2 = fullDataFile, min_uniq = min_uniq)
    } else {
      PMFresults <- PseudoMarkerFinder(groups = as.data.frame(groups), redu_data2 = file.path(location, paste0('data_processed_reclust_', filename, '.txt')), full_data2 = NULL, min_uniq = min_uniq)
    }
    if (write == TRUE) {
      write.table(PMFresults, file = file.path(location, paste0('new_PMF_results_', filename, '.txt')), sep = '\t')
    }
  }
  #Doublet Detection method 2: Pseudo_Marker_Finder
  allClusters <- unique(groups[,1])
  if (PMF == FALSE) {
    newDoubletClusters <- allClusters
  } else {
    hallmarkClusters <- as.numeric(x = unique(PMFresults[, 2]))
    newDoubletClusters <- setdiff(allClusters, hallmarkClusters)
  }
  #Doublet Detection method 1: Is_A_Doublet
  uniqueClusters <- as.character(x = unique(groups[, 2]))
  DeconCalledFreq <- as.data.frame(matrix(nrow = length(x = allClusters), ncol = 1), row.names = uniqueClusters)
  for (clus in 1:length(x = allClusters)) {
    #modified this line, was originally "clus in allClusters"
    temp1 <- subset(doubletTable$isADoublet, Group_Cluster == uniqueClusters[clus])
    if (nrow(temp1) == 0) {
      #not an original cluster, only a new doublet cluster
      DeconCalledFreq[clus, 1] <- 100
    } else {
      DeconCalledFreq[clus, 1] <- (length(x = which(temp1$isADoublet == TRUE))/nrow(x = temp1))*100
    }
  }

  #Combine to find real doublets
  if (PMF == FALSE) {
    finalDoublets <- row.names(x = doubletTable$isADoublet)[which(doubletTable$isADoublet$isADoublet == TRUE)] #this gives you the names of cells called as doublets by deconvolution
  } else {
    finalDoublets <- intersect(row.names(x = doubletTable$isADoublet)[which(doubletTable$isADoublet$isADoublet == TRUE)], row.names(x = subset(groups, groups[, 1] %in% newDoubletClusters)))
  }
  #Results
  finalDoubletCellCall <- groups[row.names(x = groups) %in% finalDoublets, ]
  finalNotDoubletCellCall <- groups[!(row.names(x = groups) %in% finalDoublets), ]
  if (write == TRUE) {
    write.table(finalDoubletCellCall, file = file.path(location, paste0('Final_doublets_groups_', filename, '.txt')), sep = '\t')
    write.table(finalNotDoubletCellCall, file = file.path(location, paste0('Final_nondoublets_groups_', filename, '.txt')), sep = '\t')
  }

  #Subset expression matrix for doublets and save
  doublets_matrix <- cbind(og_processed_data[, 1], og_processed_data[, which(colnames(x = og_processed_data) %in% row.names(x = finalDoubletCellCall))])
  if (write == TRUE) {
    write.table(doublets_matrix, file = file.path(location, paste0('Final_doublets_exp_', filename, '.txt')), sep = '\t')
  }

  #Heatmap of cells removed as doubets
  if (heatmap == TRUE) {
    message('Creating doublets heatmap')
    breaks <- seq(0, #start point of color key
                  as.numeric(quantile(data.matrix(doublets_matrix[2:nrow(x = doublets_matrix), 2:ncol(x = doublets_matrix)]), 0.99)),  #end point of color key
                  by = 0.05) #length of sub-division
    mycol <- colorpanel(n = length(x = breaks) - 1, low = 'black', high = 'yellow') #heatmap colors
    suppressWarnings(DDheatmap(data.matrix(doublets_matrix[2:nrow(x = doublets_matrix), 2:ncol(x = doublets_matrix)]), #the data matrix
                               Colv = FALSE, # No clustering of columns
                               Rowv = FALSE, #no clustering of rows
                               col = mycol, #colors used in heatmap
                               dendrogram = 'none', #turn of dendrogram generation
                               ColSideColors = as.color(Renumber(doublets_matrix[1, 2:ncol(x = doublets_matrix)]), alpha = 1, seed = 4), #column color bar
                               RowSideColors = as.color(Renumber(doublets_matrix[2:nrow(x = doublets_matrix),1]), alpha = 1, seed = 2), # row color bar
                               breaks = breaks, #color key details
                               trace = 'none', #no trace on map
                               na.rm = TRUE, #ignore missing values
                               margins = c(5,5), # size and layout of heatmap window
                               labRow = NA, #turn off gene labels
                               labCol = NA, #turn off cell labels
                               xlab = 'Samples', #x axis title
                               ylab =  'Genes', # y axis title
                               main = paste0('Doublets: ', filename))) #main title)
  }
  #Subset expression matrix for non-doublets and save
  nondoublets_matrix <- cbind(og_processed_data[, 1], og_processed_data[, which(colnames(x = og_processed_data) %in% row.names(x = finalNotDoubletCellCall))])
  if (write == TRUE) {
    write.table(nondoublets_matrix, file = file.path(location, paste0('Final_nondoublets_exp_', filename, '.txt')), sep = '\t')
  }
  #New heatmap of non-doublet cells
  if (heatmap == TRUE) {
    message('Creating non-doublets heatmap')
    breaks <- seq(0, #start point of color key
                  as.numeric(quantile(data.matrix(nondoublets_matrix[2:nrow(x = nondoublets_matrix), 2:ncol(x = nondoublets_matrix)]), 0.99)),  #end point of color key
                  by = 0.05) #length of sub-division
    mycol <- colorpanel(n = length(x = breaks) - 1, low = 'black', high = 'yellow') #heatmap colors
    suppressWarnings(DDheatmap(data.matrix(nondoublets_matrix[2:nrow(x = nondoublets_matrix), 2:ncol(x = nondoublets_matrix)]), #the data matrix
                               Colv = FALSE, # No clustering of columns
                               Rowv = FALSE, #no clustering of rows
                               col = mycol, #colors used in heatmap
                               dendrogram = 'none', #turn of dendrogram generation
                               ColSideColors = as.color(Renumber(x = nondoublets_matrix[1, 2:ncol(x = nondoublets_matrix)]), alpha = 1, seed = 4), #column color bar
                               RowSideColors = as.color(Renumber(x = nondoublets_matrix[2:nrow(x = nondoublets_matrix),1]), alpha = 1, seed = 2), # row color bar
                               breaks = breaks, #color key details
                               trace = 'none', #no trace on map
                               na.rm = TRUE, #ignore missing values
                               margins = c(5, 5), # size and layout of heatmap window
                               labRow = NA, #turn off gene labels
                               labCol = NA, #turn off cell labels
                               xlab = 'Samples', #x axis title
                               ylab =  'Genes', # y axis title
                               main = paste0('Non-Doublets: ', filename))) #main title
  }
  #last message
  message('Finished!')
  return(list(data_processed = data,
              groups_processed = groups,
              DRS_doublet_table = doubletTable$isADoublet,
              DRS_results = doubletTable$resultsreadable,
              PMF_results = PMFresults,
              Final_doublets_groups = finalDoubletCellCall,
              Final_nondoublets_groups = finalNotDoubletCellCall,
              Final_doublets_exp = doublets_matrix,
              Final_nondoublets_exp = nondoublets_matrix,
              Synth_doublet_info = doubletCellsInput2))
}

