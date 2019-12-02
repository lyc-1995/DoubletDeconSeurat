#' Pseudo MarkerFinder
#'
#' This function uses ANOVA to look for unique gene expression in each possible doublet cluster.
#'
#' @param groups Processed groups file from CleanUpInput.
#' @param redu_data2 Processed data from CleanUpInput (or RemoveCellCycle) path, automatically written.
#' @param full_data2 cleaned full expression matrix from CleanUpInput.
#' @param min_uniq minimum number of unique genes required for a cluster to be rescued
#'
#' @return new_table - non-doublet clusters, as determined by the "Remove" and "Rescue" steps.
#'
#' @keywords Marker Finder ANOVA
#'
#' @export
#'
PseudoMarkerFinder <- function(
  groups,
  redu_data2,
  full_data2,
  min_uniq = 4
) {
  #Define the expression file
  if (!is.null(full_data2)) {
    rawFilePath <- full_data2
  } else {
    rawFilePath <- redu_data2
  }

  #Define possible doublet clusters
  doublet_pre <- cbind(unique(groups[grep('even|one|two', groups[, 2]), 1]), unique(groups[grep('even|one|two', groups[, 2]), 1]))

  #Non-doublet clusters vs doublet clusters
  doublet <- as.numeric(x = unique(doublet_pre[, 1]))
  nonDoublet <- as.numeric(x = unique(groups[which(!(as.numeric(x = groups[, 1]) %in% doublet)), 1]))

  #How many genes per time are read in (length block of input file)
  #The more cells (columns), the larger the memory for same number of genes (rows)
  genesPerTime <- 1000

  #Get the cell names (first row of file assumed)
  cellNames <- str_replace(scan(rawFilePath, character(), nlines = 1, sep = "\t", quiet = TRUE)[-1], '-', '\\.')
  reclust_groups <- groups
  reclust_groups$cells <- row.names(x = reclust_groups)
  colnames(x = reclust_groups) <- c('clusterNr', 'clusterName', 'cellName')

  #Make sure only to evaluate cells output from AltAnalyze
  processedCells <- cellNames %in% reclust_groups$cellName
  reclust_groups <- reclust_groups[match(cellNames[processedCells], reclust_groups$cellName),]

  #Read the total lines (e.g. genes) in the file (needed for parallel)
  nLines <- countLines(file = rawFilePath)[1]

  #Run in parallel the code for each chunk of data read from the file
  cl <- makeCluster(detectCores())
  registerDoParallel(cl = cl)

  anovaResult <- foreach(linePos = 0:floor(x = nLines/genesPerTime), .packages = c('dplyr'), .combine = 'rbind') %dopar% {

    #Read in a block from the input file (different blocks analysed in parallel) and pre-process
    startPos <- 1 + linePos*genesPerTime
    genesPerTime <- ifelse(linePos == floor(x = nLines/genesPerTime), nLines - startPos, genesPerTime)

    #Read in the chunk of data as a bix matrix
    myChunk <- scan(rawFilePath, character(), skip = startPos, nlines = genesPerTime, sep = '\t')
    nCells <- length(x = myChunk)/genesPerTime
    myChunk <- t(sapply(0:(genesPerTime - 1), function(x) {
      myChunk[(1 + nCells*x):(nCells*x + nCells)]
    }))

    #Get the genes we're testing (first value on each line)
    chunkGenes <- unlist(myChunk[, 1])

    #Convert the rest into numeric dataframe and add the groups (from altAnalyze)
    myChunk <- t(data.matrix(myChunk[, c(FALSE, processedCells)]))
    myChunk <- as.matrix(myChunk)
    class(myChunk) <- 'numeric'
    myChunk <- as.data.frame(myChunk)
    colnames(x = myChunk) <- paste('gene', 1:genesPerTime, sep = '')
    myChunk$group <- reclust_groups$clusterNr

    #ANOVA in bulk (main optimisation)
    # By precalculating sums and sum of squares for each group in each gene for the whole chunk at once,
    # we save a lot of time downstream
    # https://www.youtube.com/watch?v=ynx04Qgqdrc
    nGroups <- myChunk  %>% count(group) %>% group_by(group)
    nCells <- nrow(x = reclust_groups)
    allSums <- myChunk %>% group_by(group) %>% arrange(group) %>% summarise_all(funs(sum))
    allSums$theCount <- nGroups$n
    nGroups <- nrow(x = nGroups)
    myChunk[, -ncol(x = myChunk)] <- myChunk[, -ncol(x = myChunk)]^2
    allSquares <- myChunk %>% group_by(group) %>% summarise_all(funs(sum))

    #Needed to perform an optimised Tukey test on every doublets vs. all non-doublet
    invertCount <- 1/allSums$theCount
    q <- qtukey(.95, nGroups, df = (nCells - nGroups))

    # This is the fast ANOVA algorithm for each gene based on pre-calculated sum and sum of squares
    do.call(rbind, lapply(1:genesPerTime, function(x) {
      Cx <- sum(allSums[, x+1])^2/nCells
      SSt <- sum(allSquares[, x+1])-Cx
      SSa <- sum(unlist(allSums[, x+1])^2/allSums$theCount) - Cx
      SSw <- SSt - SSa
      MMSa <- SSa/(nGroups - 1)
      MSSw <- SSw/(nCells - nGroups)

      Fratio <- MMSa/MSSw

      #If the F ratio is significant, we can reject null hypothesis and gene has group with significant difference
      if (is.na(Fratio) | qf(.95, df1 = (nGroups - 1), df2 = (nCells - nGroups)) >= Fratio) {
        data.frame(gene = chunkGenes[x], doublet = -1, stringsAsFactors = FALSE) #ANOVA Not significant (code -1)
      } else {

        #Perform Tukey and mean comparisons
        myMeans <- unlist(allSums[, x + 1])/allSums$theCount
        signifDoublets <- doublet[sapply(doublet, function(myDoublet) {
          #Check Tukey tests
          all(abs(myMeans[myDoublet] - myMeans[nonDoublet]) >= q*sqrt(MSSw/2*(invertCount[myDoublet] + invertCount[nonDoublet]))) &&
            all(myMeans[myDoublet] >= myMeans[nonDoublet])     #doublet mean has to be larger than all non-doublets
        })]
        if (length(x = signifDoublets) == 0) {
          #Not all Tukey tests are significant or doublet mean is not higher than all (code 0)
          data.frame(gene = chunkGenes[x], doublet = 0, stringsAsFactors = FALSE)
        } else {
          #All are for these doublets
          data.frame(gene = chunkGenes[x], doublet = signifDoublets, stringsAsFactors = FALSE)
        }
      }
    }))
  }
  stopImplicitCluster()
  stopCluster(cl = cl)

  #Table the anova results
  anovaResultRedu <- anovaResult[anovaResult[, 2] %in% doublet, ]
  summedResults <- table(anovaResultRedu[, 2])
  # cat(paste0('Unique Genes By Cluster: ', summedResults), file = log_file_name, append = TRUE, sep = "\n")
  message(paste0('Unique Genes By Cluster: ', summedResults, ' ', '\n'))

  unique_rescued_clusters <- names(summedResults[summedResults >= min_uniq]) #min genes unique
  all_rescued_clusters <- c(as.character(x = unique_rescued_clusters), as.character(x = nonDoublet))
  new_table <- as.data.frame(matrix(ncol = 2, nrow = length(x = all_rescued_clusters)))
  new_table[, 2] <- all_rescued_clusters

  return(new_table)
}
