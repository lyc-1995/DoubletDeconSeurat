% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SyntheticDoublets.R
\name{SyntheticDoublets}
\alias{SyntheticDoublets}
\title{Synthetic Doublets}
\usage{
SyntheticDoublets(
  data,
  groups,
  groupsMedoids,
  newMedoids,
  num_doubs,
  only50,
  location,
  seed
)
}
\arguments{
\item{data}{Processed data from Clean_Up_Input (or Remove_Cell_Cycle).}

\item{groups}{Processed groups file from Clean_Up_Input.}

\item{groupsMedoids}{New groups file based on blacklisted clusters from Blacklist_Groups}

\item{newMedoids}{New combined medoids from Blacklist_groups.}

\item{num_doubs}{The user defined number of doublets to make for each pair of clusters}

\item{only50}{Use only synthetic doublets created with 50\%/50\% mix of parent cells, as opposed to the extended option of 30\%/70\% and 70\%/30\%, default is TRUE.}

\item{location}{Directory where output should be stored}

\item{seed}{Set a random seed.}
}
\value{
averagesAverages = average deconvolution profiles for each combination of cell types.

doubletCellsInput2 = deconvolution profiles for synthetic doublet cells (for quality control).
}
\description{
This function creates synthetic doublets by averaging gene expression profiles from each combination of clusters to generate deconvolution profiles for each type of doublet.
}
\keyword{synthetic}
