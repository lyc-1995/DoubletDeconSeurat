% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICGS2toICGS1.R
\name{ICGS2toICGS1}
\alias{ICGS2toICGS1}
\title{ICGS2_to_ICGS1}
\usage{
ICGS2toICGS1(rawDataFile, groupsFile, log_file_name = NULL)
}
\arguments{
\item{rawDataFile}{ICGS expression or counts file (in ICGS expression file format).}

\item{groupsFile}{ICGS groups file.}

\item{log_file_name}{used for saving run notes to log file}
}
\value{
processed - data.frame of genes by samples with a row of cell clusters (column_clusters-flat) and a column of gene clusters (row_clusters-flat) when available.

groups - groups file with cell names matching the expression file.
}
\description{
This function converts the new ICGS2 expression file format to the ICGS1 format. In making this change, the cluster numbers (column 1 in the groups file) are retained and not the cluster names (column 2 in the groups file).
}
\keyword{ICGS}
