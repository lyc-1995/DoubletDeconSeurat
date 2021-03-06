% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PrepDoubletDecon.R
\name{PrepDoubletDecon}
\alias{PrepDoubletDecon}
\title{Seurat V3 Pre Process}
\usage{
PrepDoubletDecon(
  object,
  assay = NULL,
  markers = NULL,
  top.n = 50,
  group.by = NULL,
  output.dir = NULL
)
}
\arguments{
\item{object}{Seurat object which has been pre-processed with clustering and dimensional reduction}

\item{assay}{Name of assay to pull the count data from; default is 'RNA'.}

\item{markers}{Output of FindAllMarkers. If NULL, will recalculate it with default parameters.}

\item{top.n}{Number of markers for the top_n function. Default is 50.}

\item{group.by}{Name of one metadata columns to group cells by; If NULL, will use current identity.}

\item{output.dir}{If set, save the output files to .txt format. Defauly is NULL.}
}
\value{
a list with three elements:
\code{newExpressionFile} - Seurat expression file in ICGS format (ICGS genes)
\code{newFullExpressionFile} - Seurat expression file in ICGS format (all genes)
\code{newGroupsFile} - Groups file ICGS format
}
\description{
Pre-process Seurat V3 object for DoubletDecon
}
\references{
DePasquale, E. A., Schnell, D. J., Van Camp, P. J., Valiente-Alandí, Í., Blaxall, B. C., Grimes, H. L., ... & Salomonis, N. (2019). DoubletDecon: Deconvoluting Doublets from Single-Cell RNA-Sequencing Data. Cell reports, 29(6), 1718-1727.
}
\keyword{Seurat}
