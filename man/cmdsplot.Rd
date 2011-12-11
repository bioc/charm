\name{cmdsplot}
\alias{cmdsplot}
\title{
Classical multi-dimensional scaling plot of charm data.
}
\description{
Unsupervised clustering of charm data samples by classical multi-dimensional scaling.
}
\usage{
cmdsplot(labcols, expcol, rawData, p, okqc=1:nrow(p), noXorY=TRUE, outfile="./cmds_topN.pdf", topN=c(100000,1000)) 
}
\arguments{
  \item{labcols}{
vector of colors, one for each group being plotted (each unique value of pData(rawData)[,expcol])
}
  \item{expcol}{
name of the column of pData(rawData) the values of which the points are colored differently by.
}
  \item{rawData}{
the TilingFeatureSet object that p came from.
}
  \item{p}{
the matrix of percentage methylation values (scale: 0, 1) from using methp on rawData. One column per sample
}
  \item{okqc}{
vector of indices identifying which rows of p to use.
}
  \item{noXorY}{
logical value for whether or not (TRUE or FALSE) to ignore probes in chr X and Y?
}
   \item{outfile}{
file argument to pdf()
}
  \item{topN}{
Use only the  
}
}


\author{
Martin Aryee <aryee@jhu.edu>, Peter Murakami, Rafael Irizarry
}
\seealso{
\code{\link{regionPlot}}, \code{\link{dmrFinder}}, \code{\link{dmrFdr}}
}
\examples{
# See dmrFdr
}
