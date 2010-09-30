\name{dmrPlot}
\alias{dmrPlot}
\title{
Plot differentially methylated regions (DMRs)
}
\description{
Plot differentially methylated regions (DMRs) from tiling microarray data.
}
\usage{
dmrPlot(dmr, which.table=1:length(dmr$tabs), which.plot=1:30, legend.size=1, all.lines=TRUE, all.points=FALSE, colors.l, colors.p, outpath=".", plot.p=TRUE)
}
\arguments{
  \item{dmr}{
a list object as returned by dmrFinder.
}
  \item{which.table}{
a vector of indices identifying which tables in the dmr list to plot regions from.
}
  \item{which.plot}{
a vector of indices identifying which regions (rows) from each table to plot.
}
  \item{legend.size}{
cex argument for the legend (factor by which to magnify/shrink the legend).
}
  \item{all.lines}{
if TRUE, plot the smooth lines for all groups.  If FALSE, only for the 2 groups being compared.
}
  \item{all.points}{
if TRUE, plot the points for all groups.  If FALSE, only for the 2 groups being compared.
}
   \item{colors.l}{
a vector of line colors, one color for each group whose line is to be plotted (in alphabetical order).
}
  \item{colors.p}{
a vector of point colors, one color for each group whose points are to be plotted (in alphabetical order).
}
  \item{outpath}{
where to save the output pdf file.
}
  \item{plot.p}{
set to FALSE if you want to plot the methlation values (the "l" output from dmrFinder) instead of the percentage methylation values (the "p" output).  If dmrFinder was run on l instead of p, plot.p=FALSE necessarily.
}

}
\details{
This function plots the differentially methylated regions (DMRs).
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
