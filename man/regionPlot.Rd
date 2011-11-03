\name{regionPlot}
\alias{regionPlot}
\title{
Plot user-supplied genomic regions.
}
\description{
Plot any given genomic regions from tiling microarray data.
}
\usage{
regionPlot(tab, dmr, outfile, which.plot, which.groups=colnames(dmr$gm), cl=1+1:ncol(dmr$gm), legend.size=1, buffer=3000, plot.p=TRUE)
}
\arguments{
  \item{tab}{
a data frame with columns chr, start, and end identifying the regions to be plotted from the data.
}
  \item{dmr}{
a list object as returned by dmrFinder, providing the data to be plotted.
}
  \item{outfile}{
a character string giving the name of the pdf file that will be saved.  Include the full path if file is not to be saved in the current working directory.
}
  \item{which.plot}{
a vector of indices identifying which regions (rows) from tab to plot.
}
  \item{which.groups}{
a character vector of names (or a numeric vector of indices for the columns of dmr$gm)  identifying which groups to plot.
}
  \item{cl}{
a vector of line and point colors, one for each group in which.groups in alphabetical order by group name.
}
  \item{legend.size}{
cex argument for the legend (factor by which to magnify/shrink the legend).
}
  \item{buffer}{
An integer to control how many basepairs to show on either side of the plotted regions.
}
  \item{plot.p}{
set to FALSE if you want to plot the methlation values (the "l" output from dmrFinder) instead of the percentage methylation values (the "p" output).  If dmrFinder was run on l instead of p, plot.p=FALSE necessarily.
}
}
\details{
This function enables plotting of any regions, not just DMRs.
}


\author{
Martin Aryee <aryee@jhu.edu>, Peter Murakami, Rafael Irizarry
}
\seealso{
\code{\link{dmrPlot}}, \code{\link{dmrFinder}}, \code{\link{dmrFdr}}
}
\examples{
# See dmrFdr
}
