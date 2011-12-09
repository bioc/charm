\name{dmrPlot}
\alias{dmrPlot}
\title{
Plot differentially methylated regions (DMRs) found using the dmrFinder function.
}
\description{
Plot differentially methylated regions (DMRs) from tiling microarray data that were identified using the dmrFinder function.  To plot DMRs identified using the dmrFind function, use the plotDMRs function.
}
\usage{
dmrPlot(dmr, which.table=1:length(dmr$tabs), which.plot=1:30, legend.size=1, all.lines=TRUE, all.points=FALSE, colors.l, colors.p, outpath=".", cpg.islands, Genome, plotG=FALSE, dat=NULL, buffer=NULL, plot.p=TRUE)
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
  \item{cpg.islands}{
a table with columns "chr","start", and "end" for CpG islands to plot in the second panel.
}
  \item{Genome}{
the BSgenome object for the organism based upon which your array was designed.
}
  \item{plotG}{
if TRUE, a third panel of each DMR plot will show the difference between the median green channel value (after subtracting probe medians and correcting for gc content) between the 2 groups. If true, dat must be provided.
}
  \item{dat}{
if plotG=TRUE, this must be provided.  It is the TilingFeatureSet object used when estimating the matrix of percent methylation estimates used in dmrFinder when it produced dmr.
}
  \item{buffer}{
Number of base pairs to plot on either side of each DMR candidate.
}
  \item{plot.p}{
set to FALSE if you want to plot the methlation values (the "l" output from dmrFinder) instead of the percentage methylation values (the "p" output).  If dmrFinder was run on l instead of p, plot.p=FALSE necessarily.
}

}
\details{
This function plots the differentially methylated regions (DMRs).  The second panel shows the location of CpG's with ticks on the bottom (islands are colored) and the location of mcrbc recognition sites with ticks on the top.
}


\author{
Martin Aryee <aryee@jhu.edu>, Peter Murakami, Rafael Irizarry
}
\seealso{
\code{\link{plotDMRs}}, \code{\link{regionPlot}}, \code{\link{dmrFinder}}, \code{\link{dmrFdr}}
}
\examples{
# See dmrFdr
}
