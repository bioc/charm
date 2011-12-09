\name{regionPlot}
\alias{regionPlot}
\title{
Plot user-supplied genomic regions using data returned by the dmrFinder function.
}
\description{
Plot any given genomic regions from tiling microarray data using data returned by the dmrFinder function.
}
\usage{
regionPlot(tab, dmr, cpg.islands, Genome, outfile="./regions.pdf", which.plot=1:10, plot.these, cl, legend.size=1, buffer=3000, plot.p=TRUE, plotG=FALSE, dat=NULL, grs=NULL)
}
\arguments{
  \item{tab}{
a data frame with columns chr, start, and end identifying the regions to be plotted from the data.
}
  \item{dmr}{
a list object as returned by dmrFinder, providing the data to be plotted.
}
  \item{cpg.islands}{
a table with columns "chr","start", and "end" for CpG islands to plot in the second panel.
}
  \item{Genome}{
the BSgenome object for the organism based upon which your array was designed.
}
  \item{outfile}{
a character string giving the name of the pdf file that will be saved.  Include the full path if file is not to be saved in the current working directory.
}
  \item{which.plot}{
a vector of indices identifying which regions (rows) from tab to plot.
}
  \item{plot.these}{
if dmr results are from unpaired analysis, specify the groups to plot.  Default is colnames(dmr$gm).  If dmr results are from paired analysis, specify the comparisons to plot.  Default is colnames(dmr$sMD).
}
  \item{cl}{
a vector of line and point colors, one for each element of plot.these in alphabetical order by group name.
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
  \item{plotG}{
if TRUE, a third panel of each DMR plot will show the difference between the median green channel value (after subtracting probe medians and correcting for gc content) between the 2 groups. If true, dat must be provided.
}
  \item{dat}{
if plotG=TRUE, this must be provided.  It is the TilingFeatureSet object used when estimating the matrix of percent methylation estimates used in dmrFinder when it produced dmr.
}
  \item{grs}{
if plotG=TRUE, this must be provided.  It is the names of the two groups whose difference in G should be plotted in the 3rd panel.
}

}
\details{
This function enables plotting of any regions, not just DMRs, using the results of dmrFinder.  The second panel shows the location of CpG's with ticks on the bottom (islands are colored) and the location of mcrbc recognition sites with ticks on the top.
}


\author{
Martin Aryee <aryee@jhu.edu>, Peter Murakami, Rafael Irizarry
}
\seealso{
\code{\link{plotRegions}}, \code{\link{dmrPlot}}, \code{\link{dmrFinder}}, \code{\link{dmrFdr}}
}
\examples{
# See dmrFdr
}
