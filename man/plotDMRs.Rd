\name{plotDMRs}
\alias{plotDMRs}
\title{
Plot differentially methylated regions (DMRs) found using the dmrFind function.
}
\description{
Plot differentially methylated regions (DMRs) from tiling microarray data that were identified using the dmrFind function.
}
\usage{
plotDMRs(dmrs, Genome, cpg.islands, exposure, outfile, which_plot=1:50, all.points=TRUE, plot.these=NULL, ADD=3000, cols=c("black","blue","red","gray","brown","pink","orange"), legend.size=1, smoo="loess", SPAN=300, DELTA=36, point.info=FALSE, pch.groups=NULL, panel3="pvalues", G=NULL, seq=NULL)
}
\arguments{
  \item{dmrs}{
a list object as returned by dmrFind.
}
  \item{Genome}{
the BSgenome object for the organism based upon which your array was designed.
}
  \item{cpg.islands}{
a table with columns "chr","start", and "end" for CpG islands to plot in the second panel.
}
  \item{exposure}{
The covariate of interest.
}
  \item{outfile}{
the name of the file to save (including the full path) 
}
  \item{which_plot}{
numeric vector of indices identifying which DMR candidates from dmrs$dmrs to plot.
}
  \item{all.points}{
if TRUE, plot the points for all groups.  If FALSE, only for the 2 groups being compared.
}
   \item{plot.these}{
if all.points=FALSE (and covariate is not continuous), then which groups do you want to plot?
}
  \item{SPAN}{
see DELTA. Only used if smoo="loess"
}
  \item{DELTA}{
span parameter in loess smoothing will = SPAN/(DELTA * number of probes in the plotted region).  Only used if smoo="loess".
}
  \item{smoo}{
"loess" for loess smoother or "runmed" for running median smoother (runmed with k=3). This does not  need to be the same as the smoo argument to dmrFind.
}
  \item{ADD}{
Number of base pairs to plot on either side of each DMR candidate (if it is covered on the array).
}
  \item{cols}{
vector of colors to use, one for each group (if covariate is categorical)
}
  \item{point.info}{
if TRUE, function returns a table identifying which sample is plotted with which number or letter (if pch.groups=NULL, the default).
}
  \item{legend.size}{
magnification factor for the legend
}
  \item{pch.groups}{
vector whose length is equal to the number of samples.  Each unique value will be plotted with a different point type.  
}
  \item{panel3}{
if panel3="G", the third panel of each DMR plot will show the difference between the median green channel value (after subtracting probe medians and correcting for gc content) between the 2 groups (i.e., the group defined by mod[,coef] in dmrFind minus the reference group).  If panel="G", seq argument must be provided.  If panel!="G", the 3rd panel will show -log10(dmrs$pval). G 
}
  \item{G}{
matrix of green channel intensities to use for plotting in the 3rd panel if panel3="G".
}
  \item{seq}{
vector of probe sequences corresponding to the rows of G (and dmrs$cleanp) if panel3="G".
}
}

\details{
This function plots the differentially methylated regions (DMRs) that were identified using the dmrFind function.
}


\author{
Martin Aryee <aryee@jhu.edu>, Peter Murakami, Rafael Irizarry
}
\seealso{
\code{\link{plotRegions}}, \code{\link{dmrFind}}, \code{\link{qval}}
}
\examples{
# See qval
}
