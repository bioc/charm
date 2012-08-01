\name{plotRegions}
\alias{plotRegions}
\title{
Plot user-provided regions.
}
\description{
Plot user-provided regions.
}
\usage{
plotRegions(thetable, cleanp, chr, pos, seq=NULL, Genome, cpg.islands, exposure, exposure.continuous=FALSE, outfile="./myregions.pdf", which_lines=NULL, which_points=which_lines, ADD=3000, cols=c("black","red","blue","gray","green","orange","brown"), legend.size=1, smoo="loess", SPAN=300, DELTA=36, panel3="none", G=NULL, grs=NULL)
}
\arguments{
  \item{thetable}{
a table with columns for "chr","start", and "end", identifying the genomic regions to plot (if they are covered on the array).
}
  \item{cleanp}{
the matrix of percent methylation estimates to be used for plotting 
}
  \item{chr}{
vector of chromosome labels for the probes in cleanp
}
  \item{pos}{
vector of chromosomal coordinates for the probes in cleanp
}
  \item{seq}{
vector of probe sequences corresponding to the rows of G (and cleanp). Needed only if panel3="G".
}
  \item{Genome}{
the BSgenome object for the organism based upon which your array was designed.
}
  \item{cpg.islands}{
a table with columns "chr","start", and "end" for CpG islands to plot in the second panel.
}
  \item{outfile}{
the name of the file to save (including the full path) 
}
  \item{exposure}{
The covariate of interest.
}
  \item{exposure.continuous}{
set to TRUE if exposure is a continuous variable.
}
  \item{which_lines}{
vector specifying which groups (unique elements of exposure) to plot the lines for.  If NULL (the default), plots lines for all groups.  Only applies if exposure is categorical.
}
   \item{which_points}{
vector specifying which groups (unique elements of exposure) to plot the points for.  If NULL (the default), plots points for all groups.  Only applies if exposure is categorical.
}
  \item{ADD}{
Number of base pairs to plot on either side of each DMR candidate (if it is covered on the array).
}
  \item{cols}{
vector of colors to use, one for each group (if covariate is categorical)
}
  \item{legend.size}{
magnification factor for the legend
}
  \item{smoo}{
"loess" for loess smoother or "runmed" for running median smoother (runmed with k=3).
}
  \item{SPAN}{
see DELTA. Only used if smoo="loess"
}
  \item{DELTA}{
span parameter in loess smoothing will = SPAN/(DELTA * number of probes in the plotted region).  Only used if smoo="loess".
}
  \item{panel3}{
if panel3="G", the third panel of each DMR plot will show the difference between the median green channel value (after subtracting probe medians and correcting for gc content) between the 2 groups (i.e., the group defined by mod[,coef] in dmrFind minus the reference group).  If panel="G", seq argument must be provided.  If panel!="G", the 3rd panel will show -log10(dmrs$pval).
}
  \item{G}{
matrix of green channel intensities to use for plotting in the 3rd panel if panel3="G".
}
  \item{grs}{
if panel3="G", plot difference between the median green channel value (after subtracting probe medians and correcting for gc content) between these 2 groups.
}
}

\details{
This function plots user-provided regions.
}


\author{
Martin Aryee <aryee@jhu.edu>, Peter Murakami, Rafael Irizarry
}
\seealso{
\code{\link{plotDMRs}}, \code{\link{dmrFind}}, \code{\link{qval}}
}
\examples{
# See qval
}
