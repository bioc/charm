\name{controlQC}
\alias{controlQC}
\title{
Boxplots of control and non-control probes
}
\description{
Make boxplots of the non-control probes and the control probes (after spatial and background correction but before any normalization), to confirm that the control probes have a lower distribution of intensities than the non-control probes.
}
\usage{
controlQC(rawData, controlProbes=NULL, controlIndex=NULL, IDcol, expcol, ylimits=c(-6,8), outfile="./boxplots_check.pdf", height=7, width=9) 
}
\arguments{
  \item{rawData}{
a TilingFeatureSet object
}
  \item{controlProbes}{
character string of the label assigned to non-CpG control probes in the annotation file (i.e. the container column of the .ndf file) (getContainer reports the container values for the pm probes).  This is used only if controlIndex is not provided.  Either this or controlIndex must be provided.
}
  \item{controlIndex}{
a vector of non-CpG control probe indices (can be obtained using getControlIndex).  Either this or controlProbes must be provided.
}
  \item{IDcol}{
column name of pData(rawData) with which to label samples in the output plot.
}
  \item{expcol}{
name of the column of pData(rawData) the values of which the points are colored differently by.
}
  \item{ylimits}{
ylim argument to plot()
}
   \item{outfile}{
file argument to pdf()
}
  \item{height}{
height argument to pdf()
}
  \item{width}{
width argument to pdf()
}
}

\value{
A data frame with one row per sample and columns
	\describe{
                \item{non_control}{Median methylation value among the non-control probes}
		\item{control}{Median methylation value among the control probes}
		\item{diff}{non_control-control}
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
