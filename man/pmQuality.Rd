\name{pmQuality}
\alias{pmQuality}
\title{
Calculate probe quality scores
}
\description{
pmQuality calculates probe quality for each pm probe by comparing the total DNA signal (green) to the distribution of the background probe signals. 0 means lower than all background probes. 100 means higher than all background probes.
}
\usage{
pmQuality(dat, channel="channel1", verbose=FALSE, idx=NULL) 
}
\arguments{
  \item{dat}{
a TilingFeatureSet object
}
  \item{channel}{
which channel to assess
}
  \item{verbose}{
logical value for whether or not to print message to screen
}
  \item{idx}{
vector of probe indices specifying which probes to evaluate.  By default all pm probes are used.
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
