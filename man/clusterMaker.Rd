\name{clusterMaker}
\alias{clusterMaker}
\title{
Redefine array regions given chromosomal coordinates.
}
\description{
Redefine array regions given chromosomal coordinates.
}
\usage{
clusterMaker(chr,pos,order.it=TRUE,maxGap=300)
}
\arguments{
  \item{chr}{
Vector of chromosome names
}
  \item{pos}{
Vector of positions within chromosomes.
}
  \item{order.it}{
order probes if they are not already ordered.
}
  \item{maxGap}{
maximum allowable gap between probe start positions for probes to be grouped into the same region. 
}
}

\details{
Redefine array regions given chromosomal coordinates.
}


\author{
Martin Aryee <aryee@jhu.edu>, Peter Murakami, Rafael Irizarry
}
\seealso{
\code{\link{dmrFind}}, \code{\link{plotDMRs}}, \code{\link{plotRegions}}
}
\examples{
# See qval
}
