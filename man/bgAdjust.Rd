\name{bgAdjust}
\alias{bgAdjust}
\title{
Remove background
}
\description{
Estimate and remove background signal using anti-genomic background probes
}
\usage{
bgAdjust(dat, copy=TRUE)
}
\arguments{
  \item{dat}{
a TilingFeatureSet
}
  \item{copy}{
Only relevant when using disk-backed objects. If TRUE a copy will be made leaving the original object (dat) unchanged. The input object will not be preserved if copy=FALSE}
}
\details{
Background signal removal using a modified version of the RMA convolution model. The background signal level is estimated within GC-strata using anti-genomic background probes.
}
\value{
a TilingFeatureSet

}
\author{
Martin Aryee <aryee@jhu.edu>
}

\examples{
# See normalizeBetweenSamples
}
