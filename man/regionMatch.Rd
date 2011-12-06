\name{regionMatch}
\alias{regionMatch}
\title{
Given two data frames with columns "chr", "start", and "end", identify the nearest region in one to the other.
}
\description{
Given two data frames with columns "chr", "start", and "end", identify the nearest region in one (object2) to the other (object1).
}
\usage{
regionMatch(object1, object2, verbose=TRUE)
}
\arguments{
  \item{object1}{
Data frame with columns "chr", "start", and "end".  For each of the rows of this data frame, find the nearest regions in object2.
}
  \item{object2}{
Data frame with columns "chr", "start", and "end".
}
  \item{verbose}{
print progress.
}
}

\details{
Given two data frames with columns "chr", "start", and "end", identify the nearest region in one (object2) to the other (object1).
}

\value{
A data frame with one row corresponding for each row of object1, and columns:
	\describe{
                \item{dist}{bp separating the region in object1 and the nearest region in object2}
		\item{matchIndex}{row index of object2 for the region that is closest to the region in object1}
		\item{type}{"inside" if nearest region in object2 is wholly contained inside the region in object 1, "cover" if nearest region in object2 covers the whole region in object1, "disjoint" if there is no overlap between the region in object1 and the nearest region in object1, and "overlap" if the region in object1 and the nearest region in object2 overlap but one does not wholly cover the other.}
		\item{amountOverlap}{amount of overlap between the region in object1 and the nearest region in object2}
		\item{insideDist}{for type="inside" regions, the smaller (fewest bp) of 1., the end position of region in object1 - the end position of the nearest region in object2, and 2., the start position of region in object1 - the start position of the nearest region in object2. }
		\item{size1}{number of bp in region in object1}
		\item{size2}{number of bp in nearest region in object2}
	}
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
