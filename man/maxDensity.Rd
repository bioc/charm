\name{maxDensity}
\alias{maxDensity}
\title{
Find the mode of a density function
}
\description{
Calculate a density function and find the max point
}
\usage{
maxDensity(x, n.pts = 2^14, minPoints=30)
}
\arguments{
  \item{x}{
a data vector
}
  \item{n.pts}{
Number of points to use in density estimation
}
  \item{minPoints}{
Minimum number of data points to accept
}
}
\details{
This function finds the maximum of a density function. It is identical to the (unexported) max.density function in affy except that it returns an NA if the number of data points provided is less than \code{minPoints}}
\value{
a numeric value 
}
\author{
Martin Aryee <aryee@jhu.edu>
}


\examples{
x <- rnorm(1000)
maxDensity(x)
}


