\name{countGC}
\alias{countGC}
\title{
Count probe GC content 
}
\description{
Return the GC content for each probe 
}
\usage{
countGC(dat, type = "pm", idx)
}
\arguments{
  \item{dat}{
a TilingFeatureSet object 
}
  \item{type}{
pm or bg probes
}
  \item{idx}{
An optional vector of probe indices for which to return GC content. If not specified, values for all pm (or bg) probes will be returned.
}

}
\details{
This function returns the sum of #G + #C in the pm or bg probes.
}
\value{
a numeric vector
}
\author{
Martin Aryee <aryee@jhu.edu>
}

\seealso{
\code{\link{readCharm}}
}
\examples{
	if (require(charmData)) {
		phenodataDir <- system.file("extdata", package="charmData")
		pd <- read.delim(file.path(phenodataDir, "phenodata.txt"))
		pd <- subset(pd, sampleID=="441_liver")
		dataDir <- system.file("data", package="charmData")
		setwd(dataDir)
		rawData <- readCharm(files=pd$filename, sampleKey=pd)
		ngc <- countGC(rawData)
		head(ngc)
	}
}
