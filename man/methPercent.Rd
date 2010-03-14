\name{methPercent}
\alias{methPercent}
\title{
Estimate percentage DNA methylation from log-ratios
}
\description{
Estimate percentage DNA methylation from log-ratios
}
\usage{
methPercent(m, pmIndex, ngc, commonParams = TRUE)
}
\arguments{
  \item{m}{
a matrix of M-values (methylation log-ratios). One column per sample.
}
  \item{pmIndex}{
A vector of probe indices to use in the calculation. Usually set to the indices of the pm probes (excluding background and other non-specific controls) by using pmIndex=pmindex(dat)
}
  \item{ngc}{
a vector with GC-content of probes. Same length as nrow(m)
}
  \item{commonParams}{
boolean indicating whether a common set of parameters should be used for all samples when converting M-values to percentage methylation.
}
}
\details{
This function estimates percentage DNA methylation from normalized methylation log-ratios (M-values).
}
\value{
a matrix of percentage methylation estimates. Same dimensions as m
}
\author{
Martin Aryee <aryee@jhu.edu>
}

\examples{
	if (require(charmData) & require(BSgenome.Hsapiens.UCSC.hg18)) {
		phenodataDir <- system.file("extdata", package="charmData")
		pd <- read.delim(file.path(phenodataDir, "phenodata.txt"))
		pd <- subset(pd, sampleID=="441_liver")
		dataDir <- system.file("data", package="charmData")
		setwd(dataDir)
		# Read in raw data
		rawData <- readCharm(files=pd$filename, sampleKey=pd)
		# Find non-CpG control probes
		ctrlIdx <- getControlIndex(rawData, subject=Hsapiens)	
		# Get normalized methylation log-ratios
		m <- methp(rawData, controlIndex=ctrlIdx, returnM=TRUE)
		# Estimate percentage methylation
		ngc <- countGC(rawData)
		pmIndex <- pmindex(rawData)
		p <- methPercent(m, pmIndex=pmIndex, ngc=ngc)
	}
}