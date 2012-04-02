\name{dmrFdr}
\alias{dmrFdr}
\title{
Calculate FDR q-values for differentially methylated regions (DMRs)
}
\description{
Estimate false discovery rate q-values for a set of differentially methylated regions (found using the dmrFinder function) using a permutation approach.  For differentially methylated regions found using the dmrFind function, use the qval function instead.
}
\usage{
dmrFdr(dmr, compare = 1, numPerms = 1000, seed = NULL, verbose = TRUE)
}
\arguments{
  \item{dmr}{
a dmr object as returned by \code{\link{dmrFinder}}
}
  \item{compare}{
The dmr table for which to calculate DMRs. See details.
}
  \item{numPerms}{
Number of permutations
}
  \item{seed}{
Random seed (for reproducibility)
}
  \item{verbose}{
Boolean
}
}
\details{
This function estimates false discovery rate q-values for a dmr object returned by \code{\link{dmrFinder}}. dmrFinder can return a set of DMR tables with one or more pair-wise comparisons between groups. dmrFdr currently only calculated q-values for one of these at a time. The dmr table to use (if the dmr object contains more than one) is specified by the compare option.
}
\value{
a list object in the same format as the input, but with extra p-val and q-val columns for the tabs element.
}
\author{
Martin Aryee <aryee@jhu.edu>
}


\seealso{
\code{\link{qval}}, \code{\link{dmrFinder}}, \code{\link{dmrPlot}}, \code{\link{regionPlot}}
}
\examples{
	if (require(charmData) & require(BSgenome.Hsapiens.UCSC.hg18)) {
		phenodataDir <- system.file("extdata", package="charmData")
		pd <- read.delim(file.path(phenodataDir, "phenodata.txt"))
#		pd <- subset(pd, tissue \%in\% c("liver", "colon"))
		# Validate format of sample description file		
#		res <- validatePd(pd)
#		dataDir <- system.file("data", package="charmData")
#		setwd(dataDir)
		# Read in raw data
#		rawData <- readCharm(files=pd$filename, sampleKey=pd)
		# Find non-CpG control probes
#		ctrlIdx <- getControlIndex(rawData, subject=Hsapiens)
		# Estimate methylation
#		p <- methp(rawData, controlIndex=ctrlIdx)
		# Find differentially methylated regions
#		grp <- pData(rawData)$tissue
#		dmr <- dmrFinder(rawData, p=p, groups=grp, 
#                        removeIf=expression(nprobes<4 | abs(diff)<.05 | abs(maxdiff)<.05),
#			compare=c("liver", "colon"), cutoff=0.95)
#		head(dmr$tabs[[1]])
		# Estimate false discovery rate for DMRs
#		dmr <- dmrFdr(dmr, numPerms=3, seed=123) 
#		head(dmr$tabs[[1]])

                ##Not run:
                ## Plot top 10 DMRs:
                #cpg.cur = read.delim("http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-hg18.txt", as.is=TRUE)
                #dmrPlot(dmr=dmr, which.table=1, which.plot=1:5, legend.size=1, all.lines=TRUE, all.points=TRUE, colors.l=c("blue","black"), colors.p=c("blue","black"), outpath=".", cpg.islands=cpg.cur, Genome=Hsapiens)

                ## plot any given genomic regions using this data, supplying the regions in a data frame that must have columns with names "chr", "start", and "end": 
                #mytab = data.frame(chr=as.character(c(dmr$tabs[[1]]$chr[1],"chrY",dmr$tabs[[1]]$chr[-1])), start=as.numeric(c(dmr$tabs[[1]]$start[1],1,dmr$tabs[[1]]$start[-1])), end=as.numeric(c(dmr$tabs[[1]]$end[1],100,dmr$tabs[[1]]$end[-1])), stringsAsFactors=FALSE)[1:5,]
                #regionPlot(tab=mytab, dmr=dmr, cpg.islands=cpg.cur, Genome=Hsapiens, outfile="./myregions.pdf", which.plot=1:5, plot.these=c("liver","colon"), cl=c("blue","black"), legend.size=1, buffer=3000)
                ## note that region 2 is not plotted since it is not on the array.

                ## Example of paired analysis:
#                pData(rawData)$pair = c(1,2,1,2) ## fake pairing information for this example.
#                dmr2 <- dmrFinder(rawData, p=p, groups=grp, 
#                                  compare=c("colon", "liver"),
#                                  removeIf=expression(nprobes<4 | abs(diff)<.05 | abs(maxdiff)<.05),
#                                  paired=TRUE, pairs=pData(rawData)$pair, cutoff=0.95)

                #dmrPlot(dmr=dmr2, which.table=1, which.plot=c(3), legend.size=1, all.lines=TRUE, all.points=TRUE, colors.l=c("black"), colors.p=c("black"), outpath=".", cpg.islands=cpg.cur, Genome=Hsapiens)
                #regionPlot(tab=mytab, dmr=dmr2, cpg.islands=cpg.cur, Genome=Hsapiens, outfile="myregions.pdf", which.plot=1:5, plot.these=c("colon-liver"), cl=c("black"), legend.size=1, buffer=3000)
	}
}
