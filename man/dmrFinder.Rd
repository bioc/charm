\name{dmrFinder}
\alias{dmrFinder}
\title{
Find differentially methylated regions (DMRs)
}
\description{
Find differentially methylated regions (DMRs) from tiling microarray data.  If you want to adjust for covariates (including batch effects estimated using SVA) or if your covariate of interest is continuous, use the dmrFind function.
}
\usage{
dmrFinder(eset=NULL, groups, p=NULL, l=NULL, chr=NULL, pos=NULL, pns=NULL,
	sdBins=NULL, controlIndex=NULL, 
	controlProbes=c("CONTROL_PROBES", "CONTROL_REGIONS"), Indexes=NULL, 
	filter=NULL, package=NULL, ws=7, verbose=TRUE, compare="all",  
	withinSampleNorm="loess", betweenSampleNorm="quantile", 
	cutoff=0.995, sortBy="ttarea",paired=FALSE,pairs=NULL,DD=NULL,COMPS=NULL,
        removeIf=expression(nprobes<3),...)
}
\arguments{
  \item{eset}{
a TilingFeatureSet
}
  \item{groups}{
a vector of group labels for the samples in eset 
}
  \item{p}{
a matrix of percentage methylation values (scale: 0, 1). One column per sample
}
  \item{l}{
a matrix of methylation values (scale: -Inf, Inf), typically log-ratios.}
  \item{chr}{
vector of chromosome labels for the probes in eset, p or l
}
  \item{pos}{
vector of chromosomal coordinates for the probes in eset, p or l
}
  \item{pns}{
vector of region names for the probes in eset, p or l
}
  \item{sdBins}{
not currently implemented
}
  \item{controlIndex}{
vector of indices of non-CpG control probes (which can be obtained using getControlIndex), to be passed on to the methp function if neither the p nor the l arguments are provided, in which case either this or the controlProbes argument must be provided.
}
  \item{controlProbes}{
names of probe containers corresponding to control probes (e.g., see values returned by getContainer), to be passed on to the methp function if neither the p nor the l arguments are provided, in which case either this or the controlIndex argument must be provided.
}
  \item{Indexes}{
not currently used
}
  \item{filter}{
smoothing window weights. See details
}
  \item{package}{
annotation package name
}
  \item{ws}{
smoothing window size parameter. See details.
}
  \item{verbose}{
Verbose progress reporting
}
  \item{compare}{
the groups between which to find DMRs.
}
  \item{withinSampleNorm}{
within-sample normalization method. "loess" or "none"
}
  \item{betweenSampleNorm}{
between-sample normalization method. "quantile", "sqn" or "none"
}
  \item{cutoff}{
t-statistic cutoff used to identify probes as being in a DMR
}
  \item{sortBy}{
sort column for the DMR table. "area", "ttarea", "avg.diff", or "max.diff".
}
  \item{paired}{
if TRUE, do comparisons within pairs of samples.  FALSE by default.
}
  \item{pairs}{
if paired=TRUE, this must be provided.  a vector of pair identifiers for the samples in eset.  values must be the same within pairs and different between pairs.
}
  \item{DD}{
DD object returned by dmrFinder when paired=TRUE.  This argument may be ignored. 
}
  \item{COMPS}{
comps object returned by dmrFinder.  This argument may be ignored.
}
  \item{removeIf}{
expression indicating which DMRs to drop from the DMR tables that get returned.  The negation of this is used as the subset argument to the subset function when it is called on the final DMR table before it is returned.  If NULL, no DMRs will be subsetted out from the final table before it is returned.  DMR table column names to use are listed below.  E.g., to drop all DMRs with less than 4 probes, set removeIf=expression(nprobes<4).
}
  \item{\dots}{
further options to be passed to methp
}
}
\details{
This function finds differentially methylated regions (DMRs). The sortBy parameter can be used to sort the DMRs by area (# probes x average difference), t-statistic area (# probes x average t-statistic), average difference, or maximum difference.
}
\value{
A list with
\item{tabs}{A list of DMR tables, one per comparison with columns:
	\describe{
                \item{chr}{chromosome of DMR (bp)}
		\item{start}{start of DMR (bp)}
		\item{end}{end of DMR (bp)}
		\item{p1}{if paired=FALSE, and p!=NULL or l=NULL, average percentage methylation of all probes between start and end for group 1}
		\item{p2}{if paired=FALSE, and p!=NULL or l=NULL, average percentage methylation of all probes between start and end for group 2}
		\item{m1}{if paired=FALSE, p=NULL and l!=NULL, average methylation l (logit(percentage methylation) if l=NULL) of all probes between start and end for group 1}
		\item{m2}{if paired=FALSE, p=NULL and l!=NULL, average methylation l (logit(percentage methylation) if l=NULL) of all probes between start and end for group 2}
		\item{regionName}{name of the tiling region in which the DMR is found (These names come from the NDF file)}
		\item{indexStart}{index of first probe in DMR.  This indexes the output of dmrFinder, *not* the input.}
		\item{indexEnd}{index of last probe in DMR.  This indexes the output of dmrFinder, *not* the input.}
                \item{nprobes}{number of probes for the DMR, i.e., indexEnd-indexStart+1}
                \item{diff}{average percentage methylation difference within the DMR if paired=FALSE, and average l (logit(percentage) methylation if l=NULL) difference within the DMR if paired=TRUE}
                \item{maxdiff}{maximum percentage methylation difference within the DMR if paired=FALSE, and maximum l (logit(percentage) methylation if l=NULL) difference within the DMR if paired=FALSE}
		\item{area}{nprobes x average difference}	
		\item{ttarea}{nprobes x (average probe level t-statistic for between group difference)}		
	}
}
\item{p}{A matrix of percentage methylation estimates (NOTE: the probe order may differ from that of the input p matrix since probes are sorted into chromosomal order)}
\item{l}{This contains methylation log-ratios if they were passed to the function. Otherwise it contains logit-transformed percentage methylation estimates. (NOTE: the probe order may differ from that of the input l matrix since probes are sorted into chromosomal order)}
\item{chr }{a vector of chromosomes corresponding to the rows of p and l}
\item{pos }{a vector of positions corresponding to the rows of p and l}
\item{pns }{a vector of probe region names corresponding to the rows of p and l}
\item{index }{a vector identifying which subset of the input probes (i.e. which elements of the input chr, pos, and pns, and rows of the input p and/or l) were used to search for DMRs.  The output objects (chr, pos, pns, p, l, etc) are this subset of probes from the input.  Therefore, e.g., while tabs$indexStart:tabs$indexEnd indexes the elements or rows of the output objects for each DMR candidate in tabs, index[tabs$indexStart:tabs$indexEnd] indexes the elements or rows of the input objects.}
\item{gm }{if paired=FALSE, group medians of the l matrix}
\item{DD }{if paired=TRUE, a list of within-pair differences for each comparison}
\item{sMD }{if paired=TRUE, a matrix of smoothed mean within-pair differences for each comparison}
\item{groups }{a vector of group labels}
\item{args }{the DMR finder parameter vector}
\item{comps }{the vector of pairwise group comparisons}
\item{package }{the array annotation package name}
}
\author{
Martin Aryee <aryee@jhu.edu>, Peter Murakami, Rafael Irizarry
}
\seealso{
\code{\link{dmrFind}}, \code{\link{readCharm}}, \code{\link{methp}}, \code{\link{dmrFdr}}
}
\examples{
# See dmrFdr
}
