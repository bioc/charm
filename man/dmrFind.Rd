\name{dmrFind}
\alias{dmrFind}
\title{
Identify DMR candidates using a regression-based approach and correcting for batch effects.
}
\description{
Identify DMR candidates using a regression-based approach and correcting for batch effects.
}
\usage{
dmrFind(p=NULL, logitp=NULL, svs=NULL, mod, mod0, coeff, pns, chr, pos, only.cleanp=FALSE, only.dmrs=FALSE, rob=TRUE, use.limma=FALSE, smoo="weighted.loess", k=3, SPAN=300, DELTA=36, use="sbeta", Q=0.99, min.probes=3, min.value=0.075, keepXY=TRUE, sortBy="area.raw", verbose=TRUE, vfilter=NULL, nsubsets=50, ...)
}
\arguments{
  \item{p}{
matrix of methylation percentage estimates. Either this or logitp must be provided. Can be an ff_matrix object.
}
  \item{logitp}{
matrix of logit-transformed methylation percentage estimates.  Either this or p must be provided. Can be an ff_matrix object.
}
  \item{svs}{
surrogate variables whose effect will be corrected for.  This should be svaobj$sv, where svaobj is the object returned by sva().  Setting svs=0 will result in sva not being used.
}
  \item{mod}{
The mod argument provided to sva() which yielded svs.  This should be a design matrix with all the adjustment covariates and (in the rightmost column(s)) your covariate of interest.
}
  \item{mod0}{
The mod0 argument provided to sva() which yielded svs.  This should be a design matrix with just the adjustment covariates.  Thus it should be the same as mod, excluding the rightmost column(s) for the covariate of interest.
}
  \item{coeff}{
a character or numeric index for the column of mod that identifies the covariate column of interest.
}
  \item{pns}{
vector of region names for the probes corresponding to rows of p or logitp.
}
  \item{chr}{
vector of chromosomal identifiers for the probes corresponding to rows of p or logitp.
}
  \item{pos}{
vector of chromosomal coordinates for the probes corresponding to rows of p or logitp.
}
  \item{only.cleanp}{
if TRUE, only return the matrix of methylation percentage estimates after removing the batch effects (columns of sv) 
}
   \item{only.dmrs}{
if TRUE, do not return the matrix of methylation percentage estimates that have had the batch effects (columns of sv) removed (called cleanp).
}
   \item{rob}{
One of the outputs of dmrFind is cleanp, which is the input p matrix after removing batch effects identified by SVA.  By default these are the only effects removed from the p matrix.  However, if you set rob=FALSE, then the other adjustment variables in mod and mod0 (all other variables besides the covariate of interest) are also removed.  This will affect the methylation levels shown in plots using plotDMRs, plotRegions, and any other function that uses the cleanp output of dmrFind.  It does not affect the selection of DMRs, though, except in the application of the filter at the end of the function where DMRs with an average difference between the 2 groups (or the average correlation between methylation and the covariate if the covariate is continuous) of less than the min.value argument are filtered out, since this step uses the cleanp matrix to calculate those averages or that correlation.
}
  \item{use.limma}{
Use the linear modeling approach (borrowing strength across probes) of lmFit in the limma package.
}
  \item{smoo}{
which method to use for smoothing.  "weighted loess", "loess", or "runmed".
}
  \item{k}{
k argument to runmed() if smoo="runmed".
}
  \item{SPAN}{
see DELTA. Only used if smoo="loess"
}
  \item{DELTA}{
span parameter in loess smoothing will = SPAN/(DELTA * number of probes in the plotted region).  Only used if smoo="loess".
}
  \item{use}{
If "sbeta", identify DMRs by segmenting the smoothed effect estimates.  If "swald", identify DMRs by segmenting the smoothed wald statistics.
}
  \item{Q}{
Identify DMRs as the consecutive groups of probes whose smoothed effect estimate (if use="sbeta") or smoothed wald statistics (if use="swald") exceed this quantile.
}
  \item{min.probes}{
The minimum allowable number of probes in a DMR candidate.
}
  \item{min.value}{
The minimum allowable average difference in methylation percentage between the 2 groups if covariate is categorical, or the minimum average correlation between methylation and the covariate if covariate is continuous. 
}
  \item{keepXY}{
if FALSE, exclude DMRs in "chrX" and "chrY". 
}
  \item{sortBy}{
column of DMR table to sorty by. 
}
  \item{verbose}{
print progress messages if TRUE.
}
  \item{vfilter}{
vfilter argument to sva function.  The number of most variable probes to use when building SVs--must be between 100 and m, where m is the total number of probes.  vfilter=NULL by default, which means vfilter=m.  Setting this to something smaller, like 100000, will typically yield satisfactory SVs and is advisable when the size of the data is very large and memory is limited (as when the p and/or logitp arguments given to dmrFind are ff_matrix objects).
}
  \item{nsubsets}{
used if p or logitp are ff_matrix objects.  Rather than doing computations on the whole logitp matrix, break up its m rows into chunks of m/nsubsets rows.  Default is 50, but if even that uses too much memory, set this higher.  Results do not depend on the value chosen.
}
  \item{...}{
Additional arguments passed to sva()
}
}

\details{
Identify DMR candidates using a regression-based approach and correcting for batch effects.
}

\value{
If only.cleanp=TRUE, only the cleanp matrix (below) is returned.  Otherwise, the function returns a list with
\item{dmrs}{A data frame with all DMR candidates, with columns:
	\describe{
                \item{chr}{chromosome of DMR}
		\item{start}{start of DMR (bp)}
		\item{end}{end of DMR (bp)}
		\item{value}{average value of the smoothed effect estimate within the DMR if use="sbeta" (the default), or the average value of the smoothed wald statistic within the DMR if use="swald"}
		\item{area}{nprobes x value}
		\item{pns}{name of the region on the array in which the DMR candidate was identified.}
		\item{indexStart}{index of first probe in DMR.  This indexes chr, pos, pns, and cleanp}
		\item{indexEnd}{index of last probe in DMR.  This indexes chr, pos, pns, and cleanp}
                \item{nprobes}{number of probes for the DMR, i.e., indexEnd-indexStart+1}
                \item{avg}{average (across probes) percentage methylation difference within the DMR if covariate is categorical, or average (across probes) correlation between cleanp and covariate if covariate is continuous.}
                \item{max}{maximum (across probes) percentage methylation difference within the DMR if covariate is categorical, or maximum (across probes) correlation between cleanp and covariate if covariate is continuous.}
		\item{area.raw}{nprobes x avg}
	}
}
\item{pval}{a vector of p-values for the t-test at each probe (in same order as rows of cleanp)}
\item{pns}{a vector of probe region names corresponding to the rows of cleanp}
\item{chr}{a vector of chromosomes corresponding to the rows of cleanp}
\item{pos}{a vector of positions corresponding to the rows of cleanp}
\item{args}{A list containing all the arguments provided to dmrFind.  If svs was not provided, svs here will be the surrogate variables obtained from sva.}
If only.dmrs=FALSE, 
\item{cleanp}{The matrix of percentage methylation estimates, after subtracting batch effects.  If rob=FALSE, the effects of the other adjustment covariates are removed also.}
If covariate is continuous,
\item{beta}{the effect estimate at each probe.}
\item{sbeta}{the smoothed effect estimate at each probe.}
}

\author{
Martin Aryee <aryee@jhu.edu>, Peter Murakami, Rafael Irizarry
}
\seealso{
\code{\link{plotDMRs}}, \code{\link{plotRegions}}, \code{\link{qval}}
}
\examples{
# See qval
}
