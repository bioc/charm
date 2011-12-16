\name{qval}
\alias{qval}
\title{
Obtain False Discovery Rate q-values for the DMR candidates returned by dmrFind.
}
\description{
Obtain False Discovery Rate q-values for the DMR candidates returned by dmrFind.
}
\usage{
qval(p=NULL, logitp=NULL, dmr, numiter=500, seed=54256, verbose=FALSE, mc=1, return.permutations=FALSE)
}
\arguments{
  \item{p}{
The matrix of methylation percentage estimates provided to dmrFind (argument p) when it returned dmr, if p was provided. 
}
  \item{logitp}{
The matrix of logit-transformed methylation percentage estimates provided to dmrFind (argument logitp) when it returned dmr, if logitp was provided. 
}
  \item{dmr}{
The output of dmrFind with the dmr table to obtain q-values for.
}
  \item{numiter}{
The number of times to iterate the resampling procedure. 
}
  \item{seed}{
random seed.
}
  \item{verbose}{
verbose argument to dmrFind.
}
  \item{mc}{
Number of CPUs to use for parallel processing.
}
  \item{return.permutations}{
if TRUE, also return a matrix defining the resampling at each iteration.
}
}

\details{
The q-value for region i is estimated by a resampling procedure.  It is estimated as x/N, where N is the number of dmrs in the DMR table returned by dmrFind whose statistic (defined by the sortBy argument to dmrFind) is >= the statistic for region i, and x is the average (across all iterations' DMR tables) number of DMRs with statistic >= the statistic for region i.  Monotonicity is then enforced and q-values are capped at 1.  For example, to estimate the q-value for the 5th dmr in your dmr table that you ranked by area, N is 5.  If it has area=100, then x is the average number of dmrs with area>=100, averaging over all the tables produced by random-resampling iterations.  To the extent that the proportion of the genome covered by the array that is truly non-differentially methylated is less than 1, x will slightly overestimate the true number of false DMR candidates expected to be returned from the data, and hence the q-value will be slightly conservative.

The resampling procedure used to obtain null data for each iteration is as follows.  
1. Add the surrogate variables identified by SVA to the null model defined by the mod0 argument to dmrFind, and obtain fitted methylation values M, where rows correspond to probes and columns correspond to arrays.
2. Add the surrogate variables identified by SVA to the full model defined by the mod argument to dmrFind, and obtain residuals R, where rows correspond to probes and columns correspond to arrays.
Then, for iteration k, obtain a column-resampled version of R, Q (resampling with replacement), and obtain the methylation values to be used in iteration k as M+Q.
}

\value{
the dmr table in your dmr object (dmr$dmrs), with an additional column ("qvalue") for the Fdr q-value.
}

\author{
Martin Aryee <aryee@jhu.edu>, Peter Murakami, Rafael Irizarry
}
\seealso{
\code{\link{dmrFind}}, \code{\link{plotDMRs}}, \code{\link{plotRegions}}
}
\examples{
    ## See vignette.
    if (require(charmData) & require(BSgenome.Hsapiens.UCSC.hg18)) {
        dataDir <- system.file("data", package="charmData")
        phenodataDir <- system.file("extdata", package="charmData")
        pd <- read.delim(file.path(phenodataDir, "phenodata.txt"))
        res <- validatePd(pd)
        
        ## Read in raw data:
        rawData <- readCharm(files=pd$filename, path=dataDir, sampleKey=pd, sampleNames=pd$sampleID)
        ## Check quality of arrays:
        #qual <- qcReport(rawData, file="qcReport.pdf")

        ## Assess individual probe qualities:
        pmq = pmQuality(rawData)
        rmpmq = rowMeans(pmq)
        okqc = which(rmpmq>75)

        ## Identify control probes as the probes at positions surrounded by a CpG-free 600bp window: 
        ctrlIdx <- getControlIndex(rawData, subject=Hsapiens, noCpGWindow=600)

        ## Check that these control probes do indeed have lower intensities than the non-control probes (after spatial and background corrections, but no normalization, since normalization uses the control probes):
        #controlQC(rawData=rawData, controlIndex=ctrlIdx, IDcol="sampleID", expcol="tissue", ylimits=c(-6,8), outfile="boxplots_check.pdf", height=7, width=9)

        chr = pmChr(rawData)
        pns = probeNames(rawData)
        pos = pmPosition(rawData)
        seq = pmSequence(rawData)
        pd  = pData(rawData)

        ## Estimate percent methylation:
        p <- methp(rawData, controlIndex=ctrlIdx, plotDensityGroups=pd$tissue) 

        ## unsupervised clustering of samples:
        #cmdsplot(labcols=c("red","black","blue"), expcol="tissue", rawData=rawData, p=p, okqc=okqc, noXorY=TRUE, outfile="cmds_topN.pdf", topN=c(100000,1000))

        ## Do not look for DMRs among control probes or probes with average probe quality score less than or equal to 75 for example:
        Index=setdiff(which(rmpmq>75),ctrlIdx)
        Index = Index[order(chr[Index], pos[Index])]
        p = p[Index,]
        seq = seq[Index]
        chr = chr[Index]
        pos = pos[Index]
        pns = pns[Index]
        pns=clusterMaker(chr,pos) 

        ## Identify DMR candidates between colon and liver (in this example not adjusting for any other covariates (besides batch)):
        mod0 = matrix(1,nrow=nrow(pd),ncol=1)
        mod  = model.matrix(~1 +factor(pd$tissue,levels=c("liver","colon","spleen")))
        thedmrs = dmrFind(p=p, mod=mod, mod0=mod0, coeff=2, pns=pns, chr=chr, pos=pos)

        ## Obtain FDR q-values for each DMR candidate. In practice, numiter should be set much higher.
        withq = qval(p=p, dmr=thedmrs, numiter=2, verbose=FALSE, mc=1)

        #### Plotting not run:
        ## Plot DMR candidates 1,2, and 4, for example.  First have to load a table of CpG islands.
        cpg.cur = read.delim("http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-hg18.txt", as.is=TRUE)
        #plotDMRs(dmrs=thedmrs, Genome=Hsapiens, cpg.islands=cpg.cur, exposure=pd$tissue, outfile="./colon-liver.pdf", which_plot=c(1,2,4), all.points=FALSE, plot.these=c("colon","liver"), smoo="loess", ADD=3000, cols=c("black","red","blue"))

        ## Plot DMR candidates, and in the 3rd panel plot the difference in average green channel:
        dat0 = spatialAdjust(rawData, copy=FALSE)
        dat0 = bgAdjust(dat0, copy=FALSE)
        G = pm(dat0)[,,1] #from oligo
        G = G[Index,]
        #plotDMRs(dmrs=thedmrs, Genome=Hsapiens, cpg.islands=cpg.cur, exposure=pd$tissue, outfile="./colon-liver2.pdf", which_plot=c(1), all.points=FALSE, plot.these=c("colon","liver"), smoo="loess", ADD=3000, cols=c("black","red","blue"), panel3="G", G=G, seq=seq)

        ## Example if covariate of interest is continuous:
        pd$x = c(1,2,3,4,5,6)
        mod0 = matrix(1,nrow=nrow(pd),ncol=1)
        mod  = model.matrix(~1 +pd$x)
        coeff = 2
        thedmrs2 = dmrFind(p=p, mod=mod, mod0=mod0, coeff=coeff, pns=pns, chr=chr, pos=pos)

        ## If covariate of interest is continuous, you can still plot it like it is categorical by categorizing it for plotting purposes:
        groups = as.numeric(cut(mod[,coeff],c(-Inf,2,4,Inf))) #You can change these cutpoints.
        pd$groups = c("low","medium","high")[groups]
        #plotDMRs(dmrs=thedmrs2, Genome=Hsapiens, cpg.islands=cpg.cur, exposure=pd$groups, outfile="./test.pdf", which_plot=c(1), smoo="loess", ADD=3000, cols=c("black","red","blue"))

        ## Otherwise, if covariate of interest is continuous, plot will show correlation with covariate:
        #plotDMRs(dmrs=thedmrs2, Genome=Hsapiens, cpg.islands=cpg.cur, exposure=pd$x, outfile="./x.pdf", which_plot=c(1), all.points=FALSE, plot.these=c("colon","liver"), smoo="loess", ADD=3000, cols=c("black","red","blue"))

        ## Plot arbitrary regions:
        mytable = thedmrs$dmrs[,c("chr","start","end")]
        mytable[2,] = c("chr1",1,1000) #not on array
        mytable$start = as.numeric(mytable$start)
        mytable$end = as.numeric(mytable$end)
        #plotRegions(thetable=mytable[1:5,], cleanp=thedmrs$cleanp, chr=chr, pos=pos, Genome=Hsapiens, cpg.islands=cpg.cur, outfile="myregions.pdf", exposure=pd$tissue, exposure.continuous=FALSE)
    }
}
