##This file contains both dmrPlot and dmrPlot_paired functions.

dmrPlot <- function(dmr=dmr, which.table=1:length(dmr$tabs), which.plot=1:30, legend.size=1, all.lines=TRUE, all.points=FALSE, colors.l, colors.p, outpath=".") {
  
    if(is.null(dmr$DD)) {
        if(missing(colors.p)) colors.p=2:(length(unique(dmr$groups))+1)
        if(missing(colors.l)) colors.l=2:(length(unique(dmr$groups))+1)       
        dmrPlot_unpaired(dmr=dmr, which.table=which.table, which.plot=which.plot, legend.size=legend.size, all.lines=all.lines, all.points=all.points, colors.l=colors.l, colors.p=colors.p, outpath=outpath)
    } else {
        if(missing(colors.p)) colors.p = 2:(ncol(dmr$sMD)+1)
        if(missing(colors.l)) colors.l = 2:(ncol(dmr$sMD)+1)
        dmrPlot_paired(dmr=dmr, which.table=which.table, which.plot=which.plot, legend.size=legend.size, all.lines=all.lines, all.points=all.points, colors.l=colors.l, colors.p=colors.p, outpath=outpath)
    }
    
}

#This function is just results_function.R except there is no option to return the table(s) and it does not plot the gene panel.
dmrPlot_unpaired <- function(dmr=dmr, which.table=1:length(dmr$tabs), which.plot=1:30, legend.size=1, all.lines=TRUE, all.points=FALSE, colors.l, colors.p, outpath=".", buffer=NULL){

stopifnot(inherits(which.plot,c("numeric","integer")))
stopifnot(inherits(which.table,c("numeric","integer")))
if(length(dmr$tabs)<max(which.table)) stop("which.table specifies more comparisons than were made.")
if(!is.null(buffer) & !inherits(buffer,c("numeric","integer"))) stop("if provided, buffer must be numeric.")
require(RColorBrewer)

FACS=dmr$groups
if(all.points)  stopifnot(length(colors.p)==length(unique(FACS)))
if(all.lines)   stopifnot(length(colors.l)==length(unique(FACS)))
if(!all.points) stopifnot(length(colors.p)%in%c(2, length(unique(FACS))))
if(!all.lines)  stopifnot(length(colors.l)%in%c(2, length(unique(FACS))))

if(dmr$package %in% 
          c("pd.feinberg.hg18.me.hx1", 
            "pd.feinberg.hg18.me.hx1.noxy", 
            "pd.0701.hg18.tile.05.hx1.v2", 
            "pd.0701.hg18.tile.05.hx1.charmed", 
            "pd.charm.hg18.example",
            "pd.081229.hg18.promoter.medip.hx1")) {
        spec="human"; build="hg18"
       } else if(dmr$package %in% 
          c("pd.feinberg.mm8.me.hx1", 
            "pd.2006.07.18.mm8.refseq.promoter")) {
        spec="mouse"; build="mm8"
       } else if (dmr$package %in% 
          c("pd.091020.hg19.charm.hx1", 
            "pd.091117.hg19.charm.hd4.hx1",
            "pd.100428.hg19.charm.hx1")){
        spec="human"; build="hg19"
       } else if (dmr$package %in% 
          c("pd.2006.10.31.rn34.refseq.promoter",
            "pd.090519.rn34.charm.hx1")) {
        spec="rat"; build="rn4"    
       } else { 
         stop("This array design not accomodated yet.") 
       }

if(spec=="human"){
    require(paste("BSgenome.Hsapiens.UCSC.",build,sep=""), character.only=TRUE)
    if(build=="hg18"){
        con <- url("http://rafalab.jhsph.edu/CGI/CGI-Hsapiens-hg18.rda")
    } else if(build=="hg19"){
        con <- url("http://rafalab.jhsph.edu/CGI/CGI-Hsapiens-hg19.rda")
    }
}
if(spec=="mouse"){
    require(paste("BSgenome.Mmusculus.UCSC.",build,sep=""), character.only=TRUE)
    if(build=="mm8"){
        con <- url("http://rafalab.jhsph.edu/CGI/CGI-Mmusculus-mm8.rda")
    } else if(build=="mm9"){
        con <- url("http://rafalab.jhsph.edu/CGI/CGI-Mmusculus-mm9.rda")
    }
}
if(spec=="rat") {
    require(paste("BSgenome.Rnorvegicus.UCSC.",build,sep=""), character.only=TRUE)
    if(build=="rn4"){
        con <- url("http://rafalab.jhsph.edu/CGI/CGI-Rnorvegicus-rn4.rda")
    }
}
load(con)
close(con)
cpg.cur = cgi
ocpgi=data.frame(chr=I(cpg.cur[,1]), start=as.numeric(cpg.cur[,2]), end=as.numeric(cpg.cur[,3]))
chr=dmr$chr
ocpgi=ocpgi[ocpgi$chr%in%unique(as.character(chr)),]

#source("functions_supp.R")
tabs=dmr$tabs
pns=dmr$pns
Indexes=split(seq(along=pns),pns)
pos=dmr$pos
#In dmrFinder dmr$gm only gets columns for groups that are in >=1 comparison, so:
all.compared = all(unique(FACS)%in%colnames(dmr$gm))
if(is.null(dmr$p)){
    if(all.lines & !all.compared){
        gmp = get.tog(l=dmr$l,groups=FACS,compare=comp(FACS),verbose=TRUE)$lm
    } else{ gmp = dmr$gm }
} else{
    if(all.lines & !all.compared){
        l=log(dmr$p)-log(1-dmr$p)
        tog = get.tog(l=l,groups=FACS,compare=comp(FACS),verbose=TRUE)
          #FACS does id the group that each column of l belongs to.
        gmp = exp(tog$lm)/(1+exp(tog$lm))
    } else {
        gmp = exp(dmr$gm)/(1+exp(dmr$gm))
    }
}
sMM=get.smooth(gmp,Indexes=Indexes,filter=NULL,columns=1:ncol(gmp),ws=dmr$args$ws)
if(is.null(dmr$p)) M=dmr$l else M=dmr$p

for(h in 1:nrow(dmr$comps)) stopifnot(identical(paste(colnames(dmr$gm)[dmr$comps[h,]],collapse="-"), names(tabs)[h]))

#################Make output table and pdf file of plots##############
has = which(!sapply(tabs[which.table],is.null))
for(object.i in which.table[has]){
    message("Plotting table ",object.i)
    #################Make output table:###############################
    obj = tabs[[object.i]]
    grs = colnames(dmr$gm)[dmr$comps[object.i,]] #this is preferable, if we have COMPS
    if(!any(colnames(obj)=="index")) obj$index=match(obj$regionName,names(Indexes))

    ################Make pdf file of plots###############
    pdf(file.path(outpath,paste(names(tabs)[object.i],".pdf",sep="")),width=11,height=8,pointsize=10)
    palette(rev(brewer.pal(8,"Dark2")[c(2,1,3,4,8)]))
    ADD1=1500

    g1 = grs[1]
    g2 = grs[2]
    if(all.points){
        showpts = 1:ncol(M)
        th = as.numeric(factor(rank(FACS)))
        tcols = colors.p[th]
    }
    if(!all.points){
        showpts = which(FACS%in%c(g1,g2))
        if(length(colors.p)==length(unique(FACS))){
            th = as.numeric(factor(rank(FACS)))
            tcols = colors.p[th][showpts]
        } else{ #length(colors.p)==2
            th = as.numeric(factor(rank(FACS[showpts])))
            tcols = colors.p[th]
        }
    }
    if(all.lines){
        comps.l = match(sort(colnames(sMM)),colnames(sMM))
        cl = colors.l
    }
    if(!all.lines){
        comps.l = match(c(g1,g2),colnames(sMM))
        if(length(colors.l)==length(unique(FACS))){
            ao = match(c(g1,g2),sort(unique(FACS)))
            cl = colors.l[ao]
        } else{ #length(colors.l)==2
            cl = colors.l[order(c(g1,g2))]
        }
    }

    for(i in intersect(which.plot,1:nrow(obj))){
      cat(i," ")
      if(is.null(dmr$p)){
          YLIM=c(-0.5,3)
          ylabel="M"
      } else{
          YLIM=c(0, ifelse(all.lines,1+ncol(sMM)*.03*legend.size,1+.06*legend.size) )
          ylabel="p"
      }
      thechr=obj$chr[i]
      if(is.null(buffer)) ADD=max(ADD1*2,obj$end[i]-obj$start[i]+132) else ADD=buffer
      start=floor((obj$start[i]+obj$end[i])/2)-ADD
      end=start+2*ADD
    
      Index=Indexes[[obj$index[i]]]
      x=pos[Index] 
      start=max(start,min(x))
      end=min(end,max(x))
      
      Index=Index[x>=start & x <= end]
      Index=Index[order(pos[Index])]
      
      ##PLOT TISSUES
      #par(mfrow=c(3,1))
      #layout(matrix(1:3,ncol=1),heights=c(0.6,0.2,0.2))
      layout(matrix(1:2,ncol=1),heights=c(0.6,0.2))
      #par(mar=c(0,2.5,0.25,1.1),oma=c(0,0,2,0))
      par(mar=c(0,3.5,0.25,1.1),oma=c(1.1,0,2,0),mgp=c(2.5,.5,0))

      matplot(pos[Index],M[Index,showpts],col=tcols,cex=0.6,ylim=YLIM,xlab="",
              xaxt="n",ylab=ylabel,xlim=c(start,end),lty=1,las=1)
      
      matlines(pos[Index],sMM[Index,comps.l],lwd=2,lty=1,col=cl)
      
      abline(h=0,lwd=2,lty=3)
      abline(h=1,lwd=2,lty=3)
      abline(v=obj$start[i],lty=2)
      abline(v=obj$end[i],lty=2)
      
      nms=colnames(sMM)
      ast = match(c(g1,g2),nms)
      if((all.points|all.lines) & length(nms)>2) nms[ast] = paste(nms[ast],"*")
      legend("topright",nms[comps.l],col=cl,lty=1,lwd=2,cex=legend.size,bg="white")
      
      ##PLOT CPG ISLANDS
      par(mar=c(3.5,3.5,0.25,1.1))
      if(spec=="human") seq<-Hsapiens[[as.character(thechr) ]]
      if(spec=="mouse") seq<-Mmusculus[[as.character(thechr) ]]
      if(spec=="rat")   seq<-Rnorvegicus[[as.character(thechr) ]]
      Index=max(start,1):min(end,length(seq))
      #subseq<-seq[Index]  #using older version of Biostrings. For newer:
      subseq<-subseq(seq,start=Index[1],end=Index[length(Index)])
      cpgs=start(matchPattern("CG",subseq))+Index[1]-1
      ##mcrbc=sort(c(start(matchPattern("ACG",subseq))+Index[1]-1,
      ##  start(matchPattern("GCG",subseq))+Index[1]-1))
      cuts=seq(min(Index),max(Index),8)
      scpgs=cut(cpgs,cuts,include.lowest=TRUE)
      x = (cuts[1:(length(cuts)-1)]+cuts[2:(length(cuts))])/2
      y = table(scpgs)/8
      SPAN=400/diff(range(x))
      d = loess(y~x,span=SPAN,degree=1)
      plot(x,d$fitted,type="l",ylim=c(0,0.15),xlab="Location",
           ylab="CpG density",xlim=c(start,end),las=1)
      Rug(cpgs)
    ##  Rug(mcrbc,side=3)
      Index1 = which(ocpgi[,1]==as.character(thechr) &
                     ((ocpgi[,2] > start & ocpgi[,2]< end) |
                      (ocpgi[,3] > start & ocpgi[,3]< end)))
      if(length(Index1)>0)
          sapply(Index1,function(j) Rug(unlist(ocpgi[j,2:3]),col=5,lwd=3,side=1))

      mtext(paste("ID:",i,"--",as.character(thechr),":",start,"-",end,sep=""),
            side=3,cex=1.5,outer=TRUE)
    }
    cat("\n")
    dev.off()
}
cat("\nPlotting finished.\n")
}


dmrPlot_paired <- function(dmr=dmr, which.table=1:length(dmr$tabs), which.plot=1:30, legend.size=1, all.lines=TRUE, all.points=FALSE, colors.l, colors.p, outpath=".", buffer=NULL){
#NB: all.points and all.lines refer to all comparisons that were made, not that could be made!
  
if(is.null(dmr$DD)) stop("DMRs are not from paired comparisons. Use dmrPlot().")
stopifnot(inherits(which.plot,c("numeric","integer")))
stopifnot(inherits(which.table,c("numeric","integer")))
if(length(dmr$tabs)<max(which.table)) stop("which.table specifies more comparisons than were made.")
if(!is.null(buffer) & !inherits(buffer,c("numeric","integer"))) stop("if provided, buffer must be numeric.")
require(RColorBrewer)

sMD = dmr$sMD
FACS = dmr$groups
if(!all.points & length(colors.p)==1) colors.p=rep(colors.p,ncol(sMD))
stopifnot(length(colors.p)==ncol(sMD))
if(all.lines)  stopifnot(length(colors.l)==ncol(sMD))
if(!all.lines) stopifnot(length(colors.l)%in%c(1, ncol(sMD)))

if(dmr$package %in% 
          c("pd.feinberg.hg18.me.hx1", 
            "pd.feinberg.hg18.me.hx1.noxy", 
            "pd.0701.hg18.tile.05.hx1.v2", 
            "pd.0701.hg18.tile.05.hx1.charmed", 
            "pd.charm.hg18.example",
            "pd.081229.hg18.promoter.medip.hx1")) {
        spec="human"; build="hg18"
       } else if(dmr$package %in% 
          c("pd.feinberg.mm8.me.hx1", 
            "pd.2006.07.18.mm8.refseq.promoter")) {
        spec="mouse"; build="mm8"
       } else if (dmr$package %in% 
          c("pd.091020.hg19.charm.hx1", 
            "pd.091117.hg19.charm.hd4.hx1",
            "pd.100428.hg19.charm.hx1")){
        spec="human"; build="hg19"
       } else if (dmr$package %in% 
          c("pd.2006.10.31.rn34.refseq.promoter",
            "pd.090519.rn34.charm.hx1")) {
        spec="rat"; build="rn4"    
       } else { 
         stop("This array design not accomodated yet.") 
       }

if(spec=="human"){
    require(paste("BSgenome.Hsapiens.UCSC.",build,sep=""), character.only=TRUE)
    if(build=="hg18"){
        con <- url("http://rafalab.jhsph.edu/CGI/CGI-Hsapiens-hg18.rda")
    } else if(build=="hg19"){
        con <- url("http://rafalab.jhsph.edu/CGI/CGI-Hsapiens-hg19.rda")
    }
}
if(spec=="mouse"){
    require(paste("BSgenome.Mmusculus.UCSC.",build,sep=""), character.only=TRUE)
    if(build=="mm8"){
        con <- url("http://rafalab.jhsph.edu/CGI/CGI-Mmusculus-mm8.rda")
    } else if(build=="mm9"){
        con <- url("http://rafalab.jhsph.edu/CGI/CGI-Mmusculus-mm9.rda")
    }
}
if(spec=="rat") {
    require(paste("BSgenome.Rnorvegicus.UCSC.",build,sep=""), character.only=TRUE)
    if(build=="rn4"){
        con <- url("http://rafalab.jhsph.edu/CGI/CGI-Rnorvegicus-rn4.rda")
    }
}
load(con)
close(con)
cpg.cur = cgi
ocpgi=data.frame(chr=I(cpg.cur[,1]), start=as.numeric(cpg.cur[,2]), end=as.numeric(cpg.cur[,3]))
chr=dmr$chr
ocpgi=ocpgi[ocpgi$chr%in%unique(as.character(chr)),]

tabs=dmr$tabs
pns=dmr$pns
Indexes=split(seq(along=pns),pns)
pos=dmr$pos

if(is.null(dmr$p)) {
    DD = dmr$DD
    sMD = dmr$sMD
} else { #then user wanted to use p (or rather, logit(p))
    #must do this since can't back-transform to get differrence in p's after the logit transformations.
    DD = get.DD(l=dmr$p, groups=FACS, compare=dmr$args$compare, pairs=dmr$args$pairs)$DD
    sMD = get.tt.paired(DD=DD,Indexes=Indexes,filter=dmr$args$filter,ws=dmr$args$ws)$sMD
}

#################Make output table and pdf file of plots##############
has = which(!sapply(tabs[which.table],is.null))
for(object.i in which.table[has]){
    message("Plotting table ",object.i)
    #################Make output table:###############################
    obj = tabs[[object.i]]
    if(!any(colnames(obj)=="index")) obj$index=match(obj$regionName,names(Indexes))

    ################Make pdf file of plots###############
    pdf(file.path(outpath,paste(names(tabs)[object.i],".pdf",sep="")),width=11,height=8,pointsize=10)
    palette(rev(brewer.pal(8,"Dark2")[c(2,1,3,4,8)]))
    ADD1=1500

    if(all.lines){
        comps.l = has
        cl = colors.l[has]
    }
    if(!all.lines){
        comps.l = object.i
        if(length(colors.l)==ncol(sMD)){
            cl = colors.l[object.i]
        } else{ #length(colors.l)==1
            cl = colors.l
        }
    }

    for(i in intersect(which.plot,1:nrow(obj))){
      cat(i," ")
      if(is.null(dmr$p)){
          YLIM=c(-1,1)
          ylabel="difference in M"
      } else{
          YLIM=c(-1, ifelse(all.lines,1+ncol(sMD)*.03*legend.size,1+.06*legend.size) )
          ylabel="difference in p"
      }
      thechr=obj$chr[i]
      if(is.null(buffer)) ADD=max(ADD1*2,obj$end[i]-obj$start[i]+132) else ADD=buffer
      start=floor((obj$start[i]+obj$end[i])/2)-ADD
      end=start+2*ADD
    
      Index=Indexes[[obj$index[i]]]
      x=pos[Index] 
      start=max(start,min(x))
      end=min(end,max(x))
      
      Index=Index[x>=start & x <= end]
      Index=Index[order(pos[Index])]
      
      ##PLOT TISSUES
      #par(mfrow=c(3,1))
      #layout(matrix(1:3,ncol=1),heights=c(0.6,0.2,0.2))
      layout(matrix(1:2,ncol=1),heights=c(0.6,0.2))
      #par(mar=c(0,2.5,0.25,1.1),oma=c(0,0,2,0))
      par(mar=c(0,3.5,0.25,1.1),oma=c(1.1,0,2,0),mgp=c(2.5,.5,0))

      matplot(pos[Index],DD[[object.i]][Index,],type="n",ylim=YLIM,xlab="",
              xaxt="n",ylab=ylabel,xlim=c(start,end),lty=1, las=1)
      if(all.points){ plotpoints = has } else{ plotpoints = object.i }
      for(ppi in plotpoints){
          matpoints(pos[Index],DD[[ppi]][Index,],col=colors.p[ppi],cex=0.6)
      }
      
      matlines(pos[Index],sMD[Index,comps.l],lwd=2,lty=1,col=cl)
      
      abline(h=0,lwd=2,lty=3)
      abline(v=obj$start[i],lty=2)
      abline(v=obj$end[i],lty=2)
      
      nms=colnames(sMD)
      if(all.points|all.lines) nms[object.i] = paste(nms[object.i],"*")
      legend("topright",nms[comps.l],col=cl,lty=1,lwd=2,cex=legend.size)
      if(ncol(DD[[object.i]]) < 15){  #62 is max possible.
          ids = dmr$args$pairs[match(colnames(DD[[object.i]]),colnames(dmr$l))]
          symbols = c(1:9,0,letters,LETTERS)
          leg2 = paste(symbols[1:length(ids)],"   ","Pair",ids)
          legend("bottomright",leg2,text.col=colors.p[object.i])
      }
      
      ##PLOT CPG ISLANDS
      par(mar=c(3.5,3.5,0.25,1.1))
      if(spec=="human") seq<-Hsapiens[[as.character(thechr) ]]
      if(spec=="mouse") seq<-Mmusculus[[as.character(thechr) ]]
      if(spec=="rat")   seq<-Rnorvegicus[[as.character(thechr) ]]
      Index=max(start,1):min(end,length(seq))
      #subseq<-seq[Index]  #using older version of Biostrings. For newer:
      subseq<-subseq(seq,start=Index[1],end=Index[length(Index)])
      cpgs=start(matchPattern("CG",subseq))+Index[1]-1
      ##mcrbc=sort(c(start(matchPattern("ACG",subseq))+Index[1]-1,
      ##  start(matchPattern("GCG",subseq))+Index[1]-1))
      cuts=seq(min(Index),max(Index),8)
      scpgs=cut(cpgs,cuts,include.lowest=TRUE)
      x = (cuts[1:(length(cuts)-1)]+cuts[2:(length(cuts))])/2
      y = table(scpgs)/8
      SPAN=400/diff(range(x))
      d = loess(y~x,span=SPAN,degree=1)
      plot(x,d$fitted,type="l",ylim=c(0,0.15),xlab="Location",
           ylab="CpG density",xlim=c(start,end), las=1)
      Rug(cpgs)
    ##  Rug(mcrbc,side=3)
      Index1 = which(ocpgi[,1]==as.character(thechr) &
                     ((ocpgi[,2] > start & ocpgi[,2]< end) |
                      (ocpgi[,3] > start & ocpgi[,3]< end)))
      if(length(Index1)>0)
          sapply(Index1,function(j) Rug(unlist(ocpgi[j,2:3]),col=5,lwd=3,side=1))

      mtext(paste("ID:",i,"--",as.character(thechr),":",start,"-",end,sep=""),
            side=3,cex=1.5,outer=TRUE)
    }
    cat("\n")
    dev.off()
}
cat("\nPlotting finished.\n")
}

##A function to plot arbitrary regions of your Martinized charm data:
##In the plots, the left side of the region is demarcated by a dashed line and the right side by a dash-dot line.  3k on either side are also plotted.  If there is less than this in a plot, it means the charm array did not cover the missing part.  In the title, the ID number indicates what row of your input table each plot is for.  If a row of your input table is not on the Charm array, it will be missing from the plots (it's ID number will be skipped).

##Arguments:
#myregions - table with chr,start,end columns for regions you want to see
#dmrpath - dmr object from dmrFinder containing the data you want to plot (the actual dmr tables aren't necessary here though)
#Regions in myregions must be using the same genome assembly version as the dmr object in dmrpath!
#outfile - what you want the name of the output file to be.  Include the full path.
#which.plot - plot these rows from your myregions table
#which.groups - names of groups you want to plot
#cl - colors for the plotted groups in alphabetic order
#legend.size - number of times to magnify the legend
#buffer - how much to plot on either side of each region

regionPlot <- function(tab, dmr, outfile, which.plot, which.groups=colnames(dmr$gm), cl=2:(ncol(dmr$gm)+1), legend.size=1, buffer=3000) {

if(!inherits(buffer,c("numeric","integer"))) stop("buffer must be numeric.")
if(missing(tab)) stop("tab argument must be provided.")
if(!is.data.frame(tab)) stop("tab must be a data frame.")
if(!all(c("chr","start","end")%in%colnames(tab))) stop("tab data frame must have columns labeled chr,start, and end.")
if(missing(dmr)) stop("dmr argument must be provided.")
if(missing(outfile)) stop("outfile argument must be provided. outfile should be the name of the output pdf file (include .pdf extension), including the full path.")
if(!inherits(outfile,"character")) stop("invalid outfile argument")
if(missing(which.plot)) which.plot=1:min(50,nrow(tab))
if(nrow(tab)<min(which.plot)) stop("which.plot argument starts too high given that there are only ",nrow(tab)," regions in the table.")
stopifnot(inherits(which.plot,c("numeric","integer")))
stopifnot(length(cl)==length(which.groups))
stopifnot(all(which.groups%in%colnames(dmr$gm)))
require(RColorBrewer)
tab=tab[which.plot,]
tab$id = 1:nrow(tab)
pns=dmr$pns
Indexes=split(seq(along=pns),pns)
chr=dmr$chr
pos=dmr$pos

if(dmr$package %in% 
          c("pd.feinberg.hg18.me.hx1", 
            "pd.feinberg.hg18.me.hx1.noxy", 
            "pd.0701.hg18.tile.05.hx1.v2", 
            "pd.0701.hg18.tile.05.hx1.charmed", 
            "pd.charm.hg18.example",
            "pd.081229.hg18.promoter.medip.hx1")) {
        spec="human"; build="hg18"
       } else if(dmr$package %in% 
          c("pd.feinberg.mm8.me.hx1", 
            "pd.2006.07.18.mm8.refseq.promoter")) {
        spec="mouse"; build="mm8"
       } else if (dmr$package %in% 
          c("pd.091020.hg19.charm.hx1", 
            "pd.091117.hg19.charm.hd4.hx1",
            "pd.100428.hg19.charm.hx1")){
        spec="human"; build="hg19"
       } else if (dmr$package %in% 
          c("pd.2006.10.31.rn34.refseq.promoter",
            "pd.090519.rn34.charm.hx1")) {
        spec="rat"; build="rn4"    
       } else { 
         stop("This array design not accomodated yet.") 
       }

if(spec=="human"){
    require(paste("BSgenome.Hsapiens.UCSC.",build,sep=""), character.only=TRUE)
    if(build=="hg18"){
        con <- url("http://rafalab.jhsph.edu/CGI/CGI-Hsapiens-hg18.rda")
    } else if(build=="hg19"){
        con <- url("http://rafalab.jhsph.edu/CGI/CGI-Hsapiens-hg19.rda")
    }
}
if(spec=="mouse"){
    require(paste("BSgenome.Mmusculus.UCSC.",build,sep=""), character.only=TRUE)
    if(build=="mm8"){
        con <- url("http://rafalab.jhsph.edu/CGI/CGI-Mmusculus-mm8.rda")
    } else if(build=="mm9"){
        con <- url("http://rafalab.jhsph.edu/CGI/CGI-Mmusculus-mm9.rda")
    }
}
if(spec=="rat") {
    require(paste("BSgenome.Rnorvegicus.UCSC.",build,sep=""), character.only=TRUE)
    if(build=="rn4"){
        con <- url("http://rafalab.jhsph.edu/CGI/CGI-Rnorvegicus-rn4.rda")
    }
}
load(con)
close(con)
cpg.cur = cgi
ocpgi=data.frame(chr=I(cpg.cur[,1]), start=as.numeric(cpg.cur[,2]), end=as.numeric(cpg.cur[,3]))
chr=dmr$chr
ocpgi=ocpgi[ocpgi$chr%in%unique(as.character(chr)),] 

if(is.null(dmr$p)) {
    gmp = dmr$gm[,which.groups]
} else {
    gmp = exp(dmr$gm[,which.groups])/(1+exp(dmr$gm[,which.groups]))
}
sMM=get.smooth(gmp,Indexes=Indexes,filter=NULL,columns=1:ncol(gmp),ws=dmr$args$ws)
if(is.null(dmr$p)) M=dmr$l else M=dmr$p

wh = which(dmr$groups%in%which.groups)
th = as.numeric(factor(rank(dmr$groups[wh])))
comps.l = order(colnames(sMM))

#################Make output pdf file of plots##############
pdf(file=outfile, width=11, height=8, pointsize=10)

palette(rev(brewer.pal(8,"Dark2")[c(2,1,3,4,8)]))
for(i in intersect(which.plot,1:nrow(tab))) {
  cat(i," ")
  if(is.null(dmr$p)){
      YLIM=c(-0.5,3)
      ylabel="M"
  } else{
      YLIM=c(0, 1+ncol(sMM)*.03*legend.size)
      ylabel="p"
  }

  ##PLOT TISSUES:
  thechr = tab$chr[i]
  Index=which(dmr$chr== thechr) #only that one large chr11 region could cause trouble here.
  Index=Index[which(dmr$pos[Index]>=tab$start[i]-buffer & dmr$pos[Index]<=tab$end[i]+buffer)]
  test =Index[which(dmr$pos[Index]>=tab$start[i]        & dmr$pos[Index]<=tab$end[i]    )] 
  if(length(Index)==0 | length(test)==0) next
  start  = min(dmr$pos[Index])
  end    = max(dmr$pos[Index])

  #par(mfrow=c(3,1))
  #layout(matrix(1:3,ncol=1),heights=c(0.6,0.2,0.2))
  layout(matrix(1:2,ncol=1),heights=c(0.6,0.2))
  #par(mar=c(0,2.5,0.25,1.1),oma=c(0,0,2,0))
  par(mar=c(0,3.5,0.25,1.1),oma=c(1.1,0,2,0),mgp=c(2.5,.5,0))

  matplot(dmr$pos[Index], M[Index,wh], col=cl[th], lty=1, ylab=ylabel, ylim=YLIM, cex=.6, xlab="position", xlim=c(start,end), xaxt="n", las=1)
  matlines(dmr$pos[Index],sMM[Index,comps.l],lwd=2,lty=1,col=cl)

  abline(h=0,lwd=2,lty=3)
  abline(h=1,lwd=2,lty=3)
  abline(v=tab$start[i],lty=2)
  abline(v=tab$end[i], ,lty=4)
  legend("topright",colnames(sMM)[comps.l],col=cl,lty=1,lwd=2,cex=legend.size)

      ##PLOT CPG ISLANDS
      par(mar=c(3.5,3.5,0.25,1.1))
      if(spec=="human") seq<-Hsapiens[[as.character(thechr) ]]
      if(spec=="mouse") seq<-Mmusculus[[as.character(thechr) ]]
      if(spec=="rat")   seq<-Rnorvegicus[[as.character(thechr) ]]
      Index=max(start,1):min(end,length(seq))
      #subseq<-seq[Index]  #using older version of Biostrings. For newer:
      subseq<-subseq(seq,start=Index[1],end=Index[length(Index)])
      cpgs=start(matchPattern("CG",subseq))+Index[1]-1
      ##mcrbc=sort(c(start(matchPattern("ACG",subseq))+Index[1]-1,
      ##  start(matchPattern("GCG",subseq))+Index[1]-1))
      cuts=seq(min(Index),max(Index),8)
      scpgs=cut(cpgs,cuts,include.lowest=TRUE)
      x = (cuts[1:(length(cuts)-1)]+cuts[2:(length(cuts))])/2
      y = table(scpgs)/8
      SPAN=400/diff(range(x))
      d = loess(y~x,span=SPAN,degree=1)
      plot(x,d$fitted,type="l",ylim=c(0,0.15),xlab="Location",
           ylab="CpG density",xlim=c(start,end), las=1)
      Rug(cpgs)
    ##  Rug(mcrbc,side=3)
      Index1 = which(ocpgi[,1]==as.character(thechr) &
                     ((ocpgi[,2] > start & ocpgi[,2]< end) |
                      (ocpgi[,3] > start & ocpgi[,3]< end)))
      if(length(Index1)>0)
          sapply(Index1,function(j) Rug(unlist(ocpgi[j,2:3]),col=5,lwd=3,side=1))
      
        mtext(paste("ID:",tab$id[i],"--",as.character(thechr),":",tab$start[i],"-",tab$end[i],sep=""),
              side=3,cex=1.5,outer=TRUE)
    }
  dev.off()
  cat("\nPlotting finished.\n")
}


#####################################################
####Functions needed by the above plotting functions:

mypar <-function(a=1,b=1,brewer.n=8,brewer.name="Dark2",...){
  require(RColorBrewer)
  par(mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0)) 
  par(mfrow=c(a,b),...)
  palette(brewer.pal(brewer.n,brewer.name)) 
}

get.smooth <- function(lm,filter=NULL,Indexes,ws=3,verbose=TRUE,columns=1:ncol(lm)){
##this function assumes genome position of lm[Indexes[[1],column] are ordered:
  if(is.null(filter)){
    Tukey = function(x) pmax(1 - x^2,0)^2
    filter= Tukey(seq(-ws,ws)/(ws+1));filter=filter/sum(filter)
  }
  if(!isTRUE(all.equal(sum(filter),1))) stop("filter must sum to 1.")
  
  slm=matrix(NA,nrow(lm),ncol(lm))
  colnames(slm) = colnames(lm)
  if(verbose) cat("Smoothing:\n")
  if(verbose) pb = txtProgressBar(min=1,max=length(Indexes),initial=0,style=3)
  for(i in seq(along=Indexes)){
    #if(verbose) if(i%%1000==0) cat(".")
    Index=Indexes[[i]]
    for(r in columns){
        slm[Index,r]=myfilter(lm[Index,r],filter)
    }
    if(verbose) setTxtProgressBar(pb,i)
  }
  if(verbose){ close(pb); cat("Done.\n") }
  return(slm)
}

#This is the rug function from before R-2.8.0, which colors the axis between ticks:
Rug <- function (x, ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"), 
    quiet = getOption("warn") < 0, ...) {
    x <- as.vector(x)
    ok <- is.finite(x)
    x <- x[ok]
    if (!quiet) {
        u <- par("usr")
        u <- if (side%%2 == 1) {
            if (par("xlog")) 
                10^u[1:2]
            else u[1:2]
        }
        else {
            if (par("ylog")) 
                10^u[3:4]
            else u[3:4]
        }
        if (any(x < u[1] | x > u[2])) 
            warning("some values will be clipped")
    }
    Axis(side = side, at = x, labels = FALSE, lwd = lwd, col = col, 
        tck = ticksize, ...)
}
