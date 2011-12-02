##This file contains both dmrPlot and dmrPlot_paired functions.

dmrPlot <- function(dmr, which.table=1:length(dmr$tabs), which.plot=1:30, legend.size=1, all.lines=TRUE, all.points=FALSE, colors.l, colors.p, outpath=".", cpg.islands, Genome, mcrbc=TRUE, plotG=FALSE, dat=NULL, buffer=NULL, plot.p=TRUE) {
  
    if(is.null(dmr$DD)) {
        if(missing(colors.p)) colors.p=1+1:length(unique(dmr$groups))
        if(missing(colors.l)) colors.l=1+1:length(unique(dmr$groups))       
        dmrPlot_unpaired(dmr=dmr, which.table=which.table, which.plot=which.plot, legend.size=legend.size, all.lines=all.lines, all.points=all.points, colors.l=colors.l, colors.p=colors.p, outpath=outpath, cpg.islands=cpg.islands, Genome=Genome, mcrbc=mcrbc, plotG=plotG, dat=dat, buffer=buffer, plot.p=plot.p)
    } else {
        if(missing(colors.p)) colors.p = 1+1:ncol(dmr$sMD)
        if(missing(colors.l)) colors.l = 1+1:ncol(dmr$sMD)
        dmrPlot_paired(dmr=dmr, which.table=which.table, which.plot=which.plot, legend.size=legend.size, all.lines=all.lines, all.points=all.points, colors.l=colors.l, colors.p=colors.p, outpath=outpath, cpg.islands=cpg.islands, Genome=Genome, mcrbc=mcrbc, plotG=plotG, dat=dat, buffer=buffer, plot.p=plot.p)
    }
    
}

#This function is just results_function.R except there is no option to return the table(s) and it does not plot the gene panel.
dmrPlot_unpaired <- function(dmr, which.table=1:length(dmr$tabs), which.plot=1:30, legend.size=1, all.lines=TRUE, all.points=FALSE, colors.l, colors.p, outpath=".", cpg.islands, Genome, mcrbc=TRUE, plotG=FALSE, dat=NULL, buffer=NULL, plot.p=plot.p) {

#require(RColorBrewer)
#require(Biobase) #for rowMedians

if(plotG & is.null(dat)) stop("If plotG=TRUE, you must provide dat.")
stopifnot(inherits(which.plot,c("numeric","integer")))
stopifnot(inherits(which.table,c("numeric","integer")))
if(length(dmr$tabs)<max(which.table)) stop("which.table specifies more comparisons than were made.")
if(!is.null(buffer) & !inherits(buffer,c("numeric","integer"))) stop("if provided, buffer must be numeric.")

if(all.points)  stopifnot(length(colors.p)==length(unique(dmr$groups)))
if(all.lines)   stopifnot(length(colors.l)==length(unique(dmr$groups)))
if(!all.points) stopifnot(length(colors.p)%in%c(2, length(unique(dmr$groups))))
if(!all.lines)  stopifnot(length(colors.l)%in%c(2, length(unique(dmr$groups))))
FACS=dmr$groups

if(plotG) {
    #require(oligo)
    G = oligo:::pm(dat)[,,1]
    seq = pmSequence(dat)
    stopifnot(length(dmr$groups)==ncol(G))
    if(is.null(dmr$p)) { stopifnot(all(colnames(dmr$l)==sampleNames(dat)))
    } else stopifnot(all(colnames(dmr$p)==sampleNames(dat)))

    #require(Biobase) #for rowMedians
    G = log2(preprocessCore::normalize.quantiles(G))
    # remove probe effects from G:
    G = G - rowMedians(G)
    # GC sequence effect correction:
    G = remove_gc_effect(G, seq)
    G = G[dmr$index,]
}

chr=dmr$chr
ocpgi=cpg.islands[cpg.islands$chr%in%unique(as.character(chr)),]

tabs=dmr$tabs
pns=dmr$pns
Indexes=split(seq(along=pns),pns)
pos=dmr$pos
#In dmrFinder dmr$gm only gets columns for groups that are in >=1 comparison, so:
all.compared = all(unique(FACS)%in%colnames(dmr$gm))
if(is.null(dmr$p) | !plot.p){
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
if(is.null(dmr$p) | !plot.p) M=dmr$l else M=dmr$p

for(h in 1:nrow(dmr$comps)) stopifnot(identical(paste(colnames(dmr$gm)[dmr$comps[h,]],collapse="-"), names(tabs)[h]))

#################Make output table and pdf file of plots##############
has = which(!sapply(tabs[which.table],is.null))
for(object.i in which.table[has]){
    message("Plotting regions in table ",object.i)
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

    if(plotG) {
        mes = data.frame(gr1 = rowMedians(G[,which(FACS==grs[1])]), 
                         gr2 = rowMedians(G[,which(FACS==grs[2])]))
        mes = mes[,1]-mes[,2]
        stopifnot(length(mes)==length(pos)) #pos is dmr$pos, defined above
    }

    for(i in intersect(which.plot,1:nrow(obj))){
      message("Plotting region ",i)
      if(is.null(dmr$p) | !plot.p) {
          YLIM=c(-0.5,3)
          ylabel="M"
      } else {
          YLIM=c(0, ifelse(all.lines,1+ncol(sMM)*.03*legend.size,1+.06*legend.size) )
          ylabel="p"
      }
      thechr=obj$chr[i]
      if(is.null(buffer)) ADD=max(ADD1*2,obj$end[i]-obj$start[i]+132) else ADD=buffer
      start=floor((obj$start[i]+obj$end[i])/2)-ADD
      end  =start+2*ADD
    
      Index=Indexes[[obj$index[i]]]
      x=pos[Index]
      start=max(start,min(x))
      end  =min(end,max(x))
      
      Index=Index[x>=start & x <= end]
      Index=Index[order(pos[Index])]
      start=min(pos[Index])
      end  =max(pos[Index])
      
      ##PLOT METHYLATION:
      if(plotG) {
          layout(matrix(1:3,ncol=1),heights=c(0.6,0.2,0.2))
      } else {
          layout(matrix(1:2,ncol=1),heights=c(0.6,0.2))
      }
      par(mar=c(0,3.5,0.25,1.1),oma=c(0,0,2,0),mgp=c(2.5,.5,0))

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
      if(plotG) {
          plot_CpG(thechr=thechr, xx=start:end, ocpgi=ocpgi, Genome=Genome, mcrbc=mcrbc)
      } else {
          par(mar=c(3.5,3.5,0.25,1.1))
          plot_CpG(thechr=thechr, xx=start:end, ocpgi=ocpgi, Genome=Genome, mcrbc=mcrbc, xlab="Location", xaxt=NULL)
      }

      if(plotG) { ##Plot difference in G
          par(mar=c(3.5,3.5,0.25,1.1))
          plot(pos[Index], mes[Index], ylab="Difference in median G", xlab="Location",
               ylim=c(-2.5,2.5), type="b")
          abline(h=0,lty=3)
      }

      mtext(paste("ID:",i,"--",as.character(thechr),":",start,"-",end,sep=""),
            side=3,cex=1.5,outer=TRUE)
    }
    cat("\n")
    dev.off()
}
cat("\nPlotting finished.\n")
}


dmrPlot_paired <- function(dmr, which.table=1:length(dmr$tabs), which.plot=1:30, legend.size=1, all.lines=TRUE, all.points=FALSE, colors.l, colors.p, outpath=".", cpg.islands, Genome, mcrbc, plotG, dat, buffer, plot.p) {
#NB: all.points and all.lines refer to all comparisons that were made, not that could be made!
  
#require(RColorBrewer)
#require(Biobase) #for rowMedians

if(plotG & is.null(dat)) stop("If plotG=TRUE, must provide dat.")
stopifnot(inherits(which.plot,c("numeric","integer")))
stopifnot(inherits(which.table,c("numeric","integer")))
if(length(dmr$tabs)<max(which.table)) stop("which.table specifies more comparisons than were made.")
if(!is.null(buffer) & !inherits(buffer,c("numeric","integer"))) stop("if provided, buffer must be numeric.")

if(!all.points & length(colors.p)==1) colors.p=rep(colors.p,ncol(dmr$sMD))
stopifnot(length(colors.p)==ncol(dmr$sMD))
if(all.lines)  stopifnot(length(colors.l)==ncol(dmr$sMD))
if(!all.lines) stopifnot(length(colors.l)%in%c(1, ncol(dmr$sMD)))
sMD = dmr$sMD
FACS = dmr$groups

if(plotG) {
    #require(oligo)
    G = oligo:::pm(dat)[,,1]
    seq = pmSequence(dat)
    stopifnot(length(dmr$groups)==ncol(G))
    if(is.null(dmr$p)) { stopifnot(all(colnames(dmr$l)==sampleNames(dat)))
    } else stopifnot(all(colnames(dmr$p)==sampleNames(dat)))

    #require(Biobase) #for rowMedians
    G = log2(preprocessCore::normalize.quantiles(G))
    # remove probe effects from G:
    G = G - rowMedians(G)
    # GC sequence effect correction:
    G = remove_gc_effect(G, seq)
    G = G[dmr$index,]
}

chr=dmr$chr
ocpgi=cpg.islands[cpg.islands$chr%in%unique(as.character(chr)),]

tabs=dmr$tabs
pns=dmr$pns
Indexes=split(seq(along=pns),pns)
pos=dmr$pos

if(is.null(dmr$p) | !plot.p) {
    DD = dmr$DD
    sMD = dmr$sMD
} else { #then user wanted to use p (or rather, logit(p))
    #must do this since can't back-transform to get differrence in p's after the logit transformations.
    DD = get.DD(l=dmr$p, groups=FACS, compare=dmr$args$compare, pairs=dmr$args$pairs)$DD
    sMD = get.tt.paired(DD=DD,Indexes=Indexes,filter=dmr$args$filter,ws=dmr$args$ws)$sMD
}

#################Make output table and pdf file of plots##############
has = which(!sapply(tabs[which.table],is.null))
for(object.i in which.table[has]) {
    message("Plotting table ",object.i)
    #################Make output table:###############################
    obj = tabs[[object.i]]
    if(!any(colnames(obj)=="index")) obj$index=match(obj$regionName,names(Indexes))

    ################Make pdf file of plots###############
    pdf(file.path(outpath,paste(names(tabs)[object.i],"_paired.pdf",sep="")), width=11, height=8, pointsize=10)
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

    if(plotG) {
        grs = strsplit(colnames(sMD)[object.i],"-")[[1]]
        mes = data.frame(gr1 = rowMedians(G[,which(FACS==grs[1])]), 
                         gr2 = rowMedians(G[,which(FACS==grs[2])]))
        mes = mes[,1]-mes[,2]
        stopifnot(length(mes)==length(pos)) #pos is dmr$pos, defined above
    }

    for(i in intersect(which.plot,1:nrow(obj))){
      message("Plotting region ",i)
      if(is.null(dmr$p) | !plot.p){
          YLIM=c(-6,6)
          ylabel="difference in M"
      } else {
          YLIM=c(-1, ifelse(all.lines,1+ncol(sMD)*.03*legend.size,1+.06*legend.size) )
          ylabel="difference in p"
      }
      thechr=obj$chr[i]
      if(is.null(buffer)) ADD=max(ADD1*2,obj$end[i]-obj$start[i]+132) else ADD=buffer
      start=floor((obj$start[i]+obj$end[i])/2)-ADD
      end  =start+2*ADD
    
      Index=Indexes[[obj$index[i]]]
      x=pos[Index] 
      start=max(start,min(x))
      end  =min(end,max(x))
      
      Index=Index[x>=start & x <= end]
      Index=Index[order(pos[Index])]
      start=min(pos[Index])
      end  =max(pos[Index])
      
      ##PLOT METHYLATION:
      if(plotG) {
          layout(matrix(1:3,ncol=1),heights=c(0.6,0.2,0.2))
      } else {
          layout(matrix(1:2,ncol=1),heights=c(0.6,0.2))
      }
      #par(mar=c(0,2.5,0.25,1.1),oma=c(0,0,2,0))
      par(mar=c(0,3.5,0.25,1.1),oma=c(1.1,0,2,0),mgp=c(2.5,.5,0))

      matplot(pos[Index],DD[[object.i]][Index,],type="n",ylim=YLIM,xlab="",
              xaxt="n",ylab=ylabel,xlim=c(start,end), lty=1, las=1)
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
      if(plotG) {
          plot_CpG(thechr=thechr, xx=start:end, ocpgi=ocpgi, Genome=Genome, mcrbc=mcrbc)
      } else {
          par(mar=c(3.5,3.5,0.25,1.1))
          plot_CpG(thechr=thechr, xx=start:end, ocpgi=ocpgi, Genome=Genome, mcrbc=mcrbc, xlab="Location", xaxt=NULL)
      }

      if(plotG) { ##Plot difference in G
          par(mar=c(3.5,3.5,0.25,1.1))
          plot(pos[Index], mes[Index], ylab="Difference in median G", xlab="Location",
               ylim=c(-2.5,2.5), type="b")
          abline(h=0,lty=3)
      }

      mtext(paste("ID:",i,"--",as.character(thechr),":",start,"-",end,sep=""),
            side=3,cex=1.5,outer=TRUE)
    }
    cat("\n")
    dev.off()
}
cat("\nPlotting finished.\n")
}



regionPlot <- function(tab, dmr, cpg.islands, Genome, outfile, which.plot=1:10, plot.these, cl, legend.size=1, buffer=3000, mcrbc=TRUE, plot.p=TRUE, plotG=FALSE, dat=NULL, grs=NULL) {
  
    if(is.null(dmr$DD)) {
        if(missing(plot.these)) plot.these = colnames(dmr$gm)
        if(missing(cl))         cl = 1+1:ncol(dmr$gm)       
        regionPlot_unpaired(tab=tab, dmr=dmr, cpg.islands=cpg.islands, Genome=Genome, outfile=outfile, which.plot = which.plot, plot.these=plot.these, cl=cl, legend.size=legend.size, buffer=buffer, mcrbc=mcrbc, plot.p=plot.p, plotG=plotG, dat=dat, grs=grs)
    } else {
        if(missing(plot.these)) plot.these = colnames(dmr$sMD)
        if(missing(cl))         cl = 1+1:length(plot.these)
        regionPlot_paired(tab=tab, dmr=dmr, cpg.islands=cpg.islands, Genome=Genome, outfile=outfile, which.plot = which.plot, plot.these=plot.these, cl=cl, legend.size=legend.size, buffer=buffer, mcrbc=mcrbc, plot.p=plot.p, plotG=plotG, dat=dat, grs=grs)
    }
    
}
  
regionPlot_unpaired <- function(tab, dmr, cpg.islands, Genome, outfile, which.plot=1:10, plot.these=colnames(dmr$gm), cl=1+1:ncol(dmr$gm), legend.size=1, buffer=3000, mcrbc=TRUE, plot.p=TRUE, plotG=FALSE, dat=NULL, grs=NULL) {

#require(RColorBrewer)

if(plotG & (is.null(dat) | is.null(grs))) stop("If plotG=TRUE, must provide dat and grs arguments.")
if(plotG) stopifnot(all(grs%in%dmr$groups))

if(!inherits(buffer,c("numeric","integer"))) stop("buffer must be a number.")
if(!is.data.frame(tab)) stop("tab must be a data frame.")
if(!all(c("chr","start","end")%in%colnames(tab))) stop("tab data frame must have columns labeled chr,start, and end.")
if(missing(outfile)) stop("outfile argument must be provided. outfile should be the name of the output pdf file (include .pdf extension), including the full path.")
if(!inherits(outfile,"character")) stop("invalid outfile argument")
stopifnot(inherits(which.plot,c("numeric","integer")))

stopifnot(length(cl)==length(plot.these))
stopifnot(all(plot.these%in%colnames(dmr$gm)))

tab$id = 1:nrow(tab)
tab = tab[,c("chr","start","end","id")]
Indexes=split(seq(along=dmr$pns),dmr$pns)

if(plotG) {
    #require(oligo)
    G = oligo:::pm(dat)[,,1]
    seq = pmSequence(dat)
    stopifnot(length(dmr$groups)==ncol(G))
    if(is.null(dmr$p)) { stopifnot(all(colnames(dmr$l)==sampleNames(dat)))
    } else stopifnot(all(colnames(dmr$p)==sampleNames(dat)))

    #require(Biobase) #for rowMedians
    G = log2(preprocessCore::normalize.quantiles(G))
    # remove probe effects from G:
    G = G - rowMedians(G)
    # GC sequence effect correction:
    G = remove_gc_effect(G, seq)
    G = G[dmr$index,]
}

ocpgi=cpg.islands[cpg.islands$chr%in%unique(as.character(dmr$chr)),]

if(is.null(dmr$p) | !plot.p) {
    gmp = dmr$gm[,plot.these]
} else {
    gmp = exp(dmr$gm[,plot.these])/(1+exp(dmr$gm[,plot.these]))
}
sMM=get.smooth(gmp,Indexes=Indexes,filter=NULL,columns=1:ncol(gmp),ws=dmr$args$ws)
if(is.null(dmr$p) | !plot.p) M=dmr$l else M=dmr$p

wh = which(dmr$groups%in%plot.these)
th = as.numeric(factor(rank(dmr$groups[wh])))
comps.l = order(colnames(sMM))

#################Make output pdf file of plots##############

if(is.null(dmr$p) | !plot.p) {
    YLIM=c(-0.5,3)
    ylabel="M"
} else {
    YLIM=c(0, 1+ncol(sMM)*.03*legend.size)
    ylabel="p"
}
palette(rev(brewer.pal(8,"Dark2")[c(2,1,3,4,8)]))

if(plotG) {
    mes = data.frame(gr1 = rowMedians(G[,which(dmr$groups==grs[1])]), 
                     gr2 = rowMedians(G[,which(dmr$groups==grs[2])]))
    mes = mes[,1]-mes[,2]
    stopifnot(length(mes)==length(dmr$pos))
}

pdf(file=outfile, width=11, height=8, pointsize=10)
for(i in intersect(which.plot,1:nrow(tab))) {
  cat(i," ")

  ##PLOT METHYLATION:
  thechr = tab$chr[i]
  Index=which(dmr$chr== thechr)
  Index=Index[which(dmr$pos[Index]>=tab$start[i]-buffer & dmr$pos[Index]<=tab$end[i]+buffer)]
  test =Index[which(dmr$pos[Index]>=tab$start[i]        & dmr$pos[Index]<=tab$end[i]    )] 
  if(length(Index)==0 | length(test)==0) next
  start  = min(dmr$pos[Index])
  end    = max(dmr$pos[Index])

  if(plotG) {
      layout(matrix(1:3,ncol=1),heights=c(0.6,0.2,0.2))
  } else {
      layout(matrix(1:2,ncol=1),heights=c(0.6,0.2))
  }
  #par(mar=c(0,2.5,0.25,1.1),oma=c(0,0,2,0))
  par(mar=c(0,3.5,0.25,1.1),oma=c(1.1,0,2,0),mgp=c(2.5,.5,0))

  matplot(dmr$pos[Index], M[Index,wh], col=cl[th], lty=1, ylab=ylabel, ylim=YLIM, cex=0.6, xlab="position", xlim=c(start,end), xaxt="n", las=1)
  matlines(dmr$pos[Index],sMM[Index,comps.l],lwd=2,lty=1,col=cl)
  abline(h=0,lwd=2, lty=3)
  abline(h=1,lwd=2, lty=3)
  abline(v=tab$start[i], lty=2)
  abline(v=tab$end[i], lty=4)
  legend("topright",colnames(sMM)[comps.l],col=cl,lty=1,lwd=2,cex=legend.size)

  ##PLOT CPG ISLANDS
  if(plotG) {
      plot_CpG(thechr=thechr, xx=start:end, ocpgi=ocpgi, Genome=Genome, mcrbc=mcrbc)
  } else {
      par(mar=c(3.5,3.5,0.25,1.1))
      plot_CpG(thechr=thechr, xx=start:end, ocpgi=ocpgi, Genome=Genome, mcrbc=mcrbc, xlab="Location", xaxt=NULL)
  }
      
  if(plotG) { ##Plot difference in G
      par(mar=c(3.5,3.5,0.25,1.1))
      plot(dmr$pos[Index], mes[Index], ylab="Difference in median G", xlab="Location",
           ylim=c(-2.5,2.5), type="b")
      abline(h=0,lty=3)
  }

  mtext(paste("ID:",tab$id[i],"--",as.character(thechr),":",tab$start[i],"-",tab$end[i],sep=""),
        side=3,cex=1.5,outer=TRUE)
  }
  dev.off()
  cat("\nPlotting finished.\n")
}


regionPlot_paired <- function(tab, dmr, cpg.islands, Genome=Genome, outfile, which.plot = 1:10, plot.these=colnames(dmr$sMD), cl=1+1:length(plot.these), legend.size=1, buffer=3000, plot.p=TRUE, mcrbc=TRUE, plotG=FALSE, dat=NULL, grs=NULL) {

#require(RColorBrewer)

if(plotG & (is.null(dat) | is.null(grs))) stop("If plotG=TRUE, must provide dat and grs arguments.")
if(plotG) stopifnot(all(grs%in%dmr$groups))

if(!inherits(buffer,c("numeric","integer"))) stop("buffer must be a number.")
if(!is.data.frame(tab)) stop("tab must be a data frame.")
if(!all(c("chr","start","end")%in%colnames(tab))) stop("tab data frame must have columns labeled chr,start, and end.")
if(missing(outfile)) stop("outfile argument must be provided. outfile should be the name of the output pdf file (include .pdf extension), including the full path.")
if(!inherits(outfile,"character")) stop("invalid outfile argument")
stopifnot(inherits(which.plot,c("numeric","integer")))

stopifnot(length(cl)==length(plot.these))
if(is.character(plot.these)) {
    stopifnot(all(plot.these%in%colnames(dmr$sMD)))
} else if(is.numeric(plot.these)|is.integer(plot.these)) {
    stopifnot(all(plot.these<=ncol(dmr$sMD)))
} else stop("invalid plot.these argument")

tab$id = 1:nrow(tab)
tab = tab[,c("chr","start","end","id")]
Indexes=split(seq(along=dmr$pns),dmr$pns)

if(plotG) {
    #require(oligo)
    G = oligo:::pm(dat)[,,1]
    seq = pmSequence(dat)
    stopifnot(length(dmr$groups)==ncol(G))
    if(is.null(dmr$p)) { stopifnot(all(colnames(dmr$l)==sampleNames(dat)))
    } else stopifnot(all(colnames(dmr$p)==sampleNames(dat)))

    #require(Biobase) #for rowMedians
    G = log2(preprocessCore::normalize.quantiles(G))
    # remove probe effects from G:
    G = G - rowMedians(G)
    # GC sequence effect correction:
    G = remove_gc_effect(G, seq)
    G = G[dmr$index,]
}

ocpgi=cpg.islands[cpg.islands$chr%in%unique(as.character(dmr$chr)),]

if(is.null(dmr$p) | !plot.p) {
    DD = dmr$DD
    sMD = dmr$sMD
} else { #then user wanted to use p (or rather, logit(p))
    #must do this since can't back-transform to get differrence in p's after the logit transformations.
    DD = get.DD(l=dmr$p, groups=dmr$groups, compare=dmr$args$compare, pairs=dmr$args$pairs)$DD
    sMD = get.tt.paired(DD=DD,Indexes=Indexes,filter=dmr$args$filter,ws=dmr$args$ws)$sMD
}

#################Make output pdf file of plots##############

if(is.null(dmr$p) | !plot.p){
    YLIM=c(-6,6)
    ylabel="difference in M"
} else{
    YLIM=c(-1, 1+length(plot.these)*.03*legend.size)
    ylabel="difference in p"
}
palette(rev(brewer.pal(8,"Dark2")[c(2,1,3,4,8)]))

if(plotG) {
    mes = data.frame(gr1 = rowMedians(G[,which(dmr$groups==grs[1])]), 
                     gr2 = rowMedians(G[,which(dmr$groups==grs[2])]))
    mes = mes[,1]-mes[,2]
    stopifnot(length(mes)==length(dmr$pos))
}

pdf(outfile, width=11, height=8, pointsize=10)
for(i in intersect(which.plot,1:nrow(tab))) {
  cat(i," ")

  ##PLOT METHYLATION:
  thechr = tab$chr[i]
  Index=which(dmr$chr==thechr)
  Index=Index[which(dmr$pos[Index]>=tab$start[i]-buffer & dmr$pos[Index]<=tab$end[i]+buffer)]
  test =Index[which(dmr$pos[Index]>=tab$start[i]        & dmr$pos[Index]<=tab$end[i]       )] 
  if(length(Index)==0 | length(test)==0) next
  start = min(dmr$pos[Index])
  end   = max(dmr$pos[Index])

  if(plotG) {
      layout(matrix(1:3,ncol=1),heights=c(0.6,0.2,0.2))
  } else {
      layout(matrix(1:2,ncol=1),heights=c(0.6,0.2))
  }
  par(mar=c(0,3.5,0.25,1.1),oma=c(0,0,2,0),mgp=c(2.5,.5,0))

  matplot(dmr$pos[Index],DD[[1]][Index,],type="n",ylim=YLIM,xlab="",
          xaxt="n",ylab=ylabel,xlim=c(start,end),lty=1, las=1)
  if(is.character(plot.these)) wh = match(plot.these,names(DD)) else wh = plot.these
  pcl = rep(NA,length(DD)); pcl[wh] = cl
  for(ppi in wh){
      matpoints(dmr$pos[Index],DD[[ppi]][Index,],col=pcl[ppi],cex=0.6)
  }  
  matlines(dmr$pos[Index],sMD[Index,plot.these],lwd=2,lty=1,col=cl)
  abline(h=0,lwd=2,lty=3)
  abline(v=tab$start[i],lty=2)
  abline(v=tab$end[i],lty=2)
  legend("topright",colnames(sMD)[wh],col=cl,lty=1,lwd=2,cex=legend.size)

  ##PLOT CPG ISLANDS
  if(plotG) {
      plot_CpG(thechr=thechr, xx=start:end, ocpgi=ocpgi, Genome=Genome, mcrbc=mcrbc)
  } else {
      par(mar=c(3.5,3.5,0.25,1.1))
      plot_CpG(thechr=thechr, xx=start:end, ocpgi=ocpgi, Genome=Genome, mcrbc=mcrbc, xlab="Location", xaxt=NULL)
  }

  if(plotG) { ##Plot difference in G
      par(mar=c(3.5,3.5,0.25,1.1))
      plot(dmr$pos[Index], mes[Index], ylab="Difference in median G", xlab="Location",
           ylim=c(-2.5,2.5), type="b")
      abline(h=0,lty=3)
  }

  mtext(paste("ID:",i,"--",as.character(thechr),":",start,"-",end,sep=""),
        side=3, cex=1.5, outer=TRUE)
  }
dev.off()
unlink(c("tmp.txt","tmp.cod"))
cat("\n Plotting finished.\n")
}


#####################################################
####Functions needed by the above plotting functions:

## plot_CpGs is in pipeline1.R

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
  if(verbose) pb = txtProgressBar(min=1,max=length(Indexes),initial=0,style=1)
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

#####################################################
#### Additional functions that produce plots:
##Plot distribution of control and non-control probes, to check that they worked:
controlQC <- function(rawData,controlProbes=NULL,controlIndex=NULL,IDcol,expcol,ylimits=c(-6,8),outfile="./boxplots_check.pdf",height=7,width=9) {
    p2 = methp(rawData, betweenSampleNorm="none", withinSampleNorm=FALSE, scale=FALSE, returnM=TRUE, 
               controlProbes=controlProbes, controlIndex=controlIndex) 
    if(is.null(controlIndex)) {
        if(is.null(controlProbes)) stop("if is.null(controlIndex), controlProbes must be provided.")
        cont = getControlIndex(rawData, controlProbes=controlProbes)
    } else cont = controlIndex
    #stopifnot(identical(cont, which(getContainer(rawData)%in%controlProbes)))
    stopifnot(length(getContainer(rawData))==nrow(p2))
    ord = order(pData(rawData)[,expcol])
    #tmp = sapply(strsplit(pData(rawData)[,"arrayUT"],"_"), function(x) x[1])
    stopifnot(all(colnames(p2[,ord])==pData(rawData)[ord,IDcol]))

    medm = apply(p2[-cont,ord],2,median)
    medc = apply(p2[cont,ord],2,median)
    medd = medm-medc

    pdf(file=outfile, width=width, height=height)
    par(mfrow=c(3,1))
    boxplot(p2[-cont,ord], outline=FALSE, names=pData(rawData)[ord,IDcol], las=3,
            ylab="M (no normalization)", main="non-control probes", cex.lab=1.5, ylim=ylimits)
    abline(h=0,lty=3)
    boxplot(p2[cont,ord], outline=FALSE, xaxt="n", las=3, ylab="M (no normalization)",
            main="control probes", cex.lab=1.5, ylim=ylimits)
    abline(h=0,lty=3)
    plot(medd~c(1:ncol(p2[,ord])),main="median difference",xlab="",ylab="",las=3)
    abline(h=0,lty=3)
    dev.off()
}

cmdsplot <- function(labcols, expcol, rawData, p, okqc=1:nrow(p), noXorY=TRUE, outfile="./cmds_topN.pdf", topN=c(100000,1000)) {
    stopifnot(expcol%in%colnames(pData(rawData)))
    if(missing(labcols)) labcols = 1+1:length(unique(pData(rawData)[,expcol]))

    stopifnot(length(labcols)>=length(unique(pData(rawData)[,expcol])))
    stopifnot(ncol(p)==nrow(pData(rawData)))
    stopifnot(all(colnames(p)==rownames(pData(rawData))))

    require(genefilter) #for rowSds()
    thechr = pmChr(rawData)
    thepos = pmPosition(rawData)
    stopifnot(length(thechr)==length(pmindex(rawData)))
    stopifnot(length(thechr)==nrow(p))

    thechr = thechr[okqc]
    thepos = thepos[okqc]
    p = p[okqc,]
    if(noXorY) p3=p[!thechr%in%c("chrX","chrY"),] else p3=p
    
    v = rowSds(p3)
    o = order(v, decreasing=TRUE)
    lcol = labcols[as.numeric(factor(rank(pData(rawData)[,expcol])))]
    subt = paste("Does",ifelse(noXorY," not "," "),"use probes in sex chromosomes.",sep="")
    pdf(file=outfile, width=7, height=7)
    for(k in topN) {
        message(k)
        if(k=="all") k = nrow(p3) else k=as.numeric(k)
        xlabel = paste("Using the",k,"most variable probes, out of",nrow(p3),"total.")
        d = dist(t(p3[o[1:k],]))
        cmds = cmdscale(d,k=2)
        plot(cmds[,1], cmds[,2], xlab = xlabel, ylab = "", main = "cMDS", col=lcol, sub=subt)
        tmp = tapply(lcol,pData(rawData)[,expcol],unique)
        legend("topleft",pch=1,legend=names(tmp),col=tmp)
    }    
    dev.off()
}
