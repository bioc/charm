
dmrFind <- function(p=NULL, logitp=NULL, svs=NULL, mod, mod0, coeff, pns, chr, pos, only.cleanp=FALSE, only.dmrs=FALSE, rob=TRUE, use.limma=FALSE, smoo="weighted.loess", k=3, SPAN=300, DELTA=36, use="sbeta", Q=0.99, min.probes=3, min.value=0.075, keepXY=TRUE, sortBy="area.raw", verbose=TRUE) {
    if(only.cleanp & only.dmrs) stop("only.cleanp and only.dmrs cannot both be TRUE.")
    if(!sortBy%in%c("area","area.raw","avg","max")) stop("sortBy must be area, area.raw, avg, or max.")
    if(min.value<0 | min.value>1) stop("min.value must be between 0 and 1.")
    if(Q<0 | Q>1) stop("Q must be between 0 and 1.")

    for(j in 1:ncol(mod0)) if(!isTRUE(all.equal(mod0[,j],mod[,j]))) stop("This function requires that all columns of mod not in mod0 (i.e., the column(s) for your covariate of interest) come after the columns of mod0, i.e., the first n columns of mod must be the same as mod0, where n is the number of columns in mod0.") # If this function were modified to enable the covariate of interest columns to come at the start or elsewhere, e.g., changing mod[,-c(1:ncol(mod0)),drop=FALSE] below, it would require matching on column names of mod and mod0, but then the function would require column names to be the same, a new requirement and one often not met for the intercept column.

    if(is.null(p) & is.null(logitp)) stop("Either p or logitp must be provided.")
    if(!is.null(p)) {
        stopifnot(min(p)>=0 & max(p)<=1)
        stopifnot(ncol(p)==nrow(mod))
    }
    if(!is.null(logitp)) {
        stopifnot(ncol(logitp)==nrow(mod))
    } else logitp = log(p)-log(1-p)

    errmsg = "coeff must be a character or numeric index for the column of mod that is the coefficient of interest."
    if(is.character(coeff)) {
        if(!coeff%in%colnames(mod)) stop(errmsg)
        coiIndex = which(colnames(mod)==coeff)
    } else if(is.numeric(coeff)) {
        if(coeff>ncol(mod)) stop(errmsg)
        coiIndex = coeff        
    } else stop(errmsg)

    if(is.null(svs)) {
        require(sva)
        cat("Running SVA\n")
        svaobj = sva(logitp,mod=mod,mod0=mod0,method="irw")
        svs = svaobj$sv
    }
    args = list(svs=svs, mod=mod, mod0=mod0, coeff=coeff, use.limma=use.limma, smoo=smoo, SPAN=SPAN, DELTA=DELTA, use=use, Q=Q, min.probes=min.probes, min.value=min.value, keepXY=keepXY, sortBy=sortBy) #save for obtaining q-values.

    svs = as.matrix(svs)
    svs = svs[,!duplicated(svs[1,])] #remove duplicate sv's
    if(svs==0) X = mod else X = cbind(mod, svs)
    P = ncol(X)
    if(rob) no = 1:ncol(mod) else no = c(1, which(!colnames(mod)%in%colnames(mod0)))
    if(use.limma) {
        #require(limma)
        if(verbose) cat("\nRegression (limma)\n")
        fit=limma::lmFit(logitp, X)

        ##Get cleanp:
        cleanp = ilogit(logitp-fit$coef[,-no,drop=FALSE]%*%t(X[,-no,drop=FALSE]))
        if(only.cleanp) return(cleanp)
      
        beta=fit$coef[,coiIndex]
        if(verbose) cat("Obtaining estimates for ",colnames(fit$coef)[coiIndex],"\n")
        fit2=limma::ebayes(fit, proportion=0.01)
        wald=fit2$t[,coiIndex]
        pval=fit2$p.value[,coiIndex]
        ## sigma here only gets used when smoo=="weighted.loess": 
        #sigma=fit$sigma # same as sigma when !use.limma.
        sigma=sqrt(fit2$s2.post) #could do this instead, to borrow strength.
    } else {
        ##below is just bare version of lm
        if(verbose) cat("\nRegression\n")
        Hat=solve(t(X)%*%X)%*%t(X)
        beta=(Hat%*%t(logitp))
        N=ncol(logitp)
        sigma=genefilter::rowSds(logitp-t(X%*%beta))*sqrt((N-1)/(N-P))

        ##Get cleanp:
        cleanp=ilogit(logitp-t(X[,-no,drop=FALSE]%*%beta[-no,,drop=FALSE]))
        colnames(cleanp)=colnames(logitp)
        if(only.cleanp) return(cleanp)
  
        if(verbose) cat("Obtaining estimates for ",rownames(beta)[coiIndex],"\n")
        wald=beta[coiIndex,]/(sqrt(diag(solve(t(X)%*%X))[coiIndex])*sigma)
        beta=beta[coiIndex,]
        ##All one sided tests, not 2-sided:
        #pval=ifelse(wald<0, pt(wald,df=N-P), 1-pt(wald,df=N-P))
        pval=pt(-abs(wald),df=N-P)*2
    }

    ## Find DMRs:
    if(verbose) cat("Smoothing\n")
    pnsIndexes = split(seq(along = pns), pns)
    sbeta = dosmooth(stat=beta, pnsIndexes=pnsIndexes, pos=pos, smoo=smoo, verbose=verbose, k=k, SPAN=SPAN, DELTA=DELTA, sigma=sigma) #sigma only gets used if smoo=="weighted.loess"
    if(use=="sbeta") { ##Finding based on sbeta:
        swald = NULL
        odmrs = regionFinder(x=sbeta,pns,chr,pos,cutoff=quantile(sbeta,Q),verbose=verbose)
    } else if(use=="swald") { ##Finding based on swald:
        swald = dosmooth(stat=wald, pnsIndexes=pnsIndexes, pos=pos, smoo=smoo, verbose=verbose, k=k, SPAN=SPAN, DELTA=DELTA, sigma=NULL)
        odmrs = regionFinder(x=swald,pns,chr,pos,cutoff=quantile(swald,Q),verbose=verbose)
    }
    ##Remove dmrs in sex chromosomes if keepXY==FALSE:
    if(!keepXY) odmrs=odmrs[which(!odmrs$chr%in%c("chrX","chrY")),]
    ##Remove dmrs that don't have at least min.probes probes:
    odmrs=odmrs[which(odmrs$nprobes>min.probes),]
    ##Remove dmrs in which the average difference or correlation is not at least min.value:
    ##(Can't just use the 'value' column regionFinder(..., y=sbeta) since sbeta is from logit data)
    contin = FALSE
    if(all(mod[,coiIndex]%in%c(0,1))) { #covariate is categorical
        if(verbose) message("Covariate recognized as categorical.")
        base = which(rowSums(mod[,-c(1:ncol(mod0)),drop=FALSE])==0)
        grp1 = which(mod[,coiIndex]==1)
        mat = cbind(rowMeans(cleanp[,grp1,drop=FALSE]), rowMeans(cleanp[,base,drop=FALSE]))
        ## Mean diff:
        res = t(apply(odmrs[,c("indexStart","indexEnd")],1,
                      function(se) colMeans(mat[se[1]:se[2],,drop=FALSE])))
        x = res[,1]-res[,2]
        ## Max diff:
        res0 = mat[,1]-mat[,2]
        x2 = apply(odmrs[,c("indexStart","indexEnd")],1,
                       function(se) {
                           rt = res0[se[1]:se[2]]
                           rt[which.max(abs(rt))]
                       })      
    } else { #covariate is continuous
        if(verbose) message("Covariate recognized as continuous.")
        contin = TRUE
        mat = rowCor(m=cleanp,y=mod[,coiIndex])
        ## Mean correlation:
        x = apply(odmrs[,c("indexStart","indexEnd")],1,
                  function(se) median(mat[se[1]:se[2]]))
        ## Max correlation:
        x2 = apply(odmrs[,c("indexStart","indexEnd")],1,
                   function(se) {
                       rt = mat[se[1]:se[2]]
                       rt[which.max(abs(rt))]
                   })
    }
    stopifnot(length(x)==nrow(odmrs))
    odmrs$avg = x
    odmrs$max = x2
    odmrs$area.raw = abs(x*odmrs$nprobes) #unsigned area, like area column from regionFinder
    x[is.na(x)] = 1 #maximum possible so it will definitely exceed min.value
    dmrs = odmrs[abs(x)>min.value,] #If not for this there would be no difference in the dmrs returned by dmrFind when rob=TRUE or rob=FALSE.
    
    if(nrow(dmrs)==0) {
        if(verbose) cat("\nNo DMRs found (using ",use,").\n")
    } else {
        dmrs=dmrs[order(-abs(dmrs[,sortBy])),]
        if(verbose) cat("\nFound",nrow(dmrs),"potential DMRs\n")
    }
    if(!verbose) cat(".") #for qval when verbose=FALSE.
    if(only.dmrs) {
        if(!contin) return(list(dmrs=dmrs, pval=pval, pns=pns, chr=chr, pos=pos, args=args))
        if(contin)  return(list(dmrs=dmrs, pval=pval, pns=pns, chr=chr, pos=pos, args=args, beta=beta, sbeta=sbeta))        
    } else {
       if(!contin)  return(list(dmrs=dmrs, pval=pval, pns=pns, chr=chr, pos=pos, args=args, cleanp=cleanp))
       if(contin)   return(list(dmrs=dmrs, pval=pval, pns=pns, chr=chr, pos=pos, args=args, beta=beta, sbeta=sbeta, cleanp=cleanp)) #whether or not the result contains beta and sbeta is how resultMaker knows that the analysis is of continuous phenotype.
   }
}

rowCor <- function(m, y) { #faster version of apply(m,1,function(x) cor(x,y))
    #require(genefilter)
    stopifnot(is.matrix(m)|is.data.frame(m))
    stopifnot(is.vector(y))
    stopifnot(ncol(m)==length(y))
    X = cbind(1,y)
    Hat=solve(t(X)%*%X)%*%t(X)
    beta=t(Hat%*%t(m))[,"y"]
    beta * sd(y)/rowSds(m)
}

dosmooth <- function(stat, pnsIndexes, pos, smoo="runmed", k=3, SPAN=300, DELTA=36, sigma=NULL, verbose=TRUE) {
#SPAN  = 300 ##used by loess. we smooth 300 base pairs
#DELTA = 36  ##span parameter in loess smoothing will = SPAN/(DELTA* no. of probes in region)
    if(smoo=="runmed") { #assumes stat is already in order by genomic location
        sstat = stat
        if(verbose) { pb=txtProgressBar(min=0,max=length(pnsIndexes),initial=0); j=0 }
        for(ind in pnsIndexes) {
            if(length(ind)>k) {
                sstat[ind]=runmed(stat[ind],k,endrule="constant")  
            } else sstat[ind]=median(stat[ind])
            if(verbose) { j=j+1; setTxtProgressBar(pb,j) }
        }
        if(verbose) close(pb)        
    } else if(smoo=="loess") {
        sstat = stat
        if(verbose) { pb=txtProgressBar(min=0,max=length(pnsIndexes),initial=0); j=0 }
        for(ind in pnsIndexes) {
            if(length(ind)>4) {
                thespan=SPAN/(DELTA*length(pos[ind]))
                fit1=loess(stat[ind]~pos[ind],span=thespan,degree=1,family="symmetric")    
                sstat[ind]=fit1$fitted
            } else sstat[ind]=median(stat[ind])
            if(verbose) { j=j+1; setTxtProgressBar(pb,j) }
        }
        if(verbose) close(pb)        
    } else if(smoo=="weighted.loess") { #use Andrew's way
        sstat = loessPns(stat=stat, se=sigma, pnsIndexes=pnsIndexes, pos=pos, verbose=verbose)
    } else stop("Choices for smoo are runmed, loess, and weighted.loess.")
    sstat
}

loessPns <- function(stat, se=NULL, pnsIndexes, pos, numProbes = 8, verbose=TRUE) {
        ##pnsIndexes = split(seq(along = pns), pns) #takes a second, so do outside function for speed
	#require(limma)
	sstat = stat

        if(verbose) pb=txtProgressBar(min=0,max=length(pnsIndexes),initial=0)        
	if(is.null(se)) {
		for(i in seq(along=pnsIndexes)) {
			#if(i %% 1e4 == 0) cat(".")
			Index = pnsIndexes[[i]]
			ss = numProbes/length(Index) # make a span
			if(ss > 1) ss = 0.75
						
			if(length(Index) > 4) {
			    sstat[Index] = loessFit(stat[Index], pos[Index], span = ss)$fitted
			} else sstat[Index] = median(stat[Index])
                        if(verbose) setTxtProgressBar(pb,i)
		}
	} else {
		for(i in seq(along=pnsIndexes)) {
			#if(i %% 1e4 == 0) cat(".")
			Index = pnsIndexes[[i]]
			ss = numProbes/length(Index) # make a span
			if(ss > 1) ss = 0.75
						
			if(length(Index) > 4) {
			    sstat[Index] = loessFit(stat[Index], pos[Index], span = ss, weights = 1/se[Index])$fitted
			} else sstat[Index] = median(stat[Index])
                        if(verbose) setTxtProgressBar(pb,i)
		}
	}
        if(verbose) close(pb)
	return(sstat)
}

######################################################################
######################################################################

qval <- function(p=NULL, logitp=NULL, dmr, numiter=500, seed=54256, verbose=FALSE, mc=1, return.permutations=FALSE, method="direct") {
    if(!any(c("direct","pool","fwer")%in%method)) stop("method must = direct, pool, or fwer.")
    #require(parallel)
    if(is.null(p) & is.null(logitp)) stop("Either p or logitp must be provided (the same as you provided to dmrFind when it produced dmr).")
    if(!is.null(p) & "cleanp"%in%names(dmr)) {
        if(nrow(dmr$cleanp)!=nrow(p) | ncol(dmr$cleanp)!=ncol(p)) stop("p must be the same as the p that you provided to dmrFind when it produced dmr.")        
    }
    if(!is.null(logitp)) {
        if("cleanp"%in%names(dmr)) if(nrow(dmr$cleanp)!=nrow(logitp) | ncol(dmr$cleanp)!=ncol(logitp)) stop("logitp must be the same as the logitp that you provided to dmrFind when it produced dmr.")
    } else logitp = log(p)-log(1-p)
    #stopifnot(all(dmr$dmrs[,dmr$args$sortBy]>=0))

    ## Test that result is same as original:
    dmr2 = dmrFind(logitp=logitp, svs=dmr$args$svs, mod=dmr$args$mod, mod0=dmr$args$mod0, only.cleanp=FALSE, only.dmrs=TRUE, coeff=dmr$args$coeff, use.limma=dmr$args$use.limma, smoo=dmr$args$smoo, SPAN=dmr$args$SPAN, DELTA=dmr$args$DELTA, use=dmr$args$use, Q=dmr$args$Q, min.probes=dmr$args$min.probes, min.value=dmr$args$min.value, pns=dmr$pns, chr=dmr$chr, pos=dmr$pos, keepXY=dmr$args$keepXY, sortBy=dmr$args$sortBy, verbose=verbose)
    orig_tab = dmr$dmrs    
    TFvec = vector()
    for(hg in 1:ncol(dmr2$dmrs)) TFvec[hg] = isTRUE(all.equal(dmr2$dmrs[,hg],orig_tab[,hg]))
    dimsame = identical(dim(dmr2$dmrs),dim(orig_tab))
    if(any(TFvec==FALSE) | !dimsame) stop("Make sure your p or logitp argument is the same as the one you passed to dmrFinder when it produced dmr.")

    ## Now start the iterations:
    svs = as.matrix(dmr$args$svs)
    svs = svs[,!duplicated(svs[1,])] #remove duplicate sv's
    if(isTRUE(svs==0)) {
        X  = dmr$args$mod
        X0 = dmr$args$mod0
    } else {
        X  = cbind(dmr$args$mod, svs)
        X0 = cbind(dmr$args$mod0, svs)
    }
    #fit1 = limma::lmFit(logitp,X)
    #r_ij = residuals(fit1, logitp)
    #fit0 = limma::lmFit(logitp,X0)
    #n_ij = fitted(fit0)
    ## To use less memory and go faster, do this instead of lmFit:
    r_ij = lm_fit(dat=logitp,X=X,out="res")
    n_ij = lm_fit(dat=logitp,X=X0,out="fit")

    colsamp = matrix(NA, nrow=numiter, ncol=ncol(r_ij))
    set.seed(seed)
    for(j in 1:nrow(colsamp)) colsamp[j,] = sample(1:ncol(r_ij),replace=TRUE)
    fun1 <- function(h, method) {
        dmr2 = dmrFind(logitp= n_ij + r_ij[,colsamp[h,]], svs=dmr$args$svs, mod=dmr$args$mod, mod0=dmr$args$mod0, only.cleanp=FALSE, only.dmrs=TRUE, coeff=dmr$args$coeff, use.limma=dmr$args$use.limma, smoo=dmr$args$smoo, SPAN=dmr$args$SPAN, DELTA=dmr$args$DELTA, use=dmr$args$use, Q=dmr$args$Q, min.probes=dmr$args$min.probes, min.value=dmr$args$min.value, pns=dmr$pns, chr=dmr$chr, pos=dmr$pos, keepXY=dmr$args$keepXY, sortBy=dmr$args$sortBy, verbose=verbose)
        ret = vector("list",3)
        ## abs() here makes these tests 2-sided:
        if("direct"%in%method) {
            if(nrow(dmr2$dmrs)==0) add=0 else add = abs(dmr2$dmrs[,dmr$args$sortBy])
            #rm(dmr2); gc(); gc()
            ret[[1]] = counts(add, abs(orig_tab[,dmr$args$sortBy]))
        }
        if("pool"%in%method) {
            if(nrow(dmr2$dmrs)>0) ret[[2]] = abs(dmr2$dmrs[,dmr$args$sortBy])
        }
        if("fwer"%in%method) {
            if(nrow(dmr2$dmrs)==0) ret[[3]] = 0 else ret[[3]] = max(abs(dmr2$dmrs[,dmr$args$sortBy]))
        }
        ret
    }
    message("\nBeginning the iterations")
    nullnum0 = mclapply(1:nrow(colsamp), fun1, method=method, mc.cores=mc)
    message("\nFinished the iterations")
    dig=7
    if("direct"%in%method) {
        nullnum = sapply(nullnum0, function(x) x[[1]])
        ##Obtain FDR estimates:
        stopifnot(nrow(nullnum)==nrow(orig_tab))
        qv = rowMeans(nullnum)/seq(1,nrow(orig_tab))
        orig_tab$qvalue = round(monotonic(pmin(1,qv)), dig)
    }
    if("pool"%in%method) {
        nulldist = unlist(sapply(nullnum0, function(x) x[[2]]))
        Fn = ecdf(nulldist + (1e-9))
        pv = 1-Fn(abs(orig_tab[,dmr$args$sortBy]))
        orig_tab$pvalue.pool = round(pv,dig)
        orig_tab$qvalue.pool = round(edge.qvalue(pv)$qvalues,dig)
    }
    if("fwer"%in%method) {
        maxs = sapply(nullnum0, function(x) x[[3]])
        Fnn = ecdf(maxs + (1e-9))
        pvv = 1-Fnn(abs(orig_tab[,dmr$args$sortBy]))
        orig_tab$pvalue.fwer = round(pvv, 7) #max(1,nchar(numiter)-1)
    }
    
    if(return.permutations) return(list(q=orig_tab, permutations=colsamp)) else return(orig_tab)
}

lm_fit <- function(dat, X, out) {
    ## Use this instead of lmFit because this uses less memory and is faster.  Estimated betas are the same anyway--lmFit is only useful for (shrunk) se's, which we don't need for this.
    Hat=solve(t(X)%*%X)%*%t(X)
    if(out=="fit") {
        return(t(X%*% (Hat%*%t(dat)) ))
    } else if(out=="res") {
        return(dat - t(X%*% (Hat%*%t(dat)) ))
    } else stop("invalid out arg")
}
counts <- function(x1,x2) {
    if(length(x1)==0) rep(0,length(x2)) else {
        ec = ecdf(x1+(1e-9)) #add small amount to make ecdf < rather than <=.
        out = 1-ec(x2)
        stopifnot(all(out>=0 & out<=1))
        round(out*length(x1)) # round since some numbers have very small imprecision
    }

    ##This is the slower way to do this function:
    #ret = vector()
    #for(i in 1:length(x2))  ret[i] <- sum(x1>=x2[i])
    #ret
}
monotonic <- function(x,increasing=TRUE) {
    if(length(x)>1) {
        if(increasing)  for(i in 2:length(x))  x[i] = max(x[i],x[i-1])
        if(!increasing) for(i in 2:length(x))  x[i] = min(x[i],x[i-1])
    }
    x
}

# qvalue calculator that doesn't use tcltk
edge.qvalue <- function(p,lambda = seq(0, 0.9, 0.05), pi0.method = "smoother",
                         fdr.level = NULL, robust = FALSE,smooth.df = 3, smooth.log.pi0 = FALSE, ...) {

  err.func <- "edge.qvalue"
  if (min(p) < 0 || max(p) > 1) {
    err.msg(err.func,"P-values not in valid range.")
    return(invisible(1))
  }
  if (length(lambda) > 1 && length(lambda) < 4) {
    err.msg(err.func,"If length of lambda greater than 1, you need at least 4 values.")
    return(invisible(1))
  }
  if (length(lambda) > 1 && (min(lambda) < 0 || max(lambda) >= 1)) {
    err.msg(err.func,"Lambda must be in [0,1).")
    return(invisible(1))
  }
  m <- length(p)
  if (length(lambda) == 1) {
    if (lambda < 0 || lambda >= 1) {
      err.msg(err.func,"Lambda must be in [0,1).")
      return(invisible(1))
    }
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0 <- min(pi0, 1)
  } else {
    pi0 <- rep(0, length(lambda))
    for (i in 1:length(lambda)) {
      pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
    }
    if (pi0.method == "smoother") {
      if (smooth.log.pi0){ 
        pi0 <- log(pi0)
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0 <- predict(spi0, x = max(lambda))$y
      }
      if (smooth.log.pi0) {
        pi0 <- exp(pi0)
      }
      pi0 <- min(pi0, 1)
    }
    else if (pi0.method == "bootstrap") {
      minpi0 <- min(pi0)
      mse <- rep(0, length(lambda))
      pi0.boot <- rep(0, length(lambda))
      for (i in 1:100) {
        p.boot <- sample(p, size = m, replace = TRUE)
        for (i in 1:length(lambda)) {
          pi0.boot[i] <- mean(p.boot > lambda[i])/(1 - lambda[i])
        }
        mse <- mse + (pi0.boot - minpi0)^2
      }
      pi0 <- min(pi0[mse == min(mse)])
      pi0 <- min(pi0, 1)
    }
    else {
      err.msg(err.func,"'pi0.method' must be one of 'smoother' or 'bootstrap'")
      return(invisible(1))
    }
  }
  if (pi0 <= 0) {
    err.msg(err.func,"The estimated pi0 <= 0. Check that you have valid\np-values or use another lambda method.")
    return(invisible(1))
  }
  if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 1)) {
    err.msg(err.func,"'fdr.level' must be within (0,1].")
    return(invisible(1))
  }
  u <- order(p)
  qvalue.rank <- function(x) {
    idx <- sort.list(x)
    fc <- factor(x)
    nl <- length(levels(fc))
    bin <- as.integer(fc)
    tbl <- tabulate(bin)
    cs <- cumsum(tbl)
    tbl <- rep(cs, tbl)
    tbl[idx] <- tbl
    return(tbl)
  }
  v <- qvalue.rank(p)
  qvalue <- pi0 * m * p/v
  if (robust) {
    qvalue <- pi0 * m * p/(v * (1 - (1 - p)^m))
  }
  qvalue[u[m]] <- min(qvalue[u[m]], 1)
  for (i in (m - 1):1) {
    qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 1)
  }
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, pvalues = p, fdr.level = fdr.level, significant = (qvalue <= fdr.level), lambda = lambda)
  }
  else {
    retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, pvalues = p, lambda = lambda)
  }
  class(retval) <- "qvalue"
  return(retval)
}

######################################################################
######################################################################

plotDMRs <- function(dmrs, Genome, cpg.islands, exposure, outfile="./dmr_plots.pdf", which_plot=1:50, all.points=TRUE, plot.these=NULL, ADD=3000, cols=c("black","blue","red","gray","brown","pink","orange"), legend.size=1, smoo="loess", SPAN=300, DELTA=36, point.info=FALSE, pch.groups=NULL, panel3="pvalues", G=NULL, seq=NULL) {

    args = dmrs$args
    if(panel3=="G") {      
        if(is.null(G)|is.null(seq)) stop("If panel3=G, then G and seq must be provided.")
        stopifnot(length(seq)==nrow(G))
        stopifnot(length(seq)==length(dmrs$pos))
        stopifnot(all(colnames(G)==colnames(dmrs$cleanp)))
        stopifnot(ncol(G)==nrow(args$mod))        
    }
    stopifnot(length(exposure)==ncol(dmrs$cleanp))
    stopifnot(inherits(Genome,"BSgenome"))
    stopifnot(is.data.frame(cpg.islands))
    stopifnot(all(c("chr","start","end")%in%colnames(cpg.islands)))
    ocpgi = cpg.islands
    dmrs$chr = as.character(dmrs$chr)
    ocpgi = ocpgi[ocpgi$chr%in%unique(dmrs$chr),]
    
    if("beta"%in%names(dmrs)) beta = dmrs$beta else beta = NULL
    if("sbeta"%in%names(dmrs)) sbeta = dmrs$sbeta else sbeta = NULL
    ##plot as continuous if analysis was of continuous exposure and exposure not recoded for plotting.  Otherwise plot as categorical:
    if(!is.null(beta) & all(exposure==args$mod[,args$coeff])) cont=TRUE else cont=FALSE
    if(!cont) {
        groups = as.numeric(factor(exposure))
        gNames = levels(factor(exposure))
        stopifnot(length(cols)>=length(unique(exposure)))
        if(!all.points) stopifnot(all(plot.these%in%gNames))
        if(!is.null(pch.groups)) stopifnot(length(pch.groups)==length(exposure))
    }
    dmrs$dmrs$chr = as.character(dmrs$dmrs$chr)   

    if(panel3=="G") {
        #require(Biobase) #for rowMedians
        G = log2(preprocessCore::normalize.quantiles(G))
        # remove probe effects from G:
        G = G - rowMedians(G)
        # GC sequence effect correction:
        G = remove_gc_effect(G, seq)

        if(is.null(beta)) { #not !cont
            g1 = which(args$mod[,args$coeff]==1)
            rmP = 1:ncol(args$mod0)
            g0 = which(rowSums(args$mod[,-rmP,drop=FALSE])==0)
            mes = cbind(rowMedians(G[,g1]), rowMedians(G[,g0]))
            colsg = c(unique(cols[groups[g1]]), unique(cols[groups[g0]]))
            stopifnot(length(colsg)==2)
            stopifnot(nrow(mes)==length(dmrs$pos))           
        } else mes = rowCor(G,exposure)
    }

    if(!cont) {
        if(all.points) ind = 1:length(exposure) else ind = which(exposure%in%plot.these)
        gIndexes=split(1:ncol(dmrs$cleanp),groups)
        if(!is.null(pch.groups)) {
            thepch = as.numeric(factor(pch.groups))[ind]
            thepch.lb = tapply(thepch, as.character(pch.groups)[ind],unique)
        } else thepch = pch.groups #if null, keep NULL and it will plot numbers and letters.
    }
    
    ###################################################################
    cat("Making",min(length(which_plot),nrow(dmrs$dmrs)),"figures\n")
    pdf(outfile,width=11,height=8)
    for(i in intersect(which_plot,1:nrow(dmrs$dmrs))) {
      cat("Plotting DMR candidate",i,"\n")
      layout(matrix(1:3,ncol=1),heights=c(0.6,0.3,0.3))
      par(mar=c(0,3.6,1.1,1.1),oma=c(0,0,1,0),mgp=c(2.6,.5,0))

      Index = which(dmrs$chr==dmrs$dmrs$chr[i]) 
      Index = Index[which(dmrs$pos[Index]>=dmrs$dmrs$start[i]-ADD & dmrs$pos[Index]<=dmrs$dmrs$end[i]+ADD)]
      test  = Index[which(dmrs$pos[Index]>=dmrs$dmrs$start[i]     & dmrs$pos[Index]<=dmrs$dmrs$end[i]    )]
      if(length(Index)==0 | length(test)==0) stop("Region doesn't exist on array.")

      xx=dmrs$pos[Index]
      if(!cont) {
          cn = colnames(dmrs$cleanp[,ind])
          matplot(xx,dmrs$cleanp[Index,ind], col=cols[groups[ind]], cex=0.4,
                  xlim=range(xx), ylim=c(0,1), ylab="Methylation", xaxt="n", xlab="",
                  main=paste("DMR ",i," - ",dmrs$dmrs$chr[i],":",dmrs$dmrs$start[i],"-", 
                             dmrs$dmrs$end[i],sep=""),las=1, pch=thepch)
          for(j in seq(along=gIndexes)) {
              yy=rowMeans(dmrs$cleanp[Index,gIndexes[[j]],drop=FALSE])
              if(smoo=="runmed") {
                  yy2 = runmed(yy,3,endrule="constant")
              } else {
                  fit1=loess(yy~xx,degree=1,span=SPAN/(DELTA*length(xx)),family="symmetric")
                  yy2 =fit1$fitted        
              }
              lines(xx,yy2,col=cols[as.numeric(names(gIndexes))[j]],lwd=2)
          }
          legend("topleft",gNames,col=cols[seq(along=gNames)],lty=1,lwd=2,cex=legend.size)
          if(!is.null(pch.groups)) legend("topright", pch=thepch.lb, legend=names(thepch.lb))
      } else {
          pc = rowCor(dmrs$cleanp[Index,],exposure)
          plot(xx,pc,ylim=c(-1,1), xlab="", xaxt="n", xlim=range(xx), las=1,
               ylab="Pearson correlation with phenotype",
               main=paste("DMR ",i," - ",dmrs$dmrs$chr[i],":",dmrs$dmrs$start[i],"-", dmrs$dmrs$end[i],sep=""))
          #plot(xx,beta[Index],cex=0.6,ylim=range(sbeta),xlab="",
          #     xaxt="n", xlim=range(xx), las=1,
          #     ylab="difference in logit(p) associated with having phenotype 1 unit larger",
          #     main=paste("DMR ",i," - ",dmrs$dmrs$chr[i],":",dmrs$dmrs$start[i],"-", dmrs$dmrs$end[i],sep=""))
          #lines(xx,sbeta[Index],lwd=2,lty=1)
          abline(h=0,lty=3)
      }
      abline(v=dmrs$dmrs$start[i], lty=2)
      abline(v=dmrs$dmrs$end[i], lty=2)

      ############
      ##PLOT CPG ISLANDS AND MCRBC RECOGNITION SITES:
      par(mar=c(0,3.6,0.5,1.1))
      plot_CpG(thechr=dmrs$dmrs$chr[i],xx=xx,ocpgi=ocpgi,Genome=Genome,mcrbc=TRUE)

      ##PLOT DIFFERENCE IN G, GENES, OR MANHATTAN PLOT:
      if(panel3=="G") {
          par(mar=c(5.1,3.6,0.5,1.1))
          if(is.null(beta)) {
              #matplot(xx, mes[Index,], col=colsg, ylab="Median corrected G", xlab="Location",
              #        ylim=c(-2,2), xlim=range(xx), type="b", las=1, pch=21)
              matplot(xx, mes[Index,1]-mes[Index,2], ylab="Difference in median G", xlab="Location",
                      ylim=c(-2.5,2.5), xlim=range(xx), type="b", las=1, pch=21)
          } else {
              matplot(xx, mes[Index], ylab="cor(corrected G, exposure)", xlab="Location",
                      ylim=c(-1,1), xlim=range(xx), type="b", las=1, pch=21)
          }
          abline(h=0,lty=3)
      } else {
          par(mar=c(4.1,3.6,0.5,1.1), mgp=c(1.6,.5,0))
          yy = -log10(dmrs$pval[Index])
          plot(xx,yy, ylab="-log10(pvalue)", xlab="Location",
               ylim=c(0,7), xlim=range(xx), las=1)
          if(smoo=="runmed") {
              yy2 = runmed(yy,3,endrule="constant")
          } else {
              fit1=loess(yy~xx,degree=1,span=SPAN/(DELTA*length(xx)),family="symmetric")
              yy2 =fit1$fitted        
          }
          lines(xx,yy2)
          #abline(h=-log10(.05),lty=2)
        }
    }
    dev.off()
    if(point.info & is.null(pch.groups) & !cont) {
        pts = c(1:9,0,letters,LETTERS)
        pts = rep(pts, ceiling(length(cn)/length(pts)))
        ret = cbind(cn, pts[1:length(cn)])
        colnames(ret) = c("sample","symbol")
        ret
    }
}

plotRegions <- function(thetable, cleanp, chr, pos, seq=NULL, Genome, cpg.islands, exposure, exposure.continuous=FALSE, outfile="./myregions.pdf", all.points=TRUE, plot.these=NULL, ADD=3000, cols=c("black","red","blue","gray","green","orange","brown"), legend.size=1, smoo="loess", SPAN=300, DELTA=36, panel3="none", G=NULL, grs=NULL) {

    if(is.factor(exposure)|!exposure.continuous) cont=FALSE else cont=TRUE  
    if(panel3=="G") {
        if(is.null(G)|is.null(grs)|is.null(seq)) stop("If panel3=G, then G ,grs and seq must be provided.")
        stopifnot(length(seq)==nrow(G))
        stopifnot(length(seq)==length(pos))
        stopifnot(all(colnames(G)==colnames(cleanp)))
        if(!cont) {
            stopifnot(length(grs)==2)
            stopifnot(all(grs%in%exposure))
        }
    }
    stopifnot(length(exposure)==ncol(cleanp))
    stopifnot(length(chr)==length(pos))
    stopifnot(length(chr)==nrow(cleanp))
    stopifnot(is.numeric(thetable$start))
    stopifnot(is.numeric(thetable$end))  
    stopifnot(inherits(Genome,"BSgenome"))
    stopifnot(is.data.frame(cpg.islands))
    stopifnot(all(c("chr","start","end")%in%colnames(thetable)))
    stopifnot(all(c("chr","start","end")%in%colnames(cpg.islands)))
    ocpgi = cpg.islands
    chr = as.character(chr)
    ocpgi = ocpgi[ocpgi$chr%in%unique(chr),]
    thetable$chr = as.character(thetable$chr)
          
    if(!cont) {
        ### for plots, this defines the groups to be shown in color.
        groups = as.numeric(factor(exposure))
        gNames = levels(factor(exposure))
        if(!all.points & !all(plot.these%in%exposure)) stop("plot.these not all in exposure")
        if(length(cols)<length(unique(exposure))) stop("Not as many colors as groups.")
    }  

    if(panel3=="G") {
        #require(Biobase) #for rowMedians
        G = log2(preprocessCore::normalize.quantiles(G))
        # remove probe effects from G:
        G = G - rowMedians(G)
        # GC sequence effect correction:
        G = remove_gc_effect(G, seq)

        if(!cont) {
            g1 = which(exposure==grs[1])
            g0 = which(exposure==grs[2])        
            mes = cbind(rowMedians(G[,g1]), rowMedians(G[,g0]))
            colsg = c(unique(cols[groups[g1]]), unique(cols[groups[g0]]))
            stopifnot(length(colsg)==2)
            stopifnot(nrow(mes)==length(pos)) #pos is dmr$pos, defined above
        } else mes = rowCor(G,exposure)
    }

    if(!cont) {
        if(all.points) ind = 1:length(exposure) else ind = which(exposure%in%plot.these)
        gIndexes=split(1:ncol(cleanp),groups)
    }

    ###################################################################
    cat("Making",nrow(thetable),"figures\n")
    pdf(outfile,width=11,height=8)
    for(i in 1:nrow(thetable)) {
        if(panel3=="G") {
            layout(matrix(1:3,ncol=1),heights=c(0.6,0.3,0.3))
        } else layout(matrix(1:2,ncol=1),heights=c(0.6,0.3))
        par(mar=c(0,3.6,1.1,1.1),oma=c(0,0,1,0),mgp=c(2.6,.5,0))

        ##PLOT REGIONS:
        Index = which(chr==thetable$chr[i])
        Index = Index[which(pos[Index]>=thetable$start[i]-ADD & pos[Index]<=thetable$end[i]+ADD)]
        test  = Index[which(pos[Index]>=thetable$start[i]     & pos[Index]<=thetable$end[i]    )] 
        if(length(Index)==0 | length(test)==0) {
            message("Region ",i," not on the array.")
            next
        }

        xx=pos[Index]
        if(!cont) {
            matplot(xx,cleanp[Index,ind],col=cols[groups[ind]], cex=0.4,
                    xlim=range(xx), ylim=c(0,1), ylab="Methylation", xaxt="n", xlab="",
                    main=paste("Region ",i," - ",thetable$chr[i],":",thetable$start[i],"-", 
                               thetable$end[i],sep=""),las=1)
            for(j in seq(along=gIndexes)){
                yy=rowMeans(cleanp[Index,gIndexes[[j]],drop=FALSE])
                if(smoo=="runmed") {
                    yy2 = runmed(yy,3,endrule="constant")
                } else {
                    fit1=loess(yy~xx,degree=1,span=SPAN/(DELTA*length(xx)),family="symmetric")
                    yy2 =fit1$fitted        
                }
                lines(xx,yy2,col=cols[as.numeric(names(gIndexes))[j]],lwd=2)
            }
            legend("topleft",gNames,col=cols[seq(along=gNames)],lty=1,lwd=2,cex=legend.size)
        } else {
            pc = rowCor(cleanp[Index,],exposure)
            plot(xx,pc,ylim=c(-1,1), xlab="", xaxt="n", xlim=range(xx), las=1,
                 ylab="Pearson correlation with phenotype",
                 main=paste("Region ",i," - ",thetable$chr[i],":",thetable$start[i],"-", thetable$end[i],sep=""))
            #plot(xx,beta[Index],cex=0.6,ylim=range(sbeta),xlab="",
            #     xaxt="n", xlim=range(xx), las=1,
            #     ylab="difference in logit(p) associated with having phenotype 1 unit larger",
            #     main=paste("Region ",i," - ",Dmrs$chr[i],":",Dmrs$start[i],"-", Dmrs$end[i],sep=""))
            #lines(xx,sbeta[Index],lwd=2,lty=1)
            abline(h=0,lty=3)
        }      
        abline(v=thetable$start[i], lty=2)
        abline(v=thetable$end[i], lty=2)

        ############
        ##PLOT CPG ISLANDS AND MCRBC RECOGNITION SITES:
        if(panel3=="G") {
            plot_CpG(thechr=thetable$chr[i],xx=xx,ocpgi=ocpgi,Genome=Genome,mcrbc=TRUE)
        } else {
            par(mar=c(4.1,3.6,0.5,1.1))
            plot_CpG(thechr=thetable$chr[i],xx=xx,ocpgi=ocpgi,Genome=Genome,mcrbc=TRUE,
                     xaxt="s",xlab="Location")
        }
        ##PLOT DIFFERENCE IN G:
        if(panel3=="G") {
            par(mar=c(4.1,3.6,0.5,1.1))
            if(!cont) {
                #matplot(xx, mes[Index,], col=colsg, ylab="Median corrected G", xlab="Location",
                #        ylim=c(-2,2), xlim=range(xx), type="b", las=1, pch=21)
                matplot(xx, mes[Index,1]-mes[Index,2], ylab="Difference in median G", xlab="Location",
                        ylim=c(-2.5,2.5), xlim=range(xx), type="b", las=1, pch=21)
            } else {
                matplot(xx, mes[Index], ylab="cor(corrected G, exposure)", xlab="Location",
                        ylim=c(-1,1), xlim=range(xx), type="b", las=1, pch=21)
            }
            abline(h=0,lty=3)
        }
    }
    dev.off()
}

###Remove sequence effect:
remove_gc_effect <- function(M, seq) {
    #require(Biostrings)
    stopifnot(nrow(M)==length(seq))
    #dset = DNAStringSet(seq)
    freq = alphabetFrequency(seq)
    gccontent = rowSums(freq[,c("C", "G")]) / seq@ranges@width
    gc = round(gccontent*100); gc[gc<=20]<-20; gc[gc>=80]<-80
    
    sapply(1:ncol(M),function(i) {
      cat(".")
      meds=tapply(M[,i],gc,median)
      mads=tapply(M[,i],gc,mad)
      Index=match(gc,as.numeric(names(meds)))
      bias=meds[Index]
      sds=mads[Index]
      return((M[,i]-bias)/sds*mad(M[,i]))
    })
}

plot_CpG <- function(thechr,xx,ocpgi,Genome,mcrbc=TRUE,xlab="",xaxt="n") {
    #Genome should be a BSgenome object, e.g., Hsapiens, Mmusculus, Rnorvegicus, Mmulatta, etc.
      ##PLOT CPG ISLANDS AND MCRBC RECOGNITION SITES:
      seq<-Genome[[thechr]]
      Index2=c(min(xx),max(xx))
      if(Index2[1]<1 | Index2[2]>length(seq)) stop("position outside length of chromosome.")
      subseq<-subseq(seq,start=Index2[1],end=Index2[2])
      cpgs=start(matchPattern("CG",subseq))+Index2[1]-1
      ##Find the C in ACG or GCG if mcrbc=TRUE:
      if(mcrbc) Mcrbc=sort(c(start(matchPattern("ACG",subseq))+Index2[1],
                             start(matchPattern("GCG",subseq))+Index2[1]))
      cuts=seq(Index2[1],Index2[2],8)
      scpgs=cut(cpgs,cuts,include.lowest=TRUE)
      x = (cuts[1:(length(cuts)-1)]+cuts[2:(length(cuts))])/2
      y = table(scpgs)/8
      SPAN2=400/diff(range(x))
      d = loess(y~x,span=SPAN2,degree=1)
      plot(x,d$fitted,type="l",ylim=c(0,0.15),xlab=xlab,xaxt=xaxt,
           ylab="CpG density",xlim=range(xx),las=1)
      Rug(cpgs)
      if(mcrbc) Rug(Mcrbc,side=3)
      Index1 = which(ocpgi[,"chr"]==as.character(thechr) &
                   ((ocpgi[,"start"] > min(xx) & ocpgi[,"start"]< max(xx)) |
                    (ocpgi[,"end"  ] > min(xx) & ocpgi[,"end"]  < max(xx))))
      if(length(Index1)>0)
          sapply(Index1,function(j) Rug(unlist(ocpgi[j,c("start","end")]),col="red",lwd=3,side=1))
}

######################################################################
######################################################################

logit <- function(x) log(x)-log(1-x)
ilogit <- function(x) 1/(1+exp(-x))

getSegments <- function(x,factor,cutoff=quantile(abs(x),0.99),verbose=TRUE){

  Indexes=split(seq(along=x),factor)
  regionID=vector("numeric",length(x))
  LAST = 0
  
  segmentation = vector("numeric", length(x))
  type = vector("numeric", length(x))

  for (i in seq(along = Indexes)) {
    if (verbose) if (i%%1000 == 0) cat(".")
    Index = Indexes[[i]]
    y = x[Index]
    z = sign(y) * as.numeric(abs(y) > cutoff)
    w = cumsum(c(1, diff(z) != 0)) + LAST
    segmentation[Index] = w
    type[Index] = z
    LAST = max(w)
  }
  if(verbose) cat("\n")
  ##add a vector of the pns
  res=list(upIndex=split(which(type>0),segmentation[type>0]),
    dnIndex=split(which(type<0),segmentation[type<0]),
    zeroIndex=split(which(type==0),segmentation[type==0]))
  names(res[[1]])<-NULL
  names(res[[2]])<-NULL
  names(res[[3]])<-NULL
  return(res)
}


##you can pass cutoff through the ...
regionFinder<-function(x,regionNames,chr,position,y=x,
                       summary=mean,ind=seq(along=x),order=TRUE,oneTable=TRUE,
                       ...){
  if(is.unsorted(order(chr,position))) stop("regionFinder inputs must be in order by genomic location.")

  Indexes=getSegments(x[ind],regionNames[ind],...)
  res=vector("list",2)
  for(i in 1:2){
    res[[i]]=data.frame(chr=sapply(Indexes[[i]],function(Index) chr[ind[Index[1]]]),
         start=sapply(Indexes[[i]],function(Index) min(position[ind[Index]])),
         end=sapply(Indexes[[i]],function(Index) max(position[ind[Index]])),
         value=sapply(Indexes[[i]],function(Index) summary(y[ind[Index]])),
         area=sapply(Indexes[[i]],function(Index) abs(sum(y[ind[Index]]))),
         pns=sapply(Indexes[[i]],function(Index) regionNames[ind[Index]][1]),
         indexStart=sapply(Indexes[[i]],function(Index) min(ind[Index])),
         indexEnd=sapply(Indexes[[i]],function(Index) max(ind[Index])))
    res[[i]]$nprobes=res[[i]]$indexEnd-res[[i]]$indexStart+1
    #res[[i]]$length=res[[i]]$end-res[[i]]$start+0 #could add 50 to go to end of last probe
  }
  names(res)=c("up","dn")
  if(order & !oneTable){
    if(nrow(res$up)>0) res$up=res$up[order(-res$up$area),]
    if(nrow(res$dn)>0) res$dn=res$dn[order(-res$dn$area),]
  }
  if(oneTable){
    res=rbind(res$up,res$dn)
    if(order & nrow(res)>0) res=res[order(-res$area),]
  }
  return(res)
}

regionMatch <- function(object1, object2, verbose=TRUE) {
    #require(IRanges)
    m1 = object1; m2 = object2 #keep arguments named object1,object2 for backward compatibility
    stopifnot(all(c("chr","start","end")%in%colnames(m1)))
    stopifnot(all(c("chr","start","end")%in%colnames(m2)))
    if(any(is.na(m1[,c("chr","start","end")])) | any(is.na(m2[,c("chr","start","end")]))) stop("No missing values allowed in chr, start, or end columns.")
    m1$chr = as.character(m1$chr)
    m2$chr = as.character(m2$chr)
    stopifnot(is.numeric(m1$start) & is.numeric(m1$end))
    stopifnot(is.numeric(m2$start) & is.numeric(m2$end))
    stopifnot(all(m1$end>=m1$start))
    stopifnot(all(m2$end>=m2$start))
              
    ret = matrix(NA,nrow=nrow(m1),ncol=7)
    colnames(ret) = c("dist","matchIndex","type","amountOverlap","insideDist","size1","size2")
    sp1 = split(1:nrow(m1), m1$chr)
    sp2 = split(1:nrow(m2), m2$chr)          
    for(i in names(sp1)) {
        if(verbose) cat(i," ")
        inds1 = sp1[[i]]
        if(i %in% names(sp2)) {
            inds2 = sp2[[i]]
            m1IR = IRanges(start=m1$start[inds1],end=m1$end[inds1])
            m2IR = IRanges(start=m2$start[inds2],end=m2$end[inds2])
            ret[inds1,"matchIndex"] = inds2[ nearest(m1IR,m2IR) ]
        } else ret[inds1,"matchIndex"] = NA
    }
    if(verbose) cat("\n")
    inds = ret[,"matchIndex"]
    m2b = m2[inds,]
    ## if all(is.na(inds)) & length(inds)<nrow(m2), R returns a matrix (of all NAs) with nrow(m2) rows instead of length(inds), so have to do this:
    if(all(is.na(inds)) & length(inds)<nrow(m2)) m2b = m2b[1:length(inds),]
    ret = data.frame(ret, stringsAsFactors=FALSE)
    
    ret$type[(m1$start>m2b$start & m1$end<=m2b$end) |
             (m1$start>=m2b$start & m1$end<m2b$end)] = "inside"
    ret$type[m1$start<=m2b$start & m1$end>=m2b$end] = "cover"
    ret$type[m1$start>m2b$end] = "disjointR"
    ret$type[m1$end<m2b$start] = "disjointL"
    ret$type[is.na(ret$matchIndex)] = "disjoint"
    ret$type[(m1$start>m2b$start & m1$start<=m2b$end) & m1$end>m2b$end] = "overlapR"    
    ret$type[m1$start<m2b$start & (m1$end>=m2b$start & m1$end<m2b$end)] = "overlapL"

    ret$dist = 0
    ret$dist[ret$type=="disjoint"]  = NA
    ret$dist[ret$type=="disjointR"] = m2b$end[ret$type=="disjointR"] - m1$start[ret$type=="disjointR"]
    ret$dist[ret$type=="disjointL"] = m2b$start[ret$type=="disjointL"] - m1$end[ret$type=="disjointL"]
    ret$amountOverlap[ret$type=="overlapR"] = -1*(m2b$end[ret$type=="overlapR"]-m1$start[ret$type=="overlapR"]+1)
    ret$amountOverlap[ret$type=="overlapL"] = m1$end[ret$type=="overlapL"]-m2b$start[ret$type=="overlapL"]+1
    ret$type[ret$type%in%c("disjointR","disjointL")] = "disjoint"
    ret$type[ret$type%in%c("overlapR","overlapL")] = "overlap"

    ## insideDist column:
    insideIndex = ret$type=="inside" #no missing ret$type at this point
    tmp0 = cbind(m1$end[insideIndex]  -m2b$end[insideIndex],
                 m1$start[insideIndex]-m2b$start[insideIndex])
    tmp = apply(abs(tmp0),1,which.min)
    tmpinsidedist = tmp0[,1]
    tmpIndex = tmp==2
    tmpinsidedist[tmpIndex] = tmp0[tmpIndex,2]
    ret$insideDist[insideIndex] = tmpinsidedist

    ## size1 and size2 columns:
    ret$size1 = m1$end -m1$start +1
    ret$size2 = m2b$end-m2b$start+1
    
    ret
}

clusterMaker <- function(chr,pos,order.it=TRUE,maxGap=300) {
  nonaIndex=which(!is.na(chr) & !is.na(pos))
  Indexes=split(nonaIndex,chr[nonaIndex])
  clusterIDs=rep(NA,length(chr))
  LAST=0
  for(i in seq(along=Indexes)){
    Index=Indexes[[i]]
    x=pos[Index]

    if(order.it){ Index=Index[order(x)];x=pos[Index] }

    y=as.numeric(diff(x)>maxGap)
    z=cumsum(c(1,y))
    clusterIDs[Index]=z+LAST
    LAST=max(z)+LAST
  }
  clusterIDs
}  
