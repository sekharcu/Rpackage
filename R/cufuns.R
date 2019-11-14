
cvsum =function (object, ...) {
  UseMethod("cvsum")
}

cvsum.data.frame = function (object, maxsum = 7L, digits = max(3L, getOption("digits") -
                                                                 3L), ...)
{
  ncw <- function(x) {
    z <- nchar(x, type = "w")
    if (any(na <- is.na(z))) {
      z[na] <- nchar(encodeString(z[na]), "b")
    }
    z
  }
  z <- lapply(X = as.list(object), FUN = cvsum.default, maxsum = maxsum,
              digits = 12L, ...)
  nv <- length(object)
  nm <- names(object)
  lw <- numeric(nv)
  nr <- if (nv)
    max(unlist(lapply(z, NROW)))
  else 0
  for (i in seq_len(nv)) {
    sms <- z[[i]]
    if (is.matrix(sms)) {
      cn <- paste(nm[i], gsub("^ +", "", colnames(sms),
                              useBytes = TRUE), sep = ".")
      tmp <- format(sms)
      if (nrow(sms) < nr)
        tmp <- rbind(tmp, matrix("", nr - nrow(sms),
                                 ncol(sms)))
      sms <- apply(tmp, 1L, function(x) paste(x, collapse = "  "))
      wid <- sapply(tmp[1L, ], nchar, type = "w")
      blanks <- paste(character(max(wid)), collapse = " ")
      wcn <- ncw(cn)
      pad0 <- floor((wid - wcn)/2)
      pad1 <- wid - wcn - pad0
      cn <- paste0(substring(blanks, 1L, pad0), cn, substring(blanks,
                                                              1L, pad1))
      nm[i] <- paste(cn, collapse = "  ")
      z[[i]] <- sms
    }
    else {
      # sms <- format(sms, digits=digits, drop0trailing=TRUE)
      lbs <- format(names(sms))
      sms <- paste0(lbs, ":", sms, "  ")
      lw[i] <- ncw(lbs[1L])
      length(sms) <- nr
      z[[i]] <- sms
    }
  }
  if (nv) {
    z <- unlist(z, use.names = TRUE)
    dim(z) <- c(nr, nv)
    if (any(is.na(lw)))
      warning("probably wrong encoding in names(.) of column ",
              paste(which(is.na(lw)), collapse = ", "))
    blanks <- paste(character(max(lw, na.rm = TRUE) + 2L),
                    collapse = " ")
    pad <- floor(lw - ncw(nm)/2)
    nm <- paste0(substring(blanks, 1, pad), nm)
    dimnames(z) <- list(rep.int("", nr), nm)
  }
  else {
    z <- character()
    dim(z) <- c(nr, nv)
  }
  attr(z, "class") <- c("table")
  msage="Skewness, Kurtosis & Normality testing not done if SD=0 or n<6.\n"
  cat(msage)
  msage="Normality test is by Shapiro-Wilk\n"
  cat(msage)
  msage="If normality testing is done, p-value is displayed only if it is <0.1 ; a value of 0.0 indicates <0.0001\n\n"
  cat(msage)
  z
}

cvsum.default = function (object, ..., digits = max(3L, getOption("digits") -
                                                      3L))
{
  if (is.factor(object))
    return(summary.factor(object, ...))
  else if (is.matrix(object))
    return(summary.matrix(object, digits = digits, ...))
  value <- if (is.logical(object))
    c(Mode = "logical", {
      tb <- table(object, exclude = NULL)
      if (!is.null(n <- dimnames(tb)[[1L]]) && any(iN <- is.na(n))) dimnames(tb)[[1L]][iN] <- "NA's"
      tb
    })
  else if (is.numeric(object)) {
    nas <- is.na(object)
    object <- object[!nas]
    qtl <- stats::quantile(object)
    len=length(object)
    av=mean(object)
    stdev=sd(object)
    DoNorm=1
    if (len<6 | stdev==0) DoNorm=0
    cv=100*sd(object)/mean(object)
    se=sd(object)/sqrt(len)
    if (length(object)>3 & sd(object)!=0) {
      sk=skewness(object, type=2)
      kt=kurtosis(object, type=2)
      skp=2*(1-pt(abs(sk/(sqrt(6/len))), len-1))
      ktp=2*(1-pt(abs(kt/(sqrt(24/len))), len-1))
      sh=shapiro.test(object)
      sh2=sh[[2]] }
    #sh=as.numeric(shapiro.test(object))
    #cv shows as Inf if mean is zero

    qq <- c(len,signif(av,3),signif(stdev,3),round(cv,1),signif(se,2),signif(qtl[1L:5L],3))
    if (length(object)>3 & sd(object)!=0) {
      qq <- c(qq,round(sk,2),round(kt,2))
      if (skp<.1 | ktp<.1 | sh2<.1) qq <- c(qq,round(skp,3),round(ktp,3),round(sh2,3))
    }
    qq <- format(qq, digits=digits, drop0trailing=TRUE)
    nams = c("N","Mean","SD","CV%","SEM","Min.","1st Q","Median","3rd Q","Max.")
    names(qq) <- nams
    if (length(object)>3 & sd(object)!=0) {
      names(qq) <- c(nams,"Skew","Kurt.")
      if (skp<.1 | ktp<.1 | sh2<.1) names(qq) <- c(nams,"Skew","Kurt.","Skew p","Kurt p","norm ?")
    }
    if (any(nas))
      c(qq, `NA's` = sum(nas))
    else qq
  }
  else if (is.recursive(object) && !is.language(object) &&
           (n <- length(object))) {
    sumry <- array("", c(n, 3L), list(names(object), c("Length",
                                                       "Class", "Mode")))
    ll <- numeric(n)
    for (i in 1L:n) {
      ii <- object[[i]]
      ll[i] <- length(ii)
      cls <- oldClass(ii)
      sumry[i, 2L] <- if (length(cls))
        cls[1L]
      else "-none-"
      sumry[i, 3L] <- mode(ii)
    }
    sumry[, 1L] <- format(as.integer(ll))
    sumry
  }
  else c(Length = length(object), Class = class(object), Mode = mode(object))
  class(value) <- c("summaryDefault", "table")
  value
}


cv1way = function(x,group,ebars=1,dots=0,horiz=FALSE,padj="none",
                  namex=NULL, namey=NULL, namegroup=NULL, main=NULL,
                  xlab=NULL,ylab=NULL, log="", cilwd=1, dens=NULL,color="black",fill="grey",pos=0){
  if (is.null(namex)) namex=deparse(substitute(x))
  if (is.null(namegroup)) namegroup=deparse(substitute(group))
  if (is.null(main)) main=paste(namex, "at different levels of", namegroup)
  # xloc = x
  # grouploc = group
  # names(xloc) = namex
  # names(grouploc) = namegroup
  cat(paste("\n",namex,"compared across",namegroup,"groups","\n\n"))
  if (ebars !=4) {myaov = aov(x ~ group)
  print(summary(myaov))
  }
  if (ebars == 4) {myaov = kruskal.test(x ~ group)
  myaov$data.name = paste(deparse(substitute(x)),"by",deparse(substitute(group)))
  print(myaov)
  }
  cat("\n")
  desc = by(x,group,cvsum)
  print(desc)
  msage="If normality testing is done, p-value is displayed only if it is <0.1 ; a value of 0.0 indicates <0.0001\n\n"
  cat("\n", msage, "\n")

  normflag=F
  normmsg = "DATA FAIL NORMALITY TEST. LOOK FOR DATA ERRORS IN 'Min' AND 'Max' VALUES. IF DATA ARE NOT NORMAL,YOU MUST USE THE NONPARAMETRIC DUNNTEST (ebars=4)"
  shby = by(x,group,shapiro.test)
  nlev = nlevels(group)
  for (i in 1:nlev) {
    pval = shby[[i]]$p.value
    if (pval<.05) normflag = T
  }
  if (normflag && ebars != 4) cat("\n", normmsg, "\n")
  pool.sd=T
  bt=bartlett.test(x,group)
  bt$data.name = paste(namex,"across",namegroup,"groups")
  if(bt[3]<.05) pool.sd=F

  if (ebars != 4) {
    if (pool.sd) {
      pairout = pairwise.t.test(x,group,p.adjust=padj)
      g <- factor(group)
      s <- tapply(x, g, sd, na.rm = TRUE)
      n <- tapply(!is.na(x), g, sum)
      degf <- n - 1
      total.degf <- sum(degf)
      pooled.sd <- sqrt(sum(s^2 * degf)/total.degf)
      pairout$data.name = paste(namex,"compared across",namegroup,"groups")
      pairout$method = paste(pairout$method,signif(pooled.sd,3))
      print(pairout)
    }

    else {
      cat("\n","COULD NOT USE POOLED SD DUE TO UNEQUAL VARIANCES. LOOK FOR DATA ERRORS, FOCUSING ON GROUP(S) WITH LARGE SD AND UNEXPECTED MIN/MAX IN SUMMARY ABOVE","\n")
      print(bt)
      pairout = pairwise.t.test(x,group,p.adjust=padj,pool.sd=F)
      pairout$data.name = paste(namex,"compared across",namegroup,"groups")
      print(pairout)
    }

    cat("\nEach p-value compares groups above and to the left\n")
    #  cat(paste("\nPooled SD =", signif(pooled.sd,3), "\n"))
  }

  if (ebars==4) {
    pairout = kwAllPairsDunnTest(x,group,p.adjust.method="none")
    pairout$data.name = paste(namex,"compared across",namegroup,"groups")
    print(pairout)
  }
  if (is.null(xlab)) xlab=namegroup
  if (is.null(ylab)) ylab=namex
  ds = data.frame(x,group)
  if (ebars==1) {addon="mean_sd"
  main=paste(main,"\n Mean +/- S.D.")}
  if (ebars==2) {addon="mean_se"
  main=paste(main,"\n Mean +/- S.E.")}
  if (ebars==3) {addon="mean_ci"
  main=paste(main,"\n Mean and 95% CL")}
  if (ebars==4) {addon="median_iqr"
  main=paste(main,"\n Median and IQR")}
  if (dots==1) addon=c(addon, "dotplot")
  pal = "lancet"
  if (fill=="group") color = fill

  p = ggbarplot(ds,x="group",y="x",add = addon, color=color, fill=fill, xlab=namegroup, ylab=namex);
  ggpar(p, main=main, palette=pal) + theme(plot.title = element_text(hjust=0.5))
}

cv2way = function(x,group1,group2, ebars=1,dots=0, padj="none",
                  namex=NULL, namey=NULL, namegroup1=NULL, namegroup2=NULL, main=NULL,
                  log="", cilwd=1, dens=NULL,color="black",pos=0){

  if (is.null(namex)) namex=deparse(substitute(x))
  if (is.null(namegroup1)) namegroup1=deparse(substitute(group1))
  if (is.null(namegroup2)) namegroup2=deparse(substitute(group2))

  GroupPre = paste(group1, group2, sep="+")
  Group = factor(GroupPre)
  GrpNam = paste(namegroup1,"+",namegroup2,sep="")
  if (is.null(main)) main=paste(namex, "at different levels of", GrpNam)

  desc = by(x,Group,cvsum)
  print(desc)
  msage="If normality testing is done, p-value is displayed only if it is <0.1 ; a value of 0.0 indicates <0.0001\n\n"
  cat("\n", msage, "\n")

  normflag=F
  normmsg = "DATA FAIL NORMALITY TEST. LOOK FOR DATA ERRORS IN 'Min' AND 'Max' VALUES. IF DATA ARE NOT NORMAL,YOU MUST USE THE NONPARAMETRIC DUNNTEST (ebars=4)"
  shby = by(x,Group,shapiro.test)
  nlev = nlevels(Group)
  for (i in 1:nlev) {
    pval = shby[[i]]$p.value
    if (pval<.05) normflag = T
  }
  if (normflag && ebars != 4) cat("\n", normmsg, "\n")
  pool.sd=T
  bt=bartlett.test(x,Group)
  bt$data.name = paste(namex,"across",GrpNam,"groups")
  if(bt[3]<.05) pool.sd=F

  if (ebars != 4) {
    if (pool.sd) {
      pout = pairwise.t.test(x, Group, p.adjust=padj)
      pout$data.name = paste(namex,"compared across",GrpNam,"groups")
      print(pout)
    }
    else {
      cat("\n","COULD NOT USE POOLED SD DUE TO UNEQUAL VARIANCES. LOOK FOR DATA ERRORS, FOCUSING ON GROUP(S) WITH LARGE SD AND UNEXPECTED MIN/MAX IN SUMMARY ABOVE","\n")
      print(bt)
      pout = pairwise.t.test(x,Group,p.adjust=padj,pool.sd=F)
      pout$data.name = paste(namex,"compared across",GrpNam,"groups")
      print(pout)
    }
    cat("\nEach p-value compares groups above and to the left\n\n")
  }
  if (ebars==4) {
    pout = kwAllPairsDunnTest(x,Group,p.adjust.method="none")
    pout$data.name = paste(namex,"compared across",GrpNam,"groups")
    print(pout)
  }

  gp1 = nlevels(group1)
  gp2 = nlevels(group2)
  g1names = levels(group1)
  g2names = levels(group2)
  beg=gp2+gp1
  end=beg+((gp2-1)*(gp1-1)) - 1
  M = matrix(beg:end, nrow=gp2-1, ncol=gp1-1, byrow=T)
  #   print(M)
  mylm = lm(x ~ group1 + group2 + group1*group2)
  print(summary(mylm))
  for (Mrow in 0:(gp2-1)) {  # factor 1 effects at each level of factor 2
    cat(g2names[Mrow+1],"\n")
    for (j in 1:(gp1-1)) {   # each factor 1 level vs first
      vec = c(rep(0,end)); vec[j+1] = 1; if (Mrow > 0) vec[M[Mrow,j]] = 1
      est = estimable(mylm, vec)
      cat(g1names[j+1],"-",g1names[1],": ",signif(est[,1],3),"+/-",
          signif(est[,2],3),", p=",signif(est[,5],3),vec,"\n",sep=" ")
    }
    if (gp1>2) for (j in 1:(gp1-2)) {  # comparing factor 1 levels not first
      for (k in (j+1):(gp1-1)) {
        vec = c(rep(0,end)); vec[j+1] = -1; vec[k+1] = 1;
        if (Mrow > 0) vec[M[Mrow,c(j,k)]] = c(-1,1)
        est = estimable(mylm, vec)
        cat(g1names[k+1],"-",g1names[j+1],":",signif(est[,1],3),"+/-",
            signif(est[,2],3),", p=",signif(est[,5],3),vec,"\n",sep=" ")
      }
    }
    cat("\n")
  }
  for (Mrow in 1:(gp2-1)) {   # each level of factor 2 vs first
    cat(g2names[Mrow+1],"-",g2names[1],"\n")
    for (j in 1:gp1) {        # at each level of factor 1
      vec = c(rep(0,end)); vec[gp1+Mrow] = 1; if (j>1) vec[M[Mrow,j-1]] = 1
      est = estimable(mylm, vec)
      cat(g1names[j],":",signif(est[,1],3),"+/-",signif(est[,2],3),
          ", p=",signif(est[,5],3),vec,"\n",sep=" ")
    }
    for (j in 1:(gp1-1)) {   # each factor 1 level vs first
      vec = c(rep(0,end)); vec[M[Mrow,j]] = 1
      est = estimable(mylm, vec)
      cat(g1names[j+1],"-",g1names[1],": ",signif(est[,1],3),"+/-",
          signif(est[,2],3),", p=",signif(est[,5],3),vec,"\n",sep=" ")
    }
    if (gp1>2) for (j in 1:(gp1-2)) {  # comparing factor 1 levels not first
      for (k in (j+1):(gp1-1)) {
        vec = c(rep(0,end)); vec[M[Mrow,c(j,k)]] = c(-1,1)
        est = estimable(mylm, vec)
        cat(g1names[k+1],"-",g1names[j+1],":",signif(est[,1],3),"+/-",
            signif(est[,2],3),", p=",signif(est[,5],3),vec,"\n",sep=" ")
      }
    }
    cat("\n")
  }
  if (gp2>2) for (Mrow in 1:(gp2-2)) {  # comparing factor 2 levels not first
    for (j in (Mrow+1):(gp2-1)) {
      cat(g2names[j+1],"-",g2names[Mrow+1],"\n")
      for (k in 1:gp1) {        # at each level of factor 1
        vec = c(rep(0,end)); vec[gp1+Mrow] = -1; vec[gp1+j] = 1;
        if (k>1) {vec[M[Mrow,k-1]] = -1; vec[M[j,k-1]] = 1}
        est = estimable(mylm, vec)
        cat(g1names[k],":",signif(est[,1],3),"+/-",signif(est[,2],3),
            ", p=",signif(est[,5],3),vec,"\n",sep=" ")
      }
      for (k in 1:(gp1-1)) {   # each factor 1 level vs first
        vec = c(rep(0,end))
        vec[M[c(Mrow,j),k]] = c(-1,1)
        cat(g1names[k+1],"-",g1names[1],": ")
        est = estimable(mylm, vec)
        cat(" ",signif(est[,1],3),"+/-",signif(est[,2],3),
            ", p=",signif(est[,5],3),vec,"\n",sep=" ")
      }
      if (gp1>2) for (i in 1:(gp1-2)) {  # comparing factor 1 levels not first
        for (k in (i+1):(gp1-1)) {
          vec = c(rep(0,end))
          vec[M[c(Mrow,j),c(i,k)]] = c(1,-1,-1,1)
          cat(g1names[k+1],"-",g1names[j+1],": ")
          est = estimable(mylm, vec)
          cat(" ",signif(est[,1],3),"+/-",signif(est[,2],3),
              ", p=",signif(est[,5],3),vec,"\n",sep=" ")
        }
      }
      cat("\n")
    }
  }

  ds = data.frame(x,Group, group1, group2)
  if (ebars==1) {addon="mean_sd"
  main=paste(main,"\n Mean +/- S.D.")}
  if (ebars==2) {addon="mean_se"
  main=paste(main,"\n Mean +/- S.E.")}
  if (ebars==3) {addon="mean_ci"
  main=paste(main,"\n Mean and 95% CL")}
  if (ebars==4) {addon="median_iqr"
  main=paste(main,"\n Median and IQR")}
  if (dots==1) addon=c(addon, "dotplot")
  if (color=="yes") pal="lancet"
  if (color=="black") pal="grey"

  p = ggbarplot(ds,x="group1",y="x",add = addon, color="group2", fill="group2", xlab=namegroup1, ylab=namex, position = position_dodge(0.8))

  ggpar(p, main=main, legend.title=namegroup2, palette=pal) + theme(plot.title = element_text(hjust=0.5))

}

cvroc = function(LRobj, Case)
{
  DepVar = deparse(substitute(Case))
  levs = levels(Case)
  tablenote = paste("Note: In the logistic regression above, the logistic function is to predict", DepVar,"=",levs[2],"\n")
  cat(tablenote, "\n")
  mytab = table(fitted.values(LRobj) > 0.5, Case)
  print(mytab)
  tablenote = paste("\nFALSE is subjects with logistic fn<0, predicted to have",DepVar,"=",levs[1])
  cat(tablenote, "\n")
  tablenote = paste("TRUE is subjects with logistic fn>0, predicted to have",DepVar,"=",levs[2])
  cat(tablenote, "\n")
  tablenote = paste("True positive rate is the same as sensitivity, and false positive rate is 1-specificity")
  cat(tablenote, "\n")
  PredObj = prediction(fitted(LRobj), Case)
  perfAUC = performance(PredObj, measure="auc")
  AUCtitle =  paste(DepVar, "\n", "AUC = ", strtrim(perfAUC@"y.values",5))
  perfPlot = performance(PredObj, "tpr", "fpr")
  ROCR::plot(perfPlot, main=AUCtitle)
}


cwkm = function(timevar, statvar, ttmtvar, kmtype="Survival", title="Default", pvalue=T)
{
  ds = data.frame(Time=timevar,Status=statvar,Factor=ttmtvar)
  SurvTst = survdiff(Surv(Time,Status) ~ Factor, data=ds)
  facnam = deparse(substitute(ttmtvar))
  names(SurvTst$n) = c(paste(facnam,substr(names(SurvTst$n[1]),start=7,stop=30),sep=""), paste(facnam,substr(names(SurvTst$n[2]),start=7,stop=30),sep=""))
  SurvTst$call = paste("Kaplan-Meier model:",deparse(substitute(timevar)),"to",deparse(substitute(statvar)),"(event) depends on",facnam)
  print(SurvTst)
  SurvObj = survfit(Surv(Time,Status) ~ Factor, data=ds)
  cat("    by the log-rank or Mantel-Haenszel test\n")
  TtmtNam = deparse(substitute(ttmtvar))
  StatNam = deparse(substitute(statvar))
  ttl = paste("Kaplan-Meier curves for ",StatNam," at different levels of ",TtmtNam,sep="")
  if (title!="Default") ttl=title
  levs=levels(ttmtvar)
  if (kmtype=="Survival"){
    ggsurvplot(SurvObj, risk.table=T, title=ttl, data=ds, legend.title=TtmtNam, legend.labs=levs, pval=pvalue)
  }
  else {
    ggsurvplot(SurvObj, risk.table=T, fun="event", title=ttl, data=ds, legend.title=TtmtNam, legend.labs=levs, pval=pvalue)
  }
}


cvcov1way = function(depvar, covar, indvar, xs=NULL) {
  if (!is.factor(indvar) || !is.numeric(covar))
    stop("second argument should be numeric and third should be categorical")
  dvar = deparse(substitute(depvar))
  cvar = deparse(substitute(covar))
  ivar = deparse(substitute(indvar))
  #ds = data.frame(deparse(substitute(depvar)), deparse(substitute(covar)), deparse(substitute(indvar)))
  #ds = data.frame(dvar=depvar,cvar=covar,ivar=indvar)

  beg=rep(0,2);  levs=nlevels(indvar)
  levnams=levels(indvar); levnams[levs+1]=levnams[1]

  coeffs = c("(Intercept)",cvar)
  for (i in 2:levs) {
    coeffnam = paste(ivar,levnams[i],sep="")
    coeffs = c(coeffs,coeffnam)
  }

  if (is.null(xs)) {
    sets=1; zeros=levs
    fit = lm(depvar ~ covar+indvar)
    fit$call = paste("lm(formula = ",dvar," ~ ",cvar,"+",ivar,")",sep="")
  }

  else {
    sets=length(xs); zeros=levs*2
    fit = lm(depvar ~ covar+indvar+covar*indvar)
    fit$call = paste("lm(formula = ",dvar," ~ ",cvar," + ",ivar," + ",cvar,"*",ivar,sep="")
    for (i in 2:levs) {
      coeffnam = paste(cvar,":",ivar,levnams[i],sep="")
      coeffs = c(coeffs,coeffnam)
    }
  }
  names(fit$coefficients) = coeffs
  sumlm = summary(fit)
  print(sumlm)

  for (k in 1:sets) {
    for (i in 1:(levs-1)) {
      for (j in (i+1):(levs)) {
        vec=c(rep(0,zeros)); vec[i] = -1; vec[j] = 1
        contrasts = c(beg, vec[1:(levs-1)])
        if (!is.null(xs)) {
          x=xs[k]; vec[i+levs]=-x; vec[j+levs]=x
          contrasts=c(beg,vec[1:(levs-1)],vec[(levs+1):(zeros-1)])
        }
        est = estimable(fit, contrasts)
        cat("\n", levnams[j+1], "minus", levnams[i+1])
        if (!is.null(xs)) {cat(" ","at ",cvar,"=", x,sep="")}
        tval=signif(est[,3],3)
        pval=signif(est[,5],3)
        cat(" ",est[,1],est[,2],tval,pval,contrasts,sep="  ")
      }
    }
  }
  cat("\n\nFor the graph: \n")
  cat("ycalc = predict(fit) \n")
  cat("datf = data.frame(",dvar,",",ivar,",",cvar,",ycalc)\n", sep="")
  cat("ggplot(datf, aes(covar,color=",ivar,"))", sep="")
  cat(" + geom_point(aes(y=",dvar,")) + geom_line(aes(y=ycalc)) \n", sep="")

  ycalc = predict(fit)
  datf = data.frame(depvar, indvar, covar, ycalc)
  p <- ggplot(datf, aes(covar,color=indvar)) + geom_point(aes(y=depvar)) + geom_line(aes(y=ycalc))
  p + labs(x=cvar) + labs(colour=ivar) + labs(y=dvar)
}

cvrep1way = function(depvar, GroupVar, Subject) {
  dvar = deparse(substitute(depvar))
  rvar = deparse(substitute(GroupVar))
  ivar = deparse(substitute(Subject))
  beg=rep(0,1);  levs=nlevels(GroupVar)
  levnams=levels(GroupVar); levnams[levs+1]=levnams[1]
  zeros = levs
  fitrm = lme(depvar ~ GroupVar, random=~1|Subject/GroupVar)
  coeffs = "(Intercept)"
  for (i in 2:levs) {
    coeffnam = paste(rvar,levnams[i],sep="")
    coeffs = c(coeffs,coeffnam)
  }
  names(fitrm$coefficients$fixed) = coeffs
  #names(fitrm$coefficients$random) = c(rvar,ivar)
  fitrm$call$fixed[2]=dvar
  fitrm$call$fixed[3]=rvar
  sumlm = summary(fitrm)
  print(sumlm)

  for (i in 1:(levs-1)) {
    for (j in (i+1):(levs)) {
      vec=c(rep(0,zeros)); vec[i] = -1; vec[j] = 1
      contrasts = c(beg, vec[1:(levs-1)])
      est = estimable(fitrm, contrasts)
      cat("\n", levnams[j+1], "minus", levnams[i+1])
      cat(" ",signif(est[,1],3)," +/- ",signif(est[,2],3),"  p=",signif(est[,5],3),sep="")
    }
  }
  cat("\n")
}


cvrep2way = function(depvar, GroupVar, Factor, Subject) {
  dvar = deparse(substitute(depvar))
  rvar = deparse(substitute(GroupVar))
  ivar = deparse(substitute(Subject))
  fvar = deparse(substitute(Factor))
  beg=rep(0,1);  levs=nlevels(GroupVar)
  levnams=levels(GroupVar); levnams[levs+1]=levnams[1]
  zeros = levs
  mylm = lme(depvar ~ GroupVar+Factor+GroupVar*Factor, random=~1|Subject/GroupVar)
  # coeffs = "(Intercept)"
  # for (i in 2:levs) {
  #   coeffnam = paste(rvar,levnams[i],sep="")
  #   coeffs = c(coeffs,coeffnam)
  # }
  # names(fitrm$coefficients$fixed) = coeffs
  # #names(fitrm$coefficients$random) = c(rvar,ivar)
  # fitrm$call$fixed[2]=dvar
  # fitrm$call$fixed[3]=rvar
  sumlm = summary(mylm)
  print(sumlm)

  gp1 = nlevels(GroupVar)
  gp2 = nlevels(Factor)
  g1names = levels(GroupVar)
  g2names = levels(Factor)
  beg=gp2+gp1
  end=beg+((gp2-1)*(gp1-1)) - 1
  M = matrix(beg:end, nrow=gp2-1, ncol=gp1-1, byrow=T)
  #   print(M)
  cat("\nComparison of",rvar,"Effects Between Pairs of",fvar,"Groups\n    (Effect difference +/- SE, p-value)\n\n")

  for (Mrow in 0:(gp2-1)) {  # factor 1 effects at each level of factor 2
    cat(g2names[Mrow+1],"\n")
    for (j in 1:(gp1-1)) {   # each factor 1 level vs first
      vec = c(rep(0,end)); vec[j+1] = 1; if (Mrow > 0) vec[M[Mrow,j]] = 1
      est = estimable(mylm, vec)
      cat(g1names[j+1],"-",g1names[1],": ",signif(est[,1],3),"+/-",
          signif(est[,2],3),", p=",signif(est[,5],3),vec,"\n",sep=" ")
    }
    if (gp1>2) for (j in 1:(gp1-2)) {  # comparing factor 1 levels not first
      for (k in (j+1):(gp1-1)) {
        vec = c(rep(0,end)); vec[j+1] = -1; vec[k+1] = 1;
        if (Mrow > 0) vec[M[Mrow,c(j,k)]] = c(-1,1)
        est = estimable(mylm, vec)
        cat(g1names[k+1],"-",g1names[j+1],":",signif(est[,1],3),"+/-",
            signif(est[,2],3),", p=",signif(est[,5],3),vec,"\n",sep=" ")
      }
    }
    cat("\n")
  }
  for (Mrow in 1:(gp2-1)) {   # each level of factor 2 vs first
    cat(g2names[Mrow+1],"-",g2names[1],"\n")
    for (j in 1:gp1) {        # at each level of factor 1
      vec = c(rep(0,end)); vec[gp1+Mrow] = 1; if (j>1) vec[M[Mrow,j-1]] = 1
      est = estimable(mylm, vec)
      cat(g1names[j],":",signif(est[,1],3),"+/-",signif(est[,2],3),
          ", p=",signif(est[,5],3),vec,"\n",sep=" ")
    }
    for (j in 1:(gp1-1)) {   # each factor 1 level vs first
      vec = c(rep(0,end)); vec[M[Mrow,j]] = 1
      est = estimable(mylm, vec)
      cat(g1names[j+1],"-",g1names[1],": ",signif(est[,1],3),"+/-",
          signif(est[,2],3),", p=",signif(est[,5],3),vec,"\n",sep=" ")
    }
    if (gp1>2) for (j in 1:(gp1-2)) {  # comparing factor 1 levels not first
      for (k in (j+1):(gp1-1)) {
        vec = c(rep(0,end)); vec[M[Mrow,c(j,k)]] = c(-1,1)
        est = estimable(mylm, vec)
        cat(g1names[k+1],"-",g1names[j+1],":",signif(est[,1],3),"+/-",
            signif(est[,2],3),", p=",signif(est[,5],3),vec,"\n",sep=" ")
      }
    }
    cat("\n")
  }
  if (gp2>2) for (Mrow in 1:(gp2-2)) {  # comparing factor 2 levels not first
    for (j in (Mrow+1):(gp2-1)) {
      cat(g2names[j+1],"-",g2names[Mrow+1],"\n")
      for (k in 1:gp1) {        # at each level of factor 1
        vec = c(rep(0,end)); vec[gp1+Mrow] = -1; vec[gp1+j] = 1;
        if (k>1) {vec[M[Mrow,k-1]] = -1; vec[M[j,k-1]] = 1}
        est = estimable(mylm, vec)
        cat(g1names[k],":",signif(est[,1],3),"+/-",signif(est[,2],3),
            ", p=",signif(est[,5],3),vec,"\n",sep=" ")
      }
      for (k in 1:(gp1-1)) {   # each factor 1 level vs first
        vec = c(rep(0,end))
        vec[M[c(Mrow,j),k]] = c(-1,1)
        cat(g1names[k+1],"-",g1names[1],": ")
        est = estimable(mylm, vec)
        cat(" ",signif(est[,1],3),"+/-",signif(est[,2],3),
            ", p=",signif(est[,5],3),vec,"\n",sep=" ")
      }
      if (gp1>2) for (i in 1:(gp1-2)) {  # comparing factor 1 levels not first
        for (k in (i+1):(gp1-1)) {
          vec = c(rep(0,end))
          vec[M[c(Mrow,j),c(i,k)]] = c(1,-1,-1,1)
          cat(g1names[k+1],"-",g1names[j+1],": ")
          est = estimable(mylm, vec)
          cat(" ",signif(est[,1],3),"+/-",signif(est[,2],3),
              ", p=",signif(est[,5],3),vec,"\n",sep=" ")
        }
      }
      cat("\n")
    }
  }
}

