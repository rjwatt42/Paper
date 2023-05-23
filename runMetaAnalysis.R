
getMaxS<-function(S,kvals) {
  np<-length(kvals)
  if (S[1]==max(S)) {
    return(kvals[1])
  }
  if (S[np]==max(S)) {
    return(kvals[np])
  }
  return(approx(diff(S),(kvals[1:(np-1)]+kvals[2:np])/2,0)$y)
}

runMetaAnalysis<-function(metaAnalysis,metaResult){
  rs<-metaResult$result$rIV
  zs<-atanh(rs)
  ns<-metaResult$result$nval
  nkpoints<-13
  nnullpoints<-13
  niterations<-3
  reInc<-(nkpoints-1)/2/2
  
  
  if (metaAnalysis$meta_fixedAnal=="fixed") {
    kvals<-seq(-1,1,length.out=nkpoints)*0.95
    singleS<-getLogLikelihood(zs,ns,"Single",kvals,0,metaAnalysis$meta_psigAnal)
    singleKmax<-approx(diff(singleS),(kvals[1:(nkpoints-1)]+kvals[2:nkpoints])/2,0)$y
    singleSmax<-max(singleS,na.rm=TRUE)
    singleNullmax<-0
    
    bestDist<-"Single"
    bestK<-singleKmax
    bestNull<-0
    bestS<-singleSmax
    
    gaussKmax<-NA
    gaussSmax<-NA
    gaussNullmax<-NA
    expKmax<-NA
    expSmax<-NA
    expNullmax<-NA
    
  } else {
    
    # doing random effects analysis
    singleSmax<-NA
    singleKmax<-NA
    singleNullmax<-NA
    gaussSmax<-NA
    gaussKmax<-NA
    gaussNullmax<-NA
    expSmax<-NA
    expKmax<-NA
    expNullmax<-NA
    gammaSmax<-NA
    gammaKmax<-NA
    gammaNullmax<-NA
    
    if (metaAnalysis$meta_nullAnal) {
      startNullvals<-seq(0,1,length.out=nnullpoints)
    } else {
      startNullvals<-0
    }
    
    # find best Single
    if (metaAnalysis$meta_pdf=="Single" || metaAnalysis$meta_pdf=="All") {
      nullvals<-startNullvals
      kvals<-seq(-1,1,length.out=nkpoints)*0.95
      for (re in 1:niterations) {
        singleS<-getLogLikelihood(zs,ns,"Single",kvals,nullvals,metaAnalysis$meta_psigAnal)
        singleSmax<-max(singleS,na.rm=TRUE)
        use<-which(singleS==singleSmax, arr.ind = TRUE)
        singleNullmax<-nullvals[use[1,2]]
        singleKmax<-kvals[use[1,1]]
        kvals<-seq(kvals[max(1,use[1,1]-reInc)],kvals[min(use[1,1]+reInc,nkpoints)],length.out=nkpoints)
        nullvals<-seq(nullvals[max(1,use[1,2]-reInc)],nullvals[min(use[1,2]+reInc,nnullpoints)],length.out=nnullpoints)
      }
    }
    
    # find best Gauss
    if (metaAnalysis$meta_pdf=="Gauss" || metaAnalysis$meta_pdf=="All") {
      nullvals<-startNullvals
      kvals<-seq(0.01,1,length.out=nkpoints)
      for (re in 1:niterations) {
        gaussS<-getLogLikelihood(zs,ns,"Gauss",kvals,nullvals,metaAnalysis$meta_psigAnal)
        gaussSmax<-max(gaussS,na.rm=TRUE)
        use<-which(gaussS==gaussSmax, arr.ind = TRUE)
        gaussNullmax<-nullvals[use[1,2]]
        gaussKmax<-kvals[use[1,1]]
        kvals<-seq(kvals[max(1,use[1,1]-reInc)],kvals[min(use[1,1]+reInc,nkpoints)],length.out=nkpoints)
        nullvals<-seq(nullvals[max(1,use[1,2]-reInc)],nullvals[min(use[1,2]+reInc,nnullpoints)],length.out=nnullpoints)
      }
    }
    
    # find best Exp
    if (metaAnalysis$meta_pdf=="Exp" || metaAnalysis$meta_pdf=="All") {
      nullvals<-startNullvals
      kvals<-seq(0.01,1,length.out=nkpoints)
      for (re in 1:niterations) {
        expS<-getLogLikelihood(zs,ns,"Exp",kvals,nullvals,metaAnalysis$meta_psigAnal)
        expSmax<-max(expS,na.rm=TRUE)
        use<-which(expS==expSmax, arr.ind = TRUE)
        expNullmax<-nullvals[use[1,2]]
        expKmax<-kvals[use[1,1]]
        kvals<-seq(kvals[max(1,use[1,1]-reInc)],kvals[min(use[1,1]+reInc,nkpoints)],length.out=nkpoints)
        nullvals<-seq(nullvals[max(1,use[1,2]-reInc)],nullvals[min(use[1,2]+reInc,nnullpoints)],length.out=nnullpoints)
      }
    }
    
    # find best Gamma
    if (metaAnalysis$meta_pdf=="Gamma") {
      nullvals<-startNullvals
      kvals<-seq(0.01,1,length.out=nkpoints)
      for (re in 1:niterations) {
        gammaS<-getLogLikelihood(zs,ns,"Gamma",kvals,nullvals,metaAnalysis$meta_psigAnal,metaAnalysis$gamma_shape)
        gammaSmax<-max(gammaS,na.rm=TRUE)
        use<-which(gammaS==gammaSmax, arr.ind = TRUE)
        gammaNullmax<-nullvals[use[1,2]]
        gammaKmax<-kvals[use[1,1]]
        kvals<-seq(kvals[max(1,use[1,1]-reInc)],kvals[min(use[1,1]+reInc,nkpoints)],length.out=nkpoints)
        nullvals<-seq(nullvals[max(1,use[1,2]-reInc)],nullvals[min(use[1,2]+reInc,nnullpoints)],length.out=nnullpoints)
      }
    }
    
    use<-which.max(c(singleSmax,gaussSmax,expSmax,gammaSmax))
    bestDist<-c("Single","Gauss","Exp","Gamma")[use]
    bestK<-c(singleKmax,gaussKmax,expKmax,gammaKmax)[use]
    bestNull<-c(singleNullmax,gaussNullmax,expNullmax,gammaNullmax)[use]
    bestS<-c(singleSmax,gaussSmax,expSmax,gammaSmax)[use]
  }
  
  
  if (metaAnalysis$append) {
    bestDist<-c(metaResult$bestDist,bestDist)
    bestK<-c(metaResult$bestK,bestK)
    bestNull<-c(metaResult$bestNull,bestNull)
    bestS<-c(metaResult$bestS,bestS)
    singleKmax<-c(metaResult$single$kmax,singleKmax)
    singleSmax<-c(metaResult$single$Smax,singleSmax)
    singleNullmax<-c(metaResult$single$nullMax,singleNullmax)
    gaussKmax<-c(metaResult$gauss$kmax,gaussKmax)
    gaussSmax<-c(metaResult$gauss$Smax,gaussSmax)
    gaussNullmax<-c(metaResult$gauss$nullMax,gaussNullmax)
    expKmax<-c(metaResult$exp$kmax,expKmax)
    expSmax<-c(metaResult$exp$Smax,expSmax)
    expNullmax<-c(metaResult$exp$nullMax,expNullmax)
    gammaKmax<-c(metaResult$gamma$kmax,gammaKmax)
    gammaSmax<-c(metaResult$gamma$Smax,gammaSmax)
    gammaNullmax<-c(metaResult$gamma$nullMax,gammaNullmax)
  }
  
  metaResult<-list(single=list(kmax=singleKmax,Smax=singleSmax,nullMax=singleNullmax),
                   gauss=list(kmax=gaussKmax,Smax=gaussSmax,nullMax=gaussNullmax),
                   exp=list(kmax=expKmax,Smax=expSmax,nullMax=expNullmax),
                   gamma=list(kmax=gammaKmax,Smax=gammaSmax,nullMax=gammaNullmax),
                   bestDist=bestDist,
                   bestK=bestK,
                   bestNull=bestNull,
                   bestS=bestS,
                   count=length(metaResult$bestDist)+1,
                   nsims=metaResult$nsims,
                   metaAnalysis=metaAnalysis,
                   result=metaResult$result
  )
  return(metaResult)
}



