includeSingle=FALSE

getMaxLikelihood<-function(zs,ns,df1,dist,metaAnalysis) {
  nkpoints<-13
  nnullpoints<-13
  niterations<-1
  
  reInc<-(nkpoints-1)/2/3

  if (metaAnalysis$meta_nullAnal) {
    nullvals<-seq(0,1,length.out=nnullpoints)
  } else {
    nullvals<-0
  }

  if (dist=="Single") {
    kvals<-seq(-1,1,length.out=nkpoints)*0.95
  } else {
    # for kvals any smaller we hit underflow issues
    kvals<-seq(0.02,1,length.out=nkpoints)
  }
  
  # preliminary scan to find rough solution and bounds 
  for (re in 1:niterations) {
    S<-getLogLikelihood(zs,ns,df1,dist,kvals,nullvals,metaAnalysis$meta_psigAnal)
    # S of 0 indicates numerical underflow
    Smax<-max(S[S!=0],na.rm=TRUE)

    use<-which(S==Smax, arr.ind = TRUE)
    Nullmax<-nullvals[use[1,2]]
    lb2<-nullvals[max(1,use[1,2]-reInc)]
    ub2<-nullvals[min(length(nullvals),use[1,2]+reInc)]
    Kmax<-kvals[use[1,1]]
    lb1<-kvals[max(1,use[1,1]-reInc)]
    ub1<-kvals[min(length(kvals),use[1,1]+reInc)]
    
    kvals<-seq(lb1,ub1,length.out=nkpoints)
    if (metaAnalysis$meta_nullAnal) {
      nullvals<-seq(lb2,ub2,length.out=nnullpoints)
    }
  }

  # main minimization
  fun<-function(x) { -getLogLikelihood(zs,ns,df1,dist,x[1],x[2],metaAnalysis$meta_psigAnal)}
  result<-fmincon(c(Kmax,Nullmax),fun,ub=c(ub1,ub2),lb=c(lb1,lb2))
  Kmax<-result$par[1]
  Nullmax<-result$par[2]
  Smax<- -result$value
  
  # cross sections
  nullvals<-seq(0,1,length.out=65)
  kvals<-seq(0.02,1,length.out=65)
  SnullX<-getLogLikelihood(zs,ns,df1,dist,Kmax,nullvals,metaAnalysis$meta_psigAnal)
  SkX<-getLogLikelihood(zs,ns,df1,dist,kvals,Nullmax,metaAnalysis$meta_psigAnal)

  nullCI<-0
  kCI<-0
  # use<-nullvals<Nullmax
  # if (Nullmax==0) {
  #   nullCI<-c(0,
  #             approx(SnullX[!use],nullvals[!use],Smax-log(100))$y)
  # } else {
  #   nullCI<-c(approx(SnullX[use],nullvals[use],Smax-log(100))$y,
  #             approx(SnullX[!use],nullvals[!use],Smax-log(100))$y)
  # }
  # 
  # use<-kvals<Kmax
  # kCI<-c(approx(SkX[use],nullvals[use],Smax-log(100))$y,
  #           approx(SkX[!use],kvals[!use],Smax-log(100))$y)
  
  return(list(Kmax=Kmax,Nullmax=Nullmax,Smax=Smax,kCI=kCI,nullCI=nullCI,
              SX=list(kvals=kvals,SkX=as.vector(SkX),nullvals=nullvals,SnullX=as.vector(SnullX))))
}


runMetaAnalysis<-function(metaAnalysis,metaResult){

  rs<-metaResult$result$rIV
  zs<-atanh(rs)
  ns<-metaResult$result$nval
  df1<-metaResult$result$df1
  
    # doing random effects analysis
    if(metaAnalysis$meta_pdf[1]=="All"){
      metaAnalysis$meta_pdf<-c("Exp","Gauss")
    }
    
    single<-list(Kmax=NA,Nullmax=NA,Smax=NA,SX=NA)
    gauss<-list(Kmax=NA,Nullmax=NA,Smax=NA,SX=NA)
    exp<-list(Kmax=NA,Nullmax=NA,Smax=NA,SX=NA)
    gamma<-list(Kmax=NA,Nullmax=NA,Smax=NA,SX=NA)
    genexp<-list(Kmax=NA,Nullmax=NA,Smax=NA,SX=NA)
    
    for (i in 1:length(metaAnalysis$meta_pdf)) {
      d<-getMaxLikelihood(zs,ns,df1,metaAnalysis$meta_pdf[i],metaAnalysis)
      switch(metaAnalysis$meta_pdf[i],
             "Single"={single<-d},
             "Gauss"={gauss<-d},
             "Exp"={exp<-d},
             "Gamma"={gamma<-d},
             "GenExp"={genexp<-d},
      )
    }

  use<-which.max(c(single$Smax,gauss$Smax,exp$Smax,gamma$Smax,genexp$Smax))
  bestDist<-c("Single","Gauss","Exp","Gamma","GenExp")[use]
  best<-cbind(single,gauss,exp,gamma,genexp)[,use]

  metaResult<-list(single=single,
                   gauss=gauss,
                   exp=exp,
                   gamma=gamma,
                   genexp=genexp,
                   bestDist=bestDist,
                   best=best,
                   count=length(metaResult$bestDist)+1,
                   nsims=metaResult$nsims,
                   metaAnalysis=metaAnalysis,
                   result=metaResult$result
  )
  return(metaResult)
}



