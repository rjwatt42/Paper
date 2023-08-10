simData<-function(an=NULL) {
  
  if (is.null(an)) {
    metaData<-list(result=list(rIV=my_data$r_s,nval=my_data$n,df1=my_data$df1))
    metaAnal<-list(meta_fixedAnal="random",meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE,append=FALSE)
    an<-runMetaAnalysis(metaAnal,metaData)
  }
  
  nStudies=length(my_data$n)
  ns<-NA
  model="Exp"
  k=an$exp$Kmax
  pNull=an$exp$Nullmax
  actualS<-an$exp$Smax
  sigOnly=TRUE
  
  my_sim_data<-makeStudies(nStudies,ns,model,k,NA,pNull,sigOnly)
  
}

makeStudies<-function(nStudies,ns,model="Exp",k=0.325,shape=NA,pNull=0.74,sigOnly=TRUE) {
  if (sigOnly) rpts=100 
  else rpts=1
  
  minN<-10
  if (is.na(ns) || length(ns)==0) {
    ns<-round(minN+rgamma(nStudies*rpts,shape=1.0,scale=(72-minN)/1.0))
    n_rand<-TRUE
  } else {
    n_rand<-FALSE
  }
  if (length(ns)==1) {ns<-rep(ns,nStudies*rpts)}

  if (n_rand) {
    switch (model,
            "Single"={
              zp<-rep(k,nStudies*rpts)
            },
            "Exp"={
              zp<-rexp(nStudies*rpts,1/k)
            },
            "Gauss"={
              zp<-rnorm(nStudies*rpts,0,k)
            },
            "Gamma"={
              zp<-rgamma(nStudies*rpts,shape=shape,scale=k/shape)
            }
    )
    zp<-zp*(runif(nStudies*rpts)>=pNull)
    zp<-zp*sign(rnorm(nStudies*rpts))
    zs<-zp+rnorm(nStudies*rpts,0,1/sqrt(ns-3))
    
    if (sigOnly) {
      zcrit<-atanh(p2r(alpha,ns))
      ns<-ns[abs(zs)>zcrit]
      zs<-zs[abs(zs)>zcrit]
    }
    zsAll<-zs[1:nStudies]
    nsAll<-ns[1:nStudies]
    psAll<-r2p(tanh(zsAll),nsAll)
  } else {
  zsAll<-c()
  nsAll<-c()
  psAll<-c()
  uniqueN<-unique(ns)    
  for (i in 1:length(uniqueN)) {
    nsi<-sum(ns==uniqueN[i])
    nsuse<-rep(uniqueN[i],nsi*rpts)
    switch (model,
            "Single"={
              zp<-rep(k,nsi*rpts)
            },
            "Exp"={
              zp<-rexp(nsi*rpts,1/k)
            },
            "Gauss"={
              zp<-rnorm(nsi*rpts,0,k)
            },
            "Gamma"={
              zp<-rgamma(nStudies*rpts,shape=shape,scale=k/shape)
            }
    )
    zp<-zp*(runif(nsi*rpts)>=pNull)
    zp<-zp*sign(rnorm(nsi*rpts))
    zs<-zp+rnorm(nsi*rpts,0,1/sqrt(ns-3))
    
    if (sigOnly) {
      ps<-(1-pnorm(abs(zs),0,1/sqrt(ns-3)))*2
      keep<-ps<alpha
      nsuse<-nsuse[keep]
      zs<-zs[keep]
      ps<-ps[keep]
    }
    if (nsi>length(nsuse))    print(c(nsi,length(nsuse)))
    zs<-zs[1:nsi]
    nsuse<-nsuse[1:nsi]
    ps<-ps[1:nsi]
    zsAll<-c(zsAll,zs)
    nsAll<-c(nsAll,nsuse)
    psAll<-c(psAll,ps)
  }
  }
  zsAll<-zsAll[1:nStudies]
  nsAll<-nsAll[1:nStudies]
  psAll<-psAll[1:nStudies]
  return(list(r_s=tanh(zsAll),n=nsAll,p=psAll,df1=nsAll*0+1))
}



