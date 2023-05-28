
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
      zcrit<-qnorm(1-alpha/2,0,1/sqrt(ns-3))
      ns<-ns[abs(zs)>zcrit]
      zs<-zs[abs(zs)>zcrit]
    }
    zsAll<-zs[1:nStudies]
    nsAll<-ns[1:nStudies]
  } else {
  zsAll<-c()
  nsAll<-c()
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
    }
    if (nsi>length(nsuse))    print(c(nsi,length(nsuse)))
    zs<-zs[1:nsi]
    nsuse<-nsuse[1:nsi]
    zsAll<-c(zsAll,zs)
    nsAll<-c(nsAll,nsuse)
  }
  }
  return(list(r_s=tanh(zsAll),n=nsAll))
}


############################
# an example

nStudies=92877
ns<-my_data$n
ns<-NA
model="Exp"
k=0.3195403
pNull=0.725488
sigOnly=TRUE

my_sim_data<-makeStudies(nStudies,ns,model,k,NA,pNull,sigOnly)

metaAnal<-list(meta_fixedAnal="random",meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE,append=FALSE)
metaData<-list(result=list(rIV=my_sim_data$r_s,nval=my_sim_data$n))
# metaData<-list(result=list(rIV=my_data$r_s,nval=my_data$n))

an<-runMetaAnalysis(metaAnal,metaData)
showAnalysis(an)


############################
# next - 100 examples

nsims<-1000

#first we get the analysis for the real data
metaAnal<-list(meta_fixedAnal="random",meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE,append=FALSE)
metaData<-list(result=list(rIV=my_data$r_s,nval=my_data$n))

an<-runMetaAnalysis(metaAnal,metaData)

nStudies=length(my_data$n)
ns<-NA
model="Exp"
k=an$exp$Kmax
pNull=an$exp$Nullmax
actualS<-an$exp$Smax
sigOnly=TRUE


s<-c()
np<-c()
km<-c()
times<-c()
for (i in 1:nsims) {
  start<-Sys.time()
  my_data<-makeStudies(nStudies,ns,model,k,shape,pNull,sigOnly)
  metaAnal<-list(meta_fixedAnal="random",meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE,append=FALSE)
  metaData<-list(result=list(rIV=my_data$r_s,nval=my_data$n))
  
  an<-runMetaAnalysis(metaAnal,metaData)
  
  times<-c(times,Sys.time()-start)
  print(c(i,an$exp$Smax,times[i],(nsims-i)*mean(times)/60))
  s<-c(s,an$exp$Smax)
  km<-c(km,an$exp$Kmax)
  np<-c(np,an$exp$Nullmax)
}

g<-ggplot()
pts<-data.frame(s=s)
g<-g+geom_histogram(data=pts,aes(x=s, after_stat(ndensity)),color="white",fill="white")
g<-g+geom_vline(xintercept=actualS,color="red")
g<-g+geom_label(data=data.frame(x=actualS,y=1,label=paste0("Actual data = ",actualS)),aes(x=x,y=y,label=label))
g<-g+xlab("log(lk)")+ylab("Density")

g+plotTheme




