makeNetworkSample<-function(network,n,nRand,nSamples,remove=list(row=TRUE,gap=1),h,hist=TRUE) {
  
  main<-mainHist<-NULL

  use<-lower.tri(network$Stheta)
  # remove within-stage effects
  if (!is.null(remove)) {
    nstage<-length(network$stages)
    for (qs in 1:nstage) {
      qsr<-which((qs==(1:nstage) & remove$row) | abs((1:nstage)-qs)>(remove$gap+1))
      for (qr in qsr)
        use[network$stages[[qs]],network$stages[[qr]]]<-FALSE
    }
  }
  
  links<-network$fullLinks[use]
  zp<-atanh(network$Stheta[use])
  if (nSamples==0) nSamples<-length(zp)
  if (nRand) n<-getSampleSizes(nSamples,dist="Gamma",sN=n,sNRandSD=n*0.66,minN=5)
  if (length(n)<nSamples) n<-rep(n,nSamples)
  
  esd<-1/sqrt(n-3)
  err<-rnorm(nSamples,0,esd)
  
  use<-ceiling(runif(nSamples,0,1)*length(zp))
  links<-links[use]
  zp<-zp[use]*sign(runif(nSamples,-1,1))
  zs<-abs(zp+err)
  # power
  wp<-rn2w(tanh(zp),n)
  
  # find the significant results
  p<-r2p(tanh(zs),n)
  sig<-p<0.05
  zss<-zs[sig]
  zps<-zp[sig]
  wps<-rn2w(tanh(zps),n[sig])
  
  linkWhenSig<-sum(sig & (links!=0))/sum(sig)
  sigWhenLink<-sum(sig & (links!=0))/sum((links!=0))
  
  zeroWhenSig<-sum(sig & (zp==0))/sum(sig)
  sigWhenZero<-sum(sig & (zp==0))/sum((zp==0))
  
  main<-list(all=list(zp=zp,zs=zs,n=n,p=p,wp=wp,links=links),
             sig=sig,linkWhenSig=linkWhenSig,sigWhenLink=sigWhenLink,
             zeroWhenSig=zeroWhenSig,sigWhenZero=sigWhenZero,
             sigOnly=list(zp=zps[sig],zs=zs[sig],n=n[sig],p=p[sig],wp=wp[sig],links=links[sig])
  )
  
  if (hist) {
  h1a<-hist(zs[zs<max(h$breaks)],h$breaks,plot=FALSE)
  h1b<-hist(zss[zss<max(h$breaks)],h$breaks,plot=FALSE)
  
  h2a<-hist(abs(zp[abs(zp)<max(h$breaks)]),h$breaks,plot=FALSE)
  h2b<-hist(abs(zps[abs(zps)<max(h$breaks)]),h$breaks,plot=FALSE)
  
  h3a<-hist(wp,seq(0,1,length.out=101),plot=FALSE)
  h3b<-hist(wps,seq(0,1,length.out=101),plot=FALSE)

  mainHist<-list(
  h1=list(breaks=h$breaks,
          density=h1a$density,
          density1=h1b$density*sum(h1b$counts)/sum(h1a$counts)),
  h2=list(breaks=h$breaks,
          density=h2a$density,
          density1=h2b$density*sum(h2b$counts)/sum(h2a$counts)),
  h3=list(breaks=seq(0,1,length.out=101),
          density=h3a$density,
          density1=h3b$density*sum(h3b$counts)/sum(h3a$counts))
  )
  
  }
  
  return(c(main,mainHist))
  
}
