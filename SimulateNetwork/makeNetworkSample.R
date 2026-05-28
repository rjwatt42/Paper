makeNetworkSample<-function(network,n,nSamples,h,hist=TRUE) {
  
  main<-mainHist<-NULL

  # use<-(abs(network$Stheta)<1)
  use<-lower.tri(network$Stheta)
  zp<-atanh(network$Stheta[use])
  if (nSamples==0) nSamples<-length(zp)
  if (length(n)<nSamples) n<-rep(n,nSamples)
  
  esd<-1/sqrt(n-3)
  err<-rnorm(nSamples,0,esd)
  
  use<-ceiling(runif(nSamples,0,1)*length(zp))
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
  main<-list(all=list(zp=zp,zs=zs,n=n,p=p,wp=wp),
             sig=sig,
             sigOnly=list(zp=zps[sig],zs=zs[sig],n=n[sig],p=p[sig],wp=wp[sig])
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