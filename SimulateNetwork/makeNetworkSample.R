makeNetworkSample<-function(network,n,h,hist=TRUE,doReplication=FALSE,repPower) {
  
  main<-repl<-mainHist<-replHist<-NULL

  nSamples<-length(n)
  
  use<-network$Stheta<1
  zp<-atanh(network$Stheta[use])
  
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
  
  if (doReplication) {
    zpr<-zp[sig]
    
    nr<-rw2n(tanh(zss),repPower)
    esd<-1/sqrt(nr-3)
    err<-rnorm(length(nr),0,esd)
    zsr<-abs(zpr+err)
    
    pr<-r2p(tanh(zsr),nr)
    sigr<-pr<0.05
    wpr<-rn2w(tanh(zpr),nr)
    
    zssr<-zsr[sigr]
    zprs<-zpr[sigr]
    wprs<-rn2w(tanh(zprs),nr[sigr])
    
    repl=list(allrep=list(zp=zpr,zs=zsr,n=nr,p=pr,wp=wpr),
              sigrep=sigr,
              sigOnlyrep=list(zp=zprs,zs=zssr,n=nr[sigr],p=pr[sigr],wp=wprs)
    )
  } 
  
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
  
  if (doReplication) {
    h1ra<-hist(zsr[zsr<max(h$breaks)],h$breaks,plot=FALSE)
    h1rb<-hist(zssr[zssr<max(h$breaks)],h$breaks,plot=FALSE)
    
    h2ra<-hist(abs(zpr[abs(zpr)<max(h$breaks)]),h$breaks,plot=FALSE)
    h2rb<-hist(abs(zprs[abs(zprs)<max(h$breaks)]),h$breaks,plot=FALSE)
    
    h3ra<-hist(wpr,seq(0,1,length.out=101),plot=FALSE)
    h3rb<-hist(wprs,seq(0,1,length.out=101),plot=FALSE)

    replHist<-list(
      h1rep=list(breaks=h$breaks,
                 density=h1a$density,
                 density1=h1b$density*sum(h1b$counts)/sum(h1a$counts),
                 density2=h1rb$density*sum(h1rb$counts)/sum(h1a$counts)
      ),
      h2rep=list(breaks=h$breaks,
                 density=h2a$density,
                 density1=h2b$density*sum(h2b$counts)/sum(h2a$counts),
                 density2=h2rb$density*sum(h2rb$counts)/sum(h2a$counts)
      ),
      h3rep=list(breaks=seq(0,1,length.out=101),
              density=h3ra$density,
              density1=h3rb$density*sum(h3rb$counts)/sum(h3ra$counts))
    )
  } 
  } 
  
  return(c(main,repl,mainHist,replHist))
  
}