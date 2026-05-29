makeNetworkReplication<-function(samples,repPower,nIncrease,h,hist=TRUE) {
  
  repl<-replHist<-NULL

  zp<-samples$all$zp
  zs<-samples$all$zs  
  n<-samples$all$n
  wp<-samples$all$wp
  links<-samples$all$links
  sig<-samples$sig
  
  ns<-n[sig]
    zps<-zp[sig]
    zss<-zs[sig]
    wps<-wp[sig]
    linkss<-links[sig]
    
    zpr<-zps
    nr<-rw2n(tanh(zss),repPower)
    if (nIncrease) nr[nr<ns]<-ns[nr<ns]
    esd<-1/sqrt(nr-3)
    err<-rnorm(length(nr),0,esd)
    zsr<-abs(zpr+err)
    
    pr<-r2p(tanh(zsr),nr)
    sigr<-pr<0.05
    wpr<-rn2w(tanh(zpr),nr)
    
    zssr<-zsr[sigr]
    zprs<-zpr[sigr]
    wprs<-rn2w(tanh(zprs),nr[sigr])
    
    linkWhenRepl<-sum(sigr & (linkss!=0))/sum(sigr)
    replWhenLink<-sum(sigr & (linkss!=0))/sum((links!=0))
    
    zeroWhenRepl<-sum(sigr & (zpr==0))/sum(sigr)
    replWhenZero<-sum(sigr & (zpr==0))/sum(zp==0)
    # print(c(zeroWhenRepl,replWhenZero))
    
    repl=list(allrep=list(zp=zpr,zs=zsr,n=nr,p=pr,wp=wpr,no=n[sig]),
              sigrep=sigr,replWhenLink=replWhenLink,linkWhenRepl=linkWhenRepl,
              zeroWhenRepl=zeroWhenRepl,replWhenZero=replWhenZero,
              sigOnlyrep=list(zp=zprs,zs=zssr,n=nr[sigr],p=pr[sigr],wp=wprs)
    )

  if (hist) {
    h1a<-hist(zs[zs<max(h$breaks)],h$breaks,plot=FALSE)
    h1b<-hist(zss[zss<max(h$breaks)],h$breaks,plot=FALSE)
    
    h2a<-hist(abs(zp[abs(zp)<max(h$breaks)]),h$breaks,plot=FALSE)
    h2b<-hist(abs(zps[abs(zps)<max(h$breaks)]),h$breaks,plot=FALSE)

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
  
  return(c(repl,replHist))
  
}