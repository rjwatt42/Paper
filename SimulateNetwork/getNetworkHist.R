getNetworkHist<-function(fullLinks,fullModel,strengthLink,avoidZeros,h) {
  
  if (fullModel) {
    rval<-rowSums(fullLinks) # how many inputs to each node
    rval[rval>0]<-sqrt(1/rval[rval>0])*strengthLink
    rval<-matrix(rval,nrow(fullLinks),ncol(fullLinks),byrow=FALSE)
  } else rval<-strengthLink
  
  Stheta<-links2Stheta(fullLinks*rval)
  use<-Stheta<1
  if (avoidZeros) use<-use & Stheta>0
  zp<-atanh(Stheta[use])
  h1<-hist(zp[zp<max(h$breaks)],h$breaks,plot=FALSE)
  
  est<-fitdistrplus::fitdist(zp,"exp")
  
  return(list(hist=h1,zvals=zp,est=est,Stheta=Stheta,zeros=sum(Stheta==0)/sum(Stheta<1)))
}