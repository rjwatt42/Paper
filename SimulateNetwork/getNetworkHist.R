getNetworkHist<-function(network,separateZeros,h) {
  
  use<-network$Stheta<1
  if (separateZeros) use<-use & network$Stheta>0
  zp<-atanh(network$Stheta[use])
  h1<-hist(zp[zp<max(h$breaks)],h$breaks,plot=FALSE)
  
  est<-fitdistrplus::fitdist(zp,"exp")
  
  return(list(hist=h1,zvals=zp,est=est,Stheta=Stheta,zeros=sum(Stheta==0)/sum(Stheta<1)))
}