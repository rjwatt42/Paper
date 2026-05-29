getNetworkHist<-function(network,separateZeros,h) {
  
  use<-lower.tri(network$Stheta)
  if (separateZeros) use<-use & network$Stheta>0
  zp<-atanh(network$Stheta[use])
  h1<-hist(zp[zp<max(h$breaks)],h$breaks,plot=FALSE)
  est<-fitdistrplus::fitdist(zp,"exp")
  
  use<-network$fullLinks!=0
  zp<-atanh(network$fullLinks[use])
  h2<-hist(zp[zp<max(h$breaks)],h$breaks,plot=FALSE)
  
  return(list(network=network,hist=h1,hist1=h2,zvals=zp,est=est,Stheta=network$Stheta,zeros=sum(network$Stheta==0)/sum(network$Stheta<1)))
}