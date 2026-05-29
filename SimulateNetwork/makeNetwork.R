makeNetwork<-function(networkStructure=list(nStages=16,nNodesPerStage=8,rangeLink=2,probLink=0.5,strengthLink))
  {
  nvars<-networkStructure$nStages * networkStructure$nNodesPerStage
  # make a network
  fullLinks=matrix(0,nvars,nvars)
  t=1
  for (j in 1:networkStructure$nStages) {
    for (i in 1:networkStructure$nNodesPerStage) {
      if (j>1) {
        links=runif(networkStructure$rangeLink*2+1)<=networkStructure$probLink
        if (i<(networkStructure$rangeLink+1))     links[1:(networkStructure$rangeLink-i+1)]=FALSE
        if (i>(networkStructure$nNodesPerStage-networkStructure$rangeLink)) links[(networkStructure$nNodesPerStage-i+networkStructure$rangeLink*2):length(links)]=FALSE
        
        for (k in -networkStructure$rangeLink:networkStructure$rangeLink) {
          if (links[k+networkStructure$rangeLink+1] && (t+k<=nrow(fullLinks))) fullLinks[t+k,t-networkStructure$nNodesPerStage]=1
        }
      }
      t=t+1
    }
  }
  
  network<-links2Path(fullLinks,networkStructure$nStages,networkStructure$nNodesPerStage)
  
  
  if (networkStructure$fullModel) {
    rval<-rowSums(fullLinks) # how many inputs to each node
    rval[rval>0]<-sqrt(1/rval[rval>0])*networkStructure$strengthLink
    rval<-matrix(rval,nrow(fullLinks),ncol(fullLinks),byrow=FALSE)
  } else rval<-networkStructure$strengthLink
  
  fullLinks<-fullLinks*rval
  network$fullLinks<-fullLinks
  
  Stheta<-links2Stheta(fullLinks)
  d1<-matrix(diag(Stheta),nrow=nrow(Stheta),ncol=ncol(Stheta),byrow=TRUE)
  d2<-matrix(diag(Stheta),nrow=nrow(Stheta),ncol=ncol(Stheta),byrow=FALSE)
  Stheta<-Stheta/sqrt(d1)/sqrt(d2)
  network$Stheta<-Stheta
  
  return(network)
}