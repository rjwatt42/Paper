makeNetwork<-function(network=list(nStages=16,nNodesPerStage=8,rangeLink=2,probLink=0.5))
  {
  nvars<-network$nStages * network$nNodesPerStage
  # make a network
  fullLinks=matrix(0,nvars,nvars)
  t=1
  for (j in 1:network$nStages) {
    for (i in 1:network$nNodesPerStage) {
      if (j>1) {
        links=runif(network$rangeLink*2+1)<=network$probLink
        if (i<(network$rangeLink+1))     links[1:(network$rangeLink-i+1)]=FALSE
        if (i>(network$nNodesPerStage-network$rangeLink)) links[(network$nNodesPerStage-i+network$rangeLink*2):length(links)]=FALSE
        
        for (k in -network$rangeLink:network$rangeLink) {
          if (links[k+network$rangeLink+1] && (t+k<=nrow(fullLinks))) fullLinks[t+k,t-network$nNodesPerStage]=1
        }
      }
      t=t+1
    }
  }
  return(fullLinks)
}