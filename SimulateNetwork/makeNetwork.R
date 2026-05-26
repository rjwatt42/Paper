makeNetwork<-function(network=list(nStages=16,nVarsPerStage=8,rangeLink=2,probLink=0.5))
  {
  nvars<-network$nStages * network$nVarsPerStage
  # make a network
  fullLinks=matrix(0,nvars,nvars)
  t=1
  for (j in 1:network$nStages) {
    for (i in 1:network$nVarsPerStage) {
      if (j>1) {
        links=runif(network$rangeLink*2+1)<=network$probLink
        if (i<(network$rangeLink+1))     links[1:(network$rangeLink-i+1)]=FALSE
        if (i>(network$nVarsPerStage-network$rangeLink)) links[(network$nVarsPerStage-i+network$rangeLink*2):length(links)]=FALSE
        
        for (k in -network$rangeLink:network$rangeLink) {
          if (links[k+network$rangeLink+1] && (t+k<=nrow(fullLinks))) fullLinks[t+k,t-network$nVarsPerStage]=1
        }
      }
      t=t+1
    }
  }
  return(fullLinks)
}