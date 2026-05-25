makeNetwork<-function(nrows,ncols,rangeLink,probLink) {
  
  # make a network
  fullLinks=matrix(0,nrows*ncols,nrows*ncols)
  t=1
  for (j in 1:nrows) {
    for (i in 1:ncols) {
      if (j>1) {
        links=runif(rangeLink*2+1)<=probLink
        if (i<(rangeLink+1))     links[1:(rangeLink-i+1)]=FALSE
        if (i>(ncols-rangeLink)) links[(ncols-i+rangeLink*2):length(links)]=FALSE
        
        for (k in -rangeLink:rangeLink) {
          if (links[k+rangeLink+1] && (t+k<=nrow(fullLinks))) fullLinks[t+k,t-ncols]=1
        }
      }
      t=t+1
    }
  }
  return(fullLinks)
}