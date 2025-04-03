
nrows=16
ncols=8
rangeLink=2
probLink=0.5
strengthLink=0.44

ncount=100
ns=1000000


rs_show=0
for (ni in 1:ncount) {
  t=1
  fullLinks=matrix(0,nrows*ncols,nrows*ncols)

  for (j in 1:nrows) {
    for (i in 1:ncols) {
      if (j>1) {
        links=runif(rangeLink*2+1)<=probLink
        if (i<(rangeLink+1))     links[1:(rangeLink-i+1)]=FALSE
        if (i>(ncols-rangeLink)) links[(ncols-i+rangeLink*2):length(links)]=FALSE

        for (k in -rangeLink:rangeLink) {
        if (links[k+rangeLink+1]) fullLinks[t+k,t-ncols]=1
        }
      }
      
      t=t+1
    }  
  }

  
  testData=matrix(0,ns,ncols*nrows)
  for (t in seq(ncols*nrows,1,-1)) {
    use=which(fullLinks[,t]==1);
    if (length(use)==0) {
      testData[,t]=rnorm(ns);
    } else {
      usedVar=length(use)*(strengthLink^2)
      if (usedVar>1) usedVar=1
      z=0
      for (ui in 1:length(use)) 
        z=z+testData[,use[ui]]*strengthLink
      testData[,t]=z+
                sqrt(1-usedVar)*rnorm(ns)
    }
  }
  
  rs=abs(cor(testData))
  rs_show=c(rs_show, rs[upper.tri((rs))])

  zs=atan(rs_show)
  hist(zs)
  metaData<-list(result=list(rIV=rs_show,nval=ns))
  an<-runMetaAnalysis(metaAnal,metaData)
  showAnalysis(an,"")
  
}


