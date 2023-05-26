
makeStudies<-function(nStudies,model="Exp",k=0.325,pNull=0.74,sigOnly=TRUE) {
  if (sigOnly) rpts=100 
  else rpts=1
  
  ns<-round(rgamma(nStudies*rpts,shape=1.87,rate=1/30)+10)
  switch (model,
          "Exp"={
            zp<-rexp(nStudies*rpts,1/k)
          }
  )
  zp<-zp*(runif(nStudies*rpts)>=pNull)
  zp<-zp*sign(rnorm(nStudies*rpts))
  zs<-zp+rnorm(nStudies*rpts,0,1/sqrt(ns-3))
  
  if (sigOnly) {
    p<-r2p(tanh(zs),ns)
    zs<-zs[p<alpha]
    ns<-ns[p<alpha]
  }
  zs<-zs[1:nStudies]
  ns<-ns[1:nStudies]
  return(list(r_s=tanh(zs),n=ns))
}


############################
# an example

nStudies=92877
model="Exp"
k=0.325
pNull=0.72
sigOnly=TRUE

my_data<-makeStudies(nStudies,model,k,pNull,sigOnly)



