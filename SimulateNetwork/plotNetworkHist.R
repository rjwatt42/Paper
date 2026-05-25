plotNetworkHist<-function(h1,type="pop") {
  switch(type,
         "pop"={
           g1<-dataGraph(h1$hist,xlabel="z[p]",ylabel="density",
                         title=paste0("density=exp(-z[p]/",format(1/h1$est$estimate,digits=3),")",
                                      "   zeros=",format(h1$zeros*100,digits=3),"%"),
                         hist=TRUE)
           
           zi<-h1$hist$breaks
           zd<-exp(-zi*h1$est$estimate)
           zd<-zd/sum(zd)*sum(h1$hist$density)
           g1<-addG(g1,dataLine(data.frame(x=zi,y=zd),colour="#FF4400",linewidth=1.2))
         },
         "samp"={
           g1<-dataGraph(h1$hist,xlabel="z[s]",ylabel="density",
                         hist=TRUE)
         }
  )
  return(g1)
}