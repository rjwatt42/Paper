plotNetworkHist<-function(h1,type="pop") {
  switch(type,
         "pop"={
           legend<-c(#paste0("density=exp(-z[p]/",format(1/h1$est$estimate,digits=3),")"),
                     paste0("    p(z[p]>0) = ",format((1-h1$zeros)*100,digits=3),"%"),
                     paste0("    p(link) = ",format((1-mean(h1$network$fullLinks==0))*100,digits=3),"%")
           )
           h1$hist$density1=h1$hist1$density*sum(h1$hist1$counts)/sum(h1$hist$counts)
           g1<-dataGraph(h1$hist,xlabel="z[p]",ylabel="density",
                         fill<-c("white","white","#40F"),
                         legend=list(names=legend,colours=c("white","#40F")),
                         hist=TRUE)
           
           # zi<-h1$hist$breaks
           # zd<-exp(-zi*h1$est$estimate)
           # zd<-zd/sum(zd)*sum(h1$hist$density)
           # g1<-addG(g1,dataLine(data.frame(x=zi,y=zd),colour="#FF4400",linewidth=1.2))
         },
         "samp"={
           g1<-dataGraph(h1$hist,xlabel="z[s]",ylabel="density",
                         hist=TRUE)
         }
  )
  return(g1)
}