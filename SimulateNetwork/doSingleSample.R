doSingleSample<-function(network,nSamples,remove,input,gTab) {
  
  samples<-makeNetworkSample(network,input$sampleSize,input$sampleSizeRand,nSamples,remove,h,hist=TRUE)
  
  if (input$repOn) {
    replication<-makeNetworkReplication(samples,input$repPower,input$repPowerUp,h,hist=TRUE)
    
    legend<-c(
      paste0("p(links)=",format(mean(network$fullLinks!=0),digits=3)),
      paste0("p(zero)=",format(mean(network$Stheta==0),digits=3)),
      " ",
      paste0("p(sig | all)=",format(mean(samples$sig),digits=3)),
      paste0("p(repl | sig)=",format(mean(replication$sigrep),digits=3)),
      paste0(" "),
      paste0("p(repl | all)=",format(mean(samples$sig)*mean(replication$sigrep),digits=3)),
      paste0("p(repl | link)=",format(replication$replWhenLink,digits=3)),
      paste0("p(repl | zero)=",format(replication$replWhenZero,digits=3)),
      paste0(" "),
      paste0("p(link | repl)=",format(replication$linkWhenRepl,digits=3)),
      paste0("p(zero | repl)=",format(replication$zeroWhenRepl,digits=3))
    )
    g1<-dataGraph(replication$h1rep,xlabel="z[s]",ylabel="density",legend=legend,
                  hist=TRUE)
    g2<-dataGraph(replication$h2rep,
                  xlabel="z[p]",ylabel="density",legend=legend,
                  hist=TRUE)
    g3<-dataGraph(replication$h3rep,
                  xlabel="w[p]",ylabel="density",hist=TRUE)
    g4<-plotZN(samples,replication)
    
  } else {
    legend<-c(
      paste0("p(links)=",format(mean(network$fullLinks!=0),digits=3)),
      paste0("p(zero)=",format(mean(network$Stheta==0),digits=3)),
      " ",
      paste0("p(sig | all)=",format(mean(samples$sig),digits=3)),
      " ",
      paste0("p(links | sig)=",format(samples$linkWhenSig,digits=3)),
      paste0("p(zero | sig)=",format(samples$zeroWhenSig,digits=3)),
      " ",
      paste0("p(sig | links)=",format(samples$sigWhenLink,digits=3)),
      paste0("p(sig | zero)=",format(samples$sigWhenZero,digits=3))
    )
    g1<-dataGraph(samples$h1,xlabel="z[s]",ylabel="density",legend=legend,hist=TRUE)
    g2<-dataGraph(samples$h2,xlabel="z[p]",ylabel="density",legend=legend,hist=TRUE)
    g3<-dataGraph(samples$h3,xlabel="w[p]",ylabel="density",hist=TRUE)
    g4<-plotZN(samples)
  }
  
  gTab<-generate_tab("Graphs: ",titleWidth=50,
                      tabs=c("Samples","Populations","Power","z-n"),
                      tabContents=c(g1,g2,g3,g4),
                      open=1,
                      history=gTab
                      
  ) 
  return(gTab)
}