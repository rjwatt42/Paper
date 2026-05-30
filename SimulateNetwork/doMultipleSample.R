meanMinMax<-function(a,l) {
  s9<-paste0(l,format(mean(a),digits=3),
             "  {",format(min(a),digits=3),
             ",",format(max(a),digits=3),
             "}")
  
}


doMultipleSample<-function(network,newNetwork,nSamples,nCount,remove,input,gTab) {

  h1<-h1a<-h1b<-h2<-h2a<-h2b<-h3<-h3a<-0
  links<-zeros<-c()
  sigAll<-linkWhenSigAll<-sigWhenLinkAll<-zeroWhenSigAll<-sigWhenZeroAll<-c()
  
  repAll<-linkWhenReplAll<-replWhenLinkAll<-zeroWhenReplAll<-replWhenZeroAll<-c()
  
  id<-showNotification("Starting Multiple",duration=NULL)
  for (ni in 1:nCount) {
    if (newNetwork) {
      network<<-makeNetwork(networkStructure)
    }
    links<-c(links,mean(network$fullLinks!=0))
    zeros<-c(zeros,mean(network$Stheta==0))
    
    samples<-makeNetworkSample(network,input$sampleSize,input$sampleSizeRand,nSamples,remove,h,hist=TRUE)
    sigAll<-c(sigAll,mean(samples$sig))
    
    if (input$repOn) {
      replication<-makeNetworkReplication(samples,input$repPower,input$repPowerUp,h,hist=TRUE)
      
      repAll<-c(repAll,mean(replication$sigrep))
      linkWhenReplAll<-c(linkWhenReplAll,replication$linkWhenRepl)
      replWhenLinkAll<-c(replWhenLinkAll,replication$replWhenLink)
      zeroWhenReplAll<-c(zeroWhenReplAll,replication$zeroWhenRepl)
      replWhenZeroAll<-c(replWhenZeroAll,replication$replWhenZero)
      
      h1<-h1+replication$h1rep$density
      h1a<-h1a+replication$h1rep$density1
      h1b<-h1b+replication$h1rep$density2
      h2<-h2+replication$h2rep$density
      h2a<-h2a+replication$h2rep$density1
      h2b<-h2b+replication$h2rep$density2
      h3<-h3+replication$h3rep$density
      h3a<-h3a+replication$h3rep$density1
      
    } else {
      linkWhenSigAll<-c(linkWhenSigAll,samples$linkWhenSig)
      sigWhenLinkAll<-c(sigWhenLinkAll,samples$sigWhenLink)
      zeroWhenSigAll<-c(zeroWhenSigAll,samples$zeroWhenSig)
      sigWhenZeroAll<-c(sigWhenZeroAll,samples$sigWhenZero)
      
      h1<-h1+samples$h1$density
      h1a<-h1a+samples$h1$density1
      h1b<-NULL
      h2<-h2+samples$h2$density
      h2a<-h2a+samples$h2$density1
      h2b<-NULL
      h3<-h3+samples$h3$density
      h3a<-h3a+samples$h3$density1
      
    }
    if (floor(ni/(nCount/10))*(nCount/10)==ni) showNotification(paste0(ni,"/",nCount),id=id,duration=NULL)
  }
  removeNotification(id)
  
  if (input$repOn) {
    s9<-meanMinMax(links,"p(links)=")
    s10<-meanMinMax(zeros,"p(zeros)=")
    s1<-meanMinMax(sigAll,"p(sig | all)=")
    s2<-meanMinMax(repAll,"p(repl | sig)=")
    s3<-meanMinMax(sigAll*repAll,"p(repl | all)=")
    s4<-meanMinMax(replWhenLinkAll,"p(repl | link)=")
    s5<-meanMinMax(linkWhenReplAll,"p(link | repl)=")
    s7<-meanMinMax(zeroWhenReplAll,"p(zero | repl)=")
    s6<-meanMinMax(replWhenZeroAll,"p(repl | zero)=")
    legend<-c(s9,s10," ",s1,s2," ",s3,s4,s6," ",s5,s7)
  } else {
    s7<-meanMinMax(links,"p(links)=")
    s8<-meanMinMax(zeros,"p(zeros)=")
    s1<-meanMinMax(sigAll,"p(sig | all)=")
    s2<-meanMinMax(linkWhenSigAll,"p(link | sig)=")
    s3<-meanMinMax(sigWhenLinkAll,"p(sig | link)=")
    s4<-meanMinMax(sigWhenZeroAll,"p(sig | zero)=")
    s5<-meanMinMax(zeroWhenSigAll,"p(zero | sig)=")
  legend<-c(s7,s8," ",s1," ",s2,s5," ",s3,s4)
  }
  h0<-list(density=h1,breaks=h$breaks,density1=h1a,density2=h1b)
  g1<-dataGraph(h0,xlabel="z[s]",ylabel="density",legend=legend,
                hist=TRUE)
  h0<-list(density=h2,breaks=h$breaks,density1=h2a,density2=h2b)
  g2<-dataGraph(h0,xlabel="z[p]",ylabel="density",legend=legend,
                hist=TRUE)
  h0<-list(density=h3,breaks=h$breaks,density1=h3a)
  g3<-dataGraph(h0,xlabel="w[p]",ylabel="density",hist=TRUE)
  
  gTab<-generate_tab("Graphs: ",titleWidth=50,
                      tabs=c("Samples","Populations","Power"),
                      tabContents=c(g1,g2,g3),
                      open=1,
                      history=gTab
                      
  ) 
  return(gTab)
}