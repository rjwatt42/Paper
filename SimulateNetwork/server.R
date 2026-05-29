
library(pracma)
source("basicPlot.r")
source('HTMLWidget.R')
source('plotZN.R')
source("samplePower.R")
source("getSampleSizes.R")
source("doNetwork.R")
source("makeNetwork.R")
source("plotNetwork.R")
source("makeNetworkSample.R")
source("makeNetworkReplication.R")
source("getNetworkHist.R")
source("plotNetworkHist.R")

# data<-NULL
gs<-c(" "," "," ")
Stheta<-c()
h<-c()
network<-c()

server <- function(input, output) {
  initGraph("HTML",gsize=450,autoShow=FALSE,fontScale = 1)
  
  newNetwork<-TRUE
  nSamples<-100000
  nCount<-100
  # remove<-list(row=TRUE,gap=1)
  remove<-NULL
  
  observeEvent(input$actionC0, 
               {
                 if (input$actionC0==0) return()
                 
                 nSamples<-0
                 samples<-makeNetworkSample(network,input$sampleSize,input$sampleSizeRand,nSamples,remove,h,hist=TRUE)
                 replication<-makeNetworkReplication(samples,input$repPower,input$repPowerUp,h,hist=TRUE)
                 
                 legend<-c(#paste0("p(sig | all)=",format(mean(samples$sig),digits=3)),
                   paste0("p(repl | all)=",format(mean(samples$sig)*mean(replication$sigrep),digits=3)),
                   paste0("p(repl | sig)=",format(mean(replication$sigrep),digits=3)),
                   paste0("p(repl | link)=",format(replication$replWhenLink,digits=3)),
                   paste0("p(link | repl)=",format(replication$linkWhenRepl,digits=3)),
                   paste0("p(repl | zero)=",format(replication$replWhenZero,digits=3)),
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
                 g<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Samples","Populations","Power","z-n"),
                                 tabContents=c(g1,g2,g3,g4),
                                 open=1
                 ) 
                 output$mainHTML <- renderUI(HTML(g))
               }
  )
  observeEvent(input$actionC1, 
               {
                 if (input$actionC1==0) return()
                 
                 samples<-makeNetworkSample(network,input$sampleSize,input$sampleSizeRand,nSamples,remove,h,hist=TRUE)
                 replication<-makeNetworkReplication(samples,input$repPower,input$repPowerUp,h,hist=TRUE)
                 
                 legend<-c(#paste0("p(sig | all)=",format(mean(samples$sig),digits=3)),
                   paste0("p(repl | all)=",format(mean(samples$sig)*mean(replication$sigrep),digits=3)),
                   paste0("p(repl | sig)=",format(mean(replication$sigrep),digits=3)),
                   paste0("p(repl | link)=",format(replication$replWhenLink,digits=3)),
                   paste0("p(link | repl)=",format(replication$linkWhenRepl,digits=3)),
                   paste0("p(repl | zero)=",format(replication$replWhenZero,digits=3)),
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
                 g<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Samples","Populations","Power","zs-n"),
                                 tabContents=c(g1,g2,g3,g4),
                                 open=1
                 ) 
                 output$mainHTML <- renderUI(HTML(g))
               }
  )
  observeEvent(input$actionC2, 
               {
                 if (input$actionC2==0) return()
                 h1<-h1a<-h1b<-0
                 h2<-h2a<-h2b<-0
                 h3<-h3a<-0
                 sigAll1<-sigAll2<-c()
                 linkWhenReplAll<-replWhenLinkAll<-zeroWhenReplAll<-replWhenZeroAll<-c()
                 id<-showNotification("Starting Multiple",duration=NULL)
                 for (ni in 1:nCount) {
                   if (newNetwork) {
                     network<<-makeNetwork(networkStructure)
                   }
                   samples<-makeNetworkSample(network,input$sampleSize,input$sampleSizeRand,nSamples,remove,h,hist=TRUE)
                   replication<-makeNetworkReplication(samples,input$repPower,input$repPowerUp,h,hist=TRUE)
                   
                   sigAll1<-c(sigAll1,mean(samples$sig))
                   sigAll2<-c(sigAll2,mean(replication$sigrep))
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
                   
                   if (floor(ni/(nCount/10))*(nCount/10)==ni) showNotification(paste0(ni,"/",nCount),id=id,duration=NULL)
                 }
                 removeNotification(id)
                 
                 s1<-paste0("p(sig | all)=",format(mean(sigAll1),digits=3),
                            "  {",format(min(sigAll1),digits=3),
                            ",",format(max(sigAll1),digits=3),
                            "}")
                 s2<-paste0("p(repl | sig)=",format(mean(sigAll2),digits=3),
                            "  {",format(min(sigAll2),digits=3),
                            ",",format(max(sigAll2),digits=3),
                            "}")
                 s3<-paste0("p(repl | all)=",format(mean(sigAll1*sigAll2),digits=3),
                            "  {",format(min(sigAll1*sigAll2),digits=3),
                            ",",format(max(sigAll1*sigAll2),digits=3),
                            "}")
                 s4<-paste0("p(repl | link)=",format(mean(replWhenLinkAll),digits=3),
                            "  {",format(min(replWhenLinkAll),digits=3),
                            ",",format(max(replWhenLinkAll),digits=3),
                            "}")
                 s5<-paste0("p(link | repl)=",format(mean(linkWhenReplAll),digits=3),
                            "  {",format(min(linkWhenReplAll),digits=3),
                            ",",format(max(linkWhenReplAll),digits=3),
                            "}")
                 s6<-paste0("p(repl | zero)=",format(mean(replWhenZeroAll),digits=3),
                            "  {",format(min(replWhenZeroAll),digits=3),
                            ",",format(max(replWhenZeroAll),digits=3),
                            "}")
                 s7<-paste0("p(zero | repl)=",format(mean(zeroWhenReplAll),digits=3),
                            "  {",format(min(zeroWhenReplAll),digits=3),
                            ",",format(max(zeroWhenReplAll),digits=3),
                            "}")
                 
                 legend<-c(#s1,
                           s2,s3,s4,s5,s6,s7)
                 h0<-list(density=h1,breaks=h$breaks,density1=h1a,density2=h1b)
                 g1<-dataGraph(h0,xlabel="z[s]",ylabel="density",legend=legend,
                               hist=TRUE)
                 h0<-list(density=h2,breaks=h$breaks,density1=h2a,density2=h2b)
                 g2<-dataGraph(h0,xlabel="z[p]",ylabel="density",legend=legend,
                               hist=TRUE)
                 h0<-list(density=h3,breaks=h$breaks,density1=h3a)
                 g3<-dataGraph(h0,xlabel="w[p]",ylabel="density",hist=TRUE)
                 
                 g<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Samples","Populations","Power"),
                                 tabContents=c(g1,g2,g3),
                                 open=1
                 )
                 output$mainHTML <- renderUI(HTML(g))
               }
  )
  
  observeEvent(input$actionB0,
               {
                 if (input$actionB0==0) return()
                 nSamples<-0
                 samples<-makeNetworkSample(network,input$sampleSize,input$sampleSizeRand,nSamples,remove,h,hist=TRUE)

                 legend<-c(paste0("p(sig | all)=",format(mean(samples$sig),digits=3)),
                           paste0("p(links | sig)=",format(samples$linkWhenSig,digits=3)),
                           paste0("p(sig | links)=",format(samples$sigWhenLink,digits=3)),
                           paste0("p(sig | zero)=",format(samples$sigWhenZero,digits=3)),
                           paste0("p(zero | sig)=",format(samples$zeroWhenSig,digits=3))
                 )
                 g1<-dataGraph(samples$h1,xlabel="z[s]",ylabel="density",legend=legend,hist=TRUE)
                 g2<-dataGraph(samples$h2,xlabel="z[p]",ylabel="density",legend=legend,hist=TRUE)
                 g3<-dataGraph(samples$h3,xlabel="w[p]",ylabel="density",hist=TRUE)
                 g4<-plotZN(samples)
                 g<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Samples","Populations","Power","z-n"),
                                 tabContents=c(g1,g2,g3,g4),
                                 open=1
                 ) 
                 output$mainHTML <- renderUI(HTML(g))
               })
  observeEvent(input$actionB1,
               {
                 if (input$actionB1==0) return()
                 samples<-makeNetworkSample(network,input$sampleSize,input$sampleSizeRand,nSamples,remove,h,hist=TRUE)
                 
                 legend<-c(paste0("p(sig | all)=",format(mean(samples$sig),digits=3)),
                           paste0("p(links | sig)=",format(samples$linkWhenSig,digits=3)),
                           paste0("p(sig | links)=",format(samples$sigWhenLink,digits=3)),
                           paste0("p(sig | zero)=",format(samples$sigWhenZero,digits=3)),
                           paste0("p(zero | sig)=",format(samples$zeroWhenSig,digits=3))
                 )
                 g1<-dataGraph(samples$h1,xlabel="z[s]",ylabel="density",legend=legend,hist=TRUE)
                 g2<-dataGraph(samples$h2,xlabel="z[p]",ylabel="density",legend=legend,hist=TRUE)
                 g3<-dataGraph(samples$h3,xlabel="w[p]",ylabel="density",hist=TRUE)
                 g4<-plotZN(samples)
                 g<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Samples","Populations","Power","z-n"),
                                 tabContents=c(g1,g2,g3,g4),
                                 open=1
                 ) 
                 output$mainHTML <- renderUI(HTML(g))
               })
  observeEvent(input$actionB2,
               {
                 if (input$actionB2==0) return()
                 
                 h1<-h1a<-h2<-h2a<-h3<-h3a<-0
                 sigAll<-linkWhenSigAll<-sigWhenLinkAll<-zeroWhenSigAll<-sigWhenZeroAll<-c()
                 id<-showNotification("Starting Multiple",duration=NULL)
                 for (ni in 1:nCount) {
                   if (newNetwork) {
                     network<<-makeNetwork(networkStructure)
                   }
                   samples<-makeNetworkSample(network,input$sampleSize,input$sampleSizeRand,nSamples,remove,h,hist=TRUE)

                   sigAll<-c(sigAll,mean(samples$sig))
                   linkWhenSigAll<-c(linkWhenSigAll,samples$linkWhenSig)
                   sigWhenLinkAll<-c(sigWhenLinkAll,samples$sigWhenLink)
                   zeroWhenSigAll<-c(zeroWhenSigAll,samples$zeroWhenSig)
                   sigWhenZeroAll<-c(sigWhenZeroAll,samples$sigWhenZero)
                   
                   h1<-h1+samples$h1$density
                   h1a<-h1a+samples$h1$density1
                   h2<-h2+samples$h2$density
                   h2a<-h2a+samples$h2$density1
                   h3<-h3+samples$h3$density
                   h3a<-h3a+samples$h3$density1
                   if (floor(ni/(nCount/10))*(nCount/10)==ni) showNotification(paste0(ni,"/",nCount),id=id,duration=NULL)
                 }
                 removeNotification(id)
                 
                 s1<-paste0("p(sig | all)=",format(mean(sigAll),digits=3),
                            "  {",format(min(sigAll),digits=3),
                            ",",format(max(sigAll),digits=3),
                            "}")
                 s2<-paste0("p(link | sig)=",format(mean(linkWhenSigAll),digits=3),
                            "  {",format(min(linkWhenSigAll),digits=3),
                            ",",format(max(linkWhenSigAll),digits=3),
                            "}")
                 s3<-paste0("p(sig | link)=",format(mean(sigWhenLinkAll),digits=3),
                            "  {",format(min(sigWhenLinkAll),digits=3),
                            ",",format(max(sigWhenLinkAll),digits=3),
                            "}")
                 s4<-paste0("p(sig | zero)=",format(mean(sigWhenZeroAll),digits=3),
                            "  {",format(min(sigWhenZeroAll),digits=3),
                            ",",format(max(sigWhenZeroAll),digits=3),
                            "}")
                 s5<-paste0("p(zero | sig)=",format(mean(zeroWhenSigAll),digits=3),
                            "  {",format(min(zeroWhenSigAll),digits=3),
                            ",",format(max(zeroWhenSigAll),digits=3),
                            "}")
                 legend<-c(s1,s2,s3,s4,s5)
                 g1<-dataGraph(list(density=h1,breaks=h$breaks,density1=h1a),
                               xlabel="z[s]",ylabel="density",legend=legend,
                               hist=TRUE)
                 g2<-dataGraph(list(breaks=h$breaks,density=h2,density1=h2a),
                               xlabel="z[p]",ylabel="density",legend=legend,
                               hist=TRUE)
                 g3<-dataGraph(list(breaks=seq(0,1,length.out=101),density=h3,density1=h3a),
                               xlabel="w[p]",ylabel="density",hist=TRUE)
                 g<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Samples","Populations","Power"),
                                 tabContents=c(g1,g2,g3),
                                 open=1
                 ) 
                 output$mainHTML <- renderUI(HTML(g))
               })
  
  observeEvent({c(input$nStages,input$nNodesPerStage,
                  input$linkRange,input$probLink,input$strengthLink,
                  input$separateZeros,input$actionA1)}, 
               {
                 
                 networkStructure<<-list(
                   nStages=input$nStages,
                   nNodesPerStage=input$nNodesPerStage,
                   rangeLink=input$linkRange,
                   probLink=input$probLink,
                   strengthLink=input$strengthLink,
                   fullModel=TRUE
                 )
                 
                 h<<-list(breaks=seq(0,(networkStructure$strengthLink*1.5),length.out=101),density=0)
                 
                 # make a network
                 network<<-makeNetwork(networkStructure)
                 g1<-plotNetwork(network,showNames=FALSE)
                 
                 # histogram of effect sizes
                 h1<-getNetworkHist(network,input$separateZeros,h)
                 g2<-plotNetworkHist(h1)
                 gs[1]<<-g1
                 gs[2]<<-g2
                 
                 g<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Network","Effects"),
                                 tabContents=c(gs[1],gs[2]),
                                 open=1
                 ) 
                 output$mainHTML <- renderUI(HTML(g))
               })
  
  observeEvent(input$actionA2, 
               {
                 density<-0
                 ests<-zeros<-c()
                 indirect<-direct<-c()
                 id<-showNotification("Starting Multiple")
                 for (ni in 1:nCount) {
                   
                   # make a network
                   network<<-makeNetwork(networkStructure)
                   # get effect sizes & distribution
                   h1<-getNetworkHist(network,input$separateZeros,h)
                   
                   # save results
                   density<-density+h1$hist$density
                   ests<-c(ests,h1$est$estimate)
                   zeros<-c(zeros,h1$zeros)
                   if (floor(ni/(nCount/10))*(nCount/10)==ni) showNotification(paste0(ni,"/",nCount),id=id)
                 }
                 removeNotification(id)
                 
                 h1$hist$density<-density
                 h1$est$estimate<-mean(ests)
                 h1$zeros<-mean(zeros)
                 g3<-plotNetworkHist(h1)
                 gs[2]<<-g3
                 
                 g<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Network","Effects"),
                                 tabContents=c(gs[1],gs[2]),
                                 open=2
                 ) 
                 output$mainHTML <- renderUI(HTML(g))
               })
  
}
