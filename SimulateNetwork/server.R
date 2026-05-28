
library(pracma)
source("basicPlot.r")
source('HTMLWidget.R')
source("r2p.R")
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
  
  observeEvent(input$actionC0, 
               {
                 if (input$actionC0==0) return()
                 
                 n<-input$sampleSize
                 nSamples<-0
                 samples<-makeNetworkSample(network,n,nSamples,h,hist=TRUE)
                 replication<-makeNetworkReplication(samples,input$repPower,h,hist=TRUE)
                 
                 legend<-c(paste0("p[sig](original)=",format(mean(samples$sig),digits=3)),
                           paste0("p[sig](repl)=",format(mean(replication$sigrep),digits=3)),
                           paste0("p[sig](overall)=",format(mean(samples$sig)*mean(replication$sigrep),digits=3)))
                 g1<-dataGraph(replication$h1rep,xlabel="z[s]",ylabel="density",legend=legend,
                               hist=TRUE)
                 g2<-dataGraph(replication$h2rep,
                               xlabel="z[p]",ylabel="density",legend=legend,
                               hist=TRUE)
                 g3<-dataGraph(replication$h3rep,
                               xlabel="w[p]",ylabel="density",hist=TRUE)
                 g<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Samples","Populations","Power"),
                                 tabContents=c(g1,g2,g3),
                                 open=1
                 ) 
                 output$mainHTML <- renderUI(HTML(g))
               }
  )
  observeEvent(input$actionC1, 
               {
                 if (input$actionC1==0) return()
                 
                 if (input$sampleSizeV)
                   n<-getSampleSizes(nSamples,dist="Gamma",sN=input$sampleSize,sNRandSD=input$sampleSize*0.66,minN=5)
                 else n<-input$sampleSize
                 samples<-makeNetworkSample(network,n,nSamples,h,hist=TRUE)
                 replication<-makeNetworkReplication(samples,input$repPower,h,hist=TRUE)
                 
                 legend<-c(paste0("p[sig](original)=",format(mean(samples$sig),digits=3)),
                           paste0("p[sig](repl)=",format(mean(replication$sigrep),digits=3)),
                           paste0("p[sig](overall)=",format(mean(samples$sig)*mean(replication$sigrep),digits=3)))
                 g1<-dataGraph(replication$h1rep,xlabel="z[s]",ylabel="density",legend=legend,
                               hist=TRUE)
                 g2<-dataGraph(replication$h2rep,
                               xlabel="z[p]",ylabel="density",legend=legend,
                               hist=TRUE)
                 g3<-dataGraph(replication$h3rep,
                               xlabel="w[p]",ylabel="density",hist=TRUE)
                 g<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Samples","Populations","Power"),
                                 tabContents=c(g1,g2,g3),
                                 open=1
                 ) 
                 output$mainHTML <- renderUI(HTML(g))
               }
  )
  observeEvent(input$actionC2, 
               {
                 if (input$actionC2==0) return()
                 ncount<-1000
                 h1<-h1a<-h1b<-0
                 h2<-h2a<-h2b<-0
                 h3<-h3a<-0
                 sigAll1<-sigAll2<-0
                 id<-showNotification("Starting Multiple",duration=NULL)
                 for (ni in 1:ncount) {
                   if (newNetwork) {
                     network<<-makeNetwork(networkStructure)
                   }
                   if (input$sampleSizeV)
                     n<-getSampleSizes(nSamples,dist="Gamma",sN=input$sampleSize,sNRandSD=input$sampleSize*0.66,minN=5)
                   else n<-input$sampleSize
                   samples<-makeNetworkSample(network,n,nSamples,h,hist=TRUE)
                   replication<-makeNetworkReplication(samples,input$repPower,h,hist=TRUE)
                   
                   sigAll1<-sigAll1+mean(samples$sig)
                   sigAll2<-sigAll2+mean(replication$sigrep)
                   
                   h1<-h1+replication$h1rep$density
                   h1a<-h1a+replication$h1rep$density1
                   h1b<-h1b+replication$h1rep$density2
                   h2<-h2+replication$h2rep$density
                   h2a<-h2a+replication$h2rep$density1
                   h2b<-h2b+replication$h2rep$density2
                   h3<-h3+replication$h3rep$density
                   h3a<-h3a+replication$h3rep$density1
                   
                   if (floor(ni/100)*100==ni) showNotification(paste0(ni,"/",ncount),id=id,duration=NULL)
                 }
                 removeNotification(id)
                 legend<-c(paste0("p[sig](original)=",format(sigAll1/ncount,digits=3)),
                           paste0("p[sig](repl)=",format(sigAll2/ncount,digits=3)),
                           paste0("p[sig](overall)=",format(sigAll1/ncount*sigAll2/ncount,digits=3)))
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
                 n<-input$sampleSize
                 nSamples<-0
                 samples<-makeNetworkSample(network,n,nSamples,h,hist=TRUE)
                 
                 legend<-paste0("p[sig]=",format(mean(samples$sig),digits=3))
                 g1<-dataGraph(samples$h1,xlabel="z[s]",ylabel="density",legend=legend,hist=TRUE)
                 g2<-dataGraph(samples$h2,xlabel="z[p]",ylabel="density",legend=legend,hist=TRUE)
                 g3<-dataGraph(samples$h3,xlabel="w[p]",ylabel="density",hist=TRUE)
                 g<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Samples","Populations","Power"),
                                 tabContents=c(g1,g2,g3),
                                 open=1
                 ) 
                 output$mainHTML <- renderUI(HTML(g))
               })
  observeEvent(input$actionB1,
               {
                 if (input$actionB1==0) return()
                 if (input$sampleSizeV)
                   n<-getSampleSizes(nSamples,dist="Gamma",sN=input$sampleSize,sNRandSD=input$sampleSize*0.66,minN=5)
                 else n<-input$sampleSize
                 samples<-makeNetworkSample(network,n,nSamples,h,hist=TRUE)
                 
                 legend<-paste0("p[sig]=",format(mean(samples$sig),digits=3))
                 g1<-dataGraph(samples$h1,xlabel="z[s]",ylabel="density",legend=legend,hist=TRUE)
                 g2<-dataGraph(samples$h2,xlabel="z[p]",ylabel="density",legend=legend,hist=TRUE)
                 g3<-dataGraph(samples$h3,xlabel="w[p]",ylabel="density",hist=TRUE)
                 g<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Samples","Populations","Power"),
                                 tabContents=c(g1,g2,g3),
                                 open=1
                 ) 
                 output$mainHTML <- renderUI(HTML(g))
               })
  observeEvent(input$actionB2,
               {
                 if (input$actionB2==0) return()
                 
                 h1<-h1a<-h2<-h2a<-h3<-h3a<-0
                 sigAll<-0
                 id<-showNotification("Starting Multiple",duration=NULL)
                 for (ni in 1:ncount) {
                   if (newNetwork) {
                     network<<-makeNetwork(networkStructure)
                   }
                   if (input$sampleSizeV)
                     n<-getSampleSizes(nSamples,dist="Gamma",sN=input$sampleSize,sNRandSD=input$sampleSize*0.66,minN=5)
                   else n<-input$sampleSize
                   samples<-makeNetworkSample(network,n,nSamples,h,hist=TRUE)
                   
                   sigAll<-sigAll+mean(samples$sig)
                   h1<-h1+samples$h1$density
                   h1a<-h1a+samples$h1$density1
                   h2<-h2+samples$h2$density
                   h2a<-h2a+samples$h2$density1
                   h3<-h3+samples$h3$density
                   h3a<-h3a+samples$h3$density1
                   if (floor(ni/100)*100==ni) showNotification(paste0(ni,"/",ncount),id=id,duration=NULL)
                 }
                 removeNotification(id)
                 
                 legend<-paste0("p[sig]=",format(mean(sigAll/ncount),digits=3))
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
                 ncount=1000
                 ests<-zeros<-c()
                 indirect<-direct<-c()
                 id<-showNotification("Starting Multiple")
                 for (ni in 1:ncount) {
                   
                   # make a network
                   network<<-makeNetwork(networkStructure)
                   # get effect sizes & distribution
                   h1<-getNetworkHist(network,input$separateZeros,h)
                   
                   # save results
                   density<-density+h1$hist$density/ncount
                   ests<-c(ests,h1$est$estimate)
                   zeros<-c(zeros,h1$zeros)
                   if (floor(ni/100)*100==ni) showNotification(paste0(ni,"/",ncount),id=id)
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