
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
source("doSingleSample.R")
source("doMultipleSample.R")

# data<-NULL
gs<-c(" "," "," ")
gTab<-""

h<-c()
network<-c()

server <- function(input, output) {
  initGraph("HTML",gsize=450,autoShow=FALSE,fontScale = 1.1)
  
  newNetwork<-TRUE
  nSamples<-100000
  nCount<-100
  # remove<-list(row=TRUE,gap=1)
  remove<-NULL
  
  observeEvent(input$actionB0,
               {
                 if (input$actionB0==0) return()
                 gTab<<-doSingleSample(network,nSamples=0,remove,input,gTab)
                 
                 output$mainHTML <- renderUI(HTML(gTab))
               })
  observeEvent(input$actionB1,
               {
                 if (input$actionB1==0) return()
                 gTab<<-doSingleSample(network,nSamples,remove,input,gTab)
                 
                 output$mainHTML <- renderUI(HTML(gTab))
               })
  observeEvent(input$actionB2,
               {
                 if (input$actionB2==0) return()
                 gTab<<-doMultipleSample(network,newNetwork,nSamples,nCount,remove,input,gTab)
                 
                 output$mainHTML <- renderUI(HTML(gTab))
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
                 
                 gTab<<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Network","Effects"),
                                 tabContents=c(gs[1],gs[2]),
                                 open=1,
                                 history=gTab
                 ) 
                 output$mainHTML <- renderUI(HTML(gTab))
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
                 
                 gTab<<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Network","Effects"),
                                 tabContents=c(gs[1],gs[2]),
                                 open=2,
                                 history=gTab
                 ) 
                 output$mainHTML <- renderUI(HTML(gTab))
               })
  
}
