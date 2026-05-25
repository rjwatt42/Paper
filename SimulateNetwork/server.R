
library(pracma)
source("basicPlot.r")
source('HTMLWidget.R')
source("doNetwork.R")
source("makeNetwork.R")
source("plotNetwork.R")
source("getNetworkHist.R")
source("plotNetworkHist.R")

# data<-NULL
openTab<-1
startTab<-0
gs<-c(" "," "," ")
zpsMain<-c()

server <- function(input, output) {
  initGraph("HTML",gsize=450,autoShow=FALSE)
  
  openTab<<-1
  
  observeEvent({c(input$nStages,input$nVarsPerStage,
                  input$linkRange,input$probLink,input$strengthLink,
                  input$separateZeros,input$action,input$actionB)}, 
               {
                 
                 nStages=input$nStages
                 nVarsPerStage=input$nVarsPerStage
                 rangeLink=input$linkRange
                 probLink=input$probLink
                 strengthLink=input$strengthLink
                 separateZeros<-input$separateZeros
                 fullModel=TRUE
                 
                 h<-list(breaks=seq(0,(strengthLink*1.5),length.out=101),density=0)
                 
                 switch(input$action,
                        "Single"={
                          # make a network
                          fullLinks<-makeNetwork(nStages,nVarsPerStage,rangeLink,probLink)
                          g1<-plotNetwork(links2Path(fullLinks,nStages,nVarsPerStage),showNames=FALSE)
                          
                          # histogram of effect sizes
                          h1<-getNetworkHist(fullLinks,fullModel,strengthLink,separateZeros,h)
                          g2<-plotNetworkHist(h1)
                          openTab<<-1
                          gs[1]<<-g1
                          gs[2]<<-g2
                        },
                        "Multiple"={
                          density<-0
                          ncount=1000
                          ests<-zeros<-c()
                          indirect<-direct<-c()
                          id<-showNotification("Starting Multiple")
                          for (ni in 1:ncount) {
                            
                            # make a network
                            fullLinks<-makeNetwork(nStages,nVarsPerStage,rangeLink,probLink)
                            # get effect sizes & distribution
                            h1<-getNetworkHist(fullLinks,fullModel,strengthLink,separateZeros,h)
                            
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
                          gs[3]<<-g3
                          openTab<<-3
                        }
                 )
                 g<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Network","Histogram","Multiple"),
                                 tabContents=c(gs[1],gs[2],gs[3]),
                                 open=openTab
                 ) 
                 output$mainHTML <- renderUI(HTML(g))
               })
  
}