
library(pracma)
source("r2p.R")
source("samplePower.R")
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
Stheta<-c()
h<-c()

server <- function(input, output) {
  initGraph("HTML",gsize=450,autoShow=FALSE)
  
  openTab<<-1

  observeEvent({c(input$repPower,input$actionRep,input$actionC)}, 
               {
                 if (isempty(Stheta)) return()
                 
                 nSamples<-100000
                 use<-Stheta<1
                 zp<-atanh(Stheta[use])
                 
                 use<-ceiling(runif(nSamples,0,1)*length(zp))
                 zp_use<-zp[use]*sign(runif(nSamples,-1,1))
                 
                 switch(input$actionRep,
                        "Single"={
                          n<-input$sampleSize+runif(nSamples,input$sampleSize*0.5,input$sampleSize*10)
                          esd<-1/sqrt(n-3)
                          err<-rnorm(nSamples,0,esd)
                          zs<-abs(zp_use+err)
                          
                          ps<-r2p(tanh(zs),n)
                          sig<-ps<0.05
                          zss<-zs[sig]
                          zp_rep<-zp_use[sig]
                          
                          nr<-rw2n(tanh(zss),input$repPower)
                          esd<-1/sqrt(nr-3)
                          err<-rnorm(length(nr),0,esd)
                          zsr<-abs(zp_rep+err)
                          ps<-r2p(tanh(zsr),nr)
                          sig<-ps<0.05
                          zsrs<-zsr[sig]
                          zp_rep1<-zp_rep[sig]
                          
                          h1a<-hist(zsr[zsr<max(h$breaks)],h$breaks,plot=FALSE)
                          h2a<-hist(zsrs[zsrs<max(h$breaks)],h$breaks,plot=FALSE)
                          h1<-h1a$density
                          h2<-h2a$density*sum(h2a$counts)/sum(h1a$counts)

                          h3a<-hist(abs(zp_use[abs(zp_use)<max(h$breaks)]),h$breaks,plot=FALSE)
                          h4a<-hist(abs(zp_rep1[abs(zp_rep1)<max(h$breaks)]),h$breaks,plot=FALSE)
                          h3<-h3a$density
                          h4<-h4a$density*sum(h4a$counts)/sum(h3a$counts)
                          
                          h0<-list(hist=list(density=h1,breaks=h$breaks,density1=h2))
                          g1<-plotNetworkHist(h0,"samp")
                          
                          g2<-dataGraph(data.frame(x=h$breaks[2:101],y=h4a$counts/h3a$counts),
                                        xlabel="z[p]",ylabel="p(sig|rep)")
                          
                        },
                        "Multiple"={
                          ncount<-1000
                          h1<-h2<-h3<-h4<-0
                          id<-showNotification("Starting Multiple")
                          for (ni in 1:ncount) {
                            n<-input$sampleSize+runif(nSamples,input$sampleSize*0.5,input$sampleSize*10)
                            esd<-1/sqrt(n-3)
                            err<-rnorm(nSamples,0,esd)
                            zs<-abs(zp_use+err)
                            
                            ps<-r2p(tanh(zs),n)
                            sig<-ps<0.05
                            zss<-zs[sig]
                            zp_rep<-zp_use[sig]
                            
                            nr<-rw2n(tanh(zss),input$repPower)
                            esd<-1/sqrt(nr-3)
                            err<-rnorm(length(nr),0,esd)
                            zsr<-abs(zp_rep+err)
                            ps<-r2p(tanh(zsr),nr)
                            sig<-ps<0.05
                            zsrs<-zsr[sig]
                            zp_rep1<-zp_rep[sig]
                            
                            h1a<-hist(zsr[zsr<max(h$breaks)],h$breaks,plot=FALSE)
                            h2a<-hist(zsrs[zsrs<max(h$breaks)],h$breaks,plot=FALSE)
                            h1<-h1+h1a$density
                            h2<-h2+h2a$density*sum(h2a$counts)/sum(h1a$counts)
                            
                            h3a<-hist(abs(zp_use[abs(zp_use)<max(h$breaks)]),h$breaks,plot=FALSE)
                            h4a<-hist(abs(zp_rep1[abs(zp_rep1)<max(h$breaks)]),h$breaks,plot=FALSE)
                            h3<-h3+h3a$density
                            h4<-h4+h4a$density*sum(h4a$counts)/sum(h3a$counts)
                            
                            if (floor(ni/100)*100==ni) showNotification(paste0(ni,"/",ncount),id=id)
                          }
                          removeNotification(id)
                          h0<-list(hist=list(density=h1,breaks=h$breaks,density1=h2))
                          g1<-plotNetworkHist(h0,"samp")
                          
                          g2<-dataGraph(data.frame(x=h$breaks[2:101],y=h4/h3),
                                        xlabel="z[p]",ylabel="p(sig|rep)")
                        }
                 )
                 g<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Samples","Populations"),
                                 tabContents=c(g1,g2),
                                 open=1
                 ) 
                 output$mainHTML <- renderUI(HTML(g))
               }
  )
  observeEvent({c(input$sampleSize,input$actionSamp,input$actionB)}, 
               {
                 if (isempty(Stheta)) return()
                 
                 nSamples<-100000
                 use<-Stheta<1
                 zp<-atanh(Stheta[use])
                 
                 use<-ceiling(runif(nSamples,0,1)*length(zp))
                 zp_use<-zp[use]*sign(runif(nSamples,-1,1))
                 
                 switch(input$actionSamp,
                        "Single"={
                          n<-input$sampleSize+runif(nSamples,input$sampleSize*0.5,input$sampleSize*10)
                          esd<-1/sqrt(n-3)
                          err<-rnorm(nSamples,0,esd)
                          zs<-abs(zp_use+err)
                          
                          ps<-r2p(tanh(zs),n)
                          sig<-ps<0.05
                          zss<-zs[sig]
                          zps_use<-zp_use[sig]
                          
                          h1a<-hist(zs[zs<max(h$breaks)],h$breaks,plot=FALSE)
                          h2a<-hist(zss[zss<max(h$breaks)],h$breaks,plot=FALSE)
                          h1<-h1a$density
                          h2<-h2a$density*sum(h2a$counts)/sum(h1a$counts)
                          h0<-list(hist=list(density=h1,breaks=h$breaks,density1=h2))
                          g1<-plotNetworkHist(h0,"samp")
                          
                          h3a<-hist(abs(zp_use[abs(zp_use)<max(h$breaks)]),h$breaks,plot=FALSE)
                          h4a<-hist(abs(zps_use[abs(zps_use)<max(h$breaks)]),h$breaks,plot=FALSE)
                          g2<-dataGraph(data.frame(x=h$breaks[2:101],y=h4a$counts/h3a$counts),
                                        xlabel="z[p]",ylabel="p(sig)")
                          
                        },
                        "Multiple"={
                          ncount<-1000
                          h1<-h2<-h3<-h4<-0
                          id<-showNotification("Starting Multiple")
                          for (ni in 1:ncount) {
                          n<-input$sampleSize+runif(nSamples,input$sampleSize*0.5,input$sampleSize*10)
                          esd<-1/sqrt(n-3)
                          err<-rnorm(nSamples,0,esd)
                          zs<-abs(zp_use+err)
                          
                          ps<-r2p(tanh(zs),n)
                          sig<-ps<0.05
                          zss<-zs[sig]
                          
                          h1a<-hist(zs[zs<max(h$breaks)],h$breaks,plot=FALSE)
                          h2a<-hist(zss[zss<max(h$breaks)],h$breaks,plot=FALSE)
                          h1<-h1+h1a$density
                          h2<-h2+h2a$density*sum(h2a$counts)/sum(h1a$counts)
                          
                          h3a<-hist(abs(zp_use[abs(zp_use)<max(h$breaks)]),h$breaks,plot=FALSE)
                          h4a<-hist(abs(zps_use[abs(zps_use)<max(h$breaks)]),h$breaks,plot=FALSE)
                          h3<-h3+h3a$density
                          h4<-h4+h4a$density*sum(h4a$counts)/sum(h3a$counts)
                          
                          if (floor(ni/100)*100==ni) showNotification(paste0(ni,"/",ncount),id=id)
                        }
                        removeNotification(id)
                        h0<-list(hist=list(density=h1,breaks=h$breaks,density1=h2))
                        g1<-plotNetworkHist(h0,"samp")
                        
                        g2<-dataGraph(data.frame(x=h$breaks[2:101],y=h4/h3),
                                      xlabel="z[p]",ylabel="p(sig)")
                        
                        }
                 )
                 g<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Samples","Populations"),
                                 tabContents=c(g1,g2),
                                 open=1
                 ) 
                 output$mainHTML <- renderUI(HTML(g))
               }
  )
  observeEvent({c(input$nStages,input$nVarsPerStage,
                  input$linkRange,input$probLink,input$strengthLink,
                  input$separateZeros,input$actionPop,input$actionA)}, 
               {
                 
                 nStages=input$nStages
                 nVarsPerStage=input$nVarsPerStage
                 rangeLink=input$linkRange
                 probLink=input$probLink
                 strengthLink=input$strengthLink
                 separateZeros<-input$separateZeros
                 fullModel=TRUE
                 
                 h<<-list(breaks=seq(0,(strengthLink*1.5),length.out=101),density=0)
                 
                 switch(input$actionPop,
                        "Single"={
                          # make a network
                          fullLinks<-makeNetwork(nStages,nVarsPerStage,rangeLink,probLink)
                          network<-links2Path(fullLinks,nStages,nVarsPerStage)
                          g1<-plotNetwork(network,showNames=FALSE)
                          
                          # histogram of effect sizes
                          h1<-getNetworkHist(fullLinks,fullModel,strengthLink,separateZeros,h)
                          Stheta<<-h1$Stheta
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