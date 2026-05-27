
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
  
  observeEvent(input$actionC1, 
               {
                 if (input$actionC1==0) return()
                 
                 n<-getSampleSizes(nSamples,dist="Gamma",sN=input$sampleSize,sNRandSD=33,minN=5)
                 samples<-makeNetworkSample(network,n,h,hist=TRUE,doReplication=TRUE,input$repPower)
                 
                 legend<-c(paste0("p[sig](original)=",format(mean(samples$sig),digits=3)),
                           paste0("p[sig](repl)=",format(mean(samples$sigr),digits=3)),
                           paste0("p[sig](overall)=",format(mean(samples$sig)*mean(samples$sigr),digits=3)))
                 g1<-dataGraph(samples$h1r,xlabel="z[s]",ylabel="density",legend=legend,
                               hist=TRUE)
                 g2<-dataGraph(samples$h2r,
                               xlabel="z[p]",ylabel="density",legend=legend,
                               hist=TRUE)
                 g3<-dataGraph(samples$h3r,
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
                 h1<-h2<-h3<-h1p<-h2p<-h3p<-hw1<-hw2<-hw3<-0
                 sigAll1<-sigAll2<-0
                 id<-showNotification("Starting Multiple",duration=NULL)
                 for (ni in 1:ncount) {
                   if (newNetwork) {
                     # make a network
                     fullLinks<-makeNetwork(networkStructure)
                     # get effect sizes & distribution
                     hnet<-getNetworkHist(fullLinks,networkStructure,h)
                     use<-hnet$Stheta<1
                     zp<-atanh(hnet$Stheta[use])
                   }
                   n<-getSampleSizes(nSamples,dist="Gamma",sN=input$sampleSize,sNRandSD=33,minN=5)
                   esd<-1/sqrt(n-3)
                   err<-rnorm(nSamples,0,esd)
                   
                   use<-ceiling(runif(nSamples,0,1)*length(zp))
                   zp_use<-zp[use]*sign(runif(nSamples,-1,1))
                   zs<-abs(zp_use+err)
                   
                   ps<-r2p(tanh(zs),n)
                   sig<-ps<0.05
                   zss<-zs[sig]
                   zp_rep<-zp_use[sig]
                   sigAll1<-sigAll1+mean(sig)
                   
                   nr<-rw2n(tanh(zss),input$repPower)
                   esd<-1/sqrt(nr-3)
                   err<-rnorm(length(nr),0,esd)
                   zsr<-abs(zp_rep+err)
                   ps<-r2p(tanh(zsr),nr)
                   sigr<-ps<0.05
                   zsrs<-zsr[sigr]
                   zp_rep1<-zp_rep[sigr]
                   wr1<-rn2w(tanh(zp_use),n)
                   wr2<-rn2w(tanh(zp_rep),nr)
                   wr3<-rn2w(tanh(zp_rep[sigr]),nr[sigr])
                   sigAll2<-sigAll2+mean(sigr)
                   
                   h1a<-hist(zs[zs<max(h$breaks)],h$breaks,plot=FALSE)
                   h2a<-hist(zsr[zsr<max(h$breaks)],h$breaks,plot=FALSE)
                   h3a<-hist(zsrs[zsrs<max(h$breaks)],h$breaks,plot=FALSE)
                   h1<-h1+h1a$density
                   h2<-h2+h2a$density*sum(h2a$counts)/sum(h1a$counts)
                   h3<-h3+h3a$density*sum(h3a$counts)/sum(h1a$counts)
                   
                   h1pa<-hist(abs(zp_use[abs(zp_use)<max(h$breaks)]),h$breaks,plot=FALSE)
                   h2pa<-hist(abs(zp_rep[abs(zp_rep)<max(h$breaks)]),h$breaks,plot=FALSE)
                   h3pa<-hist(abs(zp_rep1[abs(zp_rep1)<max(h$breaks)]),h$breaks,plot=FALSE)
                   h1p<-h1p+h1pa$density
                   h2p<-h2p+h2pa$density*sum(h2pa$counts)/sum(h1pa$counts)
                   h3p<-h3p+h3pa$density*sum(h3pa$counts)/sum(h1pa$counts)
                   
                   hw1<-hw1+hist(wr1,seq(0,1,length.out=101),plot=FALSE)$density
                   hw2<-hw2+hist(wr2,seq(0,1,length.out=101),plot=FALSE)$density
                   hw3<-hw3+hist(wr3,seq(0,1,length.out=101),plot=FALSE)$density
                   
                   if (floor(ni/100)*100==ni) showNotification(paste0(ni,"/",ncount),id=id,duration=NULL)
                 }
                 removeNotification(id)
                 legend<-c(paste0("p[sig](original)=",format(mean(sig),digits=3)),
                           paste0("p[sig](repl)=",format(mean(sigr),digits=3)),
                           paste0("p[sig](overall)=",format(mean(sig)*mean(sigr),digits=3)))
                 h0<-list(density=h1,breaks=h$breaks,density1=h2,density2=h3)
                 g1<-dataGraph(h0,xlabel="z[s]",ylabel="density",legend=legend,
                               hist=TRUE)
                 
                 g2<-dataGraph(list(breaks=h$breaks,
                                    density=h1p,
                                    density1=h2p,
                                    density2=h3p),
                               xlabel="z[p]",ylabel="density",legend=legend,
                               hist=TRUE)
                 
                 g3<-dataGraph(list(density=hw1,density1=hw2,density2=hw3,breaks=seq(0,1,length.out=101)),
                               xlabel="w[p]",ylabel="density",hist=TRUE)
                 
                 g<-generate_tab("Graphs: ",titleWidth=50,
                                 tabs=c("Samples","Populations","Power"),
                                 tabContents=c(g1,g2,g3),
                                 open=1
                 )
                 output$mainHTML <- renderUI(HTML(g))
               }
  )
  
  observeEvent(input$actionB1,
               {
    if (input$actionB1==0) return()
    n<-getSampleSizes(nSamples,dist="Gamma",sN=input$sampleSize,sNRandSD=33,minN=5)
    samples<-makeNetworkSample(network,n,h,hist=TRUE)
    
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
                 
    h1<-h2<-h3<-h4<-hw1<-hw2<-0
    sigAll<-0
    id<-showNotification("Starting Multiple",duration=NULL)
    for (ni in 1:ncount) {
      if (newNetwork) {
        network<<-makeNetwork(networkStructure)
      }
      n<-getSampleSizes(nSamples,dist="Gamma",sN=input$sampleSize,sNRandSD=33,minN=5)
      samples<-makeNetworkSample(network,n,h,hist=TRUE)
      
      sigAll<-sigAll+mean(samples$sig)
      h1<-h1+samples$h1$density
      h2<-h2+samples$h1$density1
      h3<-h3+samples$h2$density
      h4<-h4+samples$h2$density1
      hw1<-hw1+samples$h3$density
      hw2<-hw2+samples$h3$density1
      if (floor(ni/100)*100==ni) showNotification(paste0(ni,"/",ncount),id=id,duration=NULL)
    }
    removeNotification(id)
    
    legend<-paste0("p[sig]=",format(mean(sigAll/ncount),digits=3))
    g1<-dataGraph(list(density=h1,breaks=h$breaks,density1=h2),
                  xlabel="z[s]",ylabel="density",legend=legend,
                  hist=TRUE)
    g2<-dataGraph(list(breaks=h$breaks,density=h3,density1=h4),
                  xlabel="z[p]",ylabel="density",legend=legend,
                  hist=TRUE)
    g3<-dataGraph(list(breaks=seq(0,1,length.out=101),density=hw1,density1=hw2),
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
                                 tabs=c("Network","Histogram","Multiple"),
                                 tabContents=c(gs[1],gs[2],gs[3]),
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
    gs[3]<<-g3

    g<-generate_tab("Graphs: ",titleWidth=50,
                    tabs=c("Network","Histogram","Multiple"),
                    tabContents=c(gs[1],gs[2],gs[3]),
                    open=3
    ) 
    output$mainHTML <- renderUI(HTML(g))
  })
  
}