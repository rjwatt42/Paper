######################################

library(pracma)
source("basicPlot.r")
source("SimulateNetwork/doNetwork.R")
source("SimulateNetwork/makeNetwork.R")
source("SimulateNetwork/plotNetwork.R")
source("SimulateNetwork/getNetworkHist.R")
source("SimulateNetwork/plotNetworkHist.R")

initGraph("HTML",gsize=650,autoShow=FALSE)

######################################

# network structure
nStages=16
nVarsPerStage=8
rangeLink=2
probLink=0.5
strengthLink=0.65
fullModel=TRUE

h<-list(breaks=seq(0,(strengthLink*1.5),length.out=101),density=0)
avoidZeros<-TRUE

######################################

# make a network
fullLinks<-makeNetwork(nStages,nVarsPerStage,rangeLink,probLink)
g<-plotNetwork(links2Path(fullLinks,nStages,nVarsPerStage),showNames=FALSE)
# showHTML(g)

h1<-getNetworkHist(fullLinks,fullModel,strengthLink,avoidZeros,h)

g1<-plotNetworkHist(h1)
showHTML(g1)

######################################

# make a histogram of effect sizes

density<-0
ncount=1000
ests<-zeros<-c()
indirect<-direct<-c()
for (ni in 1:ncount) {
  
  # make a network
  fullLinks<-makeNetwork(nStages,nVarsPerStage,rangeLink,probLink)
  # get effect sizes & distribution
  h1<-getNetworkHist(fullLinks,fullModel,strengthLink,avoidZeros,h)
  
  # save results
  density<-density+h1$hist$density/ncount
  ests<-c(ests,h1$est$estimate)
  zeros<-c(zeros,h1$zeros)
}


h1$hist$density<-density
h1$est$estimate<-mean(ests)
h1$zeros<-mean(zeros)
g<-plotNetworkHist(h1)
showHTML(g)

###############################
  # ns<-1000
  # testData=matrix(0,ns,nVarsPerStage*nStages)
  # for (t in seq(nVarsPerStage*nStages,1,-1)) {
  #   use=which(fullLinks[,t]==1)
  #   if (length(use)==0) {
  #     testData[,t]=rnorm(ns)
  #   } else {
  #     usedVar=length(use)*(strengthLink^2)
  #     if (usedVar>1) usedVar=1
  #     z=0
  #     for (ui in 1:length(use))
  #       z=z+testData[,use[ui]]*strengthLink
  #     testData[,t]=z+
  #               sqrt(1-usedVar)*rnorm(ns)
  #   }
  # }
  # 
  # Stheta1=abs(cor(testData))
  # zp<-atanh(Stheta1[Stheta1<1])
  # h1<-hist(zp[zp<(strengthLink*2)],seq(0,(strengthLink*2),length.out=101),plot=FALSE)

  


