############################

showResult<-function(worlds,prps,resk1,resp1,ns1,resk2,resp2,ns2) {
  xlim<-c(resk1/ns1,resk2/ns2)
  xlim<-c(min(xlim),max(xlim))
  ylim<-c(resp1/ns1,resp2/ns2)
  ylim<-c(min(ylim),max(ylim))
  g<-startPlot(xlim+c(-1,1)*diff(xlim)/10,ylim+c(-1,1)*diff(ylim)/10,
               xlabel="k",xticks=worlds,
               ylabel="pR+",yticks=prps)
  d<-data.frame(x=array(matrix(worlds,length(worlds),length(prps))),
                y=array(matrix(prps,length(worlds),length(prps),byrow=TRUE)))
  g<-addG(g,dataPoint(d,fill="black",size=1.5))
  d<-data.frame(x=array(resk1/ns1),y=array(resp1/ns1))
  g<-addG(g,dataPoint(d,fill="red"))
  d<-data.frame(x=array(resk2/ns2),y=array(resp2/ns2))
  for (w in 1:length(worlds))
    g<-addG(g,dataPath(data.frame(x=resk1[w,]/ns1[w,],y=resp1[w,]/ns1[w,]),colour="red"))
  for (p in 1:length(prps))
    g<-addG(g,dataPath(data.frame(x=resk1[,p]/ns1[,p],y=resp1[,p]/ns1[,p]),colour="red"))
  g<-addG(g,dataPoint(d,fill="blue"))
  for (w in 1:length(worlds))
    g<-addG(g,dataPath(data.frame(x=resk2[w,]/ns2[w,],y=resp2[w,]/ns2[w,]),colour="blue"))
  for (p in 1:length(prps))
    g<-addG(g,dataPath(data.frame(x=resk2[,p]/ns2[,p],y=resp2[,p]/ns2[,p]),colour="blue"))
  g<-addG(g,plotTitle(paste0("nsims = ",min(ns1))))
  print(g)
}

###################################

n<-42
nStudies<-10000
worlds<-seq(0.1,0.6,0.1)
prps<-seq(0.2,1,0.2)

setDesign(sN=n)
setDesign(getDesign("Psych"))
setEvidence(shortHand=TRUE,sigOnly=0)

mA1<-makeMetaAnalysis(analysisType="world",modelPDF="Exp",modelNulls=TRUE,sourceBias=0)
mA2<-makeMetaAnalysis(analysisType="world",modelPDF="Exp",modelNulls=TRUE,sourceBias=1)

if (!exists('ns1'))
  ns2<-ns1<-resp2<-resp1<-resk2<-resk1<-matrix(0,length(worlds),length(prps))

for (i in 1:1000) {
  for (w in 1:length(worlds)) { 
    for (p in 1:length(prps)) {
      setWorld(TRUE,"Exp","z",worlds[w],pRplus=prps[p])
      if (ns1[w,p]<i) {
        setEvidence(shortHand=TRUE,sigOnly=0)
        rs<-doMultiple(nStudies,NULL)
        m<-doMetaAnalysis(rs,mA1,keepStudies = TRUE)
        resk1[w,p]<-resk1[w,p]+m$best$PDFk
        resp1[w,p]<-resp1[w,p]+m$best$pRplus
        ns1[w,p]<-ns1[w,p]+1
      }
      
      if (ns2[w,p]<i) {
        setEvidence(shortHand=TRUE,sigOnly=1)
        rs<-doMultiple(nStudies,NULL)
        m<-doMetaAnalysis(rs,mA2,keepStudies = TRUE)
        resk2[w,p]<-resk2[w,p]+m$best$PDFk
        resp2[w,p]<-resp2[w,p]+m$best$pRplus
        ns2[w,p]<-ns2[w,p]+1
      }
      print(c(i,w,p))
    }
  }
}

showResult(worlds,prps,resk1,resp1,ns1,resk2,resp2,ns2)

file<-paste0('/Users/rogerwatt/Documents/GitHub/Paper/Sim',format(nStudies),'.RData')
save(ns1,resk1,resp1,ns2,resk2,resp2,ns2,file=file)


###################################

n<-42
nStudies<-100
worlds<-seq(0.1,0.6,0.05)
prps<-seq(0.1,1,0.1)

setDesign(sN=n)
setDesign(getDesign("Psych"))
setDesign(sN=200)
setEvidence(shortHand=TRUE,sigOnly=0)

mA1<-makeMetaAnalysis(analysisType="world",modelPDF="Exp",modelNulls=TRUE,sourceBias=0)
mA2<-makeMetaAnalysis(analysisType="world",modelPDF="Exp",modelNulls=TRUE,sourceBias=1)

if (!exists('ns1'))
  ns2<-ns1<-resp2<-resp1<-resk2<-resk1<-matrix(0,length(worlds),length(prps))

for (i in 1:1000) {
  for (w in 1:length(worlds)) { 
    for (p in 1:length(prps)) {
      setWorld(TRUE,"Exp","z",worlds[w],pRplus=prps[p])
      if (ns1[w,p]<i) {
        setEvidence(shortHand=TRUE,sigOnly=0)
        rs<-doMultiple(nStudies,NULL)
        m<-doMetaAnalysis(rs,mA1,keepStudies = TRUE)
        resk1[w,p]<-resk1[w,p]+m$best$PDFk
        resp1[w,p]<-resp1[w,p]+m$best$pRplus
        ns1[w,p]<-ns1[w,p]+1
      }
      
      if (ns2[w,p]<i) {
        setEvidence(shortHand=TRUE,sigOnly=1)
        rs<-doMultiple(nStudies,NULL)
        m<-doMetaAnalysis(rs,mA2,keepStudies = TRUE)
        resk2[w,p]<-resk2[w,p]+m$best$PDFk
        resp2[w,p]<-resp2[w,p]+m$best$pRplus
        ns2[w,p]<-ns2[w,p]+1
      }
      print(c(i,w,p))
    }
  }
}
  showResult(worlds,prps,resk1,resp1,ns1,resk2,resp2,ns2)
  
  save(ns1,resk1,resp1,ns2,resk2,resp2,ns2,file='/Users/rogerwatt/Documents/GitHub/Paper/Sim4.RData')
  
############################
