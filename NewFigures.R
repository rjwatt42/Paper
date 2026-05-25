###################################
setBrawEnv("RZ","z")
setBrawEnv("nNpoints",1001)
setggplot()
###################################
# Figure 1

zs<-0.35
rs<-tanh(zs)

setBrawEnv("RZ","z")
setDesign(sN=42)
setEvidence(sigOnly=0)

setWorld(TRUE,"Exp","z",0.31)
g<-showWorld(plotArea = c(0.05,0.25,0.45,0.5),showSingle=0.25)
g<-showWorldSampling(plotArea = c(0.55,0.25,0.45,0.5),
                     showSinglePopulation=0.25,
                     showSingleSample=rs,g=g)
print(g)
###################################
# Figure 2

zs<-0.35
rs<-tanh(zs)
n<-42

setBrawEnv("RZ","z")
setDesign(sN=n)
setEvidence(sigOnly=0)

world<-makeWorld(TRUE,"Single","z",0.18,pRplus=1)
setWorld(world)
g<-showWorldSampling(showSingleSample=rs,plotArea = c(0.01,0.55,0.3,0.4))
s<-exp(getLogLikelihood(atanh(rs),n,1,world$PDF,world$PDFk))
g<-addG(g,plotTitle(paste0("z[p]=",format(world$PDFk,digits=2),"  S=",format(s,digits=2))))

world<-makeWorld(TRUE,"Single","z",0.28,pRplus=1)
setWorld(world)
g<-showWorldSampling(showSingleSample=rs,plotArea = c(0.35,0.55,0.3,0.4),g=g)
s<-exp(getLogLikelihood(atanh(rs),n,1,world$PDF,world$PDFk))
g<-addG(g,plotTitle(paste0("z[p]=",format(world$PDFk,digits=2),"  S=",format(s,digits=2))))

world<-makeWorld(TRUE,"Single","z",0.48,pRplus=1)
setWorld(world)
g<-showWorldSampling(showSingleSample=rs,plotArea = c(0.69,0.55,0.3,0.4),g=g)
s<-exp(getLogLikelihood(atanh(rs),n,1,world$PDF,world$PDFk))
g<-addG(g,plotTitle(paste0("z[p]=",format(world$PDFk,digits=2),"  S=",format(s,digits=2))))

setWorld(TRUE,"Single")

result<-list(rIV=rs,nval=42,hypothesis=braw.def$hypothesis,design=braw.def$design,evidence=braw.def$evidence)
g<-showLikelihood(result,markRs=c(0.18,0.28,0.48),plotArea=c(0.25,0.05,0.5,0.5),g=g)

print(g)

###################################
# Figure 3

zs<-0.35
rs<-tanh(zs)
n<-42
prior<-makeWorld(TRUE,"Exp","z",0.31)

setBrawEnv("RZ","z")
setDesign(sN=n)

setWorld(TRUE,"Single")
setEvidence(sigOnly=0)
result<-list(rIV=rs,nval=n,hypothesis=braw.def$hypothesis,design=braw.def$design,evidence=braw.def$evidence)
g<-showLikelihood(result,showType="rp",prior=NULL,plotArea=c(0.01,0.25,0.4,0.4))
g<-addG(g,plotTitle(paste0("prior=","null")))

setWorld(TRUE,"Single")
setEvidence(sigOnly=0)
result<-list(rIV=rs,nval=n,hypothesis=braw.def$hypothesis,design=braw.def$design,evidence=braw.def$evidence)
g<-showLikelihood(result,showType="rp",prior=prior,plotArea=c(0.51,0.25,0.4,0.4),g=g)
g<-addG(g,plotTitle(paste0("prior=","exp")))

print(g)

###################################
# check figure 3 by simulation

setBrawEnv("RZ","z")

zs<-0.35
n<-42
sigOnly<-0

breaks<-seq(-1,1,0.01)*2
range<-seq(-1,1,length.out=1001)*2

prior<-makeWorld(TRUE,"Exp","z",0.31)
# prior<-makeWorld(TRUE,"Uniform","z",PDFspread=1.2)

setWorld(prior)
setDesign(sN=n)
setEvidence(shortHand=TRUE,sigOnly=sigOnly)

dens<-getLogLikelihood(zs,n,1,"Single",location=range,bias=0)
if (!is.null(prior)) {
  # zdens<-exp(-abs(range)/0.31)
  zdens<-zPopulationDist(range,braw.def$hypothesis$effect$world)
  dens<-dens+log(zdens)
}
dens<-exp(dens)
dens<-dens/sum(dens)*length(range)/(length(breaks)-1)

for (i in 1:100) {
  r1<-doMultiple(10000)
  use<-abs(atanh(r1$result$rIV)-zs)<0.02
  if (sigOnly>0) 
    use<-use & r1$result$pIV<0.05
  h1<-hist(atanh(r1$result$rpIV[use]),breaks=breaks,xlab="zp",main="")
  # title(main=paste0("MLE=",format(h1$mids[which.max(h1$counts)])))
  lines(range,dens*sum(use),col="red",lwd=2)
  title(main=paste0("MLE=",format(range[which.max(dens)])))
  print(c(i,r1$count))
}



###################################
# Figure 4

setBrawEnv("RZ","z")

zs<-0.35
n<-42
sigOnly<-0

setDesign(sN=n)
setEvidence(sigOnly=sigOnly)

rs<-tanh(zs)

setWorld(TRUE,"Exp","z",0.18)
g<-showWorldSampling(showSingleSample=rs,plotArea = c(0.01,0.55,0.3,0.4))
g<-addG(g,plotTitle("mean(Z[+])=0.18  S=0.52"))
setWorld(TRUE,"Exp","z",0.28)
g<-showWorldSampling(showSingleSample=rs,plotArea = c(0.35,0.55,0.3,0.4),g=g)
g<-addG(g,plotTitle("mean(Z[+])=0.28  S=0.56"))
setWorld(TRUE,"Exp","z",0.48)
g<-showWorldSampling(showSingleSample=rs,plotArea = c(0.69,0.55,0.3,0.4),g=g)
g<-addG(g,plotTitle("mean(Z[+])=0.48  S=0.51"))

setWorld(TRUE,"Exp","z",0.31)

result<-list(rIV=rs,nval=42,hypothesis=braw.def$hypothesis,design=braw.def$design,evidence=braw.def$evidence)
g<-showLikelihood(result,markRs=c(0.18,0.28,0.48),plotArea=c(0.25,0.05,0.5,0.5),g=g)

print(g)

###################################
# Check Figure 4 by simulation

setBrawEnv("RZ","z")

zs<-c(0.35,0.45,0.55,0.65)
n<-42
sigOnlySource<-0
sigOnlyAnalysis<-0
nsamps<-1000

breaks<-seq(-1,1,0.01)*4
range<-seq(-1,1,length.out=1001)*4

worlds<-seq(0.0,0.6,0.05)

setDesign(sN=n)
setEvidence(shortHand=TRUE,sigOnly=sigOnlySource)
 
zprime<-matrix(0,length(zs),length(worlds))

for (i in 1:1000) {
  zp<-worlds*0
  for (w in 1:length(worlds)) {
    setWorld(TRUE,"Exp","z",worlds[w],pRplus=1)
    r1<-doMultiple(nsamps)
    for (z in 1:length(zs)) {
      use<-abs(atanh(r1$result$rIV)-zs[z])<0.02
      zprime[z,w]<-zprime[z,w]+mean(use)
    }
  }
  print(c(i))
}

cols<-c("red","yellow","green","blue")
for (z in 1:length(zs)) {
  y<-zprime[z,]
  if (z==1)
    plot(worlds,y/max(y),'p')
  else
    lines(worlds,y/max(y),'p')
  lines(worlds,y/max(y),col=cols[z])
}
title(main=paste0("nsims=",format(i)))

###################################

n<-42
nStudies<-100
sigOnlySource<-1
sigOnlyAnalysis<-1
worlds<-seq(0.0,0.6,0.1)
prps<-seq(0,1,0.2)

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


############################

g<-startPlot(c(0,0.7),c(0,1.2),
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


###################################

setWorld(TRUE,"Uniform","z")
setDesign(sN=42)
setEvidence(shortHand=TRUE,sigOnly=0)

for (i in 1:100) {
  r<-doMultiple(10000)
  use<-abs(atanh(r$result$rIV)-0.4)<0.05
  usep<-r$result$pIV<0.05
  hist(atanh(r$result$rpIV[use & usep]),breaks=seq(-1,1,0.01)*2)
  print(c(i,r$count))
}


###################################
setBrawEnv("RZ","z")

zs<-0.35
n<-42
sigOnly<-0

breaks<-seq(-1,1,0.01)*2
range<-breaks[2:length(breaks)]-diff(breaks[1:2])/2

prior<-makeWorld(TRUE,"Exp","z",0.31)
prior<-makeWorld(TRUE,"Uniform","z")

setWorld(prior)
setDesign(sN=n)
setEvidence(shortHand=TRUE,sigOnly=sigOnly)

dens<-getLogLikelihood(zs,n,1,"Single",location=range,bias=0)
if (!is.null(prior)) {
  # zdens<-exp(-abs(range)/0.3)
  zdens<-zPopulationDist(range,braw.def$hypothesis$effect$world)
  dens<-dens+log(zdens)
}
dens<-exp(dens)
dens<-dens/sum(dens)

for (i in 1:100) {
  r1<-doMultiple(10000)
  use<-abs(atanh(r1$result$rIV)-zs)<0.02
  if (sigOnly>0) 
    use<-use & r1$result$pIV<0.05
  h1<-hist(atanh(r1$result$rpIV[use]),breaks=breaks,xlab="zp",main="")
  title(main=paste0("MLE=",format(h1$mids[which.max(h1$counts)])))
  lines(range,dens*sum(use),col="red",lwd=2)
  print(c(i,r1$count))
}


###################################

setBrawEnv("RZ","z")

setWorld(TRUE,"Exp","z",0.3)
setDesign(sN = 52,sNRand = TRUE,sNRandSD = 33.3)
setDesign(sN=42)
setEvidence(sigOnly=1)

dm<-doMultiple(10,NULL)
g<-showLikelihood(dm$result,norm=TRUE,plotArea=c(0.01,0.55,0.3,0.4))
g<-addG(g,plotTitle("n(r[s])=10"))
dm<-doMultiple(100,NULL)
g<-showLikelihood(dm$result,norm=TRUE,plotArea=c(0.35,0.55,0.3,0.4),g=g)
g<-addG(g,plotTitle("n(r[s])=100"))
dm<-doMultiple(1000,NULL)
g<-showLikelihood(dm$result,norm=TRUE,plotArea=c(0.69,0.55,0.3,0.4),g=g)
g<-addG(g,plotTitle("n(r[s])=1000"))
print(g)

####################################
