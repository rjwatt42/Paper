#### starting ################
devtools::load_all(path="~/Documents/GitHub/BrawPackDevelopment/")
setBrawEnv('RZ','z')
setBrawEnv('z_range',3)
setBrawEnv('maxN',350)
setBrawEnv('graphOrientation','horz')
setBrawEnv("graphicsType","HTML")
setBrawEnv("verbose",FALSE)

##### Main Parameters ######################################

distribution<-"Gamma"
# Exp Gauss Gamma GenExp
publicationBias<-TRUE
N<-50

resetBraw<-function(dist=distribution,sN=N,sNRand=FALSE,bias=publicationBias) {
  setWorld(TRUE,dist,"z",0.3,PDFshape=4,pRplus = 1)
  setDesign(sN=sN,sNRand=sNRand)
  setEvidence(shortHand=TRUE, sigOnly=bias)
}

#### make samples #####################

resetBraw()

nStudies<-10000

g<-showMultiple(doMultiple(nStudies),showType="rs")
showHTML(g)

# ans<-fitdist(atanh(abs(braw.res$multiple$result$rpIV)), distr = "gamma", method = "mle")
# print(
#   c(ans$estimate[["shape"]],
#     ans$estimate[["shape"]]/ans$estimate[["rate"]]
#     )
# )

zs<-atanh(braw.res$multiple$result$rIV)
ns<-braw.res$multiple$result$nval

#### basic llk function for scale ############################

resetBraw()

scales<-seq(0.1,0.6,length.out=101)

llks<-scales*0
for (i in 1:length(scales)) {
  llks[i]<-getLogLikelihood(abs(zs),ns,1,distribution,scales[i],shape=4,bias=publicationBias)
}

g<-dataGraph(data.frame(x=scales,y=llks),fill=NA,
             xlabel="scale",ylabel="d(llk)",
             title=paste0("MLE = ",scales[which.max(llks)]))
showHTML(g)

#### basic llk function for shape ############################

resetBraw()

shapes<-seq(2,6,length.out=101)

llks<-shapes*0
for (i in 1:length(shapes)) {
  llks[i]<-getLogLikelihood(abs(zs),ns,1,distribution,0.3,shape=shapes[i],bias=publicationBias)
  # setWorld(TRUE,distribution,"z",0.3,PDFshape=shapes[i],pRplus = 1)
  # q<-makeTheoryMultiple(braw.def$hypothesis,braw.def$design,braw.def$evidence,"rs")
  # llks[i]<-sum(log(approx(q$theoryVals,q$theoryDens_all,zs)$y))
}

g<-dataGraph(data.frame(x=shapes,y=llks),fill=NA,
             xlabel="shape",ylabel="d(llk)")
showHTML(g)


#### ditto scale and shape ############################

resetBraw()

shapes<-seq(2,6,length.out=11)
scales<-seq(0.1,0.5,length.out=11)
llks<-matrix(0,length(scales),length(shapes))
for (i in 1:length(shapes)) {
  for (j in 1:length(scales)) {
    llks[j,i]<-getLogLikelihood(abs(zs),ns,1,distribution,scales[j],shape=shapes[i],bias=publicationBias)
    # setWorld(TRUE,distribution,"z",scales[j],PDFshape=shapes[i],pRplus = 1)
    # q<-makeTheoryMultiple(braw.def$hypothesis,braw.def$design,braw.def$evidence,"rs")
    # q$theoryDens_all<-q$theoryDens_all/sum(q$theoryDens_all)
    # llks[j,i]<-sum(log(approx(q$theoryVals,q$theoryDens_all,zs)$y))
  }
}

Smax<- max(llks,na.rm=TRUE)
use<-which(llks==Smax, arr.ind = TRUE)
print(
  paste0("scale=",scales[use[1]],"  shape=",shapes[use[2]])
)


###########################
m1<-makeMetaAnalysis(On=TRUE, nstudies=nStudies,
                     analysisType="world",modelPDF=distribution,
                     sourceBias=publicationBias,assumeBias=publicationBias,
                     sourceAbs=TRUE
)


#### three different routes to pdf #############################

resetBraw()

q<-makeTheoryMultiple(braw.def$hypothesis,braw.def$design,braw.def$evidence,"rs")
use<-q$theoryVals>=0
g<-dataGraph(data.frame(x=q$theoryVals[use],y=q$theoryDens_all[use]),fill=NA,
             xlabel="z",ylabel="d(llk)")
showHTML(g)

z<-seq(0,3,length.out=1000)
q<-getLogLikelihood(z,braw.def$design$sN,1,distribution,scale=0.3,shape=4,
                    bias=publicationBias,returnVals=TRUE)
g<-dataGraph(data.frame(x=z,y=10^q),fill=NA,
             xlabel="z",ylabel="d(llk)")
showHTML(g)

switch(distribution,
       "Gamma"=PDF<-GammaSamplingPDF,
       "Exp"=PDF<-ExpSamplingPDF,
       "Gauss"=PDF<-GaussSamplingPDF,
       "GenExp"=PDF<-GenExpSamplingPDF
)
q<-PDF(z,0.45,sqrt(1/(braw.def$design$sN-3)),df1=1,spread=0,shape=4,bias=publicationBias)
g<-dataGraph(data.frame(x=z,y=q$pdf),fill=NA,
             xlabel="z",ylabel="d(llk)")
showHTML(g)

##########################################

scales<-seq(0.1,0.6,length.out=101)

res<-c()
for (i in 1:length(scales)) {
  q<-PDF(z,scales[i],sqrt(1/(braw.def$design$sN-3)),df1=1,spread=0,shape=4,bias=1)
  q$pdf<-q$pdf/q$sig_compensate
  print(sum(q$pdf))
  zsl<-sum(log(approx(z,q$pdf,abs(zs))$y))
  res<-c(res,zsl)
}

g<-dataGraph(data.frame(x=scales,y=res),fill=NA,
             xlabel="scale",ylabel="d(llk)",
             title=scales[which.max(res)]
)
showHTML(g)

#### a few example metaAnalyses #############################

resetBraw(sNRand=FALSE)

nStudies<-1000

m1<-makeMetaAnalysis(On=TRUE, nstudies=nStudies,
                     analysisType="world",modelPDF=distribution,
                     sourceBias=publicationBias,assumeBias=publicationBias,
                     sourceAbs=TRUE
)

for (i in 1:10) {
  # mle<-getMaxLikelihood(zs,ns,df1=1,distribution,m1,braw.def$hypothesis) 
  
  d1<-doMetaAnalysis(metaAnalysis = m1)
  print(c(d1$best$PDFshape,d1$best$PDFk))
}

#### explore metaAnalysis #############################

resetBraw(sNRand=FALSE)

nStudies<-1000

exploreType<-"PDFshape"
exploreType<-"PDFk"
# "PDFk"  "PDF"  "PDFshape"  "p(R+)" "NoStudies"
showType<-"PDFk"
if (showType=="PDF") model<-"All" else model<-distribution

m1<-makeMetaAnalysis(On=TRUE, nstudies=nStudies,
                     analysisType="world",modelPDF=model,
                     sourceBias=publicationBias,assumeBias=publicationBias,
                     sourceAbs=TRUE
)

setBrawDef("metaAnalysis",m1)

shapes<-seq(1,6)
ems<-c()
for (j in 1:length(shapes)) {
  setWorld(PDFshape=shapes[j])
  
  e1<-makeExplore(exploreType,vals=seq(0.1,0.6,0.1),exploreNPoints=6) 

  em<-NULL
  for (i in 1:50) {
    em<-doExplore(1,em,explore=e1,doingMetaAnalysis=TRUE)
    gm<-showExplore(em,showType = showType)
    print(c(j,i))
  }
  showHTML(gm)
  ems<-c(ems,em)
}
