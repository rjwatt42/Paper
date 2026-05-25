###############################

zs<-0.3
n<-42
setBrawEnv("RZ",'z')

m1<-makeMetaAnalysis(On=TRUE, nstudies=1,studies<-list(rIV=tanh(zs),nval=n,df1=1),
                     analysisType="world",modelPDF="Exp"
)
d1<-doMetaAnalysis(metaAnalysis = m1)

lambdas1<-d1$exp$Spts$PDFk
llk1<-d1$exp$Svals
g<-dataGraph(data.frame(x=lambdas1,y=llk1),fill='yellow',
             xlabel="z[λ]",ylabel="log likelihood",
             title=paste0("z[s]=",format(zs),'  n=',format(n)))

print(g)
###############################

nStudies=1000
biasedSamples<-FALSE
assumeBias<-FALSE

setWorld(TRUE,"Exp","z",0.3,pRplus = 1)
setDesign(sNRand=TRUE)
setBrawEnv("RZ",'z')

m1<-makeMetaAnalysis(On=TRUE, nstudies=nStudies, studies=NULL,
                     analysisType="world",modelPDF="Exp",
                     sourceBias=biasedSamples,assumeBias=assumeBias
)
d1<-doMetaAnalysis(NULL,metaAnalysis = m1)

lambdas1<-d1$exp$Spts$PDFk
llk1<-d1$exp$Svals
g<-dataGraph(data.frame(x=lambdas1,y=llk1),fill='yellow',
             xlabel="z[λ]",ylabel="log likelihood",
             title=paste0("nStudies=",format(nStudies)))

print(g)

###########################
nStudies=1000
biasedSamples<-FALSE
assumeBias<-FALSE

setWorld(TRUE,"Exp","z",0.3,pRplus = 1)
setDesign(sNRand=TRUE)
setBrawEnv("RZ",'z')

m1<-makeMetaAnalysis(On=TRUE, nstudies=nStudies, studies=NULL,
                     analysisType="world",modelPDF="Exp",
                     sourceBias=biasedSamples,assumeBias=assumeBias
)
m2<-makeMetaAnalysis(On=TRUE, nstudies=nStudies, studies=NULL,
                     analysisType="world",modelPDF="Gauss",
                     sourceBias=biasedSamples,assumeBias=assumeBias
)
d1<-doMetaAnalysis(NULL,metaAnalysis = m1)
d2<-doMetaAnalysis(NULL,metaAnalysis = m2)
xpts<-c(d1$exp$Spts$PDFk,d2$gauss$Spts$PDFk)
ypts<-c(d1$exp$Svals,d2$gauss$Svals)
g<-dataGraph(data.frame(x=d1$exp$Spts$PDFk,y=d1$exp$Svals),fill='yellow',
             xlim=c(min(xpts),max(xpts))+c(-1,1)*0.1*(max(xpts)-min(xpts)),
             ylim=c(min(ypts),max(ypts))+c(-1,1)*0.1*(max(ypts)-min(ypts)),
             xlabel="z[λ]",ylabel="log likelihood",
             title=paste0("nStudies=",format(nStudies)))
g<-dataGraph(data.frame(x=d2$gauss$Spts$PDFk,y=d2$gauss$Svals),fill='red',g=g)
g<-addG(g,dataLegend(data.frame(names=c("Exp","Gauss"),colours=c("yellow","red"))))
print(g)

###############################

nStudies=1000
biasedSamples<-TRUE

setWorld(TRUE,"Exp","z",0.3,pRplus = 1)
setDesign(sN=50,sNRand=TRUE)
setEvidence("sigOnly",biasedSamples)

m1<-makeMetaAnalysis(On=TRUE, nstudies=nStudies,
                     analysisType="world",modelPDF="Exp",
                     sourceBias=biasedSamples,assumeBias=FALSE
                     )
m2<-makeMetaAnalysis(On=TRUE, nstudies=nStudies,
                     analysisType="world",modelPDF="Exp",
                     sourceBias=biasedSamples,assumeBias=TRUE
                     )

d1<-doMetaAnalysis(metaAnalysis = m1)
d2<-doMetaAnalysis(metaAnalysis = m2)

g<-dataGraph(data=list(x=rbind(d1$exp$Spts$PDFk,d2$exp$Spts$PDFk),
                       y=rbind(d1$exp$Svals,d2$exp$Svals),
                       fill=c('yellow','red')
                  ),
             legend=list(names=c(paste0("assume: "," z[MLE]=",format(mleFull2,digits=3)),
                                 paste0("ignore:"," z[MLE]=",format(mleFull1,digits=3))
                               ),
                         legendTitle = "Publication bias"),
             xlabel="z[λ]",ylabel="log likelihood",
             title=paste0("no studies=",format(nStudies))
)

print(g)



