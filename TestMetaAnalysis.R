############
setEvidence(sigOnly=TRUE,absOnly=TRUE)
mA<-makeMetaAnalysis(TRUE,200,"fixed",sourceBias=TRUE,sourceAbs=TRUE)
mM<-doMetaAnalysis(NULL,mA)

showMetaSingle(mM)

############
setEvidence(sigOnly=TRUE,absOnly=TRUE)
mA<-makeMetaAnalysis(TRUE,110000,"fixed",sourceBias=TRUE,sourceAbs=TRUE)

nsims<-10
while (length(mM$best$Smax)<1000) {
mM<-doMetaMultiple(nsims,metaAnalysis=mA)
print(c(length(mM$best$Smax),mean(mM$best$Smax),std(mM$best$Smax)))
print(showMetaMultiple(mM,showType="metaSmax;metaRiv"))
}

############
setDesign(getDesign("Psych"))
setWorld(TRUE,"Gamma","z",0.3,2,pRplus=1)

mA<-makeMetaAnalysis(TRUE,200,"world",modelPDF="Gamma")
mM<-doMetaAnalysis(NULL,mA)

showMetaSingle(mM)

##############

setEvidence(sigOnly=TRUE,absOnly=TRUE)
setEffect(rIV=0.1)
mA<-makeMetaAnalysis(TRUE,200,"fixed",sourceBias=TRUE,sourceAbs=TRUE)
mM<-doMetaAnalysis(NULL,mA)

showMetaSingle(mM)


##############
