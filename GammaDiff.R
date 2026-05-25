
setEvidence(sigOnly=TRUE)
setDesign(sN=50,sNRand = TRUE)

shapes<-seq(1,5,0.1)
result<-shapes*0
for (i in 1:length(shapes)) {
setWorld(TRUE,"Gamma","z",0.3,shapes[i])
k<-makeTheoryMultiple(braw.def$hypothesis,braw.def$design,braw.def$evidence,"rs")
k$theoryDens_sig<-k$theoryDens_sig/sum(k$theoryDens_sig)

setWorld(TRUE,"Gamma","z",0.3,shapes[i]+0.1)
k1<-makeTheoryMultiple(braw.def$hypothesis,braw.def$design,braw.def$evidence,"rs")
k1$theoryDens_sig<-k1$theoryDens_sig/sum(k1$theoryDens_sig)

result[i]<-sum(abs(k1$theoryDens_sig-k$theoryDens_sig))
print(sum(abs(k1$theoryDens_sig-k$theoryDens_sig)))
}
