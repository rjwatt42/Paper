###############################

r<-list(hypothesis=makeHypothesis(),
        design=getDesign("Psych"),
        evidence=makeEvidence(),
        result=list(nval=my_data$n,
                    rIV=my_data$r_s,
                    rpIV=my_data$r_s,
                    pIV=my_data$p,
                    df1=my_data$df1)
        )

setHTML()
setBrawEnv("maxBins",201)
setBrawEnv("nPlotScale","linear")
showHTML(showMultiple(r,showType="n",orientation = "horz", showTheory = FALSE))

#######################

analyseBias<-FALSE
model<-"Exp"

r1<-list(hypothesis=makeHypothesis(),
        design=makeDesign(),
        evidence=makeEvidence(),
        result=list(nval=my_data$n,
                    rIV=my_data$r_s,
                    rpIV=my_data$r_s,
                    pIV=my_data$p,
                    df1=my_data$df1
        )
)

mA<-makeMetaAnalysis(analysisType="world",
                     modelNulls=TRUE,modelPDF=model,
                     sourceBias=TRUE,
                     analyseBias=analyseBias)

m1<-doMetaAnalysis(r1,metaAnalysis=mA,keepStudies = TRUE)

showMetaSingle(m1,autoYlim = FALSE)

#######################
# make a simulated set of data

rsim<-r1
world<-makeWorld(TRUE,"Exp","z",PDFk=m1$best$param1,pRPlus=m1$best$param2)

rs_all<-c()
n_all<-c()
for (n in unique(r1$result$nval)) {
  use<-r1$result$nval==n
  count<-sum(use)
  rp<-rRandomValue(world,count*10)$use
  rs<-tanh(atanh(rp)+rnorm(length(rp),0,1/sqrt(n-3)))
  p<-rn2p(rs,rep(n,count))
  keep<-which(p<0.05)
  keep<-keep[1:count]
  rs_all<-c(rs_all,rs[keep])
  n_all<-c(n_all,rep(n,count))
}
rsim$result$nval<-n_all
rsim$result$rIV<-rs_all

###############################
# by sample size

analyseBias<-FALSE
r<-r1
# r<-rsim

ww0<-c()
zk0<-c()
pn0<-c()
pb0<-c()
nu0<-c()
nn0<-c()

nTotal<-c()
nn<-10^seq(log10(10),log10(500),length.out=21)
nn<-c(sort(unique(r$result$nval)),1000)
for (ni in 1:(length(nn)-1)) {
  use<-which(r$result$nval>=nn[ni] & r$result$nval<nn[ni+1])
  r0<-r
  r0$result$rIV<-r0$result$rIV[use]
  r0$result$nval<-r0$result$nval[use]
  r0$result$df1<-r0$result$df1[use]

  mA<-makeMetaAnalysis(analysisType="world",
                       modelNulls=TRUE,modelPDF=model,
                       sourceBias=1,analyseBias=analyseBias)
  m<-doMetaAnalysis(r0,metaAnalysis=mA,keepStudies = TRUE)
  
  w<-ExpSamplingCDF(wn2z(0.5,n),m$best$PDFk,1/sqrt(n-3))
  w<-w*m$best$pRPlus+0.05*(1-m$best$pRPlus)
  ww0<-c(ww0,w)
  nu0<-c(nu0,nn[ni])
  zk0<-c(zk0,m$best$PDFk)
  pn0<-c(pn0,m$best$pRPlus)
  pb0<-c(pb0,m$best$sigOnly)
  nn0<-c(nn0,length(use))

  plot(log10(nu0),zk0,col='red','p',cex=1,pch=16,ylim=c(0,1))
  lines(log10(nu0),zk0,col='red')
  lines(log10(nu0),pn0,col='blue','p',cex=1,pch=16)
  lines(log10(nu0),pn0,col='blue')
  if (analyseBias) {
    lines(log10(nu0),pb0,col='green','p',cex=1,pch=16)
    lines(log10(nu0),pb0,col='green')
  }
  
  nTotal<-c(nTotal,rep(nn[ni],length(use)/w))
}

q<-fitdistrplus::fitdist(nTotal, distr = "gamma", method = "mle")

###############################
# by year
analyseBias<-FALSE

yr1<-c()
zk1<-c()
pn1<-c()
pb1<-c()
nn1<-c()
for (yr in unique(my_data$year)) {
  use<-which(my_data$year==yr)
  r<-list(hypothesis=makeHypothesis(),
          design=makeDesign(),
          evidence=makeEvidence(),
          result=list(nval=my_data$n[use],
                      rIV=my_data$r_s[use],
                      rpIV=my_data$r_s[use],
                      pIV=my_data$p[use],
                      df1=my_data$df1[use]
          )
  )
  
  mA<-makeMetaAnalysis(analysisType="world",
                       modelNulls=TRUE,modelPDF=model,
                       sourceBias=1,analyseBias=analyseBias)
  m<-doMetaAnalysis(r,metaAnalysis=mA,keepStudies = TRUE)
  
  yr1<-c(yr1,yr)
  zk1<-c(zk1,m$best$param1)
  pn1<-c(pn1,m$best$param2)
  pb1<-c(pb1,m$best$param3)
  nn1<-c(nn1,sum(use))
  print(paste(format(yr,width=8),
              format(m$best$param1,digits=2,nsmall=2),
              format(m$best$param2,digits=2,nsmall=2),
              format(m$best$param3,digits=2,nsmall=2)
  )
  )
  
plot(yr1,zk1,col='red','p',cex=1,pch=16,ylim=c(0,1))
lines(yr1,zk1,col='red')
lines(yr1,pn1,col='blue','p',cex=1,pch=16)
lines(yr1,pn1,col='blue')
if (analyseBias) {
  lines(yr1,pb1,col='green','p',cex=1,pch=16)
  lines(yr1,pb1,col='green')
}

}


################################
# by journal

jn1<-c()
zk1<-c()
pn1<-c()
pb1<-c()
nn1<-c()
for (jnl in unique(my_data$journal)) {
  use<-which(my_data$journal==jnl)
  r<-list(hypothesis=makeHypothesis(),
          design=makeDesign(),
          evidence=makeEvidence(),
          result=list(nval=my_data$n[use],
                      rIV=my_data$r_s[use],
                      rpIV=my_data$r_s[use],
                      pIV=my_data$p[use],
                      df1=my_data$df1[use]
          )
  )
  
  mA<-makeMetaAnalysis(analysisType="world",
                       modelNulls=TRUE,modelPDF=model,
                       sourceBias=1,analyseBias=TRUE)
  m<-doMetaAnalysis(r,metaAnalysis=mA,keepStudies = TRUE)
  
  jn1<-c(jn1,jnl)
  zk1<-c(zk1,m$best$param1)
  pn1<-c(pn1,m$best$param2)
  pb1<-c(pb1,m$best$param3)
  nn1<-c(nn1,sum(use))
  print(paste(format(jnl,width=8),
              format(m$best$param1,digits=2,nsmall=2),
              format(m$best$param2,digits=2,nsmall=2),
              format(m$best$param3,digits=2,nsmall=2)
              )
        )
  
}

plot(1:length(jn1),zk1,col='red','p',cex=1,pch=16,ylim=c(0,1))
lines(1:length(jn1),zk1,col='red')
lines(1:length(jn1),pn1,col='blue','p',cex=1,pch=16)
lines(1:length(jn1),pn1,col='blue')
lines(1:length(jn1),pb1,col='green','p',cex=1,pch=16)
lines(1:length(jn1),pb1,col='green')

################################
# individual sample sizes - with histograms

model<-"Exp"
nval<-60
zcrit<-wn2z(0.5,nval)

n<-which(nu0==nval)
use<-which(my_data$n==nval)

r<-list(hypothesis=makeHypothesis(),
        design=makeDesign(),
        evidence=makeEvidence(),
        result=list(nval=my_data$n[use],
                    rIV=my_data$r_s[use],
                    rpIV=my_data$r_s[use],
                    pIV=my_data$p[use],
                    df1=my_data$df1[use]
        ),
        metaAnalysis=makeMetaAnalysis(analysisType="world",modelNulls=TRUE,modelPDF=model)
)

m<-doMetaAnalysis(r,metaAnalysis=makeMetaAnalysis(analysisType="world",
                                                  modelNulls=TRUE,modelPDF=model,sourceBias=1),
                  keepStudies = TRUE)

zvals<-atanh(my_data$r_s[use])
zbreaks<-seq(zcrit,max(zvals)+0.0001,length.out=51)
zv<-zbreaks[2:length(zbreaks)]-0.5*diff(zbreaks[1:2])

rr<-m$best$param1
pn<-m$best$param2

switch(model,
       "Single"=z1<-SingleSamplingPDF(zv,rr,1/sqrt(nval-3)),
       "Exp"=z1<-ExpSamplingPDF(zv,atanh(rr),1/sqrt(nval-3)),
       "Gauss"=z1<-GaussSamplingPDF(zv,atanh(rr),1/sqrt(nval-3))
       )

z0<-SingleSamplingPDF(zv,0,1/sqrt(nval-3))

z<-z1$pdf*(1-pn)+z0$pdf*pn
zcrit<-wn2r(0.5,nval)
z[zv<zcrit]<-0
z<-z*(length(use)/sum(z))

h<-hist(zvals,breaks=zbreaks,plot=FALSE)
hist(zvals,breaks=zbreaks,
     main=paste0("S=",format(m$best$S,digits=4),"  n=",format(length(use)),"    S/n=",format(m$best$S/length(use),digits=4)))
lines(zv,z,col="red",lwd=3)
for (i in 1:length(z)) lines(c(0,0)+zv[i],c(h$counts[i],z[i]),col="black",type="o",cex=0.5,pch=16,lwd=2)

###############################
