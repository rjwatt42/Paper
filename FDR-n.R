########################################
# get w and fdr for different sample sizes

nmaxs<-round(10^seq(log10(20),log10(500),length.out=10))
metaAnal<-list(meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE)

my_data_sim<-simData(an)


nRes<-c()
kRes<-c()
sRes<-c()

w1<-c()
w0<-c()
w<-c()
fdr<-c()

w1s<-c()
w0s<-c()
ws<-c()
fdrs<-c()

nn<-c()
for (i in 2:length(nmaxs)) {
  if (i==1) {
    use<-(my_data$n<=nmaxs[i])
  } else {
    use<-(my_data$n>nmaxs[i-1] & my_data$n<=nmaxs[i])
  }
  n<-median(my_data$n[use])
  nn<-c(nn,n)
  
  metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$n[use],df1=my_data$df1[use]))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$best$Nullmax)
  kRes<-cbind(kRes,an$best$Kmax)
  sRes<-cbind(sRes,an$best$Smax)
  showAnalysis(an,paste0("n=",round(nmaxs[i])," "))

  w1t<-ExpSamplingCDF(atanh(p2r(0.05,n)),an$best$Kmax,1/sqrt(n-3))
  w0t<-alpha
  wt<-w0t*an$best$Nullmax+w1t*(1-an$best$Nullmax)
  
  w0<-c(w0,w0t)  
  w1<-c(w1,w1t)  
  w<-c(w,wt)  
  fdr<-c(fdr,w0t*an$best$Nullmax/wt)
  
  if (i==1) {
    use<-(my_data_sim$n<=nmaxs[i])
  } else {
    use<-(my_data_sim$n>nmaxs[i-1] & my_data_sim$n<=nmaxs[i])
  }
  metaData<-list(result=list(rIV=my_data_sim$r_s[use],nval=my_data_sim$n[use],df1=my_data_sim$df1[use]))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes_sim<-cbind(nRes_sim,an$best$Nullmax)
  kRes_sim<-cbind(kRes_sim,an$best$Kmax)
  sRes_sim<-cbind(sRes_sim,an$best$Smax)
  showAnalysis(an,paste0("sim:n=",round(nmaxs[i])," "))
  
  w1ts<-ExpSamplingCDF(atanh(p2r(0.05,n)),an$best$Kmax,1/sqrt(n-3))
  w0ts<-alpha
  wts<-w0ts*an$best$Nullmax+w1ts*(1-an$best$Nullmax)
  
  w0s<-c(w0s,w0ts)  
  w1s<-c(w1s,w1ts)  
  ws<-c(ws,wts)  
  fdrs<-c(fdrs,w0ts*an$best$Nullmax/wts)
  
}

pts<-data.frame(n=nn,w=w,fdr=fdr,ws=ws,fdrs=fdrs)

doublePlot(nn,'n',rbind(w,ws),'w',y2=rbind(fdr,fdrs),ylb2='fdr',
           xlog=TRUE,
           legendLabels=c("real","sim"),legendX='right')
  
#######################################