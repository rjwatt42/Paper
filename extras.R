#########################################
# run GenExp (very slow as convolution done numerically)
# 
shapes<-2^seq(-5,2,length.out=8)

resultK<-c()
resultNull<-c()
resultS<-c()
for (si in 1:length(shapes)) {
  metaAnal<-list(meta_fixedAnal="random",meta_pdf="GenExp",shape=shapes[si],meta_psigAnal=TRUE,meta_nullAnal=TRUE,append=FALSE)
  metaData<-list(result=list(rIV=my_data$r_s,nval=my_data$n,df1=my_data$df1))
  
  an<-runMetaAnalysis(metaAnal,metaData)
  resultK[si]<-an$best$Kmax
  resultNull[si]<-an$best$Nullmax
  resultS[si]<-an$best$Smax
  showAnalysis(an,paste0("GenExp(",format(shapes[si],digits=2),")"))
}

doublePlot(log2(shapes),"log2(a)",resultK,NULL,Llabel,resultNull,NULL,Plabel,resultS,Slabel)

#########################################
# number of parameters

parameters<-1:8

nRes<-rep(NA,length(parameters))
kRes<-rep(NA,length(parameters))
sRes<-rep(NA,length(parameters))
nS<-rep(NA,length(parameters))
for (i in 1:length(parameters)) {
  use<-which(my_data$df1==parameters[i])
  nS[i]<-length(use)
  metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$n[use],df1=my_data$df1[use]))
  try({
    an<-runMetaAnalysis(metaAnal,metaData)
    nRes[i]<-an$best$Nullmax
    kRes[i]<-an$best$Kmax
    sRes[i]<-an$best$Smax
    showAnalysis(an,paste0("No parameters=",parameters[i]))
  })
}

doublePlot(parameters,"parameters",kRes,NULL,Llabel,nRes,Plabel,xtick=parameters)

#########################################
# look at observed vs expected histograms

my_data_use<-my_data
if (1==2) my_data_use<-simData()

metaAnal<-list(meta_fixedAnal="random",meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE,append=FALSE)
maxZ<-6
maxN<-250

# find the common sample sizes
hn<-hist(my_data_use$n[my_data_use$n<=maxN],breaks=(0:maxN)+0.5,plot=FALSE)
useN<-which(hn$counts>=1000)

g<-ggplot()

oe<-c()
for (i in 1:length(useN)) {
  n<-useN[i]
  use<-my_data_use$n==n
  nStudies<-sum(use)
  
  metaData<-list(result=list(rIV=my_data_use$r_s[use],nval=my_data_use$n[use],df1=my_data_use$df1[use]))
  an<-runMetaAnalysis(metaAnal,metaData)
  
  lambda<-an$best$Kmax
  p_null<-an$best$Nullmax
  nsims<-10
  
  zcrit<-atanh(p2r(alpha,n))
  zbins<-seq(0,1000)*zcrit/32
  # zbins<-zbins[zbins<=maxZ]
  zs<-atanh(abs(my_data_use$r_s[use]))
  zs<-zs[zs<max(zbins)]
  h1<-hist(zs,breaks=zbins,plot=FALSE)
  zmids<-h1$mids
  
  d1<-ExpSamplingPDF(zbins,lambda,1/sqrt(n-3),remove_nonsig = TRUE)
  d1<-d1$pdf/sum(d1$pdf)/d1$sig_pdf
  d0<-SingleSamplingPDF(zbins,0,1/sqrt(n-3),remove_nonsig = TRUE)
  d0<-d0$pdf/sum(d0$pdf)/d0$sig_pdf
  zd<-d1*(1-p_null)+d0*p_null
  zd[zbins<zcrit]<-0
  zd<-(zd[1:(length(zd)-1)]+zd[2:length(zd)])/2
  zd<-zd*nStudies
  
  oe<-rbind(oe,(h1$counts-zd)/nStudies)
  s<-data.frame(z=zmids/zcrit,zm=zd,z1=h1$counts,n=rep(n,length(zmids)))
  
  # g<-g+geom_line(data=s,aes(x=z,y=(z1-zm)/nStudies,color=n),lwd=0.5)
}

s<-data.frame(x=zmids/zcrit,y=colMeans(oe),n=rep(n,length(zmids)))
g<-g+geom_line(data=s,aes(x=x,y=y),color="white",lwd=1)

# g<-g+scale_colour_gradient(low = "red", high = "yellow")
g<-g+scale_x_continuous(limits=c(0,maxZ))
g<-g + xlab(bquote(frac(z[s],z[crit]))) + ylab(bquote(frac((obs-exp),exp))) + plotTheme
g
