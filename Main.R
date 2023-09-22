#########################################
source("SetUp.R")
alpha<-0.05
nMax<-300
nMin<-5

#########################################
source("OriginalData.R")
original_my_data<-my_data

# NB - throws loads of warnings because of some unexpected contents in excel file
# may also throw 2 errors if variable my_data already exists
# ignore both

#########################################
# find siblings etc

papers<-unique(my_data$Source)
my_data$Siblings<-my_data$Source*0
my_data$Studies<-my_data$Source*0
for (i in 1:length(papers)) {
  use<-my_data$Source==papers[i]
  my_data$Studies[use]<-sum(use)
  ns<-unique(my_data$n[use])
  for (j in 1:length(ns)) {
    usen<-use & (my_data$n==ns[j])
    my_data$Siblings[usen]<-sum(usen)
  }
}

#########################################
# 2D histogram of zs vs n

df<-data.frame(x=my_data$z_s,y=my_data$n)
use<-my_data$n<=nMax & my_data$z_s<3
df<-df[use,]

g<-ggplot(df) + geom_bin_2d(aes(x=x,y=y),bins=51,show.legend=FALSE,drop=FALSE) 

ns<-seq(5,nMax)
zs<-qnorm(1-alpha/2,0,1/sqrt(ns-3))
ds<-data.frame(x=zs,y=ns)
g<-g+geom_line(data=ds,aes(x=x,y=y),colour="red")

zs<-qnorm(1-alpha/200,0,1/sqrt(ns-3))
ds<-data.frame(x=zs,y=ns)
g<-g+geom_line(data=ds,aes(x=x,y=y),colour="orange")

g<-g+ylab("n")+xlab(expression(z[s]))
g<-g + scale_x_continuous(limits=c(0,1.5))
g<-g + scale_fill_gradient(low="white",high=maincolours$windowC)
g + plotTheme + theme(panel.background=element_rect(fill="white", colour="black"),) 

##########################################
# show distribution of year

years<-unique(my_data$year)
g<-ggplot()
df<-data.frame(d=my_data$year)
g<-g+geom_histogram(data=df,aes(x=d),binwidth=1,color="black",fill="white")
g<-g + xlab("year") + ylab("no. outputs") + plotTheme
g

##########################################
# show distribution of n

use<-my_data$n<=nMax & my_data$n>=nMin

g<-ggplot()
df<-data.frame(d=my_data$n[use])
g<-g+geom_histogram(data=df,aes(x=d),binwidth=5,color="white",fill="white")
g<-g + xlab("n") + plotTheme
g


##########################################
# show distribution of zs

use<-my_data$z_s<=1.5

g<-ggplot()
df<-data.frame(d=abs(my_data$z_s[use]))
g<-g+geom_histogram(data=df,aes(x=d),binwidth=0.02,color="white",fill="white")
g<-g + xlab(bquote(z[s])) + plotTheme
g


##########################################
# analysis comparing all models

metaAnal<-list(meta_pdf="All",meta_psigAnal=TRUE,meta_nullAnal=TRUE)
metaData<-list(result=list(rIV=my_data$r_s,nval=my_data$n,df1=my_data$df1))

anMain<-runMetaAnalysis(metaAnal,metaData)

showAnalysis(anMain,"All")
drawAnalysis(anMain,metaData)

#########################################
# test Gamma by comparing shape=1 with exponential

metaData<-list(result=list(rIV=my_data$r_s[1:100],nval=my_data$n[1:100],df1=my_data$df1[1:100]))
metaAnal<-list(meta_pdf=c("Gamma","Exp"),shape=1,meta_psigAnal=TRUE,meta_nullAnal=TRUE)

an0<-runMetaAnalysis(metaAnal,metaData)

showAnalysis(an0)
drawAnalysis(an0,metaData)

#########################################
# run Gamma (very slow as convolution done numerically)
# shape<1 gives infinity at z=0

shapes<-2^seq(0,4,length.out=7)

resultK<-c()
resultNull<-c()
resultS<-c()
for (si in 1:length(shapes)) {
  metaAnal<-list(meta_pdf="Gamma",shape=shapes[si],meta_psigAnal=TRUE,meta_nullAnal=TRUE)
  metaData<-list(result=list(rIV=my_data$r_s,nval=my_data$n,df1=my_data$df1))
  
  an<-runMetaAnalysis(metaAnal,metaData)
  resultK[si]<-an$best$Kmax
  resultNull[si]<-an$best$Nullmax
  resultS[si]<-an$best$Smax
  showAnalysis(an,paste0("Gamma(",format(shapes[si],digits=2),")"))
}

doublePlot(shapes,"k",resultK,NULL,Llabel,resultNull,NULL,Plabel,resultS,Slabel)

#########################################
# number of siblings

siblings<-seq(1,10)
metaAnal<-list(meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE)

nRes<-c()
nCI<-c()
kRes<-c()
kCI<-c()
sRes<-c()
for (i in 1:length(siblings)) {
  use<-(my_data$Siblings==siblings[i])
  metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$n[use],df1=my_data$df1[use]))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$best$Nullmax)
  kRes<-cbind(kRes,an$best$Kmax)
  nCI<-cbind(nCI,matrix(an$best$nullCI,ncol=1))
  kCI<-cbind(kCI,matrix(an$best$kCI,ncol=1))
  sRes<-cbind(sRes,an$best$Smax)
  showAnalysis(an,paste0("No siblings=",siblings[i]))
}
doublePlot(siblings,"siblings",kRes,kCI,Llabel,nRes,nCI,Plabel)

#########################################
# number of studies per paper

studies<-seq(1,21,2)
metaAnal<-list(meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE)

nRes<-c()
nCI<-c()
kRes<-c()
kCI<-c()
sRes<-c()
for (i in 1:(length(studies)-1)) {
  use<-(my_data$Siblings==1 & my_data$Studies>=studies[i] & my_data$Studies<studies[i+1])
  metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$n[use],df1=my_data$df1[use]))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$best$Nullmax)
  kRes<-cbind(kRes,an$best$Kmax)
  nCI<-cbind(nCI,matrix(an$best$nullCI,ncol=1))
  kCI<-cbind(kCI,matrix(an$best$kCI,ncol=1))
  sRes<-cbind(sRes,an$best$Smax)
  showAnalysis(an,paste0("No unique studies=",studies[i]))
}
doublePlot((studies[1:10]+studies[2:11])/2,"studies",
           kRes,kCI,Llabel,nRes,nCI,Plabel)

#########################################
# split by n
# essentially a check for power analysis

nmaxs<-round(10^seq(log10(20),log10(nMax),length.out=11))

my_data_sim<-simData(anMain)

metaAnal<-list(meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE)

nRes<-c()
nCI<-c()
kRes<-c()
kCI<-c()
sRes<-c()

nRes_sim<-c()
nCI_sim<-c()
kRes_sim<-c()
kCI_sim<-c()
sRes_sim<-c()

wRes<-c()
fdrRes<-c()

wRes_sim<-c()
fdrRes_sim<-c()

for (i in 2:length(nmaxs)) {
  useN<-(my_data$n>nmaxs[i-1] & my_data$n<=nmaxs[i])
  nvals<-my_data$n[useN]
  
  metaData<-list(result=list(rIV=my_data$r_s[useN],nval=my_data$n[useN],df1=my_data$df1[useN]))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$best$Nullmax)
  kRes<-cbind(kRes,an$best$Kmax)
  nCI<-cbind(nCI,matrix(an$best$nullCI,ncol=1))
  kCI<-cbind(kCI,matrix(an$best$kCI,ncol=1))
  sRes<-cbind(sRes,an$best$Smax)
  showAnalysis(an,paste0("n=",round(nmaxs[i])," "))
  
  wAll<-0
  wF<-0
  for (ni in min(nvals):max(nvals)) {
    weight<-sum(nvals==ni)/length(nvals)
    zcrit<-qnorm(1-alpha/2,0,1/sqrt(ni-3))
    w1<-ExpSamplingCDF(zcrit,an$best$Kmax,1/sqrt(ni-3))
    w0<-alpha
    w<-w0*an$best$Nullmax + w1*(1-an$best$Nullmax)
    wAll<-wAll+w*weight
    wF<-wF+w0/w*weight
  }
  wRes<-cbind(wRes,wAll) 
  fdrRes<-cbind(fdrRes,wF) 
  
  # simulated data
  use<-(my_data_sim$n>nmaxs[i-1] & my_data_sim$n<=nmaxs[i])
  
  metaData<-list(result=list(rIV=my_data_sim$r_s[use],nval=my_data_sim$n[use],df1=my_data_sim$df1[use]))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes_sim<-cbind(nRes_sim,an$best$Nullmax)
  kRes_sim<-cbind(kRes_sim,an$best$Kmax)
  nCI_sim<-cbind(nCI_sim,matrix(an$best$nullCI,ncol=1))
  kCI_sim<-cbind(kCI_sim,matrix(an$best$kCI,ncol=1))
  sRes_sim<-cbind(sRes_sim,an$best$Smax)
  showAnalysis(an,paste0("sim:n=",round(nmaxs[i])," "))
  
  wAll<-0
  wF<-0
  for (ni in min(nvals):max(nvals)) {
    weight<-sum(nvals==ni)/length(nvals)
    zcrit<-qnorm(1-alpha/2,0,1/sqrt(ni-3))
    w1<-ExpSamplingCDF(zcrit,an$best$Kmax,1/sqrt(ni-3))
    w0<-alpha
    w<-w0*an$best$Nullmax + w1*(1-an$best$Nullmax)
    wAll<-wAll+w*weight
    wF<-wF+w0/w*weight
  }
  wRes_sim<-cbind(wRes_sim,wAll) 
  fdrRes_sim<-cbind(fdrRes_sim,wF)
  
}

# save these results for later
nResN<-nRes
kResN<-kRes


# doublePlot(nmaxs[2:length(nmaxs)],expression(n[max]),rbind(kRes,kRes_sim),Llabel,rbind(nRes,nRes_sim),Plabel,xlog=TRUE)
doublePlot(nmaxs[2:length(nmaxs)],expression(n[max]),rbind(kRes,kRes_sim),rbind(kCI,kCI_sim),Llabel,rbind(nRes,nRes_sim),rbind(nCI,nCI_sim),Plabel,xlog=TRUE)
doublePlot(nmaxs[2:length(nmaxs)],expression(n[max]),rbind(wRes,wRes_sim),NULL,"w",rbind(fdrRes,fdrRes_sim),NULL,"FDR",xlog=TRUE)

#########################################
# split by year

years<-seq(1985,2015,2)
metaAnal<-list(meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE)

nRes<-c()
nCI<-c()
kRes<-c()
kCI<-c()
sRes<-c()

nRes_expect<-c()
kRes_expect<-c()

for (i in 1:(length(years)-1)) {
  useY<-(my_data$year>=years[i]) & (my_data$year<years[i+1]) & (my_data$n<=300)
  metaData<-list(result=list(rIV=my_data$r_s[useY],nval=my_data$n[useY],df1=my_data$df1[useY]))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$best$Nullmax)
  kRes<-cbind(kRes,an$best$Kmax)
  nCI<-cbind(nCI,matrix(an$best$nullCI,ncol=1))
  kCI<-cbind(kCI,matrix(an$best$kCI,ncol=1))
  sRes<-cbind(sRes,an$best$Smax)
  showAnalysis(an,paste0("year=",years[i]))


  ne<-0
  ke<-0
  nc<-0
  for (ni in 2:length(nmaxs)) {
    useN<-my_data$n>nmaxs[ni-1] & my_data$n<=nmaxs[ni]
    count<-sum(useN & useY)
    ne<-ne+count*nResN[ni-1]
    ke<-ke+count*kResN[ni-1]
    nc<-nc+count
  }
  nRes_expect<-c(nRes_expect,ne/nc)
  kRes_expect<-c(kRes_expect,ke/nc)
}
# doublePlot(years[1:(length(years)-1)],"year",rbind(kRes),kCI,Llabel,
#            rbind(nRes),nCI,Plabel)
doublePlot(years[1:(length(years)-1)],"year",rbind(kRes,kRes_expect),rbind(kCI),Llabel,
           rbind(nRes,nRes_expect),rbind(nCI),Plabel,legendLabels=c("actual","expected"))

#########################################
# split by journal

journals<-unique(my_data$journal)
metaAnal<-list(meta_pdf="All",meta_psigAnal=TRUE,meta_nullAnal=TRUE)

nRes<-c()
nCI<-c()
kRes<-c()
kCI<-c()
sRes<-c()
nmean<-c()
wRes<-c()
fdrRes<-c()

nRes_exp<-c()
kRes_exp<-c()

for (i in 1:(length(journals))) {
  useJ<-(my_data$journal==journals[i]) 
  metaData<-list(result=list(rIV=my_data$r_s[useJ],nval=my_data$n[useJ],df1=my_data$df1[useJ]))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$best$Nullmax)
  kRes<-cbind(kRes,an$best$Kmax)
  nCI<-cbind(nCI,matrix(an$best$nullCI,ncol=1))
  kCI<-cbind(kCI,matrix(an$best$kCI,ncol=1))
  sRes<-cbind(sRes,an$best$Smax)
  nmean<-c(nmean,median(my_data$n[useJ]))
  showAnalysis(an,paste0(journals[i],"(n=",nmean[i],")"),"Best")
  
  ne<-0
  ke<-0
  for (ni in 2:length(nmaxs)) {
    useN<-my_data$n>nmaxs[ni-1] & my_data$n<=nmaxs[ni]
    count<-sum(useN & useJ)
    ne<-ne+count*nResN[ni-1]
    ke<-ke+count*kResN[ni-1]
  }
  nRes_exp<-c(nRes_exp,ne/sum(useJ))
  kRes_exp<-c(kRes_exp,ke/sum(useJ))
  
  # probability of significant given it is a non-null
  w<-1-ExpSamplingCDF(qnorm(1-alpha/2)/sqrt(my_data$n[use]-3),an$best$Kmax,1/sqrt(my_data$n[use]-3))
  w<-mean(w)
  # probability of significant non-null or null
  wAll<-w*(1-an$best$Nullmax)+alpha*an$best$Nullmax
  wRes<-c(wRes,wAll)
  # fdr
  fdr<-alpha*an$best$Nullmax/wAll
  fdrRes<-c(fdrRes,fdr)
}

useOrder<-order(nmean)
doublePlot(journals[useOrder],"",rbind(kRes,kRes_exp)[,useOrder],NULL,Llabel,rbind(nRes,nRes_exp)[,useOrder],NULL,Plabel)
doublePlot(journals[useOrder],"",wRes[useOrder],NULL,"w",fdrRes[useOrder],NULL,"FDR")

#########################################
# split by impact factor

factors<-seq(0,1,0.2)
metaAnal<-list(meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE)

nRes<-c()
nCI<-c()
kRes<-c()
kCI<-c()
sRes<-c()
for (i in 1:(length(factors)-1)) {
  use<-(my_data$APAfactor>=factors[i] & my_data$APAfactor<factors[i+1]) 
  metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$n[use],df1=my_data$df1[use]))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$best$Nullmax)
  kRes<-cbind(kRes,an$best$Kmax)
  nCI<-cbind(nCI,matrix(an$best$nullCI,ncol=1))
  kCI<-cbind(kCI,matrix(an$best$kCI,ncol=1))
  sRes<-cbind(sRes,an$best$Smax)
  showAnalysis(an,factors[i])
}
doublePlot(factors,"Impact",kRes,kCI,Llabel,nRes,nCI,Plabel)


#########################################
# split by test type

statistics<-unique(my_data$Statistic)
statistics<-c("t","r")

metaAnal<-list(meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE)

nRes<-c()
nCI<-c()
kRes<-c()
kCI<-c()
sRes<-c()
for (i in 1:(length(statistics))) {
  use<-(my_data$Statistic==statistics[i]) 
  metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$n[use],df1=my_data$df1[use]))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$best$Nullmax)
  kRes<-cbind(kRes,an$best$Kmax)
  nCI<-cbind(nCI,matrix(an$best$nullCI,ncol=1))
  kCI<-cbind(kCI,matrix(an$best$kCI,ncol=1))
  sRes<-cbind(sRes,an$best$Smax)
  showAnalysis(an,statistics[i])
}
doublePlot(statistics,"Test",kRes,kCI,Llabel,nRes,nCI,Plabel)

#########################################
# split by 1-tail vs 2
# 
cases<-c("1 @ 0.1","1 @ 0.05","2 @ 0.05")
metaAnal<-list(meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE)

tails1<-my_data$OneTail | my_data$OneTailedInTxt

nRes<-c()
nCI<-c()
kRes<-c()
kCI<-c()
sRes<-c()
for (i in 1:(length(cases))) {
  if (i==1) {
    alpha<-0.1
    r<-c(my_1t_data$r_s,my_data$r_s[tails1])
    n<-c(my_1t_data$n,my_data$n[tails1])
    df1<-c(my_1t_data$df1,my_data$df1[tails1])

    metaData<-list(result=list(rIV=r,nval=n,df1=df1))
  }
  if (i==2) {
    alpha<-0.05
    metaData<-list(result=list(rIV=my_data$r_s[tails1],nval=my_data$n[tails1],df1=my_data$df1[tails1]))
  }
  if (i==3) {
    alpha<-0.05
    metaData<-list(result=list(rIV=my_data$r_s[!tails1],nval=my_data$n[!tails1],df1=my_data$df1[!tails1]))
  }
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$best$Nullmax)
  kRes<-cbind(kRes,an$best$Kmax)
  nCI<-cbind(nCI,matrix(an$best$nullCI,ncol=1))
  kCI<-cbind(kCI,matrix(an$best$kCI,ncol=1))
  sRes<-cbind(sRes,an$best$Smax)
  showAnalysis(an,cases[i])
}
doublePlot(cases,"Tails",kRes,kCI,Llabel,nRes,nCI,Plabel)

#########################################
# changing alpha

alpha<-0.05
alphas<-0.05*10^seq(-2,0,length.out=9)

my_data_sim<-simData(anMain)

nRes<-c()
kRes<-c()
nResCI<-c()
kResCI<-c()
sRes<-c()

nRes_expect<-c()
kRes_expect<-c()

nRes_sim<-c()
kRes_sim<-c()
nRes_simCI<-c()
kRes_simCI<-c()
sRes_sim<-c()

wRes<-c()
fdrRes<-c()

wRes_sim<-c()
fdrRes_sim<-c()
use_n<-my_data$n<=800 & my_data$n>=10
use_n_sim<-my_data_sim$n<=800 & my_data_sim$n>=10

metaAnal<-list(meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE)

for (i in 1:length(alphas)) {
  alpha<-alphas[i]
  
  useP<-my_data$p<alpha & use_n
  metaData<-list(result=list(rIV=my_data$r_s[useP],nval=my_data$n[useP],df1=my_data$df1[useP]))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$best$Nullmax)
  kRes<-cbind(kRes,an$best$Kmax)
  nResCI<-cbind(nResCI,an$best$nullCI)
  kResCI<-cbind(kResCI,an$best$kCI)
  sRes<-cbind(sRes,an$best$Smax)
  showAnalysis(an,paste0("alpha=",alphas[i]))
  
  ne<-0
  ke<-0
  for (ni in 2:length(nmaxs)) {
    useN<-my_data$n>nmaxs[ni-1] & my_data$n<=nmaxs[ni]
    count<-sum(useN & useP)
    ne<-ne+count*nResN[ni-1]
    ke<-ke+count*kResN[ni-1]
  }
  nRes_expect<-c(nRes_expect,ne/sum(useP))
  kRes_expect<-c(kRes_expect,ke/sum(useP))
  
  # probability of significant given it is a non-null
  wPlus<-1-ExpSamplingCDF(qnorm(1-alpha/2)/sqrt(my_data$n[useP]-3),an$best$Kmax,1/sqrt(my_data$n[useP]-3))
  wPlus<-mean(wPlus)
  wNull<-alpha
  # probability of significant non-null or null
  wAll<-wPlus*(1-an$best$Nullmax)+wNull*an$best$Nullmax
  wRes<-c(wRes,wAll)
  # fdr
  fdr<-wNull*an$best$Nullmax/wAll
  fdrRes<-c(fdrRes,fdr)
  
  useP_sim<-my_data_sim$p<alpha & use_n
  metaData<-list(result=list(rIV=my_data_sim$r_s[useP_sim],nval=my_data_sim$n[useP_sim],df1=my_data_sim$df1[useP_sim]))
  an_sim<-runMetaAnalysis(metaAnal,metaData)
  nRes_sim<-cbind(nRes_sim,an_sim$best$Nullmax)
  kRes_sim<-cbind(kRes_sim,an_sim$best$Kmax)
  nRes_simCI<-cbind(nRes_simCI,an_sim$best$nullCI)
  kRes_simCI<-cbind(kRes_simCI,an_sim$best$kCI)
  sRes_sim<-cbind(sRes_sim,an_sim$best$Smax)
  showAnalysis(an_sim,paste0("sim:alpha=",alphas[i]))
  
  # probability of significant given it is a non-null
  wPlus<-1-ExpSamplingCDF(qnorm(1-alpha/2)/sqrt(my_data_sim$n[useP_sim]-3),an_sim$best$Kmax,1/sqrt(my_data_sim$n[useP_sim]-3))
  wPlus<-mean(wPlus)
  wNull<-alpha
  # probability of significant non-null or null
  wAll<-wPlus*(1-an_sim$best$Nullmax)+wNull*an_sim$best$Nullmax
  wRes_sim<-c(wRes_sim,wAll)
  # fdr
  fdr<-wNull*an_sim$best$Nullmax/wAll
  fdrRes_sim<-c(fdrRes_sim,fdr)
  
}
alpha<-0.05

doublePlot(alphas,"alpha",rbind(kRes,kRes_expect),rbind(kResCI,kRes_simCI),Llabel,
           rbind(nRes,nRes_expect),rbind(nResCI,nRes_simCI),Plabel,xlog=TRUE)

doublePlot(alphas,"alpha",rbind(wRes,wRes_sim),NULL,"w",rbind(fdrRes,fdrRes_sim),NULL,"FDR",xlog=TRUE)

#########################################
# siblings by journal

journals<-unique(my_data$journal)
siblings<-seq(1,10)
metaAnal<-list(meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE)

nn<-array(0,dim=c(2,length(journals)))
tnRes<-c()
tkRes<-c()
rnRes<-c()
rkRes<-c()
nmean<-c()
psum<-c()
for (i in 1:(length(journals))) {
  usej<-(my_data$journal==journals[i]) 
  nmean<-c(nmean,median(my_data$n[usej]))
  psum<-c(psum,sum(usej))
  
  uses<-(my_data$Siblings<=3) 
  use<-usej & uses
  nn[1,i]<-sum(use)
  metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$n[use],df1=my_data$df1[use]))
  an<-runMetaAnalysis(metaAnal,metaData)
  showAnalysis(an,paste0("journal=",journals[i]))
  tnRes<-cbind(tnRes,an$best$Nullmax)
  tkRes<-cbind(tkRes,an$best$Kmax)
  
  uses<-(my_data$Siblings>=7) 
  use<-usej & uses
  nn[2,i]<-sum(use)
  metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$n[use],df1=my_data$df1[use]))
  an<-runMetaAnalysis(metaAnal,metaData)
  showAnalysis(an,paste0("journal=",journals[i]))
  rnRes<-cbind(rnRes,an$best$Nullmax)
  rkRes<-cbind(rkRes,an$best$Kmax)
}

dispOrder<-order(nmean)
doublePlot(journals[dispOrder],"",rbind(tkRes[,dispOrder],rkRes[,dispOrder]),NULL,Llabel,
                                  rbind(tnRes[,dispOrder],rnRes[,dispOrder]),NULL,Plabel,
           legendLabels=c("siblings<4","siblings>6"))

#########################################
# test-type split by journal

journals<-unique(my_data$journal)
statistics<-c("t","r")
metaAnal<-list(meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE)

nn<-array(0,dim=c(2,length(journals)))
tnRes<-c()
tkRes<-c()
rnRes<-c()
rkRes<-c()
tnResCI<-c()
tkResCI<-c()
rnResCI<-c()
rkResCI<-c()
nmean<-c()
psum<-c()
for (i in 1:(length(journals))) {
  usej<-(my_data$journal==journals[i]) 
  nmean<-c(nmean,median(my_data$n[usej]))
  psum<-c(psum,sum(usej))
  
    uses<-(my_data$Statistic=="t") 
    use<-usej & uses
    nn[1,i]<-sum(use)
    metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$n[use],df1=my_data$df1[use]))
    an<-runMetaAnalysis(metaAnal,metaData)
    showAnalysis(an,paste0("journal=",journals[i]))
    tnRes<-cbind(tnRes,an$best$Nullmax)
    tkRes<-cbind(tkRes,an$best$Kmax)
    tnResCI<-cbind(tnResCI,an$best$nullCI)
    tkResCI<-cbind(tkResCI,an$best$kCI)
    
    uses<-(my_data$Statistic=="r") 
    use<-usej & uses
    nn[2,i]<-sum(use)
    metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$n[use],df1=my_data$df1[use]))
    an<-runMetaAnalysis(metaAnal,metaData)
    showAnalysis(an,paste0("journal=",journals[i]))
    rnRes<-cbind(rnRes,an$best$Nullmax)
    rkRes<-cbind(rkRes,an$best$Kmax)
    rnResCI<-cbind(rnResCI,an$best$nullCI)
    rkResCI<-cbind(rkResCI,an$best$kCI)
}

dispOrder<-order(nmean)
doublePlot(journals[dispOrder],"",
           rbind(tkRes,rkRes)[,dispOrder],rbind(tkResCI,rkResCI)[,dispOrder],Llabel,
           rbind(tnRes,rnRes)[,dispOrder],rbind(tnResCI,rnResCI)[,dispOrder],Plabel,
           legendLabels=c("t","r"))

#######################################
