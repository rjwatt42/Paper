#########################################
source("SetUp.R")

#########################################
source("OriginalData.R")

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
# 2D histogram of rs vs n

df<-data.frame(x=abs(my_data$r_s),y=my_data$n)
use<-df$y<=250 & df$x<3
df<-df[use,]

g<-ggplot(df) + geom_bin_2d(aes(x=x,y=y),bins=51,show.legend=FALSE,drop=FALSE) 

ns<-seq(5,250)
rs<-tanh(qnorm(1-0.05/2,0,1/sqrt(ns-3)))
ds<-data.frame(x=rs,y=ns)
g<-g+geom_line(data=ds,aes(x=x,y=y),colour="red")

rs<-tanh(qnorm(1-0.005/2,0,1/sqrt(ns-3)))
ds<-data.frame(x=rs,y=ns)
g<-g+geom_line(data=ds,aes(x=x,y=y),colour="orange")

g<-g+ylab("n")+xlab(expression(z[s]))
g<-g + scale_fill_gradient(low="white",high=maincolours$windowC)
g + plotTheme + theme(panel.background=element_rect(fill="white", colour="black"),) 

##########################################
# show distribution of year

years<-unique(my_data$year)
g<-ggplot()
df<-data.frame(d=my_data$year)
g<-g+geom_histogram(data=df,aes(x=d),binwidth=2,color="black",fill="yellow")
g<-g + xlab("year") + plotTheme
g

##########################################
# show distribution of n

maxn<-250


# h<-hist(my_data$n[use],breaks=seq(0,250,5),warn.unused=FALSE,main="",xlab="n",plot=FALSE)
use<-my_data$n<=maxn & my_data$n>=10

e<-egamma(my_data$n[use])
y<-dgamma(seq(10,maxn),e$parameters["shape"],1/e$parameters["scale"])*e$sample.size

g<-ggplot()
de<-data.frame(x=seq(10,maxn),y=y)
df<-data.frame(d=my_data$n[use])
g<-g+geom_histogram(data=df,aes(x=d),binwidth=2,color="white",fill="white")
# g<-g+geom_line(data=de,aes(x=x,y=y*8),color="red")
g<-g + xlab("n") + plotTheme
g


##########################################
# all the results

metaAnal<-list(meta_fixedAnal="random",meta_pdf="All",meta_psigAnal=TRUE,meta_nullAnal=TRUE,append=FALSE)
metaData<-list(result=list(rIV=my_data$r_s,nval=my_data$n))

an<-runMetaAnalysis(metaAnal,metaData)
showAnalysis(an,"All")

# CI 
drawAnalysis(an,metaData)

#########################################
# test Gamma bby comparing shape=1 with exponential
metaAnal<-list(meta_fixedAnal="random",meta_pdf="Gamma",gamma_shape=1,meta_psigAnal=TRUE,meta_nullAnal=TRUE,append=FALSE)
metaData<-list(result=list(rIV=my_data$r_s[1:100],nval=my_data$n[1:100]))

an<-runMetaAnalysis(metaAnal,metaData)

showAnalysis(an,"Gamma")

metaAnal<-list(meta_fixedAnal="random",meta_pdf="Exp",gamma_shape=1,meta_psigAnal=TRUE,meta_nullAnal=TRUE,append=FALSE)
metaData<-list(result=list(rIV=my_data$r_s[1:100],nval=my_data$n[1:100]))

an<-runMetaAnalysis(metaAnal,metaData)

showAnalysis(an,"Exp ")
#########################################
# run Gamma
# shape<1 gives infinity at z=0
shapes<-2^seq(0,4,length.out=7)

resultK<-c()
resultNull<-c()
resultS<-c()
for (si in 1:length(shapes)) {
  metaAnal<-list(meta_fixedAnal="random",meta_pdf="Gamma",gamma_shape=shapes[si],meta_psigAnal=TRUE,meta_nullAnal=TRUE,append=FALSE)
  metaData<-list(result=list(rIV=my_data$r_s,nval=my_data$n))
  
  an<-runMetaAnalysis(metaAnal,metaData)
  resultK[si]<-an$bestK
  resultNull[si]<-an$bestNull
  resultS[si]<-an$bestS
  showAnalysis(an,paste0("Gamma(",format(shapes[si],digits=2),")"))
}

doublePlot(shapes,"k",resultK,expression(lambda),resultNull,expression(p[null]),resultS/1000,expression(log(lk)))
# doublePlot(log10(shapes),"log(Shape)",resultK,"k",resultNull,"p(null)",resultS,"S")

#########################################
# results split by unique or not

metaAnal<-list(meta_fixedAnal="random",meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE,append=FALSE)

conditions=c("Unique","Duplicated")

a<-my_data[c("Source","n")]
use<-duplicated(a)

nRes<-c()
kRes<-c()
sRes<-c()
for (i in 1:length(conditions)) {
  if (i==1) {
    use<-duplicated(a)
  } else {
    use<-!duplicated(a)
  }
  metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$n[use]))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$bestNull)
  kRes<-cbind(kRes,an$bestK)
  sRes<-cbind(sRes,an$bestS)
  showAnalysis(an,paste0("Condition=",conditions[i]))
}

doublePlot(c("Unique","Duplicated"),"",kRes,"k",nRes,"p(null)")

#########################################
# number of siblings

siblings<-seq(1,10)

nRes<-c()
kRes<-c()
sRes<-c()
for (i in 1:length(siblings)) {
  use<-(my_data$Siblings==siblings[i])
  metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$df2[use]+2))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$bestNull)
  kRes<-cbind(kRes,an$bestK)
  sRes<-cbind(sRes,an$bestS)
  showAnalysis(an,paste0("No siblings=",siblings[i]))
}
doublePlot(siblings,"siblings",kRes,"k",nRes,"p(null)")

#########################################
# number of studies

studies<-seq(1,20)

nRes<-c()
kRes<-c()
sRes<-c()
for (i in 1:length(studies)) {
  use<-(my_data$Siblings==1 & my_data$Studies==i)
  metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$df2[use]+2))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$bestNull)
  kRes<-cbind(kRes,an$bestK)
  sRes<-cbind(sRes,an$bestS)
  showAnalysis(an,paste0("No unique studies=",studies[i]))
}
doublePlot(studies,"studies",kRes,"k",nRes,"p(null)")

#########################################
# essentially a check for power analysis
# looks at different ranges of n

nmaxs<-10^seq(log10(20),log10(1000),length.out=6)

nRes<-c()
kRes<-c()
sRes<-c()
for (i in 1:length(nmaxs)) {
  if (i==1) {
    use<-(my_data$df2+2<=nmaxs[i])
  } else {
    use<-(my_data$df2+2>nmaxs[i-1] & my_data$df2+2<=nmaxs[i])
  }
  metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$df2[use]+2))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$bestNull)
  kRes<-cbind(kRes,an$bestK)
  sRes<-cbind(sRes,an$bestS)
  showAnalysis(an,paste0("n=",round(nmaxs[i])))
}
doublePlot(nmaxs,expression(n[max]),kRes,"k",nRes,"p(null)")

#########################################
# split by year
metaAnal<-list(meta_fixedAnal="random",meta_pdf="All",meta_psigAnal=TRUE,meta_nullAnal=TRUE,append=FALSE)
metaData<-list(result=list(rIV=my_data$r_s,nval=my_data$df2+2))

years<-seq(1985,2015,5)

nRes<-c()
kRes<-c()
sRes<-c()
for (i in 1:(length(years)-1)) {
  use<-(my_data$year>=years[i]) & (my_data$year<years[i+1])
  metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$n))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$bestNull)
  kRes<-cbind(kRes,an$bestK)
  sRes<-cbind(sRes,an$bestS)
  showAnalysis(an,paste0("year=",years[i]))
}
doublePlot(years,"year",kRes,"k",nRes,"p(null)")

#########################################
# split by journal

journals<-unique(my_data$journal)

nRes<-c()
kRes<-c()
sRes<-c()
for (i in 1:(length(journals))) {
  use<-(my_data$journal==journals[i]) 
  metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$df2[use]+2))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$bestNull)
  kRes<-cbind(kRes,an$bestK)
  sRes<-cbind(sRes,an$bestS)
  showAnalysis(an,journals[i])
}

doublePlot(journals,"",kRes,"k",nRes,"p(null)")

#########################################
# split by impact factor

factors<-seq(0,1,0.2)

nRes<-c()
kRes<-c()
sRes<-c()
for (i in 1:(length(factors)-1)) {
  use<-(my_data$APAfactor>=factors[i] & my_data$APAfactor<factors[i+1]) 
  metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$n[use]))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$bestNull)
  kRes<-cbind(kRes,an$bestK)
  sRes<-cbind(sRes,an$bestS)
  showAnalysis(an,factors[i])
}
doublePlot(factors,"Impact",kRes,"k",nRes,"p(null)")


#########################################
# split by test type

statistics<-unique(my_data$Statistic)

nRes<-c()
kRes<-c()
sRes<-c()
for (i in 1:(length(statistics))) {
  use<-(my_data$Statistic==statistics[i]) 
  metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$df2[use]+2))
  an<-runMetaAnalysis(metaAnal,metaData)
  nRes<-cbind(nRes,an$bestNull)
  kRes<-cbind(kRes,an$bestK)
  sRes<-cbind(sRes,an$bestS)
  showAnalysis(an,statistics[i])
}
doublePlot(statistics,"Test",kRes,"k",nRes,"p(null)")

