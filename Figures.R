##########################################
# example distribution of sample sizes

nv<-seq(5,250)
nd<-dgamma(nv-4,shape=1.4,scale=45)
nd<-nd/sum(nd)

g<-ggplot()
n<-data.frame(x=c(nv,rev(nv)),y=c(nd,nd*0))
g<-g+geom_polygon(data=n,aes(x=x,y=y),fill="white")
g<-g+xlab("n")+ylab("frequency")+scale_y_continuous(breaks=c())

g+plotTheme


##########################################
# probability densities vs lambda

nullval<-0.5
n<-42
zs<-0.25

# choose one of the following two
kvals<-seq(-0.1,0.7,0.1)
kvals<-seq(-1,1,0.01)
prior<-exp(-abs(kvals/0.3))

g<-ggplot()
plotCol<-"white"

S<-dnorm(zs,0,1/sqrt(n-3))*(1-nullval)+dnorm(zs,atanh(kvals),1/sqrt(n-3))*nullval
S<-S*prior
pts<-data.frame(x=kvals,y=S)
g<-g+geom_line(data=pts,aes(x=x,y=y),col=plotCol,lwd=1)
if (length(kvals)<10)
  g<-g+geom_point(data=pts,aes(x=x,y=y),col=plotCol,size=4)

maxS<-max(S)
maxK<-kvals[which.max(S)]
g<-g+geom_vline(xintercept=maxK,color=plotCol)
ls<-data.frame(x=maxK,y=maxS,label=paste0("max=",format(maxK,digits=3)))
g<-g+geom_text(data=ls,aes(x=x,y=y,label=label),color=plotCol,vjust=-0.5,hjust=-0.1)

g<-g+scale_y_continuous(limits=c(0,2))
g<-g+xlab(Llabel)+ylab("probability density")
g+plotTheme

##########################################
# a likelihood version - with and without publication bias

z<-seq(-2,2,0.01)
nullP<-0.45

nonnull<-ExpSamplingPDF(z,lambda=0.3,sigma=1/sqrt(42-3),shape=NA,remove_nonsig=TRUE) 
null<-SingleSamplingPDF(z,0,sigma=1/sqrt(42-3),shape=NA,remove_nonsig=TRUE) 

zcrit<-atanh(p2r(alpha,42))

# nonnull$pdf[abs(z)<zcrit]<-0
# null$pdf[abs(z)<zcrit]<-0

gain<-sum(nonnull$pdf*(1-nullP)+null$pdf*nullP)*0.01

g<-ggplot()
    pts<-data.frame(x=z,y1=nonnull$pdf*(1-nullP),y2=null$pdf*nullP,y3=nonnull$pdf*(1-nullP)+null$pdf*nullP)
    g<-g+geom_line(data=pts,aes(x=x,y=y3/gain),color="white",lwd=2,lty=1)
    g<-g+geom_line(data=pts,aes(x=x,y=y1/gain),color="green",lwd=0.6,lty=1)
    g<-g+geom_line(data=pts,aes(x=x,y=y2/gain),color="red",lwd=0.6,lty=1)

    g<-g+xlab(expression(z[s]))+ylab("Probability Density")
    g<-g+scale_y_continuous(limits=c(0,2))
    g+plotTheme    
      

##########################################
# a likelihood version

sourceDistribution="Single"
analysisDistribution<-c("Single","Gauss","Exp")
# analysisDistribution<-c("Gauss","Exp")
# analysisDistribution<-"Single"
psig<-FALSE
nullval<-0.5

n<-42
zs<-0.25
nsamp<-1

kTrans<-function(k) {k}

g<-ggplot()
rymax<-c()
rymin<-c()
if (nsamp>1 && nsamp<21) {
  for (z in zs) {
    expS<-getLogLikelihood(z,n,df1=1,analysisDistribution,kvals,nullval,psig)
    rymax<-c(rymax,max(expS))
    rymin<-c(rymin,min(expS))
    pts<-data.frame(x=kTrans(kvals),y=expS)
    g<-g+geom_line(data=pts,aes(x=x,y=y),lwd=0.4,lty=1)
  }
}

maxS<- -Inf
for (ai in 1:length(analysisDistribution)) {
  if (analysisDistribution[ai]=="Single") {
    kvals<-seq(-1,1,length.out=101)
  } else {
    kvals<-seq(0.05,1,length.out=101)
  }
  
  expS<-getLogLikelihood(zs,n,df1=1,analysisDistribution[ai],kvals,nullval,psig)
  # expS<-expS/length(zs)
  rymax<-c(rymax,max(expS))
  rymin<-c(rymin,min(expS))
  
  pts<-data.frame(x=kTrans(kvals),y=expS)
  switch(analysisDistribution[ai],
         "Single"={ccol<-"white"},
         "Gauss"={ccol<-"yellow"},
         "Exp"={ccol<-"red"},
  )
  g<-g+geom_line(data=pts,aes(x=x,y=y),col=ccol,lwd=1)
  if (length(kvals)<10)
    g<-g+geom_point(data=pts,aes(x=x,y=y),col=ccol,size=4)
  if(max(expS)>maxS) {
    mcol<-ccol
    maxS<-max(expS)
    kmax<-kvals[which.max(expS)]
    d<-analysisDistribution[ai]
  }
}

g<-g+geom_vline(xintercept=kTrans(kmax),color=mcol)
ls<-data.frame(x=kTrans(kmax),y=maxS,label=paste0("max=",d,"(",format(kmax,digits=3),")"))
g<-g+geom_text(data=ls,aes(x=x,y=y,label=label),color=mcol,vjust=-0.5,hjust=-0.1)

g<-g+xlab(Llabel)+ylab("log(likelihood)")

ylim<-c(min(rymin),max(rymax))
plotRange<-diff(ylim)
g<-g+scale_y_continuous(limits=ylim+c(-1,1)*plotRange/10)
g+plotTheme


##########################################
# plot S as function of lambda

sourceDistribution="Single"
analysisDistribution<-"Single"

sourceDistribution="Exp"
analysisDistribution<-c("Exp","Gauss")
cols<-c("yellow","red")
k<-0.325
null<-0.5
psig<-FALSE

nsamp=2000
switch (sourceDistribution,
        "Exp"={
          rp<-tanh(rexp(nsamp,1/k))
        },
        "Single"={
          rp<-rep(tanh(k),nsamp)
        })
rp<-rp*(runif(nsamp)>null)

n<-175
zs<-rnorm(nsamp,atanh(rp),1/sqrt(n-3))

kvals<-seq(0.05,1,length.out=101)

g<-ggplot()
rymax<-c()
rymin<-c()
if (nsamp<21) {
  for (z in zs) {
    expS<-getLogLikelihood(z,n,df1=1,analysisDistribution[1],atanh(kvals),null,psig)
    rymax<-c(rymax,max(expS))
    rymin<-c(rymin,min(expS))
    pts<-data.frame(x=kvals,y=expS)
    g<-g+geom_line(data=pts,aes(x=x,y=y),lwd=0.4,lty=1)
  }
}
maxD<- -Inf
for (ai in 1:length(analysisDistribution)) {
  expS<-getLogLikelihood(zs,n,df1=1,analysisDistribution[ai],atanh(kvals),null,psig)
  rymax<-c(rymax,max(expS))
  rymin<-c(rymin,min(expS))

  if (max(expS)>maxD) {
    maxD<-max(expS)
    maxDist<-ai
  }
pts<-data.frame(x=kvals,y=expS)
g<-g+geom_line(data=pts,aes(x=x,y=y),col=cols[ai],lwd=1)
}

if (length(analysisDistribution)==1) {
  nmax<-kvals[which.max(expS)]
  g<-g+geom_vline(xintercept=nmax,color="white")
  ls<-data.frame(x=nmax,y=max(expS),label=paste0("max=",format(nmax,digits=3)))
  g<-g+geom_text(data=ls,aes(x=x,y=y,label=label),color="white",vjust=-0.5,hjust=-0.1)
} else {
  nmax<-kvals[which.max(expS)]
  ls<-data.frame(x=nmax,y=maxD,label=paste0("max=",analysisDistribution[maxDist]))
  g<-g+geom_text(data=ls,aes(x=x,y=y,label=label),color=cols[maxDist],vjust=-0.5,hjust=-0.1)
}

g<-g+xlab(Llabel)+ylab("log(likelihood)")
# g<-g+scale_x_continuous(breaks=c(-2,-1,0),labels=c(0.001,0.1,1))
if (length(analysisDistribution)>1) {
  rymin<- -200
}
g<-g+scale_y_continuous(limits=c(min(rymin),max(rymax))+c(-1,1)*(rymax-rymin)/10)
g+plotTheme


##########################################
# 2D density plot of S vs k and null

sourceDistribution="Exp"
analysisDistribution<-"Exp"
psig<-FALSE
null<-0.45
k<-atanh(0.25)

showLegend<-FALSE

nsamp=100
switch (sourceDistribution,
        "Exp"={
          rp<-tanh(rexp(nsamp,1/k))
        },
        "Single"={
          rp<-rep(tanh(k),nsamp)
        })
rp<-rp*(runif(nsamp)>null)

n<-175
zs<-rnorm(nsamp,atanh(rp),1/sqrt(n-3))

kvals<-seq(0.05,1,length.out=51)
plusvals<-seq(0,0.999,length.out=101)

expS<-getLogLikelihood(zs,n,df1=1,analysisDistribution,kvals,1-plusvals,psig)
threshold<-quantile(expS,0.25)
expS[expS < threshold]<-threshold

df<-melt(expS)
g<-ggplot(df)+geom_tile(aes(x=(X2-1)/100*(plusvals[101]-plusvals[1])+plusvals[1],y=(X1-1)*(kvals[51]-kvals[1])/50+kvals[1],fill = value),show.legend=showLegend)
g<-g+stat_contour(data=df,aes(x=(X2-1)/100*(plusvals[101]-plusvals[1])+plusvals[1],y=(X1-1)*(kvals[51]-kvals[1])/50+kvals[1],z = value),breaks=max(expS)-log(10),color="red")

use<-which(expS==max(expS), arr.ind = TRUE)
g<-g+geom_vline(xintercept=plusvals[use[1,2]],color="red")
g<-g+geom_hline(yintercept=kvals[use[1,1]],color="red")

label1<-bquote(.(Llabel)==.(round(100*kvals[use[1,1]])/100))
labelLoc<-data.frame(x=plusvals[use[1,2]],y=kvals[use[1,1]])
g<-g+geom_text(data=labelLoc,aes(x=x,y=y),label=deparse(label1),color="white",vjust=-0.5,hjust=1.1,parse=TRUE)
label2<-bquote(.(Plabel)==.(round(100*plusvals[use[1,2]])/100))
labelLoc<-data.frame(x=plusvals[use[1,2]],y=kvals[use[1,1]])
g<-g+geom_text(data=labelLoc,aes(x=x,y=y),label=deparse(label2),color="white",vjust=-0.5,hjust=-0.1,parse=TRUE)

g<-g+xlab(Plabel)+ylab(Llabel)
g<-g + plotTheme + scale_fill_gradient(low="white",high=maincolours$windowC)
g+theme(panel.background = element_rect(fill="white", colour="black"))

##########################################
# show the various shapes for gamma distribution

shapes=c(1,2,4,8,16)

z<-seq(0,1.5,length.out=501)

g<-ggplot()

for (si in 1:length(shapes)) {
  d<-dgamma(z,shape=shapes[si],scale=0.3/8)
  pts<-data.frame(x=z,y=d,shape=log2(shapes[si]))
  g<-g+geom_line(data=pts,aes(x=x,y=y,col=shape))
}

g<-g+xlab(expression(z[p]))+ylab("Probability Density")
g<-g+scale_color_gradient2(name="log2(k)",low="white", high="#FF2222", 
                           limits = c(0, 4))
g<-g+theme(
  legend.background = element_rect(fill="#666666"),
  legend.key = element_rect(fill="#666666"),
  legend.text=element_text(color="white"),
  legend.title = element_text(colour="white", size=12, face="bold",hjust=0.5,vjust=4),
  legend.position = c(.95, .95),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(0, 0, 0, 0),
  legend.spacing.x = unit(2, "mm"),
  legend.spacing.y = unit(0, "mm")
)

g + plotTheme

##########################################
# show the various shapes for gen exp distribution

shapes=c(0.25,0.5,1,2,4)
scale=0.3

r<-seq(0,0.99,length.out=501)
z<-seq(0,1.5,length.out=501)

g<-ggplot()

for (si in 1:length(shapes)) {
  d<- 1-(1 - exp(-z*shapes[si]/scale))^shapes[si]
  rate<-(sum(z*d)/sum(d))/scale/scale
  d<- 1-(1 - exp(-z*shapes[si]*rate))^shapes[si]
  d<-d/sum(d)
  pts<-data.frame(x=z,y=d,shape=log2(shapes[si]))
  g<-g+geom_line(data=pts,aes(x=x,y=y,col=shape))
  # switch(si, # this is nuts
  #   g<-g+geom_line(data=pts,aes(x=x,y=y,col="0.25")),
  #   g<-g+geom_line(data=pts,aes(x=x,y=y,col="0.5")),
  #   g<-g+geom_line(data=pts,aes(x=x,y=y,col="1")),
  #   g<-g+geom_line(data=pts,aes(x=x,y=y,col="2")),
  #   g<-g+geom_line(data=pts,aes(x=x,y=y,col="4"))
  # )
  print(sum(z*d)/sum(d))
}
g<-g+xlab(expression(z[p]))+ylab("Probability Density")
# g<-g+scale_color_continuous()
g<-g+scale_color_gradient2(name="log2(a)",low="#8888FF", mid="white", high="#FF2222", 
                           limits = c(-2,2))
# g<-g+scale_color_manual(name = "a") #, values = c("0.25" = "red", "0.5" = "orange","1" = "yellow","2" = "green","4" = "blue") )
g<-g+theme(
  legend.background = element_rect(fill="#666666"),
  legend.key = element_rect(fill="#666666"),
  legend.text=element_text(color="white"),
  legend.title = element_text(colour="white", size=12, face="bold",hjust=0.5,vjust=4),
  legend.position = c(.95, .95),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(0, 0, 0, 0),
  legend.spacing.x = unit(2, "mm"),
  legend.spacing.y = unit(0, "mm")
)

g + plotTheme

##########################################

shapes=c(0.25,1)
scale=c(0.341,0.325)
pnull<-c(0.6,0.74)
types<-c("GenExp","Exp")

doCDF=TRUE
# doCDF=FALSE

r<-seq(0,0.99,length.out=5001)
z<-seq(0,15,length.out=5001)

pts<-data.frame(x=z[1:501])

for (i in 1:length(shapes)) {
  d<-GenExpSamplingPDF(z,scale[i],0,shapes[i])
  d<-d/(sum(d)*(z[2]-z[1]))
  if (doCDF) d<-cumsum(d)/sum(d)*(1-pnull[i])+pnull[i]
  pts<-cbind(pts,data.frame(d=d[1:501]))
  names(pts)[i+1]=paste0("y",i)
}

if (doCDF) pts<-rbind(rep(0,ncol(pts)),pts)

g<-ggplot(pts,aes(x=x))
varnames <- names(pts)[2:ncol(pts)]
add_lines <- lapply(varnames, function(i) geom_line(aes_q(y = as.name(i), col = i),lwd=1))
g<-g+add_lines

g<-g+xlab(expression(z[p]))+ylab("Probability Density")
g<-g+scale_color_discrete(name="",labels=types)

g<-g+theme(
  legend.background = element_rect(fill="#666666"),
  legend.key = element_rect(fill="#666666"),
  legend.text=element_text(color="white"),
  legend.title = element_text(colour="white", size=12, face="bold",hjust=0.5,vjust=4),
  legend.position = c(.95, .75),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(0, 0, 0, 0),
  legend.spacing.x = unit(2, "mm"),
  legend.spacing.y = unit(0, "mm")
)

g + plotTheme

#########################################
# n vs year

df<-data.frame(x=my_data$year,y=log10(my_data$n))
use<-my_data$n<=250 
df<-df[use,]

g<-ggplot(df) + geom_bin_2d(aes(x=x,y=y),bins=30,show.legend=FALSE,drop=FALSE) 

g<-g+ylab("n")+xlab('year')
g<-g + scale_fill_gradient(low="white",high="red")
g + plotTheme + theme(panel.background=element_rect(fill="white", colour="black"),) 

years<-1985:2013
nmean<-c()
for (year in years) {
  use<-my_data$year==year
  nmean<-c(nmean,median(my_data$n[use]))
}
df<-data.frame(x=years,y=(nmean))
g<-ggplot()+geom_line(data=df,aes(x=x,y=y),color='yellow',lwd=1)
g<-g+scale_y_continuous(limits=c(0,100))
g<-g+xlab("year")+ylab("median(n)")
g + plotTheme


#########################################
#  n vs journal

journals<-unique(my_data$journal)
years<-1985:2013


nmean<-c()
for (journal in journals) {
  use<-my_data$journal==journal
  nmean<-c(nmean,median(my_data$n[use]))
}

dispOrder<-order(nmean)

df<-data.frame(x=1:length(dispOrder),y=nmean[dispOrder])
g<-ggplot()+geom_point(data=df,aes(x=x,y=y),color='yellow',size=4)
g<-g+scale_y_continuous(limits=c(0,110))
g<-g+scale_x_continuous(breaks=1:length(dispOrder),labels=journals[dispOrder])
g<-g+xlab("journal")+ylab("median(n)")
g + plotTheme

#########################################
#  t/r vs journal

journals<-unique(my_data$journal)

propRT<-c()
for (journal in journals) {
  ts<-my_data$journal==journal & my_data$Statistic=="t"
  rs<-my_data$journal==journal & my_data$Statistic=="r"
  propRT<-c(propRT,sum(ts)/(sum(ts)+sum(rs)))
}

df<-data.frame(x=1:length(dispOrder),y=propRT[dispOrder])
g<-ggplot()+geom_point(data=df,aes(x=x,y=y),color='yellow',size=4)
g<-g+scale_y_continuous(limits=c(0.5,1))
g<-g+scale_x_continuous(breaks=1:length(dispOrder),labels=journals[dispOrder])
g<-g+xlab("journal")+ylab("t/r")
g + plotTheme


#########################################
#  n vs year by journal

journals<-unique(my_data$journal)

years<-seq(1985,2015,5)

nmean<-matrix(NA,nrow=length(journals),ncol=length(years)-1)
for (j in 1:length(journals)) {
  j_use<-my_data$journal==journals[j]
  for (y in 1:(length(years)-1)) {
    y_use<-(my_data$year>=years[y]) & (my_data$year<years[y+1])
    nmean[j,y]<-median(my_data$n[j_use & y_use])
  }
}

pyears<-(years[1:(length(years)-1)]+years[2:length(years)])/2
g<-ggplot()
for (j in 1:length(journals)) {
df<-data.frame(x=pyears,y=nmean[j,],color=factor(journals[j]))
g<-g+geom_line(data=df,aes(x=x,y=y,color=color))
g<-g+geom_point(data=df,aes(x=x,y=y,color=color),size=4)
}
g<-g+scale_color_discrete()
g<-g+scale_y_continuous(limits=c(0,max(nmean,na.rm=TRUE)+5))
g<-g+xlab("year")+ylab("median(n)")
g + plotTheme

#########################################
nbins<-101

use_r<-my_data$r_s<=0.95
use_r_sim<-my_data_sim$r_s<=0.95

use_n<-my_data$n<=800 & my_data$n>=10 & log10(my_data$p)> -6
use_n_sim<-my_data_sim$n<=800 & my_data_sim$n>=10 & log10(my_data_sim$p)> -6

h2<-hist(log10(my_data_sim$p[use_n_sim&use_r_sim]),breaks=seq(-6,0,length.out=nbins))
h1<-hist(log10(my_data$p[use_n&use_r]),breaks=seq(-6,0,length.out=nbins))

df<-data.frame(x=h1$mids,y1=h1$density,y2=h2$density)
g<-ggplot()
g<-g+geom_point(data=df,aes(x=x,y=y1-y2),color='yellow',size=1)
g<-g+geom_line(data=df,aes(x=x,y=y1-y2))
g<-g+xlab("log10(p)")+ylab("real-sim")
g + plotTheme


#########################################
# rs x n density of real vs simulated 

nbins<-c(201,51)
k=0.3195403
pNull=0.725488

use_r<-my_data$r_s<=0.95
use_r_sim<-my_data_sim$r_s<=0.95

use_n<-my_data$n<=800 & my_data$n>=10
use_n_sim<-my_data_sim$n<=800 & my_data_sim$n>=10

h1<-hist2d(my_data$z_s[use_n&use_r],log10(my_data$n[use_n&use_r]),nbins=nbins,show=FALSE)

  my_data_sim$z_s<-atanh(my_data_sim$r_s)
  h2<-hist2d(my_data_sim$z_s[use_n_sim&use_r_sim],log10(my_data_sim$n[use_n_sim&use_r_sim]),nbins=nbins,show=FALSE)
  h2<-h2$counts
  for (i in 1:length(h1$y)) {
    h2[,i]<-h2[,i]/sum(h2[,i])*sum(h1$counts[,i])
  }
  
  h3<-c()
  for (i in 2:length(h1$y.breaks)) {
    expected<-0
    for (n in ceil(10^h1$y.breaks[i-1]:floor(10^h1$y.breaks[i]))) {
      zcrit<-qnorm(1-alpha/2,0,1/sqrt(n-3))
      dens1<-1-ExpSamplingCDF(h1$x.breaks,k,1/sqrt(n-3))
      dens0<-pnorm(h1$x.breaks,0,1/sqrt(n-3))
      dens<-dens1*(1-pNull)+dens0*pNull
      dens<-diff(dens)
      dens[h1$x<=zcrit]<-0
      expected<-expected+dens
    }
    expected<-expected/sum(expected)*sum(h1$counts[,i-1])
    h3<-cbind(h3,expected)
  }


# h2<-hist2d(my_data_sim$p[use_n_sim&use_r_sim],log10(my_data_sim$n[use_n_sim&use_r_sim]),nbins=nbins)
# h1<-hist2d(my_data$p[use_n&use_r],log10(my_data$n[use_n&use_r]),nbins=nbins)

xy<-meshgrid(h1$y,h1$x)
df<-data.frame(y=as.vector(xy$X),x=as.vector(xy$Y),
               f=as.vector(h1$counts/sum(h1$counts)-h3/sum(h1$counts))*1000)

g<-ggplot(df) + geom_raster(aes(x=x,y=y,fill=(abs(f)^0.5)*sign(f)))
g<-ggplot(df) + geom_raster(aes(x=x,y=y,fill=f))

ns<-h1$y
zs<-qnorm(1-alpha/2,0,1/sqrt(10^ns-3))
ds<-data.frame(x=zs,y=ns)
g<-g+geom_line(data=ds,aes(x=x,y=y),colour="white")
dst<-data.frame(x=max(zs),y=min(ns))
g<-g+geom_text(data=dst,aes(x=x,y=y,label="p=0.05"),hjust=1.2,colour="white")

zs<-qnorm(1-alpha/2/100,0,1/sqrt(10^ns-3))
ds<-data.frame(x=zs,y=ns)
g<-g+geom_line(data=ds,aes(x=x,y=y),colour="white")
dst<-data.frame(x=max(zs),y=min(ns))
g<-g+geom_text(data=dst,aes(x=x,y=y,label="p=0.0005"),hjust=1.2,colour="white")


intense<-scales::trans_new("intense", transform=function(x){abs(x)^0.5*sign(x)},inverse=function(x){abs(x)^(1/0.5)*sign(x)})
g<-g+ylab(bquote(bold(log[10](n))))+xlab(bquote(bold(z[s])))
g<-g+scale_fill_gradient2(high = "green", mid="#666666", low = "red",name="(actual-expected)/1000",trans=intense)
g + plotTheme

#########################################
