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
g<-g+xlab(expression(lambda))+ylab("probability density")
g+plotTheme

##########################################
# a likelihood version - with and without publication bias

z<-seq(-2,2,0.01)
nullP<-0.2

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

    g<-g+ylab("Probability Density")+xlab(expression(z[s]))
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
    expS<-getLogLikelihood(z,n,analysisDistribution,kvals,nullval,psig)
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
  
  expS<-getLogLikelihood(zs,n,analysisDistribution[ai],kvals,nullval,psig)
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

g<-g+xlab(expression(lambda))+ylab("log(likelihood)")

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
    expS<-getLogLikelihood(z,n,analysisDistribution[1],atanh(kvals),null,psig)
    rymax<-c(rymax,max(expS))
    rymin<-c(rymin,min(expS))
    pts<-data.frame(x=kvals,y=expS)
    g<-g+geom_line(data=pts,aes(x=x,y=y),lwd=0.4,lty=1)
  }
}
maxD<- -Inf
for (ai in 1:length(analysisDistribution)) {
  expS<-getLogLikelihood(zs,n,analysisDistribution[ai],atanh(kvals),null,psig)
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

g<-g+xlab(expression(lambda))+ylab("log(likelihood)")
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
null<-0.5
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
nullvals<-seq(0,0.999,length.out=101)

expS<-getLogLikelihood(zs,n,analysisDistribution,kvals,nullvals,psig)
threshold<-quantile(expS,0.25)
expS[expS < threshold]<-threshold

df<-melt(expS)
g<-ggplot(df)+geom_tile(aes(x=(X2-1)/100*(nullvals[101]-nullvals[1])+nullvals[1],y=(X1-1)*(kvals[51]-kvals[1])/50+kvals[1],fill = value),show.legend=showLegend)
g<-g+stat_contour(data=df,aes(x=(X2-1)/100*(nullvals[101]-nullvals[1])+nullvals[1],y=(X1-1)*(kvals[51]-kvals[1])/50+kvals[1],z = value),breaks=max(expS)-log(10),color="red")

use<-which(expS==max(expS), arr.ind = TRUE)
g<-g+geom_vline(xintercept=nullvals[use[1,2]],color="red")
g<-g+geom_hline(yintercept=kvals[use[1,1]],color="red")
g<-g+geom_text(data=data.frame(x=nullvals[use[1,2]],y=kvals[use[1,1]],
                               label=paste("\u03BB","=",format(kvals[use[1,1]],digits=2),"\n p(null)=",format(nullvals[use[1,2]],digits=2))),
               aes(x=x,y=y,label=label),color="white",vjust=-0.5)

g<-g+xlab(bquote(p[null]))+ylab(expression(lambda))
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

g<-g+xlab(expression(z[p]))+ylab("probability density")
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
g<-g+xlab(expression(z[p]))+ylab("probability density")
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

g<-g+xlab(expression(z[p]))+ylab("log10(PDF)")
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

