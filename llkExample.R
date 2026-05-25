########################

n_samples=1000
real_lambda=0.3

lambdas<-seq(0.1,0.8,0.025)
llk<-lambdas*0

setEvidence(shortHand = TRUE)
for (j in 1:n_samples) {
  setWorld(TRUE,"Exp","z",real_lambda,pRplus=1)
  s<-doSingle()
  setPossible(targetSample=s$rIV,targetSampleN=s$nval,UsePrior="world",axisType = "z")
  for (i in 1:length(lambdas)) {
    setWorld(TRUE,"Exp","z",lambdas[i],pRplus=1)
    q<-doPossible()
    llk[i]<-llk[i]+log(q$mleTotal)
  }
}

g<-startPlot(xlim=c(0,1),ylim=c(min(llk),max(llk)),xtick=seq(0,1,0.1),ytick=seq(0,1,0.1),
             xlabel="z[λ]",ylabel="likelihood")
g<-addG(g,dataPath(data.frame(x=lambdas,y=llk)),dataPoint(data.frame(x=lambdas,y=llk)))
print(g)

###########################
# with publication bias

n_samples=1000
real_lambda=0.2
doPublicationBias=TRUE
analysePublicationBias=1

setDesign("Psych")

lambdas<-seq(0.1,0.8,0.025)
llk<-lambdas*0

setEvidence(shortHand = TRUE)
count_samples<-0
for (j in 1:n_samples) {
  setWorld(TRUE,"Exp","z",real_lambda,pRplus=1)
  s<-doSingle()
  zs<-abs(s$rIV)
  n<-s$nval
  if (!doPublicationBias || s$pIV<0.05) {
  for (i in 1:length(lambdas)) {
    llk[i]<-llk[i]+getLogLikelihood(zs,n,1,"Exp",lambdas[i],bias=analysePublicationBias)
  }
    count_samples<-count_samples+1
  }
}

mle<-lambdas[which.max(llk)]

llk<-llk-max(llk)
g<-startPlot(xlim=c(0,1),ylim=c(-count_samples,count_samples*0.1),
             xtick=seq(0,1,0.1),ytick=seq(-count_samples,0,count_samples*0.2),
             xlabel="z[λ]",ylabel="log likelihood",
             top=TRUE)
g<-addG(g,plotTitle(paste0("n_samples=",format(count_samples)),
                    "right"))
g<-addG(g,plotTitle(paste0("z[MLE]=",format(mle,digits=2)),
                    "left"))
g<-addG(g,dataPath(data.frame(x=lambdas,y=llk)),dataPoint(data.frame(x=lambdas,y=llk)))
print(g)


###########################


n_samples=1000
real_lambda=0.2
doPublicationBias=TRUE
analysePublicationBias=1

setDesign("Psych")

lambdas<-seq(0.1,0.8,0.025)
llk2<-llk1<-lambdas*0
count_samples<-0
zss<-ns<-zps<-c()
for (j in 1:n_samples) {
  n<-nDistrRand(1,braw.def$design)
  zp<-rexp(1,1/real_lambda)*sign(runif(1,-1,1))
  zs<-abs(zp+rnorm(1,0,1/sqrt(n-3)))
  zcrit<-wn2z(0.5,n)
  if (!doPublicationBias || zs>=zcrit) {
    for (i in 1:length(lambdas)) {
      llk1[i]<-llk1[i]+getLogLikelihood(zs,n,1,"Exp",lambdas[i],bias=FALSE)
      llk2[i]<-llk2[i]+getLogLikelihood(zs,n,1,"Exp",lambdas[i],bias=TRUE)
    }
    zss<-c(zss,zs)
    ns<-c(ns,n)
    count_samples<-count_samples+1
  }
  zps<-c(zps,zp)
}

mle1<-lambdas[which.max(llk1)]
llfun1<-function(x) -(getLogLikelihood(zss,ns,1,"Exp",location=x[1],prplus=1,bias=FALSE))
result<-optimize(llfun1,interval=mle1+c(-1,1)*0.1)
mleFull1<-result$minimum

mle2<-lambdas[which.max(llk2)]
llfun2<-function(x) -(getLogLikelihood(zss,ns,1,"Exp",location=x[1],prplus=1,bias=TRUE))
result<-optimize(llfun2,interval=mle2+c(-1,1)*0.1)
mleFull2<-result$minimum

# res<-fminsearch(llfun,mle,method='Hooke-Jeeves',lower=lb[np],upper=ub[np])


# llk1<-llk1-max(llk1)
# llk2<-llk2-max(llk1)
llk<-c(llk1,llk2)
g<-startPlot(xlim=c(0,1),ylim=c(min(llk),max(llk))+c(-1,1)*(max(llk)-min(llk)),
             xtick=seq(0,1,0.1),
             xlabel="z[λ]",ylabel="log likelihood",
             top=TRUE)
g<-addG(g,plotTitle(paste0("n_samples=",format(count_samples)),
                    "right"))
g<-addG(g,plotTitle(paste0("z[MLE]=",format(mleFull2,digits=2)),
                    "left"))
g<-addG(g,dataPath(data.frame(x=lambdas,y=llk1)),dataPoint(data.frame(x=lambdas,y=llk1),fill="red"))
g<-addG(g,dataPath(data.frame(x=lambdas,y=llk2)),dataPoint(data.frame(x=lambdas,y=llk2),fill="yellow"))

g<-addG(g,dataLegend(data.frame(colours=c("red","yellow"),names=c("-","+")),title="Publication Bias"))
print(g)

#############################



n_samples=1000
real_lambda=0.2
doPublicationBias=TRUE
analysePublicationBias=1

setDesign("Psych")

lambdas<-seq(0.1,0.8,0.025)
llk2<-llk1<-lambdas*0
count_samples<-0
zss<-ns<-zps<-c()
for (j in 1:n_samples) {
  n<-nDistrRand(1,braw.def$design)
  zp<-rexp(1,1/real_lambda)*sign(runif(1,-1,1))
  zs<-abs(zp+rnorm(1,0,1/sqrt(n-3)))
  zcrit<-wn2z(0.5,n)
  if (!doPublicationBias || zs>=zcrit) {
    for (i in 1:length(lambdas)) {
      llk1[i]<-llk1[i]+getLogLikelihood(zs,n,1,"Exp",lambdas[i],bias=FALSE)
      llk2[i]<-llk2[i]+getLogLikelihood(zs,n,1,"Exp",lambdas[i],bias=TRUE)
    }
    zss<-c(zss,zs)
    ns<-c(ns,n)
    count_samples<-count_samples+1
  }
  zps<-c(zps,zp)
}

mle1<-lambdas[which.max(llk1)]
llfun1<-function(x) -(getLogLikelihood(zss,ns,1,"Exp",location=x[1],prplus=1,bias=FALSE))
result<-optimize(llfun1,interval=mle1+c(-1,1)*0.1)
mleFull1<-result$minimum

mle2<-lambdas[which.max(llk2)]
llfun2<-function(x) -(getLogLikelihood(zss,ns,1,"Exp",location=x[1],prplus=1,bias=TRUE))
result<-optimize(llfun2,interval=mle2+c(-1,1)*0.1)
mleFull2<-result$minimum

# res<-fminsearch(llfun,mle,method='Hooke-Jeeves',lower=lb[np],upper=ub[np])


# llk1<-llk1-max(llk1)
# llk2<-llk2-max(llk1)
llk<-c(llk1,llk2)
g<-startPlot(xlim=c(0,1),ylim=c(min(llk),max(llk))+c(-1,1)*(max(llk)-min(llk)),
             xtick=seq(0,1,0.1),
             xlabel="z[λ]",ylabel="log likelihood",
             top=TRUE)
g<-addG(g,plotTitle(paste0("n_samples=",format(count_samples)),
                    "right"))
g<-addG(g,plotTitle(paste0("z[MLE]=",format(mleFull2,digits=2)),
                    "left"))
g<-addG(g,dataPath(data.frame(x=lambdas,y=llk1)),dataPoint(data.frame(x=lambdas,y=llk1),fill="red"))
g<-addG(g,dataPath(data.frame(x=lambdas,y=llk2)),dataPoint(data.frame(x=lambdas,y=llk2),fill="yellow"))

g<-addG(g,dataLegend(data.frame(colours=c("red","yellow"),names=c("-","+")),title="Publication Bias"))
print(g)

