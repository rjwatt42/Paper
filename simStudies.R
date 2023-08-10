############################
# an example

nStudies=92877
ns<-my_data$n
ns<-NA
model="Exp"
k=0.3195403
pNull=0.725488
sigOnly=TRUE

my_sim_data<-makeStudies(nStudies,ns,model,k,NA,pNull,sigOnly)

metaAnal<-list(meta_fixedAnal="random",meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE,append=FALSE)
metaData<-list(result=list(rIV=my_sim_data$r_s,nval=my_sim_data$n))
# metaData<-list(result=list(rIV=my_data$r_s,nval=my_data$n))

an<-runMetaAnalysis(metaAnal,metaData)
showAnalysis(an)


############################
# next - 100 examples

nsims<-100

#first we get the analysis for the real data
source("OriginalData.R")
metaAnal<-list(meta_fixedAnal="random",meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE,append=FALSE)
metaData<-list(result=list(rIV=my_data$r_s,nval=my_data$n))

an<-runMetaAnalysis(metaAnal,metaData)

nStudies=length(my_data$n)
ns<-NA
model="Exp"
k=an$exp$Kmax
pNull=an$exp$Nullmax
actualS<-an$exp$Smax
sigOnly=TRUE


s<-c()
np<-c()
km<-c()
times<-c()
for (i in 1:nsims) {
  start<-Sys.time()
  my_data_sim<-makeStudies(nStudies,ns,model,k,shape,pNull,sigOnly)
  metaAnal<-list(meta_fixedAnal="random",meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE,append=FALSE)
  metaData<-list(result=list(rIV=my_data_sim$r_s,nval=my_data_sim$n))
  
  an<-runMetaAnalysis(metaAnal,metaData)
  
  times<-c(times,Sys.time()-start)
  print(c(i,an$exp$Smax,times[i],(nsims-i)*mean(times)/60))
  s<-c(s,an$exp$Smax)
  km<-c(km,an$exp$Kmax)
  np<-c(np,an$exp$Nullmax)
}

g<-ggplot()
pts<-data.frame(s=s)
g<-g+geom_histogram(data=pts,aes(x=s, after_stat(ndensity)),color="white",fill="white")
g<-g+geom_vline(xintercept=actualS,color="red")
g<-g+geom_label(data=data.frame(x=actualS,y=1,label=paste0("Actual data = ",actualS)),aes(x=x,y=y,label=label))
g<-g+xlab("log(lk)")+ylab("Density")

g+plotTheme



