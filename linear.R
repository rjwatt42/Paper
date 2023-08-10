################################

journals<-unique(my_data$journal)
statistics<-c("t","F","r")
# nmaxs<-round(10^seq(log10(20),log10(200),length.out=10))
ncells<-length(journals)*length(statistics)*length(nmaxs)

bigResult<-c()
metaAnal<-list(meta_pdf="Exp",meta_psigAnal=TRUE,meta_nullAnal=TRUE)

j<-0
for (i1 in 1:length(journals)) {
  usej<-my_data$journal==journals[i1]
  for (i2 in 1:length(statistics)) {
    uses<-my_data$Statistic==statistics[i2]
    # for (i3 in 1:(length(nmaxs)-1)) {
    #   usen<-my_data$n>=nmaxs[i3] & my_data$n<nmaxs[i3+1]
      
      use<-usej & uses # & usen
      if (sum(use)>0) {
        metaData<-list(result=list(rIV=my_data$r_s[use],nval=my_data$n[use],df1=my_data$df1[use]))
        an<-runMetaAnalysis(metaAnal,metaData)
        nRes<-an$best$Nullmax
        kRes<-an$best$Kmax
        sRes<-an$best$Smax
        j<-j+1
        bigResult<-rbind(bigResult,c(i1,i2,nmaxs[i3],nRes,kRes,sRes))
      # }
    }
  }
}

################################

resultRawData<-data.frame(journal=factor(bigResult[,1]),
                          statistic=factor(bigResult[,2]),
                          nmax=bigResult[,3],
                          nRes=bigResult[,4],
                          kRes=bigResult[,5]
)

response<-"nRes"
formula<-paste0(response,"~journal+statistic")
lmRaw<-lm(formula=as.formula(formula),data=resultRawData,contrasts=list(journal=contr.sum,statistic=contr.sum))

rJournals<-sd(c(0,lmRaw$coefficients[2:8]))/sd(lmRaw$coefficients)
rStatistics<-sd(c(0,lmRaw$coefficients[9:10]))/sd(lmRaw$coefficients)
rNmax<-sd(lmRaw$coefficients[11]*lmRaw$model$nmax)/sd(lmRaw$coefficients)
direct<-c(rJournals,rStatistics,rNmax)

anv<-anova(lmRaw)
unique<-anv$`Sum Sq`[1:3]/sum(anv$`Sum Sq`)

lmJ<-lm(formula=as.formula(paste0(response,"~journal")),data=resultRawData,contrasts=list(journal=contr.sum))
lmS<-lm(formula=as.formula(paste0(response,"~statistic")),data=resultRawData,contrasts=list(statistic=contr.sum))
lmN<-lm(formula=as.formula(paste0(response,"~nmax")),data=resultRawData)

anvJ<-anova(lmJ)
anvS<-anova(lmS)
anvN<-anova(lmN)
total<-c(anvJ$`Sum Sq`[1]/sum(anvJ$`Sum Sq`),anvS$`Sum Sq`[1]/sum(anvS$`Sum Sq`),anvN$`Sum Sq`[1]/sum(anvN$`Sum Sq`))

print(direct)
print(unique)
print(total)

################################

