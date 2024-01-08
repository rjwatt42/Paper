
##########################################
# calculate distibution of sample sizes as tested

# first we calculate the empirical distribution
# using the overall power for each sample size
lambda<-0.332
p0<-0.72
nmax<-300

z<-seq(-2,2,length.out=2001)

n_tested<-rep(0,nmax)
for (n in 1:nmax) {
  zs1<-ExpSamplingPDF(seq(-2,2,length.out=2001),lambda,1/sqrt(n-3))
  zs0<-SingleSamplingPDF(seq(-2,2,length.out=2001),0,1/sqrt(n-3))
  zs<-zs1$pdf*(1-p0)+zs0$pdf*p0
  
  critZ<-qnorm(1-alpha/2,0,1/sqrt(n-3))
  w<-sum(zs[abs(z)<critZ])/sum(zs)
  
  n_count<-sum(my_data$n==n)
  n_orig<-round(n_count/w)
  
  n_tested[n]<-n_orig
}

plot(1:nmax,n_tested)

##########################################

n_offset<-4

nsamples<-c()
for (n in 1:nmax) {
  nsamples<-c(nsamples,rep(n,n_tested[n]))
}

fit <- fitdist(nsamples[nsamples>n_offset]-n_offset, distr = "gamma", method = "mle")
print(c(1/fit$estimate[2],fit$estimate[1]))
print(fit$loglik)

plot(fit)

