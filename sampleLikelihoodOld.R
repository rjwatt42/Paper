
ZeroSamplingPDF<-function(z,sigma) {
  exp(-0.5*(z^2/sigma^2))/sqrt(2*pi*sigma^2)
}
ZeroSamplingCDF<-function(zcrit,sigma) {
  1-(pnorm(zcrit,0,sigma)-pnorm(-zcrit,0,sigma))
}


SingleSamplingPDF<-function(z,lambda,sigma) {
  exp(-0.5*((z-lambda)^2/sigma^2))/sqrt(2*pi*sigma^2)
}
SingleSamplingCDF<-function(zcrit,lambda,sigma) {
  1-(pnorm(zcrit,lambda,sigma)-pnorm(-zcrit,lambda,sigma))
}


UniformSamplingPDF<-function(z,lambda,sigma) {
  (1-tanh(z)^2)
}
UniformSamplingCDF<-function(zcrit,lambda,sigma) {
  1-(tanh(zcrit)-tanh(-zcrit))/2
}


GaussSamplingPDFE<-function(z,lambda,sigma) {
  ni<-501
  zi<-seq(-1,1,length.out=ni)*32
  zd<-exp(-0.5*(zi/lambda)^2)
  zd<-zd/sum(zd)
  if (length(sigma)==1) sigma<-rep(sigma,length(z))
  res<-z*0
  for (i in 1:length(z)) {
    res[i]<-sum(zd*dnorm(zi,z[i],sigma[i]))
  }
  res
}
GaussSamplingCDFE<-function(zcrit,lambda,sigma) {
  ni<-101
  zi<-seq(-64,-zcrit,length.out=ni)
  zd<-GaussSamplingPDFE(zi,lambda,sigma)
  p1<-sum((zd[1:(ni-1)]+zd[2:ni])/2)*(zi[2]-zi[1])
  p1*2
}


d_zi=0.05
d_max=16
GammaSamplingPDF<-function(z,lambda,sigma,gamma_shape=1) {
  zi<-seq(-d_max,d_max,d_zi)
  zd<-dgamma(abs(zi),shape=gamma_shape,scale=lambda/gamma_shape)/2
  # zd<-zd/sum(zd)
  if (length(sigma)==1) sigma=rep(sigma,length(z))
  
  res<-z*0
  for (i in 1:length(z)) {
    res[i]<-sum(zd*dnorm(zi,z[i],sigma[i]))
  }
  res*(zi[2]-zi[1])
}
GammaSamplingCDF<-function(zcrit,lambda,sigma,gamma_shape=1) {
  res<-zcrit*0
  zcritUnique<-unique(zcrit)
  for (i in 1:length(zcritUnique)) {
    use<-which(zcrit==zcritUnique[i])
    zi<-seq(-d_max,-zcritUnique[i],d_zi)
    zi<-c(zi,-zcritUnique[i])
    zd<-GammaSamplingPDF(zi,lambda,sigma[use[1]],gamma_shape)
    areas<-(zd[1:(length(zi)-1)]+zd[2:length(zi)])/2*diff(zi)
    p1<-sum( areas )
    res[use]<-p1*2
  }
  res
}


ExpSamplingPDF<-function(z,lambda,sigma) {
  lambda<-1/lambda
  0.25*(lambda*exp(-lambda*(z-sigma^2*lambda/2))*(1+erf((z-sigma^2*lambda)/sqrt(2)/sigma)) +
          lambda*exp(-lambda*(-z-sigma^2*lambda/2))*(1+erf((-z-sigma^2*lambda)/sqrt(2)/sigma)))
}
ExpSamplingCDF<-function(zcrit,lambda,sigma) {
  lambda<-1/lambda
  z<-zcrit
  p1<-0.25*(
    exp((lambda*sigma/sqrt(2))^2)*exp(z*lambda) * erfc(lambda*sigma/sqrt(2) + z/sigma/sqrt(2))
    - exp((lambda*sigma/sqrt(2))^2)/exp(z*lambda) * erfc(lambda*sigma/sqrt(2) - z/sigma/sqrt(2))
    + 2*erf(z/sigma/sqrt(2))
  )
  1-2*p1
}

GaussSamplingPDF<-function(z,lambda,sigma) {
  sigma<-sqrt(lambda^2+sigma^2)
  exp(-0.5*(z/sigma)^2)/sqrt(2*pi)/sigma
}
GaussSamplingCDF<-function(zcrit,lambda,sigma) {
  sigma<-sqrt(lambda^2+sigma^2)
  1-(pnorm(zcrit,0,sigma)-pnorm(-zcrit,0,sigma))
}



getLogLikelihood<-function(z,n,worldDistr,worldDistK,worldDistNullP=0,p_sig=FALSE,gamma_shape=NULL) {
  sigma<-1/sqrt(n-3)
  
  # get nulls ready first
  if (any(worldDistNullP>0)) {
    nullLikelihoods<-SingleSamplingPDF(z,0,sigma)
    if (p_sig) {
      zcrit<-qnorm(1-alpha/2,0,sigma)
      gainNull<-alpha
    } else {
      zcrit<-0
      gainNull<-1
    }
  } else {
    nullLikelihoods<-0
    gainNull<-0
    zcrit<-0
  } 
  
  gainMain<-1
  res<-matrix(0,nrow=length(worldDistK),ncol=length(worldDistNullP))
  switch(worldDistr,
         "Single"={
           for (i in 1:length(worldDistK)) {
             lambda<-worldDistK[i]
             if (p_sig) {
               gainMain<-SingleSamplingCDF(zcrit,lambda,sigma)
             }
             mainLikelihoods<-SingleSamplingPDF(z,lambda,sigma)
             for (j in 1:length(worldDistNullP)) {
               likelihoods<-(mainLikelihoods*(1-worldDistNullP[j])+nullLikelihoods*worldDistNullP[j])/(gainMain*(1-worldDistNullP[j])+gainNull*worldDistNullP[j])
               res[i,j]<-sum(log(likelihoods[likelihoods>1e-300]),na.rm=TRUE)
             }
           }
         },
         "Gauss"={
           for (i in 1:length(worldDistK)) {
             lambda<-worldDistK[i]
             if (p_sig) {
               gainMain<-GaussSamplingCDF(zcrit,lambda,sigma)
             }
             mainLikelihoods<-GaussSamplingPDF(z,lambda,sigma)
             for (j in 1:length(worldDistNullP)) {
               likelihoods<-(mainLikelihoods*(1-worldDistNullP[j])+nullLikelihoods*worldDistNullP[j])/(gainMain*(1-worldDistNullP[j])+gainNull*worldDistNullP[j])
               res[i,j]<-sum(log(likelihoods[likelihoods>1e-300]),na.rm=TRUE)
             }
           }
         },
         "Exp"={
           for (i in 1:length(worldDistK)) {
             lambda<-worldDistK[i]
             z<-abs(z)
             if (p_sig) {
               gainMain<-ExpSamplingCDF(zcrit,lambda,sigma)
             } else {
               gainMain<-1
             }
             mainLikelihoods<-ExpSamplingPDF(z,lambda,sigma)
             for (j in 1:length(worldDistNullP)) {
               likelihoods<-(mainLikelihoods*(1-worldDistNullP[j])+nullLikelihoods*worldDistNullP[j])/(gainMain*(1-worldDistNullP[j])+gainNull*worldDistNullP[j])
               res[i,j]<-sum(log(likelihoods[likelihoods>1e-300]),na.rm=TRUE)
             }
           }
         },
         "Gamma"={
           for (i in 1:length(worldDistK)) {
             lambda<-worldDistK[i]
             z<-abs(z)
             if (p_sig) {
               gainMain<-GammaSamplingCDF(zcrit,lambda,sigma,gamma_shape)
             } else {
               gainMain<-1
             }
             mainLikelihoods<-GammaSamplingPDF(z,lambda,sigma,gamma_shape)
             for (j in 1:length(worldDistNullP)) {
               likelihoods<-(mainLikelihoods*(1-worldDistNullP[j])+nullLikelihoods*worldDistNullP[j])/(gainMain*(1-worldDistNullP[j])+gainNull*worldDistNullP[j])
               res[i,j]<-sum(log(likelihoods[likelihoods>1e-300]),na.rm=TRUE)
             }
           }
         }
  )
  res
}

