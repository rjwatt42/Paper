d_zi=0.05
d_max=4

SingleSamplingPDF<-function(z,lambda,sigma,shape,remove_nonsig=FALSE) {
  d1<-exp(-0.5*((z-lambda)^2/sigma^2))/sqrt(2*pi*sigma^2)
  if (remove_nonsig) {
    zcrit<-qnorm(1-alpha/2,0,sigma)
    d0<-1-(pnorm(zcrit,lambda,sigma)-pnorm(-zcrit,lambda,sigma))
  } else {
    d0<-1
  }
  return(list(pdf=d1,sig_pdf=d0))
}
SingleSamplingCDF<-function(zcrit,lambda,sigma,shape) {
  1-(pnorm(zcrit,lambda,sigma)-pnorm(-zcrit,lambda,sigma))
}


UniformSamplingPDF<-function(z,lambda,sigma,shape) {
  (1-tanh(z)^2)
}
UniformSamplingCDF<-function(zcrit,lambda,sigma,shape) {
  1-(tanh(zcrit)-tanh(-zcrit))/2
}


GaussSamplingPDF<-function(z,lambda,sigma,shape=NA,remove_nonsig=FALSE) {
  sigma2<-sqrt(lambda^2+sigma^2)
  d1<-exp(-0.5*z^2/sigma2^2)/sqrt(2*pi*sigma2^2)
  
  if (remove_nonsig) {
    zcrit<-qnorm(1-alpha/2,0,sigma)
    d0<-GaussSamplingCDF(zcrit,lambda,sigma)
  } else {
    d0<-1
  }
  return(list(pdf=d1,sig_pdf=d0))
}
GaussSamplingCDF<-function(zcrit,lambda,sigma) {
  sigma2<-sqrt(lambda^2+sigma^2)
  1-(pnorm(zcrit,0,sigma2)-pnorm(-zcrit,0,sigma2))
}


ExpSamplingPDF<-function(z,lambda,sigma,shape=NA,remove_nonsig=FALSE) {
  lambda1<-1/lambda
  d1<-0.25*(lambda1*exp(-lambda1*(z-sigma^2*lambda1/2))*(1+erf((z-sigma^2*lambda1)/sqrt(2)/sigma)) +
          lambda1*exp(-lambda1*(-z-sigma^2*lambda1/2))*(1+erf((-z-sigma^2*lambda1)/sqrt(2)/sigma)))
  
  if (remove_nonsig) {
    zcrit<-qnorm(1-alpha/2,0,sigma)
    d0<-ExpSamplingCDF(zcrit,lambda,sigma)
  } else {
    d0<-1
  }
  return(list(pdf=d1,sig_pdf=d0))
}
ExpSamplingCDF<-function(zcrit,lambda,sigma) {
  lambda1<-1/lambda
  z<-zcrit
  p1<-0.25*(
    exp((lambda1*sigma/sqrt(2))^2)*exp(z*lambda1) * erfc(lambda1*sigma/sqrt(2) + z/sigma/sqrt(2)) -
      exp((lambda1*sigma/sqrt(2))^2)/exp(z*lambda1) * erfc(lambda1*sigma/sqrt(2) - z/sigma/sqrt(2)) +
        2*erf(z/sigma/sqrt(2))
  )
  z<--zcrit
  p2<-0.25*(
      exp((lambda1*sigma/sqrt(2))^2)*exp(z*lambda1) * erfc(lambda1*sigma/sqrt(2) + z/sigma/sqrt(2))
    - exp((lambda1*sigma/sqrt(2))^2)/exp(z*lambda1) * erfc(lambda1*sigma/sqrt(2) - z/sigma/sqrt(2))
    + 2*erf(z/sigma/sqrt(2))
  )
  1-(p1-p2)
}

convolveWith<-function(zi,zpd,z,sigma) {
  d1<-z*0
  for (i in 1:length(z)) {
    zs<-zpd*dnorm(zi,z[i],sigma[i])
    d1[i]<-sum(zs)*d_zi
  }
  return(d1)
}

removeNonSig<-function(zi,zpd,sigma) {
  zcrit<-qnorm(1-alpha/2,0,sigma)
  # d2<-GammaSamplingCDF(zcrit,lambda,sigma,gamma_shape)
  d2<-zcrit*0
  zcritUnique<-unique(zcrit)
  for (i in 1:length(zcritUnique)) {
    use<-which(zcrit==zcritUnique[i])
    zi1<-seq(-d_max,-zcritUnique[i],d_zi)
    zi1<-c(zi1,-zcritUnique[i])
    
    d0<-zi1*0
    for (j in 1:length(zi1)) {
      zs<-zpd*dnorm(zi,zi1[j],sigma[use[1]])
      d0[j]<-sum(zs)*d_zi
    }
    areas<-(d0[1:(length(zi1)-1)]+d0[2:length(zi1)])/2*diff(zi1)
    d2[use]<-sum(areas)*2
  }
  return(d2)
}


GammaSamplingPDF<-function(z,lambda,sigma,gamma_shape=1,remove_nonsig=FALSE) {
  if (length(sigma)==1) {sigma<-rep(sigma,length(z))}
  
  zi<-seq(-d_max,d_max,d_zi)
  zpd<-dgamma(abs(zi),shape=gamma_shape,scale=lambda/gamma_shape)
  zpd<-zpd/(sum(zpd)*d_zi)
  
  d1<-convolveWith(zi,zpd,z,sigma)

  if (remove_nonsig) {
    d2<-removeNonSig(zi,zpd,sigma)
  } else {
    d2<-1
  }
  return(list(pdf=d1,sig_pdf=d2))
  
}
GammaSamplingCDF<-function(zcrit,lambda,sigma,gamma_shape=1) {
  res<-zcrit*0
  zcritUnique<-unique(zcrit)
  for (i in 1:length(zcritUnique)) {
    use<-which(zcrit==zcritUnique[i])
    zi<-seq(-d_max,-zcritUnique[i],d_zi)
    zi<-c(zi,-zcritUnique[i])
    zd<-GammaSamplingPDF(zi,lambda,sigma[use[1]],gamma_shape)
    areas<-(zd$pdf[1:(length(zi)-1)]+zd$pdf[2:length(zi)])/2*diff(zi)
    p1<-sum( areas )
    res[use]<-p1*2
  }
  res
}

GenExpSamplingPDF<-function(z,lambda,sigma,genexp_shape=1,remove_nonsig=FALSE) {
  
  if (all(sigma==0)) {
    zd<-1-(1-exp(-abs(z)/lambda))^genexp_shape
    zd<-zd/(sum(zd)*(z[2]-z[1]))
    return(zd)
  }

  if (length(sigma)==1) {sigma<-rep(sigma,length(z))}

  zi<-seq(-d_max,d_max,d_zi)
  zpd<-1-(1-exp(-abs(zi)/lambda))^genexp_shape
  zpd<-zpd/(sum(zpd)*d_zi)

  d1<-convolveWith(zi,zpd,z,sigma)
  
  if (remove_nonsig) {
    d2<-removeNonSig(zi,zpd,sigma)
  } else {
    d2<-1
  }
  return(list(pdf=d1,sig_pdf=d2))
  
}
GenExpSamplingCDF<-function(zcrit,lambda,sigma,shape=1) {
  res<-zcrit*0
  zcritUnique<-unique(zcrit)
  for (i in 1:length(zcritUnique)) {
    use<-which(zcrit==zcritUnique[i])
    zi<-seq(-d_max,-zcritUnique[i],d_zi)
    zi<-c(zi,-zcritUnique[i])
    zd<-GenExpSamplingPDF(zi,lambda,sigma[use[1]],shape)
    areas<-(zd[1:(length(zi)-1)]+zd[2:length(zi)])/2*diff(zi)
    p1<-sum( areas )
    res[use]<-p1*2
  }
  res
}



getLogLikelihood<-function(z,n,worldDistr,worldDistK,worldDistNullP=0,remove_nonsig=FALSE) {
  sigma<-1/sqrt(n-3)
  zcrit<-qnorm(1-alpha/2,0,sigma)
  
  # get nulls ready first
  if (any(worldDistNullP>0)) {
    nullPDF<-SingleSamplingPDF(z,0,sigma,NA,remove_nonsig)
  } else {
    nullPDF<-list(pdf=0,sig_pdf=1)
    zcrit<-0
  } 
  shape<-NA
  res<-matrix(0,nrow=length(worldDistK),ncol=length(worldDistNullP))
  switch(worldDistr,
         "Single"={
           PDF<-SingleSamplingPDF
         },
         "Gauss"={
           PDF<-GaussSamplingPDF
         },
         "Exp"={
           PDF<-ExpSamplingPDF
         },
         "Gamma"={
           PDF<-GammaSamplingPDF
           shape<-metaAnal$shape
         },
         "GenExp"={
           PDF<-GenExpSamplingPDF
           shape<-metaAnal$shape
         }
  )
  for (i in 1:length(worldDistK)) {
    lambda<-worldDistK[i]
    mainPDF<-PDF(z,lambda,sigma,shape,remove_nonsig)
    for (j in 1:length(worldDistNullP)) {
      nullP<-worldDistNullP[j]
      # make the whole source first
      sourcePDF<-mainPDF$pdf*(1-nullP)+nullPDF$pdf*nullP
      # now normalize for the non-sig
      likelihoods<-sourcePDF/(mainPDF$sig_pdf*(1-nullP)+nullPDF$sig_pdf*nullP)
      res[i,j]<-sum(log(likelihoods[likelihoods>1e-300]),na.rm=TRUE)
    }
  }
  res
}

