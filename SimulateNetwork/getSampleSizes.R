getSampleSizes<-function(nSamples,dist="Gamma",sN=50,sNRandSD=33,minN=5) {
  
  
    switch(dist,
           "Gamma"={
             return(
               round(minN+rgamma(nSamples,
                                          shape=(sN-minN)/sNRandSD,
                                          scale=sNRandSD)
               )
             )
           },
           "Gauss"={
             return(round(minN+abs(rnorm(nSamples,0,sNRandSD))))
           },
           "Exp"={
             return(round(minN+abs(rexp(nSamples,1/sNRandSD))))
           },
           "Uniform"={
             return(round(runif(nSamples,minN,minN+2*sN)))
           }
    )

}