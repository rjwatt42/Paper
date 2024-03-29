
r2p<-function(r,n,df1=1){
  if (!is.numeric(r) || !is.numeric(n)) {return(1)}
  if (any(abs(r)>1)) {
    print(paste("r2p r-exception",format(max(abs(r)),digits=3)))
    r[r>1]<-1
    r[r < -1]<- -1
  }
  if (any(abs(n)<3)) {
    print("r2p n-exception")
    n[n<3]<-4
  }
  df2<-n-(df1+1)
  
    Fvals<-r^2/(1-r^2)*df2/df1
    p<-(1-pf(Fvals,df1,df2))
    return(p)
    
  #   t_vals<-r/r2se(r,n)
  #   p<-(1-pt(abs(t_vals),df2))*2
}

p2r<-function(p,n,df1=1) {
  if (any(abs(n)<3)) {
    print("p2r n-exception")
    n[n<3]<-4
  }
  df2<-n-(df1+1)
  
  Fvals <- qf(1-p,df1,df2)
  r <- sqrt(Fvals*df1)/sqrt(Fvals*df1+df2)
  return(r)
  
  t_vals <- qt(p/2,n-2)
  r_vals <- t_vals/sqrt(t_vals^2+(n-2))
  r_vals
}

r2se<-function(r,n){
  if (any(abs(r)>1)) {
    print(paste("r2se r-exception",format(max(abs(r)),digits=3)))
    r[r>1]<-1
    r[r < -1]<- -1
  }
  if (any(abs(n)<3)) {
    print("r2se n-exception")
    n[n<3]<-4
  }
  sqrt((1-r^2)/(as.vector(n)-2))
}



