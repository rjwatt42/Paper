alpha<-0.05



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
  df2<-n-(sum(df1)+1)
  if (any(df1>1)) {
    Fvals<-r^2/(1-r^2)*df2/df1
    (1-pf(Fvals,df1,df2))
  } else {
    t_vals<-r/r2se(r,n)
    (1-pt(abs(t_vals),df2))*2
  }
  
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



