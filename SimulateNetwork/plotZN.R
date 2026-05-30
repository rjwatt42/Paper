plotZN<-function(samples,replication=NULL,field1="zs",field2="n") {

  np<-1000
  alpha<-0.3
  
  if (!is.null(replication)) {
    col1<-"#FA4"
    sig<-replication$sigrep
    useData<-replication$allrep
  } else {
    col1<-"#F40"
    sig<-samples$sig
    useData<-samples$all
  }
    
  switch(field1,
         "zs"={x<-useData$zs; xlabel<-"z[s]"; xlim=c(0,1.5)},
         "zp"={x<-useData$zp; xlabel<-"z[p]"; xlim=c(0,1.5)},
         "n"={x<-log10(useData$n); xlabel<-"log[10]n"; xlim=c(0.5,3)},
         "no"={x<-log10(useData$no); xlabel<-"log[10]no"; xlim=c(0.5,3)}
  )
  switch(field2,
         "zs"={y<-useData$zs; ylabel<-"z[s]"; ylim=c(0,1.5)},
         "zp"={y<-useData$zp; ylabel<-"z[sp]"; ylim=c(0,1.5)},
         "n"={y<-log10(useData$n); ylabel<-"log[10]n"; ylim=c(0.5,3)},
         "no"={y<-log10(useData$no); ylabel<-"log[10]no"; ylim=c(0.5,3)}
  )
  
  data<-data.frame(x=x[!sig],y=y[!sig])
  if (length(data$x)>np) data<-data[1:np,]
  g<-dataGraph(data,colour=NA,fill=col1,alpha=alpha,
               xlim=xlim,xlabel=xlabel,xticks="auto",
               ylim=ylim,ylabel=ylabel,yticks="auto"
  )
  data<-data.frame(x=x[sig],y=y[sig])
  if (length(data$x)>np) data<-data[1:np,]
  g<-addG(g,dataPoint(data,fill="#4F4",alpha=alpha))
  

return(g)
}
