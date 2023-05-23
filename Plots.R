maincolours<-list(windowC="#002D40",panelC="#005E86",graphC="#BFECFF")

mainTheme=theme(panel.background = element_rect(fill="#666666", colour="black"),
                panel.grid.major = element_line(linetype="blank"),panel.grid.minor = element_line(linetype="blank"),
                plot.background = element_rect(fill=maincolours$graphC, colour=maincolours$graphC),
                plot.margin = margin(2,2,1,1,"cm"))
SMplotTheme=theme(plot.title=element_text(size=16,face="bold"),axis.title=element_text(size=16),axis.text.x=element_text(size=10),axis.text.y=element_text(size=10))
plotTheme=mainTheme+SMplotTheme

plotBlankTheme=theme(panel.background = element_rect(fill=maincolours$graphC, colour=maincolours$graphC),
                     panel.grid.major = element_line(linetype="blank"),panel.grid.minor = element_line(linetype="blank"),
                     plot.background = element_rect(fill=maincolours$graphC, colour=maincolours$graphC),
                     axis.title=element_text(size=16,face="bold")
)
gridTheme=theme(plot.margin=margin(0,0,0,0,"cm"))



##########################################
showAnalysis<-function(an,param=NULL) {
  
  if (is.null(param)) {
    s<-""
  } else {
    if (param=="All") {
      s<-paste0(param,": ",an$bestDist)
    } else {
      s<-paste0(param,": ")
    }
  }
    print(paste0(s,
               "(",
               "H0=",format(an$bestNull*100,digits=2),"%",
               "  k=",format(an$bestK,digits=3),
               "  S=",format(an$bestS,digits=3),
               "  n=",format(length(an$result$rIV),digits=2),
               ")"))
}

##########################################
drawAnalysis<-function(an1,metaData1) {
  
  ymin<- -100
  gain<-1000
  g<-ggplot()+plotBlankTheme+theme(plot.margin=margin(0,-0.2,0,0,"cm"))
  g<-g+scale_x_continuous(limits = c(0,10),labels=NULL,breaks=NULL)+scale_y_continuous(limits = c(0,10),labels=NULL,breaks=NULL)
  
  g1<-ggplot()
  g2<-ggplot()
  
  k<-seq(0.05,1,0.01)
  nullP<-seq(0.0,0.95,0.01)
  z<-atanh(metaData1$result$rIV)
  n<-metaData1$result$nval
  
  if (any(an1$metaAnalysis$meta_pdf=="Exp") || an1$metaAnalysis$meta_pdf=="All") {
    Sk1<-getLogLikelihood(z,n,"Exp",k,an1$exp$nullMax,p_sig=TRUE)
    Sk1<-Sk1/gain
    
    pts1<-data.frame(x=k,y=Sk1)
    g1<-g1+geom_line(data=pts1,aes(x=x,y=y,col="Exp"),lwd=1)
    
    Sn1<-getLogLikelihood(z,n,"Exp",an1$exp$kmax,nullP,p_sig=TRUE)
    Sn1<-Sn1[1,]/gain
    
    pts1<-data.frame(x=nullP,y=Sn1)
    g2<-g2+geom_line(data=pts1,aes(x=x,y=y,col="Exp"),lwd=1)
  }
  
  if (any(an1$metaAnalysis$meta_pdf=="Gauss") || an1$metaAnalysis$meta_pdf=="All") {
    Sk1<-getLogLikelihood(z,n,"Gauss",k,an1$exp$nullMax,p_sig=TRUE)
    Sk1<-Sk1/gain
    
    pts1<-data.frame(x=k,y=Sk1)
    g1<-g1+geom_line(data=pts1,aes(x=x,y=y,col="Gauss"),lwd=1)
    
    Sn1<-getLogLikelihood(z,n,"Gauss",an1$exp$kmax,nullP,p_sig=TRUE)
    Sn1<-Sn1[1,]/gain
    
    pts<-data.frame(x=nullP,y=Sn1)
    g2<-g2+geom_line(data=pts,aes(x=x,y=y,col="Gauss"),lwd=1)
  }
  g1<-g1+scale_color_manual(name = NULL, values = c("Exp" = "red", "Gauss" = "yellow"))
  g2<-g2+scale_color_manual(name = NULL, values = c("Exp" = "red", "Gauss" = "yellow"))
  
  if (an1$metaAnalysis$meta_pdf=="Gamma") {
    Sk1<-getLogLikelihood(z,n,"Gamma",k,an1$exp$nullMax,p_sig=TRUE)
    Sk1<-Sk1/gain
    
    pts1<-data.frame(x=k,y=Sk1)
    g1<-g1+geom_line(data=pts,aes(x=x,y=y,col="green"),lwd=0.5)
  }

  g1<-g1+scale_y_continuous(limits=c(ymin,0))
  g1<-g1+scale_x_continuous(limits=c(-0.05,1.05),breaks=c(0,0.25,0.5,0.75,1),labels = c("0","0.25","0.5","0.75","1"))
  g1<-g1+xlab(expression(lambda))+ylab(paste0("likelihood/",gain))
  g1<-g1+theme(
    legend.background = element_rect(fill="#666666"),
    legend.key = element_rect(fill="#666666"),
    legend.text=element_text(color="white"),
    legend.title=element_blank(),
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(0, "mm"),
    legend.spacing.y = unit(0, "mm")
  )
  g1<-g1+plotTheme
  
  g2<-g2+scale_y_continuous(limits=c(ymin,0),breaks=c())
  g2<-g2+scale_x_continuous(limits=c(-0.05,1.05),breaks=c(0,0.25,0.5,0.75,1),labels = c("0","0.25","0.5","0.75","1"))
  g2<-g2+xlab("p(null)")+ylab(NULL)
  g2<-g2+theme(
    legend.background = element_rect(fill="#666666"),
    legend.key = element_rect(fill="#666666"),
    legend.text=element_text(color="white"),
    legend.title=element_blank(),
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(0, "mm"),
    legend.spacing.y = unit(0, "mm")
  )
  g2<-g2+plotTheme
  
  g<-g+annotation_custom(grob=ggplotGrob(g1+gridTheme),xmin=0.5,xmax=5.5,ymin=0.5,ymax=9.5)
  g<-g+annotation_custom(grob=ggplotGrob(g2+gridTheme),xmin=5.6,xmax=9.5,ymin=0.5,ymax=9.5)
  g
}
##########################################
doublePlot<-function(x1,xlb,y1,ylb1,y2=NULL,ylb2=NULL,y3=NULL,ylb3=NULL) {
  xtick<-NULL
  if (!is.numeric(x1)) {
    xtick<-x1
    x1<-1:length(x1)
  }
  if (!is.vector(y1)) y1<-as.vector(y1)
  if (!is.vector(y2)) y2<-as.vector(y2)
  if (length(x1)-length(y1)==1) {
    x1<-(x1[1:(length(x1)-1)]+x1[2:length(x1)])/2
  }
  
  g<-ggplot()+plotBlankTheme+theme(plot.margin=margin(0,-0.2,0,-1,"cm"))
  g<-g+scale_x_continuous(limits = c(0,10),labels=NULL,breaks=NULL)+scale_y_continuous(limits = c(0,10),labels=NULL,breaks=NULL)
  
  ylim1<-c(min(min(y1),0.1),max(max(y1)+0.1,0.6))
  g1<-ggplot()+scale_y_continuous(limits=ylim1)
  if (!is.null(xtick)) {
    g1<-g1+scale_x_continuous(limits=c(0,length(x1)+1),breaks=1:length(x1),labels=xtick)
    if (max(nchar(xtick))>2) {
      g1<-g1+theme(axis.text.x=element_text(angle = 90, hjust = 1))
    }
  } else {
    g1<-g1+scale_x_continuous(limits=c(min(x1),max(x1))+c(-1,1)*(max(x1)-min(x1))/4)
  }
  pts<-data.frame(x1=x1,y1=y1)
  if (is.null(xtick)) {
    g1<-g1+geom_line(data=pts,aes(x=x1,y=y1))
  } else {
    dx<-x1[2]-x1[1]
    for (i in 1:length(x1)) {
      ptsa<-data.frame(x1=x1[i]+dx*c(-1,1)*0.65,y1=y1[i]+c(0,0))
      g1<-g1+geom_line(data=ptsa,aes(x=x1,y=y1))
    }
    pts<-data.frame(x1=x1,y1=y1)
  }
  g1<-g1+geom_point(data=pts,aes(x=x1,y=y1),fill="yellow",size=4,shape=21)
  g1<-g1+xlab(xlb)+ylab(ylb1)
  nplot<-1
  
  if (!is.null(y2)) {
  ylim2<-c(min(min(y2),0.1),max(max(y2)+0.1,0.8))
  g2<-ggplot()+scale_y_continuous(limits=ylim2)
  if (!is.null(xtick)) {
    g2<-g2+scale_x_continuous(limits=c(0,length(x1)+1),breaks=1:length(x1),labels=xtick)
    if (max(nchar(xtick))>2) {
      g2<-g2+theme(axis.text.x=element_text(angle = 90, hjust = 1))
    }
  } else {
    g2<-g2+scale_x_continuous(limits=c(min(x1),max(x1))+c(-1,1)*(max(x1)-min(x1))/4)
  }
  pts<-data.frame(x2=x1,y2=y2)
  if (is.null(xtick)) {
    g2<-g2+geom_line(data=pts,aes(x=x2,y=y2))
  } else {
    dx<-x1[2]-x1[1]
    for (i in 1:length(x1)) {
      ptsa<-data.frame(x2=x1[i]+dx*c(-1,1)*0.65,y2=y2[i]+c(0,0))
      g2<-g2+geom_line(data=ptsa,aes(x=x2,y=y2))
    }
    pts<-data.frame(x2=x1,y2=y2)
  }
  g2<-g2+geom_point(data=pts,aes(x=x2,y=y2),fill="yellow",size=4,shape=21)
  g2<-g2+xlab(xlb)+ylab(ylb2)
  nplot<-2
  }
  
  if (!is.null(y3)) {
    ylim3<-c(min(min(y3),0.1),max(max(y3)+0.1,0.8))
    g3<-ggplot()+scale_y_continuous(limits=ylim3)
    if (!is.null(xtick)) {
      g3<-g3+scale_x_continuous(breaks=1:length(x1),labels=xtick)
      if (max(nchar(xtick))>2) {
        g3<-g3+theme(axis.text.x=element_text(angle = 90, hjust = 1))
      }
    }
    pts<-data.frame(x3=x1,y3=y3)
    g3<-g3+geom_line(data=pts,aes(x=x3,y=y3))
    g3<-g3+geom_point(data=pts,aes(x=x3,y=y3),fill="yellow",size=4,shape=21)
    g3<-g3+xlab(xlb)+ylab(ylb3)
    nplot<-3
  }
  
  if (nplot==1) {
    g<-g+annotation_custom(grob=ggplotGrob(g1+plotTheme+gridTheme),xmin=0.5,xmax=9.8,ymin=0.5,ymax=9.5)
  } 
  if (nplot==2) {
    g<-g+annotation_custom(grob=ggplotGrob(g1+plotTheme+gridTheme),xmin=0.5,xmax=4.9,ymin=0.5,ymax=9.5)
    g<-g+annotation_custom(grob=ggplotGrob(g2+plotTheme+gridTheme),xmin=5.5,xmax=9.9,ymin=0.5,ymax=9.5)
  } 
  if (nplot==3) {
    g<-g+annotation_custom(grob=ggplotGrob(g1+plotTheme+gridTheme),xmin=0.5,xmax=3.2,ymin=0.5,ymax=9.5)
    g<-g+annotation_custom(grob=ggplotGrob(g2+plotTheme+gridTheme),xmin=3.6,xmax=6.5,ymin=0.5,ymax=9.5)
    g<-g+annotation_custom(grob=ggplotGrob(g3+plotTheme+gridTheme),xmin=6.9,xmax=9.8,ymin=0.5,ymax=9.5)
  }
  g
}

