maincolours<-list(windowC="#002D40",panelC="#005E86",graphC="#BFECFF")

mainTheme=theme(panel.background = element_rect(fill="#666666", colour="black"),
                panel.grid.major = element_line(linetype="blank"),panel.grid.minor = element_line(linetype="blank"),
                plot.background = element_rect(fill=maincolours$graphC, colour=maincolours$graphC),
                plot.margin = margin(2,2,1,1,"cm"))
SMplotTheme=theme(plot.title=element_text(size=16,face="bold"),
                  axis.title=element_text(size=16,face="bold"),
                  axis.text=element_text(size=10,face="bold"))
plotTheme=mainTheme+SMplotTheme

plotBlankTheme=theme(panel.background = element_rect(fill=maincolours$graphC, colour=maincolours$graphC),
                     panel.grid.major = element_line(linetype="blank"),panel.grid.minor = element_line(linetype="blank"),
                     plot.background = element_rect(fill=maincolours$graphC, colour=maincolours$graphC),
                     axis.title=element_text(size=16,face="bold")
)
gridTheme=theme(plot.margin=margin(0,0,0,0,"cm"))



##########################################
showAnalysis<-function(an,label=NULL,param=NULL) {
  
  if (is.null(param)) {
    param<-an$metaAnalysis$meta_pdf
    if (param[1]=="All") {
      param<-c("Best")
    }
    for (i in 1:length(param)) {
      s<-param[i]
      switch(param[i],
             "Single"={
               d<-an$single
             },
             "Exp"={
               d<-an$exp
             },
             "Gauss"={
               d<-an$gauss
             },
             "Gamma"={
               d<-an$gamma
             },
             "Best"={
               s<-an$bestDist
               d<-an$best
             }
      )
      if (Pplus) {
        print(paste0(label,s,
                     "(",
                     "P+","=",format((1-d$Nullmax)*100,digits=2),"%",
                     "  k=",format(d$Kmax,digits=3),
                     "  S=",format(d$Smax,digits=3),
                     "  n=",format(length(an$result$rIV),digits=2),
                     ")"))
      } else {
        print(paste0(label,s,
                     "(",
                     "P-=",format(d$Nullmax*100,digits=2),"%",
                     "  k=",format(d$Kmax,digits=3),
                     "  S=",format(d$Smax,digits=3),
                     "  n=",format(length(an$result$rIV),digits=2),
                     ")"))
      }
      
    }
  } else {
    s<-paste0(param,"=")
    if (param=="All") {
      s<-paste0(label,": ",s,an$bestDist)
      
    }
    if (param=="Best") {
      s<-paste0(label,": ",s,an$bestDist)
    }
    
    if (Pplus) {
      print(paste0(s,
                   "(",
                   "P+","=",format((1-an$best$Nullmax)*100,digits=2),"%",
                   "  k=",format(an$best$Kmax,digits=3),
                   "  S=",format(an$best$Smax,digits=3),
                   "  n=",format(length(an$result$rIV),digits=2),
                   ")"))
    } else {
      print(paste0(s,
                   "(",
                   "P-=",format(an$best$Nullmax*100,digits=2),"%",
                   "  k=",format(an$best$Kmax,digits=3),
                   "  S=",format(an$best$Smax,digits=3),
                   "  n=",format(length(an$result$rIV),digits=2),
                   ")"))
    }
  }
}

##########################################
drawAnalysis<-function(an1,metaData1) {
  
  gain<-1000
  nStudies<-length(metaData1$result$rIV)
  g<-ggplot()+plotBlankTheme+theme(plot.margin=margin(0,-0.2,0,0,"cm"))
  g<-g+scale_x_continuous(limits = c(0,10),labels=NULL,breaks=NULL)
  g<-g+scale_y_continuous(limits = c(0,10),labels=NULL,breaks=NULL)
  
  g1<-ggplot()
  g2<-ggplot()
  
  ymax<- -Inf
  
  drawType<-an1$metaAnalysis$meta_pdf
  if (drawType[1]=="All")
    drawType<-c("Exp","Gauss")
  
  for (i in 1:length(drawType)) {
    switch(drawType[i],
           "Single"={
             d<-an1$single
           },
           "Exp"={
             d<-an1$exp
           },
           "Gauss"={
             d<-an1$gauss
           },
           "Gamma"={
             d<-an1$gamma
           }
    )
    
    k<-d$SX$kvals
    Sk1<-d$SX$SkX
    pts1<-data.frame(x=k,y=Sk1)
    
    nullP<-d$SX$nullvals
    Sn1<-d$SX$SnullX
    if (Pplus) {
      pts2<-data.frame(x=1-nullP,y=Sn1)
    } else {
      pts2<-data.frame(x=nullP,y=Sn1)
    }
    
    switch(drawType[i], # serious pain in the neck here: whats wrong with col=drawtype[i]?
           "Single"={
             g1<-g1+geom_line(data=pts1,aes(x=x,y=y,col="Single"),lwd=1)
             g2<-g2+geom_line(data=pts2,aes(x=x,y=y,col="Single"),lwd=1)
           },
           "Exp"={
             g1<-g1+geom_line(data=pts1,aes(x=x,y=y,col="Exp"),lwd=1)
             g2<-g2+geom_line(data=pts2,aes(x=x,y=y,col="Exp"),lwd=1)
           },
           "Gauss"={
             g1<-g1+geom_line(data=pts1,aes(x=x,y=y,col="Gauss"),lwd=1)
             g2<-g2+geom_line(data=pts2,aes(x=x,y=y,col="Gauss"),lwd=1)
           },
           "Gamma"={
             g1<-g1+geom_line(data=pts1,aes(x=x,y=y,col="Gamma"),lwd=1)
             g2<-g2+geom_line(data=pts2,aes(x=x,y=y,col="Gamma"),lwd=1)
           }
    )
    ymax<- max(c(Sk1,ymax))
  }
  
  g1<-g1+scale_color_manual(name = NULL, values = c("Single"="green","Exp" = "red", "Gauss" = "yellow","Gamma"="blue"))
  g2<-g2+scale_color_manual(name = NULL, values = c("Single"="green","Exp" = "red", "Gauss" = "yellow","Gamma"="blue"))
  
  ymin<-ymax-nStudies
  ymax<-ymax+nStudies/5
  
  g1<-g1+scale_y_continuous(limits=c(ymin,ymax))
  g1<-g1+scale_x_continuous(limits=c(-0.05,1.05),breaks=c(0,0.25,0.5,0.75,1),labels = c("0","0.25","0.5","0.75","1"))
  g1<-g1+xlab(Llabel)+ylab(Slabel)
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
  
  g2<-g2+scale_y_continuous(limits=c(ymin,ymax),breaks=c())
  g2<-g2+scale_x_continuous(limits=c(-0.05,1.05),breaks=c(0,0.25,0.5,0.75,1),labels = c("0","0.25","0.5","0.75","1"))
  g2<-g2+xlab(Plabel)+ylab(NULL)
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
doublePlot<-function(x1,xlb,y1,ylb1,y2=NULL,ylb2=NULL,y3=NULL,ylb3=NULL,xtick=NULL,xlog=FALSE,legendLabels=c("real","sim"),legendX="left",legendY="bottom") {
  
  if (!is.numeric(x1)) {
    xtick<-x1
    x1<-1:length(x1)
  }
  if (length(x1)-length(y1)==1) {
    x1<-(x1[1:(length(x1)-1)]+x1[2:length(x1)])/2
  }
  
  g<-ggplot()+plotBlankTheme+theme(plot.margin=margin(0,-0.2,0,-1,"cm"))
  g<-g+scale_x_continuous(limits = c(0,10),labels=NULL,breaks=NULL)+scale_y_continuous(limits = c(0,10),labels=NULL,breaks=NULL)
  
  ylim1<-c(min(y1),max(y1))+c(-1,1)*(max(y1)-min(y1))*0.1
  if (ylb1==Llabel) {
    ylim1<-c(0,max(max(y1,na.rm=TRUE)+0.1,0.6))
  }
  # if (ylb1=="w") {
  #   ylim1<-c(0,1)
  # }
  # if (ylb1=="fdr") {
  #   ylim1<-c(0,1)
  # }
  g1<-ggplot()+scale_y_continuous(limits=ylim1)
  xlim<-c(min(x1),max(x1))
  if (!is.null(xtick)) {
    if (xlog) {
      xl<-log10(x1)
      g1<-g1+scale_x_log10(limits=10^(c(min(xl),max(xl))+c(-1,1)*(max(xl)-min(xl))/4))
    } else {
      g1<-g1+scale_x_continuous(limits=c(0,length(x1)+1),breaks=1:length(x1),labels=xtick)
    }
    if (!isempty(xtick) && max(nchar(xtick))>2) {
      g1<-g1+theme(axis.text.x=element_text(angle = 90, hjust = 1))
    }
  } else {
    if (xlog) {
      xl<-log10(x1)
      g1<-g1+scale_x_log10(limits=10^(c(min(xl),max(xl))+c(-1,1)*(max(xl)-min(xl))/4))
    } else {
      g1<-g1+scale_x_continuous(limits=c(min(x1),max(x1))+c(-1,1)*(max(x1)-min(x1))/4)
    }
  }
  cols<-c("yellow","red")
  for (j in 1:nrow(y1)) {
    pts<-data.frame(x1=x1,y1=y1[j,])
    if (is.numeric(x1)) {
      g1<-g1+geom_line(data=pts,aes(x=x1,y=y1))
    } else {
      dx<-x1[2]-x1[1]
      for (i in 1:length(x1)) {
        ptsa<-data.frame(x1=x1[i]+dx*c(-1,1)*0.65,y1=y1[j,i]+c(0,0))
        g1<-g1+geom_line(data=ptsa,aes(x=x1,y=y1))
      }
    }
    g1<-g1+geom_point(data=pts,aes(x=x1,y=y1),fill=cols[j],size=4,shape=21)
  }
  
  if (nrow(y1)>1) {
    if (legendX=="left") {
      if (xlog) {
        xp<-10^(log10(xlim[1])+diff(log10(xlim))*0.01)
        xp1<-10^(log10(xlim[1])+diff(log10(xlim))*0.05)
      } else {
        xp<-xlim[1]+diff(xlim)*0.01
        xp1<-xlim[1]+diff(xlim)*0.05
      }
    } else {
      if (xlog) {
        xp<-10^(log10(xlim[2])-diff(log10(xlim))*0.01)
        xp1<-10^(log10(xlim[2])-diff(log10(xlim))*0.05)
      } else {
        xp<-xlim[2]-diff(xlim)*0.01
        xp1<-xlim[2]-diff(xlim)*0.05
      }
    }
    if (legendY=="bottom") {
      yp<-ylim1[1]+diff(ylim1)*0.2
      yp1<-ylim1[1]+diff(ylim1)*0.1
    } else {
      yp<-ylim1[2]-diff(ylim1)*0.2
      yp1<-ylim1[2]-diff(ylim1)*0.1
    }
    legend<-data.frame(x=xp,y=yp)
    g1<-g1+geom_point(data=legend,aes(x=x,y=y),fill=cols[1],size=4,shape=21)
    legendText<-data.frame(x=xp1,y=yp,label=legendLabels[1])
    g1<-g1+geom_text(data=legendText,aes(x=x,y=y,label=label),color="white",hjust=-1)
    
    legend<-data.frame(x=xp,y=yp1)
    g1<-g1+geom_point(data=legend,aes(x=x,y=y),fill=cols[2],size=4,shape=21)
    legendText<-data.frame(x=xp1,y=yp1,label=legendLabels[2])
    g1<-g1+geom_text(data=legendText,aes(x=x,y=y,label=label),color="white",hjust=-1)
  }
  g1<-g1+xlab(xlb)+ylab(ylb1)
  nplot<-1
  
  if (!is.null(y2)) {
    if (Pplus && ylb2==Plabel) {
      y2<-1-y2
    }
    ylim2<-c(min(y2),max(y2))+c(-1,1)*(max(y2)-min(y2))*0.1
    # ylim2<-c(0,1)
    g2<-ggplot()+scale_y_continuous(limits=ylim2)
    if (!is.null(xtick)) {
      if (xlog) {
        xl<-log10(x1)
        g2<-g2+scale_x_log10(limits=10^(c(min(xl),max(xl))+c(-1,1)*(max(xl)-min(xl))/4))
      } else {
        g2<-g2+scale_x_continuous(limits=c(0,length(x1)+1),breaks=1:length(x1),labels=xtick)
      }
      if (max(nchar(xtick))>2) {
        g2<-g2+theme(axis.text.x=element_text(angle = 90, hjust = 1))
      }
      if (max(nchar(xtick))>2) {
        g2<-g2+theme(axis.text.x=element_text(angle = 90, hjust = 1))
      }
    } else {
      if (xlog) {
        xl<-log10(x1)
        g2<-g2+scale_x_log10(limits=10^(c(min(xl),max(xl))+c(-1,1)*(max(xl)-min(xl))/4))
      } else {
        g2<-g2+scale_x_continuous(limits=c(min(x1),max(x1))+c(-1,1)*(max(x1)-min(x1))/4)
      }
    }
    for (j in 1:nrow(y2)) {
      pts<-data.frame(x2=x1,y2=y2[j,])
      if (is.numeric(x1)) {
        g2<-g2+geom_line(data=pts,aes(x=x2,y=y2))
      } else {
        dx<-x1[2]-x1[1]
        for (i in 1:length(x1)) {
          ptsa<-data.frame(x2=x1[i]+dx*c(-1,1)*0.65,y2=y2[j,i]+c(0,0))
          g2<-g2+geom_line(data=ptsa,aes(x=x2,y=y2))
        }
      }
      g2<-g2+geom_point(data=pts,aes(x=x2,y=y2),fill=cols[j],size=4,shape=21)
    }
    if (nrow(y2)>1) {
      legend<-data.frame(x=xp,y=ylim2[1]+diff(ylim2)*0.2)
      g2<-g2+geom_point(data=legend,aes(x=x,y=y),fill=cols[1],size=4,shape=21)
      legendText<-data.frame(x=xp1,y=ylim2[1]+diff(ylim2)*0.2,label=legendLabels[1])
      g2<-g2+geom_text(data=legendText,aes(x=x,y=y,label=label),color="white",hjust=-1)
      
      legend<-data.frame(x=xp,y=ylim2[1]+diff(ylim2)*0.1)
      g2<-g2+geom_point(data=legend,aes(x=x,y=y),fill=cols[2],size=4,shape=21)
      legendText<-data.frame(x=xp1,y=ylim2[1]+diff(ylim2)*0.1,label=legendLabels[2])
      g2<-g2+geom_text(data=legendText,aes(x=x,y=y,label=label),color="white",hjust=-1)
      
    }
    g2<-g2+xlab(xlb)+ylab(ylb2)
    nplot<-2
  }
  
  if (!is.null(y3)) {
    ylim3<-c(min(y3),max(y3))+c(-1,1)*(max(y3)-min(y3))/5
    g3<-ggplot()+scale_y_continuous(limits=ylim3)
    if (!is.null(xtick)) {
      g3<-g3+scale_x_continuous(limits=c(0,length(x1)+1),breaks=1:length(x1),labels=xtick)
      if (max(nchar(xtick))>2) {
        g3<-g3+theme(axis.text.x=element_text(angle = 90, hjust = 1))
      }
    } else {
      g3<-g3+scale_x_continuous(limits=c(min(x1),max(x1))+c(-1,1)*(max(x1)-min(x1))/4)
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

