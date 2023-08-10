
library(mnormt)      # pmnorm for logistic
library(lme4)        # lmer (mixed models)
library(MuMIn)       # r-squared for mixed models
library(readxl)      # excel
library(writexl)     # x excel
library(car)         # Anova type 3 correct
library(stringr)     # for str_* functions
library(clipr)       # for clipboard functions
library(SuppDists)   # for Johnson distributions
library(e1071)       # for skewness and kurtosis
library(pracma)      # for meshgridlibrary(EnvStats)
library(EnvStats)    # for egamma

library(ggplot2)
library(ggimage)
library(plotly)
library(reshape)


source("runMetaAnalysis.R")
source("sampleLikelihood.R")
source("r2p.R")
source("Plots.R")
source("makeStudies.R")

LabelUD<-"D"
Pplus<-FALSE
Pchar<-"P" 
# Pchar<-'\u03A9'
if (Pplus) {Ptypechar<-'+' } else {Ptypechar<-0 } #'\u2013'

Lchar<-'\u03BB'
Ltypechar<-"+"
Llabel<-Lchar

switch (LabelUD, 
        "U"={
          Plabel<-bquote(.(Pchar)^.(Ptypechar))
          Llabel<-bquote(.(Lchar)^.(Ltypechar))
        },
        "D"={
          Plabel<-bquote(.(Pchar)[.(Ptypechar)])
          Llabel<-bquote(.(Lchar)[.(Ltypechar)])
        })


Slabel<-"log(likelihood)"
