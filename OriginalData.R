####################################

# my_data <- read.table(pipe("pbpaste"), sep="\t", header = TRUE)
my_data_orig <- read_excel("fullfile.xlsx")
my_data<-my_data_orig
print(paste("ALL:",nrow(my_data)))

# remove errors
before<-nrow(my_data)
waste<-my_data$Error | is.na(my_data$Error)
my_data<-my_data[!waste,]
print(paste("Errors:",before-nrow(my_data),"=",nrow(my_data)))

# inspect df2
# remove values that are missing
before<-nrow(my_data)
waste<-grepl("[^0-9.]+",my_data$df2)
my_data<-my_data[!waste,]
my_data$df2<-as.numeric(my_data$df2)
waste<- is.na(my_data$df2)
my_data<-my_data[!waste,]
print(paste("df2 not available:",before-nrow(my_data),"=",nrow(my_data)))

# remove non-integer df2
before<-nrow(my_data)
waste<-my_data$df2!=round(my_data$df2)
my_data<-my_data[!waste,]
print(paste("df2 not integer",before-nrow(my_data),"=",nrow(my_data)))

# waste<-my_data$df2<2
# my_data<-my_data[!waste,]
# print(paste("df2<2",nrow(my_data)))

# inspect df1
# remove any that are missing
before<-nrow(my_data)
waste<-grepl("[^0-9]+",my_data$df1)
my_data<-my_data[!waste,]
my_data$df1<-as.numeric(my_data$df1)
waste<-is.na(my_data$df1)
my_data<-my_data[!waste,]
waste<-my_data$df1<1
my_data<-my_data[!waste,]
print(paste("df1 not available",before-nrow(my_data),"=",nrow(my_data)))

# remove any with df1>2
# these have a different sampling error distribution
before<-nrow(my_data)
waste<-my_data$df1>1
my_data<-my_data[!waste,]
print(paste("df1>1",before-nrow(my_data),"=",nrow(my_data)))

# remove any with ch2>df2
# these have an unknown effect size
before<-nrow(my_data)
waste<-(my_data$Statistic=="Chi2") & (my_data$Value>=my_data$df2)
my_data<-my_data[!waste,]
print(paste("chi2>df2",before-nrow(my_data),"=",nrow(my_data)))

# waste<-my_data$OneTail | my_data$OneTailedInTxt
# my_data<-my_data[!waste,]
# print(paste("1-tail",nrow(my_data)))

# calculate effect size
r_s<-rep(NA,nrow(my_data))
# r
use<-my_data$Statistic=="r"
r_s[use]<-my_data$Value[use]
# t
use<-my_data$Statistic=="t"
r_s[use]<-my_data$Value[use]/sqrt(my_data$Value[use]^2+my_data$df2[use])
# F
use<-my_data$Statistic=="F"
r_s[use]<-sqrt(my_data$Value[use]*my_data$df1[use]/(my_data$Value[use]*my_data$df1[use]+my_data$df2[use]))
# chi2
use<-my_data$Statistic=="Chi2" & my_data$df1==1
r_s[use]<-sqrt(my_data$Value[use]/(my_data$df2[use]))
my_data$df2[use]<-my_data$df2[use]+my_data$df1[use]+1

my_data$r_s<-r_s

# remove any r values that are not valid
before<-nrow(my_data)
waste<-is.na(my_data$r_s)
my_data<-my_data[!waste,]
print(paste("r not available",before-nrow(my_data),"=",nrow(my_data)))
# r==1
before<-nrow(my_data)
waste<-abs(my_data$r_s)>=1
my_data<-my_data[!waste,]
print(paste("r==1",before-nrow(my_data),"=",nrow(my_data)))

my_data$n<-my_data$df2+my_data$df1+1
my_data$p<-r2p(my_data$r_s,my_data$n,my_data$df1)

before<-nrow(my_data)
waste<-my_data$n<=3
my_data<-my_data[!waste,]
print(paste("n<4",before-nrow(my_data),"=",nrow(my_data)))

# separate out the 1-tails that have p between 0.05 and 0.1 (2-tailed p)
keep<-(my_data$OneTail | my_data$OneTailedInTxt) & (my_data$p>=0.05 & my_data$p<0.1)
my_1t_data<-my_data[keep,]

# now remove non-sig results
before<-nrow(my_data)
waste<-(my_data$p>alpha)
my_data<-my_data[!waste,]
print(paste("p>0.05",before-nrow(my_data),"=",nrow(my_data)))

my_data$z_s<-atanh(my_data$r_s)


