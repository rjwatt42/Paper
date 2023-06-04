####################################
# copy data to clipboard first

# my_data <- read.table(pipe("pbpaste"), sep="\t", header = TRUE)
my_data <- read_excel("fullfile.xlsx")
print(paste("ALL",nrow(my_data)))

waste<-my_data$Error | is.na(my_data$Error)
my_data<-my_data[!waste,]
print(paste("Errors",nrow(my_data)))

waste<-grepl("[^0-9]+",my_data$df2)
my_data<-my_data[!waste,]
my_data$df2<-as.numeric(my_data$df2)
waste<-my_data$df2!=round(my_data$df2)  | is.na(my_data$df2)
my_data<-my_data[!waste,]
print(paste("df2",nrow(my_data)))

waste<-my_data$df2<8
my_data<-my_data[!waste,]
print(paste("df2<8",nrow(my_data)))

waste<-my_data$df1>1  | is.na(my_data$df1)
my_data<-my_data[!waste,]
print(paste("df1>1",nrow(my_data)))

waste<-my_data$OneTail | my_data$OneTailedInTxt
my_data<-my_data[!waste,]
print(paste("1-tail",nrow(my_data)))

# calculate effect size
my_data$r_s<-NA
use<-my_data$Statistic=="r"
my_data$r_s[use]<-my_data$Value[use]
use<-my_data$Statistic=="t"
my_data$r_s[use]<-my_data$Value[use]/sqrt(my_data$Value[use]^2+my_data$df2[use])
use<-my_data$Statistic=="F"
my_data$r_s[use]<-sqrt(my_data$Value[use]/(my_data$Value[use]+my_data$df2[use]))
use<-my_data$Statistic=="Chi2"
my_data$r_s[use]<-sqrt(my_data$Value[use]/(my_data$df2[use]+2))

waste<-abs(my_data$r_s)>=1  | is.na(my_data$r_s)
my_data<-my_data[!waste,]
print(paste("r error",nrow(my_data)))

my_data$n<-my_data$df2+2

# remove non-sig results
my_data$p<-r2p(my_data$r_s,my_data$n)
waste<-(my_data$p>alpha)
my_data<-my_data[!waste,]
print(paste("p>0.05",nrow(my_data)))



