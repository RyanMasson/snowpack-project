#USACE_2 regression analysis

setwd('~/Documents/psu/ESM 566/snowpack-project')

library(tidyr)
library(lubridate)
library(plyr)
library(zoo)
library(car)
library(MASS)
library(leaps)

##########################################
#import cleaned data 
dta<-read.csv('DATA/USACE_2_02082023')

###############################################
#SNOTEL Data
snotel<-read.csv('MarionForks.csv') #import SNOtel data from Marion Forks 

#clean the dates of the snotel data 
snotel$date<-parse_date_time(snotel$Date, orders=c("Ymd HMS","mdY HM")) 
snotel$year<-as.numeric(format(snotel$date, "%Y"))
snotel$month<-as.numeric(format(snotel$date, "%m"))
snotel$day<-as.numeric(format(snotel$date, "%d"))
snotel$hour<-as.numeric(format(snotel$date, "%H"))
snotel$hour<-ifelse(snotel$hour>=10,snotel$hour,paste0("0",snotel$hour))
snotel$doy<-yday(snotel$date)
snotel$dowy<-ifelse(snotel$doy>247,snotel$doy-247,snotel$doy+118) #October 1 is first day of water year 
snotel<-unite(snotel, col = "dowyh", sep = ".",c('dowy','hour'),remove=FALSE) #combine day of water yea rand hour into one date column

#select only columns needed for analysis from SNOTEL Data 
projsnotel<-snotel[,c(8,2,3)]
#add a precipitation and SWE increment column 
delta_precip<-vector(mode='numeric',length=length(projsnotel$Precipitation.Accumulation..in.)) #set up vector to store results
for (i in 2:length(projsnotel$Precipitation.Accumulation..in.)){
  delta_precip[i]<-projsnotel[i,3]-projsnotel[i-1,3]}
delta_swe<-vector(mode='numeric',length=length(projsnotel$Snow.Water.Equivalent..in.))
for (i in 2:length(projsnotel$Snow.Water.Equivalent..in.)){
  delta_swe[i]<-projsnotel[i,3]-projsnotel[i-1,3]}
projsnotel$delt_precip<- delta_precip #add to main dataset


#select only columns We need for ths analysis 
colnames(dta)
dta1<-dta[,c(1,46,8:11,26,32,34,36,52:60)] 

colnames(dta1)

#turn day of water year and hour decimal into number 
dta1$dowyh<-as.numeric(dta1$dowyh)
projsnotel$dowyh<-as.numeric(projsnotel$dowyh)

#change DOWYH to "date"
colnames(dta1)[2]<-"date"
colnames(projsnotel)[1]<-"date"

#add snotel data to the project data frame 
projdta<-merge(dta1,projsnotel, by= 'date')
colnames(projdta)

#clean net SW data 
projdta$netsw<-na.locf(ifelse(projdta$netsw>=0, projdta$netsw, NA )) 
summary(projdta$netsw)

#########################################################
#create 12 hour binned data 
bins<-cut(projdta$date, breaks = seq(75.15,157.12, 0.12))
dta.12hr<-aggregate(projdta , by=list(bins), mean)

########################################################
#make 3 data frames with variables that will be used in multiple linear regression analysis 
reg12<-dta.12hr[,c(16,4,5,8,12,23)]

###############################################
#Cor Matrix 
source("cor.matrix.r")
cor.matrix(reg12)

###########################################
#subset for melt periods 
#melt12<-subset(reg12, depth_inc < 0 & depth_inc > -0.01)
melt12<-subset(reg12, depth_inc < 0 & depth_inc > -0.03) #changed 3/6


############################################
#cor Matrix of melt periods 
#cor.matrix(reg12)
cor.matrix(melt12)

#################################################################
###################################################################
#use box cox to transform data 
#use box-cox to figure out best transformation for original data 

#transform data to make all variables positive to use the box cox function 
melt12bc<-melt12
colnames(melt12bc)

melt12bc$depth_inc<-melt12bc$depth_inc+1 #make response variable positive 
melt12bc$Temp_Avg<-melt12bc$Temp_Avg+20
melt12bc$RH<-melt12bc$RH
melt12bc$delt_precip<-melt12bc$delt_precip+1
melt12bc$netsw<-melt12bc$netsw+150
melt12bc$WS_ms<-melt12bc$WS_ms+1

#perform box cox on all variables 
par(mfrow=c(1,1))
b<-boxcox(lm(melt12bc$depth_inc ~ 1), lambda = seq(-0,2000,.1))
# lambda = 500
b<-boxcox(lm(melt12bc$Temp_Avg ~ 1), lambda = seq(1,4,.1))
#lambda = 1.5
b<-boxcox(lm(melt12bc$RH ~ 1), lambda = seq(0,5,.1))
#lambda=1.5
b<-boxcox(lm(melt12bc$WS_ms ~ 1), lambda = seq(0,2,0.1)) 
#lambda = 0.5
b<-boxcox(lm(melt12bc$delt_precip ~ 1), lambda = seq(-100,-2,0.1))
#lambda = -30 
b<-boxcox(lm(melt12bc$netsw ~ 1), lambda=seq(-10,10,.1))
#lambda = -8

#apply box cox transform from previous charts to dataframe  
melt12bc$depth_inc<-melt12bc$depth_inc^500 #make response variable positive 
melt12bc$Temp_Avg<-melt12bc$Temp_Avg^1.5
melt12bc$RH<-melt12bc$RH^1.5
melt12bc$WS_ms<-melt12bc$WS_ms^0.5
melt12bc$delt_precip<-melt12bc$delt_precip^-30
melt12bc$netsw<-melt12bc$netsw^-8

#######################################################
#create full model with transformed data 
cor.matrix(melt12bc)
mod12.t<-lm(depth_inc~., data=melt12bc)
summary(mod12.t) # no significant predictors

