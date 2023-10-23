#USACE analysis

#USACE_1 regression analysis 

setwd('~/Documents/psu/ESM 566/snowpack-project')

#libraries
library(tidyr)
library(lubridate)
library(plyr)
library(zoo)
library(car)
library(MASS)
library(leaps)
##########################################
#import cleaned data 
dta<-read.csv('DATA/USACE_1_02082023')

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
# lamba = 500
b<-boxcox(lm(melt12bc$Temp_Avg ~ 1), lambda = seq(1,4,.1))
#lambds = 1.5
b<-boxcox(lm(melt12bc$RH ~ 1), lambda = seq(0,5,.1))
#lamba=1.5
b<-boxcox(lm(melt12bc$WS_ms ~ 1), lambda = seq(0,2,0.1)) 
#lambda = 0.5
b<-boxcox(lm(melt12bc$delt_precip ~ 1), lambda = seq(-100,-2,0.1))
#lambda = -30 
b<-boxcox(lm(melt12bc$netsw ~ 1), lambda=seq(-10,10,.1))
#lambda = -8

#apply box cox transform from previous charts to datafram 
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
summary(mod12.t)
# Call:
#   lm(formula = depth_inc ~ ., data = melt12bc)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.58599 -0.15640 -0.01343  0.16280  0.54522 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.223e+00  1.955e-01   6.256 5.28e-09 ***
#   Temp_Avg    -1.414e-03  8.169e-04  -1.731 0.085768 .  
# RH          -3.113e-04  8.723e-05  -3.568 0.000504 ***
#   WS_ms       -2.302e-01  6.505e-02  -3.538 0.000560 ***
#   netsw       -2.617e+16  1.759e+16  -1.487 0.139348    
# delt_precip  7.732e-02  7.875e-02   0.982 0.328004    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2421 on 130 degrees of freedom
# Multiple R-squared:  0.2755,	Adjusted R-squared:  0.2476 
# F-statistic: 9.887 on 5 and 130 DF,  p-value: 4.972e-08
par(mfrow=c(2,2))
plot(mod12.t) #view results 
shapiro.test(residuals(mod12.t)) #p=0.54
vif(mod12.t) # does not suggest multicollinearity

#check equal variance 
m<-median(predict(mod12.t)) #find a median for fitted values (Y-hat)
g1<-residuals(mod12.t)[predict(mod12.t)>m]
g2<-residuals(mod12.t)[predict(mod12.t)<m]
var.test(g1,g2) #p value = 0.2271 - good 

#passes all test ok to reduce the model 

############################################
#reduce transformed model 
#p to reject = 0.01
colnames(melt12bc)

#remove delta precip (p=0.32)
mod12.t2<-lm(depth_inc~. , data=melt12bc[,-6])
summary(mod12.t2)
#remove netsw (p=0.15)
mod12.t3<-lm(depth_inc~. , data=melt12bc[,-c(6,5)])
summary(mod12.t3)
#remove temperature (p=0.0342)
mod12.t4<-lm(depth_inc~. , data=melt12bc[,-c(2,6,5)])
summary(mod12.t4) # RH and windspeed are significant 

# Call:
#   lm(formula = depth_inc ~ ., data = melt12bc[, -c(2, 6, 5)])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.58775 -0.19123 -0.00217  0.18057  0.55545 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.193e+00  1.121e-01  10.641  < 2e-16 ***
#   RH          -3.169e-04  6.685e-05  -4.740 5.42e-06 ***
#   WS_ms       -2.914e-01  6.147e-02  -4.741 5.39e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2463 on 133 degrees of freedom
# Multiple R-squared:  0.2329,	Adjusted R-squared:  0.2213 
# F-statistic: 20.19 on 2 and 133 DF,  p-value: 2.208e-08

par(mfrow=c(2,2))
plot(mod12.t4) #view results 
shapiro.test(residuals(mod12.t4)) #(p=0.625 - good) 
vif(mod12.t4)
m<-median(predict(mod12.t4)) #find a median for fitted values (Y-hat)
(g1<-residuals(mod12.t4)[predict(mod12.t4)>m]) 
(g2<-residuals(mod12.t4)[predict(mod12.t4)<m]) 
var.test(g1,g2) #p value = 0.3026 - good 

#############################################
#perform anova to see if reduced model is as good as original model 
anova(mod12.t,mod12.t4)
# Model 1: depth_inc ~ Temp_Avg + RH + WS_ms + netsw + delt_precip
# Model 2: depth_inc ~ RH + WS_ms
# Res.Df    RSS Df Sum of Sq      F  Pr(>F)  
# 1    130 7.6180                              
# 2    133 8.0663 -3  -0.44828 2.5499 0.05855 .
#p>0.05 indicates that the models are not statistically different 

#############################################
#reduce using different method 
rst<-summary(regsubsets(depth_inc~., data=melt12bc))
rst
# Selection Algorithm: exhaustive
# Temp_Avg RH  WS_ms netsw delt_precip
# 1  ( 1 ) "*"      " " " "   " "   " "        
# 2  ( 1 ) "*"      "*" " "   " "   " "        
# 3  ( 1 ) "*"      "*" "*"   " "   " "        
# 4  ( 1 ) "*"      "*" "*"   "*"   " "        
# 5  ( 1 ) "*"      "*" "*"   "*"   "*" 

#plot 
par(mfrow=c(2,2))
plot(1:(ncol(melt12bc)-1),rst$cp, xlab="# of predictors in the model", ylab="Mallows Cp")
abline(v=3,col="red", lwd=2)
plot(1:(ncol(melt12bc)-1),rst$adjr2, xlab="# of predictors in the model", ylab="Adjusted R2")
abline(v=3,col="red", lwd=2)
plot(1:(ncol(melt12bc)-1),rst$bic, xlab="# of predictors in the model", ylab="BIC")
abline(v=3,col="red", lwd=2)

##########################################
#check relationship between Temp & RH 
colnames(melt12bc)
cor.matrix(melt12bc[,c(3,2)])
mod.trh<-lm(RH~Temp_Avg, data=melt12bc)
summary(mod.trh) #significant 
# Call:
#   lm(formula = RH ~ Temp_Avg, data = melt12bc)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -817.01 -190.88   49.83  202.10  553.01 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 1071.6702    78.7123  13.615  < 2e-16 ***
#   Temp_Avg      -4.5729     0.7847  -5.827 3.99e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 286.1 on 134 degrees of freedom
# Multiple R-squared:  0.2022,	Adjusted R-squared:  0.1962 
# F-statistic: 33.96 on 1 and 134 DF,  p-value: 3.985e-08
