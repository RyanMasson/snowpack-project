# NOT UPDATED FOR RYAN'S MACHINE, DOES NOT HAVE VALIDATION, DOES NOT HAVE NEW BINS

#USACE ANalysis 

#USACE_2 regression analysis 

#libraries

library(tidyr)
library(lubridate)
library(plyr)
library(zoo)
library(car)
library(MASS)
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

#add snotel data to the project data fram 
projdta<-merge(dta1,projsnotel, by= 'date')
colnames(projdta)

#clean net SW data 
projdta$netsw<-na.locf(ifelse(projdta$netsw>=0, projdta$netsw, NA )) 
summary(projdta$netsw)

#########################################################
#create 24 hour binned data 
bins<-cut(projdta$date, breaks = seq(75.15,157.12, 1))
dta.24hr<-aggregate(projdta , by=list(bins), mean)

# sample 80 % of data for holdout validation
spl = sample.split(dta.24hr$depth_inc, SplitRatio = 0.8)
train = subset(dta.24hr, spl==TRUE)
test = subset(dta.24hr, spl==FALSE)
print(dim(train)); print(dim(test))
dta.24hr <- train  # assigning dta.12hr to train to not have to change variable name below

########################################################
#make some plots 
colnames(dta.24hr)
plot(dta.24hr$date, dta.24hr$Temp_Avg, type='l', col='blue',xlab='Day of Water Year', ylab='Temp (deg C)')
abline(h=0)
par(new=TRUE)
plot(dta.24hr$date, dta.24hr$RH, col='cyan', ylab="")
axis(side=4)
mtext("Relative Humidy (%)", side = 4)

colnames(dta.24hr)
plot(dta.24hr$depth_inc, dta.24hr$WS_ms,col='blue',xlab='Depth Change (m)', ylab='Wind Speed (m/s)', xlim=c(-0.005,0.0075))
abline(h=0)
par(new=TRUE)
plot(dta.24hr$date, dta.24hr$RH, col='cyan', ylab="")
axis(side=4)
mtext("Relative Humidy (%)", side = 4)

########################################################
#make 3 data frames with variables that will be used in multiple linear regression analysis 
reg24<-dta.24hr[,c(16,4,5,8,12,23)]

###############################################
#Cor Matrix 
source("cor.matrix.r")
cor.matrix(reg24)

###########################################
#subset for melt periods 
melt24<-subset(reg24, depth_inc <0 & depth_inc > -0.03)

############################################
#cor Matrix of melt periods 
cor.matrix(melt24)

#################################################################
###################################################################
#use box cox to transform data 
#use box-cox to figure out best transformation for original data 

#transform data to make all variables postitive to use the box cox function 
melt24bc<-melt24
colnames(melt24bc)

melt24bc$depth_inc<-melt24bc$depth_inc+1 #make response variable positive 
melt24bc$Temp_Avg<-melt24bc$Temp_Avg+20
melt24bc$RH<-melt24bc$RH
melt24bc$delt_precip<-melt24bc$delt_precip+1
melt24bc$netsw<-melt24bc$netsw+150
melt24bc$WS_ms<-melt24bc$WS_ms+1

#perform box cox on all variables 
par(mfrow=c(1,1))
b<-boxcox(lm(melt24bc$depth_inc ~ 1), lambda = seq(-0,2000,.1))
# lamba = 500
b<-boxcox(lm(melt24bc$Temp_Avg ~ 1), lambda = seq(-4,8,.1))
#lambds = 1.5
b<-boxcox(lm(melt24bc$RH ~ 1), lambda = seq(0,5,.1))
#lamba=1
b<-boxcox(lm(melt24bc$WS_ms ~ 1), lambda = seq(0,2,0.1)) 
#lambda = 1
b<-boxcox(lm(melt24bc$delt_precip ~ 1), lambda = seq(-200,-2,0.1))
#lambda = -80 
b<-boxcox(lm(melt24bc$netsw ~ 1), lambda=seq(-10,10,.1))
#lambda = -5

#apply box cox transform from previous charts to datafram 
melt24bc$depth_inc<-melt24bc$depth_inc^500 #make response variable positive 
melt24bc$Temp_Avg<-melt24bc$Temp_Avg^1.5
melt24bc$RH<-melt24bc$RH
melt24bc$WS_ms<-melt24bc$WS_ms
melt24bc$delt_precip<-melt24bc$delt_precip^-80
melt24bc$netsw<-melt24bc$netsw^-5

#######################################################
#create full model with transformed data 
cor.matrix(melt24bc)
mod24.t<-lm(depth_inc~., data=melt24bc)
summary(mod24.t)
# Call:
#   lm(formula = depth_inc ~ ., data = melt24bc)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.30722 -0.05782  0.00648  0.13102  0.33304 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.663e+00  3.791e-01   4.388 0.000317 ***
#   Temp_Avg    -1.919e-03  2.103e-03  -0.913 0.372906    
# RH           2.813e-04  2.571e-03   0.109 0.914009    
# WS_ms       -1.735e-01  7.635e-02  -2.273 0.034842 *  
#   netsw       -3.696e+10  3.126e+10  -1.182 0.251659    
# delt_precip -2.784e-01  1.874e-01  -1.486 0.153680    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1943 on 19 degrees of freedom
# Multiple R-squared:  0.4035,	Adjusted R-squared:  0.2465 
# F-statistic:  2.57 on 5 and 19 DF,  p-value: 0.0614

par(mfrow=c(2,2))
plot(mod24.t) #view results 
shapiro.test(residuals(mod24.t)) #0.1753
vif(mod24.t) #generally ok - nesw score is 3.11 

#check equal variance 
m<-median(predict(mod24.t)) #find a median for fitted values (Y-hat)
(g1<-residuals(mod24.t)[predict(mod24.t)>m]) 
(g2<-residuals(mod24.t)[predict(mod24.t)<m]) 
var.test(g1,g2) #p-value = 0.9182 - good 

#passes all test ok to reduce the model 

############################################
#reduce transformed model 
#p to reject = 0.05
colnames(melt24bc)

#remove RH precip (p=0.914)
mod24.t2<-lm(depth_inc~. , data=melt24bc[,-3])
summary(mod24.t2)
#remove temp (p=0.354)
mod24.t3<-lm(depth_inc~. , data=melt24bc[,-c(2,3)])
summary(mod24.t3)
#remove netsw (p=0.138)
mod24.t4<-lm(depth_inc~. , data=melt24bc[,-c(2,3,5)])
summary(mod24.t4) 
#remove delt_precip (p=0.346)
mod24.t5<-lm(depth_inc~. , data=melt24bc[,-c(2,3,5,6)])
summary(mod24.t5) 
# 
# Call:
#   lm(formula = depth_inc ~ ., data = melt24bc[, -c(2, 3, 5, 6)])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.37751 -0.11753  0.01371  0.13040  0.28601 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.94299    0.11303   8.343 2.07e-08 ***
#   WS_ms       -0.15420    0.05199  -2.966  0.00692 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1944 on 23 degrees of freedom
# Multiple R-squared:  0.2766,	Adjusted R-squared:  0.2452 
# F-statistic: 8.796 on 1 and 23 DF,  p-value: 0.006924
#equation = depth change = 0.9429 -0.15* Wind Speed

par(mfrow=c(2,2))
plot(mod24.t5) #view results 
shapiro.test(residuals(mod24.t5)) #(p=0.3307 - good) 
m<-median(predict(mod24.t5)) #find a median for fitted values (Y-hat)
(g1<-residuals(mod24.t5)[predict(mod24.t5)>m]) 
(g2<-residuals(mod24.t5)[predict(mod24.t5)<m]) 
var.test(g1,g2) #p value = 0.4206 - good 

#############################################
#preform anova to see if reduced model is as good as original model 
anova(mod24.t,mod24.t5)
# Model 1: depth_inc ~ Temp_Avg + RH + WS_ms + netsw + delt_precip
# Model 2: depth_inc ~ WS_ms
# Res.Df     RSS Df Sum of Sq    F Pr(>F)
# 1     19 0.71704                         
# 2     23 0.86950 -4  -0.15247 1.01  0.427
#p>0.05 indicates that the models are not statisticaly different 

##### holdout validation after model fitting
# final reduced model is mod12.t4
# predictions on the test set
prediction <- predict.lm(mod24.t4, newdata = test, se.fit = TRUE)
mean(prediction$se.fit) # mean standard error of the prediction

#############################################
#reduce using different method 
rst<-summary(regsubsets(depth_inc~., data=melt24bc))
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
plot(1:(ncol(melt24bc)-1),rst$cp, xlab="# of predictors in the model", ylab="Mallows Cp")
abline(v=1,col="red", lwd=2)
plot(1:(ncol(melt24bc)-1),rst$adjr2, xlab="# of predictors in the model", ylab="Adjusted R2")
abline(v=2,col="red", lwd=2)
plot(1:(ncol(melt24bc)-1),rst$bic, xlab="# of predictors in the model", ylab="BIC")
abline(v=1,col="red", lwd=2)

##########################################
#check if model with just temp is sig 
mod.24.temp<-lm(depth_inc~Temp_Avg, data=melt24bc)
summary(mod.24.temp)
# Call:
#   lm(formula = depth_inc ~ Temp_Avg, data = melt24bc)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.43456 -0.10774  0.04408  0.10862  0.27909 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.044398   0.188856   5.530 1.26e-05 ***
#   Temp_Avg    -0.004026   0.001782  -2.258   0.0337 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2068 on 23 degrees of freedom
# Multiple R-squared:  0.1815,	Adjusted R-squared:  0.1459 
# F-statistic: 5.101 on 1 and 23 DF,  p-value: 0.03371

#also significant but not as significant 

cor.matrix(melt24bc[,c(2:3)])
par(mfrow=c(1,1))
plot(melt24bc$Temp_Avg, melt24bc$RH)
