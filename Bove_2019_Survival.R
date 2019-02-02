#####################################################################################################
############## Net calcification code for Bove et al. 2019 written by Colleen Bove and ##############
#####################################################################################################

### Libraries required for this code
library(survival)
library(coxme)
library(car)
library(tidyr)
library(dplyr)


#################
## Setting up the data
#################

sur <- read.csv("survival_20June17.csv") 
# set temp, tank, and pCO2 as factors
sur$tank <- as.factor(sur$tank)
sur$pCO2 <- as.factor(sur$pCO2)
sur$temp2 <- as.factor(sur$ftemp)
# create unique colony ID based on species and reef environment
sur$colony <- paste(sur$species, sur$colony, sep='_')
sur$colony <- paste(sur$colony, sur$reefzone, sep='')
sur$colony <- as.factor(sur$colony)
head(sur)
str(sur)


#################
## Set some standards for analysis
#################

treat = c("400_28","400_31","700_28", "700_31","2800_28","2800_31","280_28","280_31") # for treatment ID/order
color <- c("#74a9cf","#fb6a4a","#2b8cbe","#de2d26","#045a8d","#a50f15","#bdc9e1","#fcae91") # for colours plotting

# special for UTEN (missing the 2800_31 treat due to 0% survival)
treat.mod = c("400_28","400_31","700_28", "700_31","2800_28","280_28","280_31") # for treatment ID/order
color.mod <- c("#74a9cf","#fb6a4a","#2b8cbe","#de2d26","#045a8d","#bdc9e1","#fcae91") # for colours plotting


#################
## S. siderea
#################
 
sid<-subset(sur, species=='S') # here we are only working on SSID data, so subset

### perform the survival analysis
sid.r <- Surv(sid$tod,sid$dead.NA) # creates a survival object for use in survfit() function below
sur.sid <- survfit(sid.r ~ pCO2+temp2, data=sid) # Computes an estimate of a survival curve based on pCO2 and temperature using Kaplan-Meier; this will also be used for plotting
sur.sid # view output to see mortality event and sample size summary stats

phsidt <- coxme(sid.r~pco2*ftemp+(1|tank)+(1|colony), sid) # calculate Cox proportional hazards ratios with random effects of tank and colony
summary(phsidt) # view summary of model output
anova(phsidt, test='Chisq') # test significance with Chi sqs for effects of pCO2 and temperature on survival


### plot survival curve with SSID
plot(sur.sid, col=color, lwd=c(5,3,5,3,5,3,5,3),lty=1,xaxt="n", xlab="Timepoint", ylab="Fraction Surviving", main="SSID Survival by Treatment") # lwd= line width of each curve, altered to make all lines visable
  axis(side=1, at=c(0,30,60,90)) # add x axis label for time points
  legend("bottomleft", col=color, lty=1,legend=treat, lwd=3, bty="n") 


#################
## P. strigosa
#################

dip<-subset(sur, species=='P') # here we are only working on PSTR data, so subset

### perform the survival analysis
dip.r <- Surv(dip$tod,dip$dead.NA) # creates a survival object for use in survfit() function below
sur.dip <- survfit(dip.r ~ pCO2+temp2, data=dip) # Computes an estimate of a survival curve based on pCO2 and temperature using Kaplan-Meier; this will also be used for plotting
sur.dip # view output to see mortality event and sample size summary stats
  
phdipt <- coxme(dip.r~pco2*ftemp+(1|tank)+(1|colony), dip) # calculate Cox proportional hazards ratios with random effects of tank and colony
summary(phdipt) # view summary of model output
anova(phdipt, test='Chisq') # test significance with Chi sqs for effects of pCO2 and temperature on survival


### plot survival curve with PSTR
plot(sur.dip, col=color, lwd=c(5,5,5,4,3,4,5,5),lty=1,xaxt="n", xlab="Timepoint", ylab="Fraction Surviving", main="PSTR Survival by Treatment") # lwd= line width of each curve, altered to make all lines visable
  axis(side=1, at=c(0,30,60,90)) # add x axis label for time points
  legend("bottomleft", col=color, lty=1,legend=treat, lwd=3, bty="n") 
  
  
#################
## P. astreoides
#################

por<-subset(sur, species=='A') # here we are only working on PAST data, so subset

### perform the survival analysis
por.r <- Surv(por$tod,por$dead.NA) # creates a survival object for use in survfit() function below
sur.por <- survfit(por.r ~ pCO2+temp2, data=por) # Computes an estimate of a survival curve based on pCO2 and temperature using Kaplan-Meier; this will also be used for plotting
sur.por # view output to see mortality event and sample size summary stats

phport <- coxme(por.r~pco2*ftemp+(1|tank)+(1|colony), por) # calculate Cox proportional hazards ratios with random effects of tank and colony
summary(phport) # view summary of model output
anova(phport, test='Chisq') # test significance with Chi sqs for effects of pCO2 and temperature on survival


### plot survival curve with PAST
plot(sur.por, col=color, lwd=c(5,5,5,5,5,5,3,5),lty=1,xaxt="n", xlab="Timepoint", ylab="Fraction Surviving", main="PAST Survival by Treatment") # lwd= line width of each curve, altered to make all lines visable
  axis(side=1, at=c(0,30,60,90)) # add x axis label for time points
  legend("bottomleft", col=color, lty=1,legend=treat, lwd=3, bty="n") 

  
#################
## U. tenuifolia
#################

ten<-subset(sur, species=='T') # here we are only working on UTEN data, so subset

### perform the survival analysis
ten.r <- Surv(ten$tod,ten$dead.NA) # creates a survival object for use in survfit() function below
sur.ten <- survfit(ten.r ~ pCO2+temp2, data=ten) # Computes an estimate of a survival curve based on pCO2 and temperature using Kaplan-Meier; this will also be used for plotting
sur.ten # view output to see mortality event and sample size summary stats

phtent <- coxme(ten.r~pco2*ftemp+(1|tank)+(1|colony), ten) # calculate Cox protentional hazards ratios with random effects of tank and colony
summary(phtent) # view summary of model output
anova(phtent, test='Chisq') # test significance with Chi sqs for effects of pCO2 and temperature on survival
  

### plot survival curve with UTEN
plot(sur.ten, col=color.mod, lwd=c(5,5,5,5,5,5,4),lty=1,xaxt="n", xlab="Timepoint", ylab="Fraction Surviving", main="UTEN Survival by Treatment") # lwd= line width of each curve, altered to make all lines visable
  axis(side=1, at=c(0,30,60,90)) # add x axis label for time points
  legend("bottomleft", col=color.mod, lty=1,legend=treat.mod, lwd=3, bty="n") 



  