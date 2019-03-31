###############################################################################################################################################################
############## Net Calcification and Linear Extension by Reef Environment code for Bove et al. 2019 written by Colleen Bove and James Umbanhowar ##############
###############################################################################################################################################################

setwd("/Users/Colleen/Dropbox/Boston_Experiment/Manuscript1_calcification/R_Code/Paper_code_13Feb19")

### Libraries required for this code
library(ggplot2)
library(lme4)
library(Rmisc)
library(tidyr)
library(dplyr)


#################
## Setting the Dataframe for Net Calcification by Reef Environment
#################

### Read in data and manipulate dataframe
coral <- read.csv("rate_22May17.csv")
coral[1] <- NULL # this is an arbitrary column that can be removed
coral$tp <- as.factor(coral$tp) # set tp (time point) as a factor
coral$species <- factor(coral$species, levels = c("S", "P", "A", "T")) # reorder species levels (S = S. siderea; P = P. strigosa; A = P. astreoides; U = U. tenuifolia)
coral$pco2 <- factor(coral$pco2, levels = c("pre", "cur", "eoc", "ext")) # reorder pCO2 levels
coral$ftemp<-as.factor(coral$temp) # create factor temperature column 

# add column for individual colony ID based on species and reef environment
coral$colony <- paste(coral$species, coral$colony, sep = "") 
coral$colony <- paste(coral$colony, coral$reefzone, sep = "_")

# subset dataframe for just final T90 data (total)
rate.rz <- subset(coral, tp == "total")

# this function remove rows with any missing values
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
rate.rz <- completeFun(rate.rz, "rate") # run function to remove any missing values from "rate"


#################
## Model net calcification rate by reef environment 
#################

# mixed effect model run using lmer()
rate.rzmodel <- lmer(rate ~ species * (pco2 + ftemp) + reefzone + (1 | colony), data = rate.rz)

summary(rate.rzmodel) # view summary of model

### Performing the parametric bootstrapping of the model:

bootnum = 1500 # set number of iterations (we used 1500) between 999 and 9999
seed = 4 #seed to make results replicatable (our seed was 4)

out <- simulate(rate.rzmodel,nsim=bootnum,seed=seed,re.form=NULL) # simulate your model set number of times in a dataframe (samples using random effects)
boots <- apply(out,2,function(x) predict(lmer(x ~ species * (pco2 + ftemp) + reefzone + (1 | colony), data = rate.rz),re.form=NA)) # applies the predict (does not use random effects) FUNCTION to the columns of the 'out' dataframe

boots <- cbind(predict(rate.rzmodel,re.form=NA), boots) #combines boots matrix created above with the predicted values based on the model into a single matrix
rate.rz.a <- (cbind(rate.rz, as.data.frame(t(apply(boots, 1, function(x) c(mean(x), quantile(x, c(0.025, 0.975)))))))) #estimates mean and 95% confidence intercals for the prediction values and adds them to your new dataframe

head(rate.rz.a) # view your new dataframe so you can rename the mean and CI columns
colnames(rate.rz.a)[17:19] <- c("mean", "lowerci", "upperci") # rename mean/CI columns

write.csv(rate.rz.a,file="NetCalc_RZ.csv") # save this dataframe as a CSV to keep bootstrapped mean and CI without having to rerun


#################
## Plot net calcification by reef environment
#################

rate.rz.a<-read.csv("NetCalc_RZ.csv") # read in the CSV or ignore this is not saving the dataframe in previous step

rate.rz.a $temp2<-as.factor(rate.rz.a $temp) #set temp as factor
theme_set(theme_bw()) # set ggplot theme to bw
dodge<-position_dodge(0.4) # establish 'dodge' width
rate.rz.a$species <- factor(rate.rz.a $species, levels = c("S", "P", "A", "T")) # reorder species again

# calculate the mean treatment pCO2 for plotting
mtreat<-rate.rz.a %>% 
  group_by(treat) %>% 
  summarize(mpco2 = mean(pCO2, na.rm=TRUE))

# add the treatment mean pCO2 as a new column (mpco2)
rate.rz.a<-rate.rz.a %>% 
  left_join(mtreat, by=c("treat" = "treat"))


# plot net calcificaiton raw rates (points) behind 95% CI in response to pCO2 and temperature for all four coral species by reef environment 
ggplot()+
  geom_hline(yintercept=0, linetype=2, colour="#999999")+ # add horizontal line at '0'
  geom_point(data=rate.rz.a,aes(x=mpco2, y=rate,colour=reefzone,shape=reefzone), size=1.3, position=dodge)+ # plots points of raw net calcification rates
  geom_errorbar(data=rate.rz.a, aes(x=mpco2, ymin=lowerci, ymax=upperci, width=0,colour=reefzone), size=0.8, position=dodge)+ # plots lmer 95% confidence intervals as bars over the raw rate points 
  scale_shape_manual("Reef Environment", values=c(1,2))+ # determine shape by temperature (factor)
  scale_colour_manual("Reef Environment", labels=c("Offshore","Innshore"),values=c("black", "grey"))+ # determine colour by temperature (factor) and rename values
  scale_x_continuous(trans='log2', breaks=c(280,400,700,2800))+ # plot x axis as log scale for easy of data visualization
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())+ # remove plot grid/background
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="bottom")+ # angles x axis text and but legend on the bottom
  ylab(bquote('Net calcification rate (mg/cm2/day)'))+
  xlab("pCO2 (uatm)")+
  facet_grid(species~temp, scales="free_y") # facet wrap by species and temperature with y axis different per species


sum<-summarySE(three.a, measurevar="mean", groupvars=c("species","ftemp","pco2", "reefzone"), na.rm =T)
sum$pco2 <- factor(sum$pco2, levels = c("pre", "cur", "eoc", "ext")) # reorder pCO2 levels
sum




#################
## Modify the Dataframe for Linear Extension by Reef Environment
#################

# update dataframe for linear extension (just SSID and PAST)
coral <- subset(coral, species %in% c("S","A")) # subset for only S and A (P and T not measured for linear extension)
coral$species <- factor(coral$species, levels = c("S", "A")) # reorder species levels (S = S. siderea; A = P. astreoides)
le.rz <- subset(coral, tp == "total")

# this function remove rows with any missing values
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
le.rz <- completeFun(le.rz, "le") # run function to remove any missing values from "le"


#################
## Model linear extension  by reef environment 
#################

# mixed effect model run using lmer()
le.rzmodel <- lmer(le ~ species * (pco2 + ftemp) + reefzone + (1 | colony), data = le.rz)

summary(le.rzmodel) # view summary of model

### Performing the parametric bootstrapping of the model:

bootnum = 1500 # set number of iterations (we used 1500) between 999 and 9999
seed = 4 #seed to make results replicatable (our seed was 4)

out <- simulate(le.rzmodel,nsim=bootnum,seed=seed,re.form=NULL) # simulate your model set number of times in a dataframe (samples using random effects)
boots <- apply(out,2,function(x) predict(lmer(x ~ species * (pco2 + ftemp) + reefzone + (1 | colony), data = le.rz),re.form=NA)) # applies the predict (does not use random effects) FUNCTION to the columns of the 'out' dataframe

boots <- cbind(predict(le.rzmodel,re.form=NA), boots) #combines boots matrix created above with the predicted values based on the model into a single matrix
le.rz.a <- (cbind(le.rz, as.data.frame(t(apply(boots, 1, function(x) c(mean(x), quantile(x, c(0.025, 0.975)))))))) #estimates mean and 95% confidence intercals for the prediction values and adds them to your new dataframe

head(le.rz.a) # view your new dataframe so you can rename the mean and CI columns
colnames(le.rz.a)[17:19] <- c("mean", "lowerci", "upperci") # rename mean/CI columns

write.csv(le.rz.a,file="LE_RZ.csv") # save this dataframe as a CSV to keep bootstrapped mean and CI without having to rerun


#################
## Plot linear extension by reef environment
#################

le.rz.a<-read.csv("LE_RZ.csv") # read in the CSV or ignore this is not saving the dataframe in previous step

le.rz.a $temp2<-as.factor(le.rz.a $temp) #set temp as factor
theme_set(theme_bw()) # set ggplot theme to bw
dodge<-position_dodge(0.4) # establish 'dodge' width
le.rz.a$species <- factor(le.rz.a $species, levels = c("S", "P", "A", "T")) # reorder species again

# calculate the mean treatment pCO2 for plotting
mtreat<-le.rz.a %>% 
  group_by(treat) %>% 
  summarize(mpco2 = mean(pCO2, na.rm=TRUE))

# add the treatment mean pCO2 as a new column (mpco2)
le.rz.a<-le.rz.a %>% 
  left_join(mtreat, by=c("treat" = "treat"))


# plot linear extension raw rates (points) behind 95% CI in response to pCO2 and temperature for all four coral species by reef environment
ggplot()+
  geom_hline(yintercept=0, linetype=2, colour="#999999")+ # add horizontal line at '0'
  geom_point(data=le.rz.a,aes(x=mpco2, y=le,colour=reefzone,shape=reefzone), size=1.3, position=dodge)+ # plots points of raw linear extension rates
  geom_errorbar(data=le.rz.a, aes(x=mpco2, ymin=lowerci, ymax=upperci, width=0,colour=reefzone), size=0.8, position=dodge)+ # plots lmer 95% confidence intervals as bars over the raw le points 
  scale_shape_manual("Reef Environment", values=c(1,2))+ # determine shape by temperature (factor)
  scale_colour_manual("Reef Environment", labels=c("Offshore","Innshore"),values=c("black", "grey"))+ # determine colour by temperature (factor) and rename values
  scale_x_continuous(trans='log2', breaks=c(280,400,700,2800))+ # plot x axis as log scale for easy of data visualization
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())+ # remove plot grid/background
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="bottom")+ # angles x axis text and but legend on the bottom
  ylab(bquote('linear extension le (mg/cm2/day)'))+
  xlab("pCO2 (uatm)")+
  facet_grid(species~temp, scales="free_y") # facet wrap by species and temperature with y axis different per species


sum<-summarySE(three.a, measurevar="mean", groupvars=c("species","ftemp","pco2", "reefzone"), na.rm =T)
sum$pco2 <- factor(sum$pco2, levels = c("pre", "cur", "eoc", "ext")) # reorder pCO2 levels
sum
