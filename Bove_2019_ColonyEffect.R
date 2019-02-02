########################################################################################################################
############## Colony-level effect code for Bove et al. 2019 written by Colleen Bove and James Umbanhowar ##############
########################################################################################################################

### Libraries required for this code
library(ggplot2)
library(lme4)
library(Rmisc)
library(tidyr)
library(dplyr)
library(data.table)
library(stringr)
library(cowplot)


#################
## Setting the Dataframe
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
ceffect <- subset(coral, tp == "total")

# this function remove rows with any missing values
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
ceffect <- completeFun(ceffect, "rate") # run function to remove any missing values from "rate"


#################
## Correlation coefficients between colony intercept and slope
#################

# mixed effect model run using lmer()
colonymodel <- lmer(rate ~ species*(pco2+ftemp)+(1+species | tank)+(1 + (pco2+ftemp) | colony), data = ceffect)

ceffect$resid<-resid(colonymodel) # add model residuals to dataframe

### Extracting the bootstrapped correlation coefficients of colony from the model:

bootnum = 15 #00 # set number of iterations (we used 1500) between 999 and 9999
seed = 4 #seed to make results replicatable (our seed was 4)

out <- simulate(colonymodel,nsim=bootnum,seed=seed,re.form=NULL) # simulate your model set number of times in a dataframe (samples using random effects)
boots <- apply(out,2,function(x) predict(lmer(x ~  species * (pco2 + ftemp) + (1 + species | tank) + (1 + (pco2 + ftemp) | colony), data = ceffect, REML=T),re.form=NA)) # applies the predict (does not use random effects) FUNCTION to the columns of the 'out' dataframe to ...

boot2 <- apply(out,2,function(x) as.data.frame(VarCorr(lmer(x ~  species * (pco2 + ftemp) + (1 + species | tank) + (1 + (pco2 + ftemp) | colony), data = ceffect, REML=T)))[6:9,5]) # Fit the same model to all the simulated data and extract the 4 correlations we want (correlations between colony intercept and slopes are rows 6-9 and column 5)

corrframe <- cbind(as.data.frame(VarCorr(colonymodel))[6:9,c(2,3,5)],t(apply(boot2,1,function(x) quantile(x, c(0.025, 0.975))))) # estimate the 95% confidence intervals for the 4 correlations
colnames(corrframe)[4:5] <- c("lowerci2", "upperci2") # rename bounds of CI

# combining the dataframes                            
boots <- cbind(predict(colonymodel,re.form=NA), boots) #combine the predicted values with the simulated values into one dataframe 
ceffect.a <- cbind(ceffect, predict(colonymodel,re.form=NA), as.data.frame(t(apply(boots, 1, function(x) quantile(x, c(0.025, 0.975)))))) # estimate the 95% confidence intervals for the prediction and the 95% bootstrap
colnames(ceffect.a)[18:20] <- c("mean", "lowerci", "upperci") # rename mean and 95% CI bounds
head(ceffect.a)

## simple plot of mean correlation coefficients and 95% CI per stressor
ggplot(corrframe, aes(x=var2, y=sdcor))+
  geom_point()+
  geom_errorbar(aes(ymin= lowerci2, ymax=upperci2))+
  geom_hline(yintercept=0, linetype=2, colour="#999999")+
  xlab("Treatment")+
  ylab("Correlation coefficient")

write.csv(ceffect.a,file="colony_effect.csv") # save this dataframe as a CSV to keep bootstrapped mean and CI without having to rerun


#################
## Plot random effects of colony
#################

## create the colony random effect dataframe
ranef.colony<-ranef(colonymodel)$colony # create dataframe of the random effects of colony from the model
setDT(ranef.colony, keep.rownames = TRUE)[] # convert row names to column
ranef.colony$species<-c(rep("A", times=12), rep("P", times=12),rep("S", times=12), rep("T", times=10)) # add a species column
colnames(ranef.colony)[1:6] <- c("colony","base", "cur", "eoc","ext","temp") # rename columns


## plot random effect of colony on each stress treatment (current-day pCO2, end-of-century pCO2, extreme pCO2, 31C temperature) against 'base treatment' (pre-indusctrial pCO2 at 28C)

theme_set(theme_classic()) # set the theme
# plot of current-day pCO2
cur.p<-
  ggplot(ranef.colony, aes(x=base, y=cur, colour=species,shape=species))+
  geom_point(size=2)+
  scale_colour_manual(labels=c("PAST","PSTR","SSID","UTEN"), values=c("#009E73", "#0072B2", "#D55E00", "#CC79A7"))+
  #theme(legend.position="none",axis.text.y=element_blank())+
  ylab("colony random effect for
  present-day pCO2 at 28C")+ # this spaciing allows ylab title to be two lines
  xlab("colony random effect for
  pre-industrial pCO2 at 28C")+ # this spaciing allows ylab title to be two lines
  scale_shape_manual(labels=c("PAST","PSTR","SSID","UTEN"), values=c(15,16,17,18)) # set shape of points based on species

# plot of end-of-century pCO2
eoc.p<-ggplot(ranef.colony, aes(x=base, y=eoc, colour=species,shape=species))+
  geom_point(size=2)+
  scale_colour_manual(labels=c("PAST","PSTR","SSID","UTEN"), values=c("#009E73", "#0072B2", "#D55E00", "#CC79A7"))+
  #theme(legend.position="none",axis.text.y=element_blank())+
  ylab("colony random effect for
  end-of-century pCO2 at 28C")+ # this spaciing allows ylab title to be two lines
  xlab("colony random effect for
  pre-industrial pCO2 at 28C")+ # this spaciing allows ylab title to be two lines
  scale_shape_manual(labels=c("PAST","PSTR","SSID","UTEN"), values=c(15,16,17,18)) # set shape of points based on species

# plot of extreme pCO2
ext.p<-ggplot(ranef.colony, aes(x=base, y=ext, colour=species,shape=species))+
  geom_point(size=2)+
  scale_colour_manual(labels=c("PAST","PSTR","SSID","UTEN"), values=c("#009E73", "#0072B2", "#D55E00", "#CC79A7"))+
  #theme(legend.position="none",axis.text.y=element_blank())+
  ylab("colony random effect for
  extreme pCO2 at 28C")+ # this spaciing allows ylab title to be two lines
  xlab("colony random effect for
  pre-industrial pCO2 at 28C")+ # this spaciing allows ylab title to be two lines
  scale_shape_manual(labels=c("PAST","PSTR","SSID","UTEN"), values=c(15,16,17,18)) # set shape of points based on species

# plot of 31C temperature
temp.p<-ggplot(ranef.colony, aes(x=base, y=temp, colour=species,shape=species))+
  geom_point(size=2)+
  scale_colour_manual(labels=c("PAST","PSTR","SSID","UTEN"), values=c("#009E73", "#0072B2", "#D55E00", "#CC79A7"))+
  #theme(legend.position="none",axis.text.y=element_blank())+
  ylab("colony random effect for
  pre-industrial pCO2 at 31C")+ # this spaciing allows ylab title to be two lines
  xlab("colony random effect for
  pre-industrial pCO2 at 28C")+ # this spaciing allows ylab title to be two lines
  scale_shape_manual(labels=c("PAST","PSTR","SSID","UTEN"), values=c(15,16,17,18)) # set shape of points based on species


## plot all four figures
plot_grid(cur.p, eoc.p, ext.p, temp.p, ncol=2, nrow=2)

## plot only current pCO2, EOC pCO2, and 31C as one figure like in publication
plot_grid(cur.p, eoc.p, temp.p, ncol=3, nrow=1)