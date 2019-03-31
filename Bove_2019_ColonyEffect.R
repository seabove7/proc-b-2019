######################################################################################################################
################ Colony Effect code for Bove et al. 2019 written by Colleen Bove and James Umbanhowar ################
######################################################################################################################

### Libraries required for this code
library(ggplot2)
library(lme4)
library(Rmisc)
library(tidyr)
library(dplyr)
library(data.table)
library(stringr)
library(cowplot)
library(brms)
library(ggridges)


#################
## Setting the Dataframe
#################

### Read in data and manipulate dataframe
coral <- read.csv("rate05Feb19.csv")
coral[1] <- NULL # this is an arbitrary column that can be removed
coral[1] <- NULL # this is an arbitrary column that can be removed
coral$tp <- as.factor(coral$tp) # set tp (time point) as a factor
coral$species <- factor(coral$species, levels = c("S", "P", "A", "T")) # reorder species levels (S = S. siderea; P = P. strigosa; A = P. astreoides; U = U. tenuifolia)
coral$pco2 <- factor(coral$pco2, levels = c("pre", "cur", "eoc", "ext")) # reorder pCO2 levels
coral$ftemp<-as.factor(coral$temp) # create factor temperature column 

# add column for individual colony ID based on species and reef environment
coral$colony <- paste(coral$species, coral$colony, sep = "") 
coral$colony <- paste(coral$colony, coral$reefzone, sep = "_")

# subset dataframe for just final T90 data (total)
ceffect <- subset(coral, tp == "total") # & pco2 %in% c("pre", "cur", "eoc") & species %in% c("S","A"))

# this function remove rows with any missing values
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
ceffect <- completeFun(ceffect, "rate") # run function to remove any missing values from "rate"


#################
## Bayesian Model for Colony Effect
#################

stanncmod <- brm(rate ~ species*(pco2+ftemp) + (1+ftemp+pco2| colony),family=gaussian(), data = ceffect,cores=4)

# examine the credible intervals
summary(stanncmod) # show 95% credible
summary(stanncmod,prob=.75) # show 75% credible ...
ranef(stanncmod) # show all random effects

# Make the three plots of the mean random effects against each other.

#make df of colony random effects
colony <- data.frame(ranef(stanncmod)$colony[,1,])
colony.error <- data.frame(ranef(stanncmod)$colony[,2,])
colnames(colony.error)[1:5] <- c("base", "temp31", "cur", "eoc", "ext") # rename columns for error
colony.error$species<-c(rep("A", times=12), rep("P", times=12),rep("S", times=12), rep("T", times=10)) # add a species column
head(colony.error)

both <- cbind(colony, colony.error)
setDT(both, keep.rownames = TRUE)[]
colnames(both)[1] <- c("colony")
head(both)

#write.csv(both, "Bayesian_output19Feb19.csv")
both <- read.csv("Bayesian_output19Feb19.csv", header=T)

#Intercept vs. current
cur.fig<- ggplot(data= both, aes(Intercept, pco2cur))+
  geom_errorbarh(aes(xmin = Intercept-base, xmax = Intercept+base), colour="grey")+
  geom_errorbar(aes(ymin = pco2cur-cur, ymax = pco2cur+cur), colour="grey")+
  geom_point(size=2, aes(colour=species, shape=species))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())+ 
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="bottom")+
  scale_colour_manual(labels=c("PAST","PSTR","SSID","UTEN"), values=c("#009E73", "#0072B2", "#D55E00", "#CC79A7"))+
  scale_shape_manual(labels=c("PAST","PSTR","SSID","UTEN"), values=c(15,16,17,18))+
  ylim(-0.7, 0.7)


#Intercept vs eoc
eoc.fig<- ggplot(data= both, aes(Intercept, pco2eoc))+
  geom_errorbarh(aes(xmin = Intercept-base, xmax = Intercept+base), colour="grey")+
  geom_errorbar(aes(ymin = pco2eoc-eoc, ymax = pco2eoc+eoc), colour="grey")+
  geom_point(size=2, aes(colour=species, shape=species))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())+ 
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="bottom")+
  scale_colour_manual(labels=c("PAST","PSTR","SSID","UTEN"), values=c("#009E73", "#0072B2", "#D55E00", "#CC79A7"))+
  scale_shape_manual(labels=c("PAST","PSTR","SSID","UTEN"), values=c(15,16,17,18))+
  ylim(-0.7, 0.7)


#Intercept vs extreme
ext.fig<- ggplot(data= both, aes(Intercept, pco2ext))+
  geom_errorbarh(aes(xmin = Intercept-base, xmax = Intercept+base), colour="grey")+
  geom_errorbar(aes(ymin = pco2ext-ext, ymax = pco2ext+ext), colour="grey")+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())+ 
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="bottom")+
  scale_colour_manual(labels=c("PAST","PSTR","SSID","UTEN"), values=c("#009E73", "#0072B2", "#D55E00", "#CC79A7"))+
  scale_shape_manual(labels=c("PAST","PSTR","SSID","UTEN"), values=c(15,16,17,18))+
  geom_point(size=2, aes(colour=species, shape=species))+
  ylim(-0.7, 0.7)


#Intercept vs 31C
temp.fig<- ggplot(data= both, aes(Intercept, ftemp31))+
  geom_errorbarh(aes(xmin = Intercept-base, xmax = Intercept+base), colour="grey")+
  geom_errorbar(aes(ymin = ftemp31-temp31, ymax = ftemp31+temp31), colour="grey")+
  geom_point(size=2, aes(colour=species, shape=species))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())+ 
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="bottom")+
  scale_colour_manual(labels=c("PAST","PSTR","SSID","UTEN"), values=c("#009E73", "#0072B2", "#D55E00", "#CC79A7"))+
  scale_shape_manual(labels=c("PAST","PSTR","SSID","UTEN"), values=c(15,16,17,18))+
  ylim(-0.7, 0.7)

# use cowplot to plot the three desired figures
REplot <- plot_grid(cur.fig, eoc.fig, temp.fig, labels = c("A", "B", "C"), align = 'h', ncol = 3)
REplot


#################
## Plotting violin plots of the correlation coefficients
#################

library(magrittr)
library(dplyr)
library(forcats)
library(tidyr)
library(modelr)
library(tidybayes)
library(ggplot2)
library(ggstance)
library(ggridges)
library(cowplot)
library(rstan)
library(brms)
library(ggrepel)
theme_set(theme_tidybayes() + panel_border() + background_grid()) # set the theme

# these help it run faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Make violin plot of the correlation coefficients from brms model with 75 and 95% credible intervals!
stanncmod %>% 
  spread_draws(b_Intercept, c(cor_colony__Intercept__pco2cur, cor_colony__Intercept__pco2eoc, cor_colony__Intercept__pco2ext, cor_colony__Intercept__ftemp31)) %>%
  gather(treat, correlation, cor_colony__Intercept__ftemp31:cor_colony__Intercept__pco2ext) %>%
  ggplot(aes(y= treat, x= correlation)) +
  geom_vline(xintercept=0, linetype=2, colour="#999999") +
  geom_halfeyeh(.width = c(0.75,0.95)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())


  
