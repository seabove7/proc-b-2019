########################################################################################################################
############## Gross calcification code for Bove et al. 2019 written by Colleen Bove and James Umbanhowar ##############
########################################################################################################################

### Libraries required for this code
library(ggplot2)
library(lme4)
library(Rmisc)
library(tidyr)
library(dplyr)


#################
## Setting the Dataframe
#################

### Read in data and manipulate dataframe
coral <- read.csv("gross_20Jan18.csv")
coral[1] <- NULL # this is an arbitrary column that can be removed
coral$tp <- as.factor(coral$tp) # set tp (time point) as a factor
coral$species <- factor(coral$species, levels = c("S", "P", "A", "T")) # reorder species levels (S = S. siderea; P = P. strigosa; A = P. astreoides; U = U. tenuifolia)
coral$pco2 <- factor(coral$pco2, levels = c("pre", "cur", "eoc", "ext")) # reorder pCO2 levels
coral$ftemp<-as.factor(coral$temp) # create factor temperature column 

# add column for individual colony ID based on species and reef environment
coral$colony <- paste(coral$species, coral$colony, sep = "") 
coral$colony <- paste(coral$colony, coral$reefzone, sep = "_")

# subset dataframe for just final T90 data (total)
gcalc <- subset(coral, tp == "total")

# this function remove rows with any missing values
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
gcalc <- completeFun(gcalc, "gross") # run function to remove any missing values from "gross"


#################
## Model gross calcification
#################

# mixed effect model run using lmer()
gcalcmodel <- lmer(gross ~ species*(pco2+ftemp)+(1 | colony), data = gcalc)

summary(gcalcmodel) # view summary of model

### Performing the parametric bootstrapping of the model:

bootnum = 1500 # set number of iterations (we used 1500) between 999 and 9999
seed = 4 #seed to make results replicatable (our seed was 4)

out <- simulate(gcalcmodel,nsim=bootnum,seed=seed,re.form=NULL) # simulate your model set number of times in a dataframe (samples using random effects)
boots <- apply(out,2,function(x) predict(lmer(x ~  species * (pco2 + ftemp) + (1 | colony), data = gcalc),re.form=NA)) # applies the predict (does not use random effects) FUNCTION to the columns of the 'out' dataframe

boots <- cbind(predict(gcalcmodel,re.form=NA), boots) #combines boots matrix created above with the predicted values based on the model into a single matrix
gcalc.a <- (cbind(gcalc, as.data.frame(t(apply(boots, 1, function(x) c(mean(x), quantile(x, c(0.025, 0.975)))))))) #estimates mean and 95% confidence intercals for the prediction values and adds them to your new dataframe

head(gcalc.a) # view your new dataframe so you can rename the mean and CI columns
colnames(gcalc.a)[18:20] <- c("mean", "lowerci", "upperci") # rename mean/CI columns

write.csv(gcalc.a,file="gross_calc.csv") # save this dataframe as a CSV to keep bootstrapped mean and CI without having to rerun


#################
## Plot gross calcification
#################

gcalc.a<-read.csv("gross_calc.csv") # read in the CSV or ignore this is not saving the dataframe in previous step

gcalc.a $temp2<-as.factor(gcalc.a $temp) #set temp as factor
theme_set(theme_bw()) # set ggplot theme to bw
dodge<-position_dodge(0.4) # establish 'dodge' width
gcalc.a$species <- factor(gcalc.a $species, levels = c("S", "P", "A", "T")) # reorder species again

# calculate the mean treatment pCO2 for plotting
mtreat<-gcalc.a %>% 
  group_by(treat) %>% 
  summarize(mpco2 = mean(pCO2, na.rm=TRUE))

# add the treatment mean pCO2 as a new column (mpco2)
gcalc.a<-gcalc.a %>% 
  left_join(mtreat, by=c("treat" = "treat"))


# plot gross calcification raw rates (points) behind 95% CI in response to pCO2 and temperature for all four coral species
ggplot()+
  geom_hline(yintercept=0, linetype=2, colour="#999999")+ # add horizontal line at '0'
  geom_point(data=gcalc.a,aes(x=pCO2, y=gross,colour=temp2,shape=temp2), size=1.3)+ # plots points of raw gross calcification rates
  geom_errorbar(data=gcalc.a, aes(x=mpco2, ymin=lowerci, ymax=upperci, width=0,colour=factor(temp)), size=0.8, position=dodge)+ # plots lmer 95% confidence intervals as bars over the raw rate points 
  scale_shape_manual("Temperature", values=c(1,2))+ # determine shape by temperature (factor)
  scale_colour_manual("Temperature", labels=c("28C Estimate","31C Estimate"),values=c("#56B4E9", "#D55E00"))+ # determine colour by temperature (factor) and rename values
  scale_x_continuous(trans='log2', breaks=c(280,400,700,2800))+ # plot x axis as log scale for easy of data visualization
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())+ # remove plot grid/background
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="bottom")+ # angles x axis text and but legend on the bottom
  ylab(bquote('Gross calcification rate (mg/cm2/day)'))+
  xlab("pCO2 (uatm)")+
  facet_wrap(~species, scales="free_y") # facet wrap by species with y axis different per species

sum<-summarySE(gcalc.a, measurevar="mean", groupvars=c("species","temp2","pco2"), na.rm =T)
sum$pco2 <- factor(sum$pco2, levels = c("pre", "cur", "eoc", "ext")) # reorder pCO2 levels
sum

