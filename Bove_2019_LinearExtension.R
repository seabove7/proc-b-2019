######################################################################################################################
############## Linear Extension code for Bove et al. 2019 written by Colleen Bove and James Umbanhowar ##############
######################################################################################################################


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
coral <- read.csv("rate_22May17.csv")
coral[1] <- NULL # this is an arbitrary column that can be removed
coral <- subset(coral, species %in% c("S","A")) # subset for only S and A (P and T not measured for linear extension)
coral$tp <- as.factor(coral$tp) # set tp (time point) as a factor
coral$species <- factor(coral$species, levels = c("S", "A")) # reorder species levels (S = S. siderea; A = P. astreoides)
coral$pco2 <- factor(coral$pco2, levels = c("pre", "cur", "eoc", "ext")) # reorder pCO2 levels
coral$ftemp<-as.factor(coral$temp) # create factor temperature column 

# add column for individual colony ID based on species and reef environment
coral$colony <- paste(coral$species, coral$colony, sep = "") 
coral$colony <- paste(coral$colony, coral$reefzone, sep = "_")

# subset dataframe for just final T90 data (total)
lerate <- subset(coral, tp == "total")

# this function remove rows with any missing values
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
lerate <- completeFun(lerate, "le") # run function to remove any missing values from "rate"


#################
## Model linear extension
#################

# mixed effect model run using lmer()
lemodel <- lmer(le ~ species*(pco2+ftemp)+(1 | colony), data = lerate)

summary(lemodel) # view summary of model

### Performing the parametric bootstrapping of the model:

bootnum = 1500 # set number of iterations (we used 1500) between 999 and 9999
seed = 4 #seed to make results replicatable (our seed was 4)

out <- simulate(lemodel,nsim=bootnum,seed=seed,re.form=NULL) # simulate your model set number of times in a dataframe (samples using random effects)
boots <- apply(out,2,function(x) predict(lmer(x ~ species*(pco2+ftemp)+(1 | colony), data = lerate),re.form=NA)) # applies the predict (does not use random effects) FUNCTION to the columns of the 'out' dataframe

boots <- cbind(predict(lemodel,re.form=NA), boots) #combines boots matrix created above with the predicted values based on the model into a single matrix
lerate.a <- (cbind(lerate, as.data.frame(t(apply(boots, 1, function(x) c(mean(x), quantile(x, c(0.025, 0.975)))))))) #estimates mean and 95% confidence intercals for the prediction values and adds them to your new dataframe

head(lerate.a) # view your new dataframe so you can rename the mean and CI columns
colnames(lerate.a)[17:19] <- c("mean", "lowerci", "upperci") # rename mean/CI columns

write.csv(lerate.a,file="extension.csv") # save this dataframe as a CSV to keep bootstrapped mean and CI without having to rerun


#################
## Plot linear extension
#################

lerate.a<-read.csv("extension.csv") # read in the CSV or ignore this is not saving the dataframe in previous step

lerate.a $temp2<-as.factor(lerate.a $temp) #set temp as factor
theme_set(theme_bw()) # set ggplot theme to bw
dodge<-position_dodge(0.4) # establish 'dodge' width
lerate.a$species <- factor(lerate.a $species, levels = c("S", "A")) # reorder species again

# calculate the mean treatment pCO2 for plotting
mtreat<-lerate.a %>% 
  group_by(treat) %>% 
  summarize(mpco2 = mean(pCO2, na.rm=TRUE))

# add the treatment mean pCO2 as a new column (mpco2)
lerate.a<-lerate.a %>% 
  left_join(mtreat, by=c("treat" = "treat"))


# plot linear extension raw rates (points) behind 95% CI in response to pCO2 and temperature for all four coral species
ggplot()+
  geom_hline(yintercept=0, linetype=2, colour="#999999")+ # add horizontal line at '0'
  geom_point(data=lerate.a,aes(x=pCO2, y=le,colour=temp2,shape=temp2), size=1.3)+ # plots points of raw linear extension rates
  geom_errorbar(data=lerate.a, aes(x=mpco2, ymin=lowerci, ymax=upperci, width=0,colour=factor(temp)), size=0.8, position=dodge)+ # plots lmer 95% confidence intervals as bars over the raw rate points 
  scale_shape_manual("Temperature", values=c(1,2))+ # determine shape by temperature (factor)
  scale_colour_manual("Temperature", labels=c("28C Estimate","31C Estimate"),values=c("#56B4E9", "#D55E00"))+ # determine colour by temperature (factor) and rename values
  scale_x_continuous(trans='log2', breaks=c(280,400,700,2800))+ # plot x axis as log scale for easy of data visualization
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())+ # remove plot grid/background
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="bottom")+ # angles x axis text and but legend on the bottom
  ylab(bquote('Linear extension (mm/day)'))+
  xlab("pCO2 (uatm)")+
  facet_wrap(~species) # facet wrap by species with y axis different per species

sum<-summarySE(lerate.a, measurevar="lowerci", groupvars=c("species","temp2","pco2"), na.rm =T)
sum$pco2 <- factor(sum$pco2, levels = c("pre", "cur", "eoc", "ext")) # reorder pCO2 levels
sum
