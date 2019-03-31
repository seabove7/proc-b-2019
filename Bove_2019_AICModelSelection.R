################################################################################################################
#################### AIC Model Selection code for Bove et al. 2019 written by Colleen Bove  ####################
################################################################################################################

setwd("/Users/Colleen/Dropbox/Boston_Experiment/Manuscript1_calcification/R_Code/Paper_code_22Oct18")

### Libraries required for this code
library(tidyr)
library(lme4)


#################
## Setting the Dataframe
#################

### Read in data and manipulate dataframe
coral <- read.csv("rate_22May17.csv")
coral[1] <- NULL # this is an arbitrary column that can be removed
coral$tp <- as.factor(coral$tp) # set tp (time point) as a factor
coral$species <- factor(coral$species, levels = c("S", "P", "A", "T")) # reorder species levels (S = S. siderea; A = P. astreoides)
coral$pco2 <- factor(coral$pco2, levels = c("pre", "cur", "eoc", "ext")) # reorder pCO2 levels
coral$ftemp<-as.factor(coral$temp) # create factor temperature column 

# add column for individual colony ID based on species and reef environment
coral$colony <- paste(coral$species, coral$colony, sep = "") 
coral$colony <- paste(coral$colony, coral$reefzone, sep = "_")

# subset dataframe for just final T90 data (total)
all <- subset(coral, tp == "total")

# this function remove rows with any missing values
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
all <- completeFun(all, "rate") # run function to remove any missing values from "rate"


#################
## Fitting models
#################

#full model
full.model <- lmer(rate ~ species * pco2 * ftemp * reefzone + (1 | colony), REML=F, data = all)

#removal of single parameters (all combinations)
dropRZ.model <- lmer(rate ~ species * pco2 * ftemp + (1 | colony), REML=F, data = all)
dropT.model <- lmer(rate ~ species * pco2 * reefzone + (1 | colony), REML=F, data = all)
dropP.model <- lmer(rate ~ species * ftemp * reefzone + (1 | colony), REML=F, data = all)
dropS.model <- lmer(rate ~ reefzone * pco2 * ftemp + (1 | colony), REML=F, data = all)

#removal of 2 fixed effects
dropRZT.model <- lmer(rate ~ species * pco2 + (1 | colony), REML=F, data = all)
dropRZP.model <- lmer(rate ~ species * ftemp + (1 | colony), REML=F, data = all)
dropPT.model <- lmer(rate ~ species * reefzone + (1 | colony), REML=F, data = all)
dropRZS.model <- lmer(rate ~ ftemp * pco2 + (1 | colony), REML=F, data = all)
dropSP.model <- lmer(rate ~ ftemp * reefzone + (1 | colony), REML=F, data = all)
dropST.model <- lmer(rate ~ pco2 * reefzone + (1 | colony), REML=F, data = all)

#single fixed effect
onlyS.model <- lmer(rate ~ species + (1 | colony), REML=F, data = all)
onlyT.model <- lmer(rate ~ ftemp + (1 | colony), REML=F, data = all)
onlyP.model <- lmer(rate ~ pco2 + (1 | colony), REML=F, data = all)
onlyRZ.model <- lmer(rate ~ reefzone + (1 | colony), REML=F, data = all)

#dropping interactions
interaction1 <- lmer(rate ~ species * pco2 * ftemp + reefzone + (1 | colony), REML=F, data = all)
interaction2 <- lmer(rate ~ species * pco2 * reefzone + ftemp + (1 | colony), REML=F, data = all)
interaction3 <- lmer(rate ~ species * ftemp * reefzone + pco2 + (1 | colony), REML=F, data = all)
interaction4 <- lmer(rate ~ species * ftemp * pco2 + reefzone + (1 | colony), REML=F, data = all)
interaction5 <- lmer(rate ~ pco2 * ftemp * reefzone + species + (1 | colony), REML=F, data = all)
interaction6 <- lmer(rate ~ species * pco2 + ftemp + reefzone + (1 | colony), REML=F, data = all)
interaction7 <- lmer(rate ~ species + pco2 + ftemp * reefzone + (1 | colony), REML=F, data = all)
interaction8 <- lmer(rate ~ species + pco2 * ftemp + reefzone + (1 | colony), REML=F, data = all)
interaction9 <- lmer(rate ~ species * ftemp + reefzone + pco2  +(1 | colony), REML=F, data = all) #best AIC
interaction10 <- lmer(rate ~ species * reefzone + pco2 + ftemp +(1 | colony), REML=F, data = all)
interaction11 <- lmer(rate ~ pco2 * reefzone + species + ftemp +(1 | colony), REML=F, data = all)
interaction12 <- lmer(rate ~ species + pco2 + ftemp + reefzone +(1 | colony), REML=F, data = all)

#after removing reefzone 
noRZ1.model <- lmer(rate ~ species * ftemp + pco2 + (1 | colony), REML=F, data = all)
noRZ2.model <- lmer(rate ~ species + pco2 + ftemp + (1 | colony), REML=F, data = all) #second best
noRZ3.model <- lmer(rate ~ pco2 * ftemp + species + (1 | colony), REML=F, data = all)
noRZ4.model <- lmer(rate ~ species * pco2 + ftemp + (1 | colony), REML=F, data = all)
final.model <- lmer(rate ~ species * (pco2 + ftemp) + (1 | colony), REML=F, data = all)



 # view log likelihood, degrees of freedom, and AIC for all models
model.summary <- cbind(sapply(list(full.model,dropRZ.model,dropT.model,dropP.model,dropS.model,dropRZT.model,dropRZP.model,dropPT.model,dropRZS.model,dropSP.model,dropST.model,onlyS.model,onlyT.model,onlyP.model,onlyRZ.model,interaction1,interaction2,interaction3,interaction4,interaction5,interaction6,interaction7,interaction8,interaction9,interaction10,interaction11,interaction12,noRZ1.model,noRZ2.model,noRZ3.model,noRZ4.model,final.model
), logLik), AIC(full.model,dropRZ.model,dropT.model,dropP.model,dropS.model,dropRZT.model,dropRZP.model,dropPT.model,dropRZS.model,dropSP.model,dropST.model,onlyS.model,onlyT.model,onlyP.model,onlyRZ.model,interaction1,interaction2,interaction3,interaction4,interaction5,interaction6,interaction7,interaction8,interaction9,interaction10,interaction11,interaction12,noRZ1.model,noRZ2.model,noRZ3.model,noRZ4.model,final.model
))
colnames(model.summary)[1] <- 'LL'
model.summary
