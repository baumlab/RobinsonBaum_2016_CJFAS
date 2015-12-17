
### Script to calculate trophic positions and run trophic position ~ body mass linear models
## Accompanying Robinson & Baum (2016) CJFAS


### Functions for trophic position calculations


### hussey_tp : scaled trophic position estimates. From Hussey et al. (2013, Ecology Letters, pg. 241. Correction published May2nd 2014)

### Arguments: NTP = d15N estimates (default = 9); Nbase = baseline 15N; TPbase = baseline trophic position; Nlim and k from Hussey et al. (2013)
### Returns: scaled trophic positions

 Nlim<--5.92/-0.27
 
 k<--log((5.92 - Nlim)/-Nlim)


hussey_tp<-function(Nbase=9, NTP, TPbase=2.5){


 TP<-(log(Nlim - Nbase) - log(Nlim - NTP))/k + TPbase

return(TP)

}

# orig_tp - additive fractionation - trophic position estimates.

## Arguments: Nbase = baseline 15N (default = 9); NTP = d15N estimates; TPbase = baseline trophic position
## Returns: additive trophic positions

orig_tp<-function(Nbase=9, NTP, TPbase=2, frac=3.4) {


	TP_orig<-TPbase + (NTP - Nbase)/frac

	return(TP_orig)
}


############################################################
############ Trophic position estimates ########################
############################################################

## read in isotope data
load("isotopes_MD_JD_combine.Rdata")
## only concerned with low productivity sites (= LP) on east coast

## Find baseline 15N values
iso[iso$log2==-3.5 & iso$TROPHIC=="Carnivore" & iso$prod=="LP",]
# baseline carnivore = 10.26 dN

iso[iso$log2==2.5 & iso$TROPHIC=="Herbivore" & iso$prod=="LP",]
## multiple individuals for smallest herbivore
mean(iso$dN15[iso$log2==2.5 & iso$TROPHIC=="Herbivore" & iso$prod=="LP"])
# baseline herbivore  dN = 12.21

bases<-data.frame(Nbase = c(10.26, 12.21),
 TROPHIC=c("Carnivore", "Herbivore"),
 prod=c("LP", "LP"))

## calculate trophic position

iso$TP[iso$TROPHIC=="Carnivore" & iso$prod=="LP"]<-hussey_tp(Nbase=bases$Nbase[bases$TROPHIC=="Carnivore" & bases$prod=="LP"],
 iso$dN15[iso$TROPHIC=="Carnivore" & iso$prod=="LP"], TPbase=3)

iso$TP[iso$TROPHIC=="Herbivore" & iso$prod=="LP"]<-orig_tp(Nbase=bases$Nbase[bases$TROPHIC=="Herbivore" & bases$prod=="LP"],
 iso$dN15[iso$TROPHIC=="Herbivore" & iso$prod=="LP"], TPbase=2, frac=4.778)

############################################################
###### Run mixed models for TP ~ mass relationships ############
############################################################

require(nlme)
require(MuMIn)

####################################
### Cross-species approach ############
####################################

## Create species-level data: mean TP for each species

species<-aggregate(TP ~ sp.code, iso_lp, mean)

## Add maximum mass of each species (observed)
mass<-aggregate(mass ~ sp.code, iso_lp, max)
mass_log2<-aggregate(log2 ~ sp.code, iso_lp, max)

species$mass<-mass$mass[match(species$sp.code, mass$sp.code)]
species$log2<-mass_log2$log2[match(species$sp.code, mass_log2$sp.code)]
species$TROPHIC<-iso_lp$TROPHIC[match(species$sp.code, iso_lp$sp.code)]
species$log2_mass<-log2(species$mass)

## Need family for random structure
species$family<-iso_lp$family[match(species$sp.code, iso_lp$sp.code)]

#### model selection for random effect structure
 mod_int_gls<-gls(TP ~ log2_mass*TROPHIC, species)
## interaction, random intercept = species
mod_int_sp<-lme(TP ~ log2_mass*TROPHIC , random = ~1 | sp.code, species)
### interaction, random = family
mod_int_fam<-lme(TP ~ log2_mass*TROPHIC , random = ~1  | family, species)
## interaction, random intercept, nested species in family
mod_int_nest<-lme(TP ~ log2_mass*TROPHIC , random = ~1 | family/sp.code, species)


### random slope - not enough observations

AICc(mod_int_gls, mod_int_fam,mod_int_sp, mod_int_nest)

### Random effect of family is best model

### test fixed effects structure

mod1<-lme(TP ~ log2_mass, random= ~ 1 | family, species, method="ML")
mod2<-lme(TP ~ log2_mass*TROPHIC, random= ~ 1 | family, species, method="ML")

AICc(mod1,mod2)

## trophic interaction important. 
#  Refit with REML for parameter coefficients (recommendation by Zuur)

mod2_cross<-lme(TP ~ log2_mass*TROPHIC, random= ~ 1 | family, species, method="REML")
mod1_cross<-lme(TP ~ log2_mass, random= ~ 1 | family, species, method="REML")
summary(mod1_cross)
summary(mod2_cross)

AICc(mod1_cross, mod2_cross)

r.squaredGLMM(mod2_cross)
r.squaredGLMM(mod2_size)


############################################################
############ All-individuals approach ########################
############################################################

### Model selection for best random structure

mod_int_gls<-gls(TP ~ log2*TROPHIC, iso_lp)
## interaction, random intercept species
mod_int<-lme(TP ~ log2*TROPHIC , random = ~1 | species, iso_lp)
## intearction, random intercept family
mod_int_fam<-lme(TP ~ log2*TROPHIC , random = ~1 | family, iso_lp)
## interaction, random intercept, nested species in family
mod_int_nest<-lme(TP ~ log2*TROPHIC , random = ~1 | family/species, iso_lp)
## interaction, random slope + intercept - species
mod_int2<-lme(TP ~ log2*TROPHIC , random = ~1 +log2 | species, iso_lp)
## interaction, random slope + intercept - family
mod_int2_fam<-lme(TP ~ log2*TROPHIC , random = ~1 +log2 | family, iso_lp)


## interaction, random slope + intercept, nested species in family
mod_int2_nest<-lme(TP ~ log2*TROPHIC , random = ~1 + log2 | family/species, iso_lp)

AICc(mod_int_gls, mod_int,mod_int_fam, mod_int_nest,  mod_int2,mod_int2_fam,  mod_int2_nest)

## Random slope + intercept, with family + species 2 way random effect, is best model.


### Model selection on the fixed structure
ctrl <- lmeControl(opt='optim')
## no interaction, random intercept + slope
mod1<-lme(TP ~ log2 , random = ~1 + log2 | family/species, iso_lp, control=ctrl, method="ML")
## interaction, random intercept + slope
mod2<-lme(TP ~ log2*TROPHIC , random = ~1 + log2 | family/species, iso_lp, method="ML", control=ctrl)

AICc(mod1, mod2)

## now refit with REML for parameter coefficients (Zuur recommendation)
mod1_size<-lme(TP ~ log2 , random = ~1 + log2 | family/species, iso_lp, control=ctrl, method="REML")
mod2_size<-lme(TP ~ log2*TROPHIC , random = ~1 + log2 | family/species, iso_lp, method="REML", control=ctrl)



r.squaredGLMM(mod1_size)
r.squaredGLMM(mod2_size)
summary(mod1_size)
summary(mod2_size)


######## End of trophic position ~ body mass analyses ########
