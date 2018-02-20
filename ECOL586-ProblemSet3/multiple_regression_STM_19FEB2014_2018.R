# Global NPP: Bivariate and multivariate regression analyses in R.
# Prepared by sean.michaletz@gmail.com, Sept 2013
# Updated by sean.michaletz@gmail.com, Feb 2014

# Required: data file global_npp.txt, packages ggplot2, smatr, car, lmSupport
# Note: These are versions of the data and code used in Michaletz et al. (2014) Nature 512:39.

# Note (Feb 2017): Since writing this code, we've learned that the NPP data from the source
# "Luo (1996); Ni et al. (2001)" was estimated from biomass and age, meaning data from this
# source should not be used in regressions of NPP on biomass and age, as is done in this file.
# Thus while these data are still useful for exploring bivariate versus multiple
# regression analyses, the results produced in this file are not strictly valid, and have
# been updated. For more information, see: 
# Michaletz ST, Cheng D, Kerkhoff AJ, Enquist BJ. 2016. Corrigendum: Convergence of 
# terrestrial plant production across global climate gradients. Nature 537:432.


#################################################################
### code chunk number 1: set working directory and load data 
#################################################################
setwd("E:/")

df <- read.csv("global_npp.csv", header=TRUE)


# Look at dataframe and explain variables!


#################################################################
### code chunk number 2: libraries and options
#################################################################
library(ggplot2)
theme_set(theme_bw()) # Turn off grey background in ggplot2
library(smatr)
library(car)
library(lmSupport)


#################################################################
### code chunk number 3: Climate space diagram (Whittaker)  
#################################################################
qplot(T_ma, ppt_mm, data=df, color=age_class, 
      xlab=expression('Mean annual temperature'~(degree*C)), 
      ylab="Mean annual precipitation (mm)", xlim=c(30, -20), ylim=c(0,5000))


#################################################################
### code chunk number 4: bivariate regressions  
#################################################################
## It is well-established that NPP varies with mean annual temperature (MAT)
## and mean annual precipitation (MAP).  Consequently, MAT and MAP are
## considered primary drivers of variation in NPP.

## These relationships are apparent in our data.  For example, we can
## regress NPP over MAT:
ModelMAT <- lm(anpp_gm2y ~ T_ma, df)
summary(ModelMAT)
qplot(T_ma, anpp_gm2y, data=df, geom = c("point", "smooth"), 
      method="lm", xlab=expression('Mean annual temperature'~(degree*C)), 
      ylab=expression('NPP '~(g~m^-2~yr^-1)))

## This produces a highly significant relationship (p < 2.2x10^-16) with r2 = 0.28.

## We can also regress NPP over MAP
modelMAP <- lm(anpp_gm2y ~ ppt_mm, df)
summary(modelMAP)
qplot(ppt_mm, anpp_gm2y, data=df, geom = c("point", "smooth"), 
      method="lm", xlab="Mean annual precipitation (mm)", 
      ylab=expression('NPP '~(g~m^-2~yr^-1)))

## This again yields a highly significant relationship (p < 2.2x10^-16) w/ r2 = 0.28

## Finally, we can regress NPP on total stand biomass by age class
## Model II SMA regression
sma_age <- sma(log10(df$anpp_gm2y)~log10(df$biomass_gm2)*df$age_class)
summary(sma_age)
plot(sma_age)
## test for common slope
slope.com(log10(df$anpp_gm2y), log10(df$biomass_gm2), df$age_class, method='SMA')

## NPP is significantly influenced by total biomass and decreases with age.
## Slopes are significantly heterogeneous (p = 0.028) so a unique slope required
## for each age class.

## Thus, we would conclude that NPP is governed by MAT, MAP, biomass, and age.
## However, bivariate regressions using individual regressors can be misleading
## (because spurious relationships can result from confounding variables)!  
## Here we extend recent metabolic scaling theory to link abiotic and biotic
## variables to ecosystem NPP.  Multiple regression is then used to fit this theory
## to data.  This multivariate approach properly accounts for covariation among 
## regressors (i.e. potentially confounding variables).  

## For example, MAT is also strongly correlated with growing season length
model_l_T <- lm(lgs_mo ~ T_ma, df)
summary(model_l_T)
qplot(T_ma, lgs_mo, data=df, geom = c("point", "smooth"), 
      method="lm", xlab=expression('Mean annual temperature'~(degree*C)), 
      ylab="Growing season length (mo)")

## Now look at new theory in slides and Michaletz et al. (2014)


#################################################################
### code chunk number 5: multiple regression  
#################################################################

## For simplicity, we consider growing season estimates of 
## temperautre and precipitation, as these are more relevant
## for plant metabolism and physiology.

modelNPP <- lm(log(anpp_gm2y)~GSinvBoltzT+log(GSppt_mm)+log(lgs_mo)+
                 log(biomass_gm2)+log(age_yr), df)
summary(modelNPP)
modelEffectSizes(modelNPP)  # provides partial estimates of r2 (i.e., pEta-sqr) and p

confint(modelNPP) # confidence intervals!
vif(modelNPP) # Check for collinearity; all variance-inflation factors <<< 10, so OK

# Now, produce "partial regression" plots that account for covariation in 
# included independent variables to preserve and express the proper relationship 
# (slope and variance) between predictor and each regressor.
avPlots(modelNPP, col="gray", col.lines="black", main="", 
        ylab=expression('Adjusted NPP'~(g ~m^-2 ~yr^-1)), grid=FALSE)

## After accounting for covariates, NPP is effectively invariant with temperature
## and precipitation! Total stand biomass and age emerge as best predictors.

## Next, in order to evaluate hypothesized relationships between independent 
## variables and NPP (E, alpha, and l_gs), we use multiple regression to estimate 
## model parameters.  However, as MST considers instantaneous metabolic rates, we 
## first define a monthly net primary productivity NPP/l_gs (see theory)

modelNPP_l <- lm(log(anpp_gm2y/lgs_mo)~GSinvBoltzT+log(GSppt_mm)+
                   log(biomass_gm2)+log(age_yr), df)
summary(modelNPP_l)
modelEffectSizes(modelNPP_l)
vif(modelNPP_l) # Check for collinearity; all variance-inflation factors <<<<< 10, so OK
avPlots(modelNPP_l, col="gray", col.lines="black", main="", 
        ylab=expression('Adjusted NPP'~(g ~m^-2 ~yr^-1)), grid=FALSE)

## We can see that E = -0.11603 eV (opposite in direction to hypothesized 0.32 eV)
## alpha = 0.620 (very close to theorized 3/5 = 0.6)


#################################################################
### code chunk number 6: multiple regression (SIMPLER MODEL)
#################################################################
## Since temperature and precipitation are relatively unimportant, 
## we compete fit of full model for NPP with fit of simpler model containing
## only biomass and age.

## First, fit simpler model to obtain age- and mass-scaling exponents
modelNPP_simple <- lm(log(anpp_gm2y)~log(biomass_gm2)+log(age_yr), df)
summary(modelNPP_simple)
modelEffectSizes(modelNPP_simple)  # provides partial estimates of r2 (pEta-sqr) and p
vif(modelNPP_simple)
## Second, use these fitted scaling exponents to rescale age and biomass data:
qplot(log10((biomass_gm2^0.81978)*(age_yr^-0.69687)), log10(anpp_gm2y), 
      data=df, geom = c("point", "smooth"), method="lm", 
      xlab=expression('log'~(M^alpha ~a^alpha_a)), ylab="log(NPP)")
simple <- lm(log10(df$anpp_gm2y)~log10((df$biomass_gm2^0.81978)*(df$age_yr^-0.69687)))
summary(simple)

## This model is more parsimonius than the full model yet explains essentially
## the same amount of variation in NPP.  Consequently, stand age and 
## total stand biomass are the most important predicotrs global 
## variation in NPP.
