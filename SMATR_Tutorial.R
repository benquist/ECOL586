
############################################################### 
# SMATR tutorial
# Brian J. Enquist, Jan 2018
#
# This turoial uses and builds on several references listed below using
# the bigtrees dataset assessing the allometry of stem diameter 
# and height, a brain size and body size dataset for some mammals,
# and a dataset with looking at leaf longevity and leaf mass per 
# unit area
##############################################################

# The key functions in this package are the ability to conduct ANCOVAs using Model II regression (sma and ma), which will fit SMA and MA respectively, and construct confidence intervals or test hypotheses about slope or elevation in one or several samples, depending on how the arguments are used.

# This is a handy package for testing for differences or similarities between groups, treatments, places etc.

##################################################
# Overview of SMATR functions
#################################################
#1 
# sma(y~x) 
# will fit a SMA (for y against x) and return confi-dence intervals for the slope and elevation (Pitman 1939; Warton et al., 2006).

#  2 
# sma(y~x*groups) 
# will test for common slope (Warton & Weber 2002) amongst several SMAs (for y against x), fitted separately for each level of the factor groups.

# 3 
# sma(y~x+groups) 
# will test for common elevation (Warton et al., 2006) amongst several SMAs (for y against x), fitted with common slope but with separate elevations for each level of the factor groups.

#####Additional arguments can be specified to perform some additional tasks: 

# 1 sma(y~x, slope.test=B) 
# will test the hypothesis that the SMA (for y against x) has slope B (Pitman 1939).

#3 The argument 
# multcomp=T
# when used in comparing multiple lines, will return pairwise comparisons of slopes (or elevations or locations along common-slope SMAs), and multcompmethod¼“adjust” will use adjusted P-values (via the ‘Sidak adjust-ment’, Westfall & Young, 1993) to control family-wise error rate in a conservative way.

#4 The argument 
# intercept=F
# used in combination with most of the above functions, will force lines through the origin. This is necessary, for example, when analysing phylogenetically independent contrasts (Felsenstein 1985).

#5 The argument 
# log=“xy” 
# will log 10-transform variables prior to analysis.

# Checking assumptions

#Using the 
# plot 
#function can also be used to produce residual plots, to check the critical assumptions of linearity and equal variance at all fitted values (via a residuals vs. fitted values plot) and the assumption of normally distributed residuals (via a normal quantile plot). Normality can be important for inference when sample size is small. Figure 1b,c was generated using the following commands:

# plot(ft, which=“residual”) 
# plot(ft, which=“qq”) 
# summary


#######################################################################
# Intro Tutorial using smatr

install.packages("smatr")
library(smatr)

#####  Using the Big Trees Dataset
## Note, here we are using an updated version where Class, Family, Genus, and species are broken out so that we can assess groups. PLease make sure to pull this new version of the dataset.
bigtrees.data <- read.csv(file="/Users/brianjenquist/GitHub/R/ECOL586_Biological_Scaling/bigtrees_2_Taxa.csv",header=T)

names(bigtrees.data)

sma(Height.in. ~ Diameter..in., log = "xy", data = bigtrees.data) 
sma(Height.in. ~ Diameter..in., log = "xy", data = bigtrees.data, slope.test=0.666667)
sma(Height.in. ~ Diameter..in.*Class, log = "xy", data = bigtrees.data, slope.test=0.666667) 


ftComSlope = sma(Height.in. ~ Diameter..in.*bigtrees.data$Class, log = "xy", data = bigtrees.data)
summary(ftComSlope)
plot(ftComSlope)


# Example: does brain size scale as the 2/3 power of body size? (Thus brain surface area would scale against mass)
# https://cran.r-project.org/web/packages/MASS/MASS.pdf

library(MASS)
plot(Animals, log = "xy")

ftBrainBody = lm(log(Animals$brain) ~ log(Animals$body))
confint(ftBrainBody)
ftBodyBrain = lm(log(Animals$body) ~ log(Animals$brain))
confint(ftBodyBrain)

ft = sma(brain ~ body, data = Animals, log = "xy")
summary(ft)

# what happens if you reverse the axes? Do you get the same answer?
ft = sma(body ~ brain, data = Animals, log = "xy")
summary(ft)

#now use major axis regression 
ma(brain ~ body, data = Animals, log = "xy")
# TESTING FOR EVIDENCE AGAINST A SLOPE (of 2/3)
ft = sma(brain ~ body, data = Animals, log = "xy", slope.test = 2/3)

# USEFUL PLOTS
ft = ma(brain ~ body, data = Animals, log = "xy")
plot(ft)  #plots the data with the line
plot(ft, which = "residual")  #plots residuals
abline(a = 0, b = 0)
qqnorm(residuals(ft))

# hmm - outliers… do they have undue influence on the fit?
# You can fit SMATR lines that are less sensitive to outliers. These methods were developed by Sara Taskinen (U Jyvaskyla, Finland) using Huber's M estimation (fancy method which downweights outliers).

ftRobust = sma(brain ~ body, data = Animals, log = "xy", robust = T)
summary(ftRobust)

# note that while the interpretation hasn't changed here (2/3 still in interval) the fitted line is now much steeper and the CI is much narrower - it is getting thrown less by the outliers so gives a better estimate of the line

plot(ftRobust)

# SMATR TO COMPARE LINES Sometimes you have multiple lines that you want to compare. It might be of interest to compare slopes, or to assume common slope and to compare elevations, or compare locations along a common axis. A cool graph of this, courtesy of Dan Falster (Macquarie):

# Consider leaf longevity vs Leaf Mass per Area of a bunch of species at four sites with varying soil nutrients/rainfall:

data(leaflife)
ftComSlope = sma(longev ~ lma * site, log = "xy", data = leaflife)
summary(ftComSlope)
plot(ftComSlope)

#Is there evidence that the slope of the lma-longevity SMA varies across sites?
#Multiple comparisons - which sites differ from which in SMA slope? (using Sidak adjustment)

sma(longev ~ lma * site, log = "xy", data = leaflife, multcomp = T, multcompmethod = "adjust")
#Subset to sites 2 and 4 to compare elevation at these low soil P sites:
leaf.low.soilp <- subset(leaflife, soilp == "low")
ftElev <- sma(longev ~ lma + rain, log = "xy", data = leaf.low.soilp)
summary(ftElev)

#To test for evidence of shift along the SMA between sites:

ftShift <- sma(longev ~ lma + rain, log = "xy", type = "shift", data = leaf.low.soilp)
summary(ftShift)

#### Assignment - datasources - assess differences in scaling slopes and intercepts comparing differing groups.
### scaling metabolism, brain size, 



##############   REFERENCES  ############################

# Warton et al. (2012) SMATR 3 – an R package for estimation and inference about allometric lines. Methods in Ecology and Evolution 2012, 3, 257–259

#Warton et al (2006) Bivariate line-fitting methods for allometry, Biological Reviews http://onlinelibrary.wiley.com/doi/10.1017/S1464793106007007/abstract

#Smith (2009) Use and misuse of the reduced major axis for line-fitting, American Journal of Physical Anthropology http://onlinelibrary.wiley.com/doi/10.1002/ajpa.21090/abstract

#Warton et al (2012) smatr 3 - an R package for estimation and testing of allometric lines, Methods in Ecology and Evolution http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2011.00153.x/full

#Taskinen & Warton (2013) Robust tests for one or more allometric lines, Journal of Theoretical Biology http://www.sciencedirect.com/science/article/pii/S0022519313002257

# see also http://rstudio-pubs-static.s3.amazonaws.com/13754_c928a9debf464a0095a2cd1f5f988c0a.html


