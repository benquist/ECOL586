
########################################
#
# ECOL 586 Big Trees analysis
#
#  B.J. Enquist January 2018
#########################################

bigtrees.data <- read.csv(file="/Users/brianjenquist/GitHub/R/ECOL586_Biological_Scaling/bigtrees.csv",header=T)

names(bigtrees.data)
hist(bigtrees.data $Diameter..in.)
hist(log(bigtrees.data $Diameter..in.))
hist(log10(bigtrees.data $Diameter..in.))

plot((Height.feet)~(Diameter..in.),bigtrees.data)
plot(log(Height.feet)~log(Diameter.in),bigtrees.data)
plot(log(Height.feet)~log(Diameter..in.),bigtrees.data)

summary(lm(log(Height.feet)~log(Diameter..in.),bigtrees.data))
plot((log(Height.feet))~(log(Diameter..in.)),bigtrees.data)
lm.HD_Line <-lm(lm((log(Height.feet))~(log(Diameter..in.)),bigtrees.data))
abline(lm.HD_Line)

HD_ratio <-((bigtrees.data$Height.feet)/(bigtrees.data$Diameter..in.))
plot(HD_ratio~(log(Diameter..in.)), bigtrees.data)
plot(log10(HD_ratio)~(log10(Diameter..in.)), bigtrees.data)
plot(log10(bigtrees.data$Diameter..in.),log10(bigtrees.data$Height.feet),
     xlab="log(diameter (in))",ylab="log(height (in))")


source("/Users/brianjenquist/GitHub/R/ECOL586_Biological_Scaling/slopes.R")
slopes


log_height <-(log10(bigtrees.data$Height.feet))
log_diameter <-(log10(bigtrees.data$Diameter..in.))
btreg <- slopes(log10(bigtrees.data$Diameter..in.),log10(bigtrees.data$Height.feet))
options(max.print=999999)
btreg

#Let’s put the OLS (Y |X) red line on the graph
plot(log10(bigtrees.data$Diameter..in.),log10(bigtrees.data$Height.feet),xlab="log(diameter (in))",ylab="log(height (in))")
abline(btreg$parameters[1,1],btreg$parameters[1,3],lty=1, col=2)

#Let’s put the OLS (X |Y) line on the graph
abline(btreg$parameters[2,1],btreg$parameters[2,3],lty=1)

#Let’s put the RMA line (dashed line) on the graph
abline(btreg$parameters[5,1],btreg$parameters[5,3],lty=2)

#Here is the OLS bisector line (dashed blue line) on the graph
abline(btreg$parameters[3,1],btreg$parameters[3,3],lty=2, col=4)

names(btreg)

# let's extract out the parameters for the intercept (A) and the slope (B) for each of the regression models
btreg$parameters
# the correlation 
btreg$correlation



