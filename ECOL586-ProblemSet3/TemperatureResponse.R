#=========================================================================================
# ECOL586: Plotting and fitting temperature response curves with Boltzmann-Arrhenius and
# second-order polynomials
#=========================================================================================
# Prepared by sean.michaletz@gmail.com



#=========================================================================================
# Contents
#=========================================================================================

#  Part 1:  Plotting temperature response curves with ggplot2
#  Part 2:  Modified Arrhenius plots and Boltzmann-Arrhenius model fitting


#=========================================================================================
# Introductory code: Load packages/libraries, set paths, load datasets
#=========================================================================================

#--Clear memory
rm(list=ls(all=T))

# << DIRECTORIES >>
#--Set working directory
setwd("E:/")

# << DATASETS >>
#--Load temperature response data.
df <- read.csv("temperatureResponse.csv", header=T)

# << PACKAGES >>
library(ggplot2)
library(polynom)
library(scales)


#=========================================================================================
# Part 1: Plotting temperature response curves with ggplot2
#=========================================================================================

#--Plot all curves on a single set of axes.
ggplot(df, aes(x=T_C, y=Rate)) + 
  geom_point() + 
  xlab(expression('Temperature'~(degree*C))) + 
  ylab(expression("Rate")) +
  theme_bw(base_size=12) + 
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_colour_manual(values=rep("gray60", times=197))

#--Plot with fitted second-order polynomials
ggplot(df, aes(x=T_C, y=Rate)) + 
  stat_smooth(method="lm", se=TRUE, fill=NA, formula=y ~ poly(x,2,raw=TRUE), 
              aes(colour=factor(curveID)), data=df) +
  geom_point() + 
  xlab(expression('Temperature'~(degree*C))) + 
  ylab(expression("Rate")) +
  theme_bw(base_size=12) + 
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_colour_manual(values=rep("black", times=197))

#--Plot a single curve (e.g. curve 1).
ggplot(subset(df, df$curveID=="1"), aes(x=T_C, y=Rate)) + 
  stat_smooth(method="lm", se=TRUE, fill=NA, formula=y ~ poly(x,2,raw=TRUE),aes(colour=factor(curveID))) +
  geom_point(aes(color=as.factor(curveID))) + 
  xlab(expression('Temperature'~(degree*C))) + 
  ylab(expression("Rate")) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))


# Fit second-order polynomial to a curve (e.g. curveID 1)
model1 <- lm(Rate ~ poly(T_C, 2, raw=TRUE), data=subset(df, df$curveID=="1"))
summary(model1)
 

#=========================================================================================
# Part 2: Modified Arrhenius plots and Boltzmann-Arrhenius model fitting
#=========================================================================================

#--Calculate the inverse temperature 1/kT - note T is converted from C to K.
df$invT <- 1/(8.617e-5*(df$T_C + 273.15))

#--Plot all curves on a single set of axes.
ggplot(df, aes(x=invT, y=Rate)) +
  geom_point(color="black") +                           
  xlab(expression(paste('Temperature ', '1/',italic('kT'), ' (',  eV^{-1}, ')'))) + 
  ylab(expression("Rate")) +
  scale_x_continuous(sec.axis = sec_axis(trans = ~ (1/(.*0.00008617))-273.15 , name = expression(paste("Temperature (", degree, C, ")")))) +
  scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3), 
                     labels = trans_format("log", math_format(e^.x))) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))


#--Plot a single curve (e.g., curveID = 2).
ggplot(subset(df, df$curveID=="2"), aes(invT, Rate, group=factor(curveID))) +
  geom_point(color="black") +                           
  xlab(expression(paste('Temperature ', '1/',italic('kT'), ' (',  eV^{-1}, ')'))) + 
  ylab(expression("Rate")) +
  scale_x_continuous(sec.axis = sec_axis(trans = ~ (1/(.*0.00008617))-273.15 , name = expression(paste("Temperature (", degree, C, ")")))) +
  scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3), 
                     labels = trans_format("log", math_format(e^.x))) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

#--Plot curves 1 to 2 on individual axes.
for(i in seq_along(1:2)) {
  print(
    ggplot(subset(df, df$curveID==i), aes(invT, Rate, group=factor(curveID))) +
      geom_point(color="black") +                           
      xlab(expression(paste('Temperature ', '1/',italic('kT'), ' (',  eV^{-1}, ')'))) + 
      ylab(expression("Rate")) +
      scale_x_continuous(sec.axis = sec_axis(trans = ~ (1/(.*0.00008617))-273.15 , name = expression(paste("Temperature (", degree, C, ")")))) +
      scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3), 
                         labels = trans_format("log", math_format(e^.x))) +
      theme_bw(base_size=12) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"))
  )
}

# Fitting the Boltzmann-Arrhenius model to a single curve (e.g. curveID 1) to estimate 
# the apparent activation energy E (eV). Note that since B = B_0*e^(-E/kT), the correct 
# value of E requires reversing the sign on the invT coefficient
model2 <- lm(log(Rate) ~ invT,  data=subset(df, df$curveID=="1"))
summary(model2)
confint(model2)




