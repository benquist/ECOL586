
############################################################################################################
##  Slopes routine - Brian J. Enquist Jan 23 2018 with additions from Brian Maitner
##
## Originally written in S+ by Andrew P. Allen at the University of New Mexico, borrowing heavily from a 
## Fortran program written by Isobe et al. This program computes estimates and standard deviations for the 
## slopes and intercepts of six different regression models. Details and methodological considerations are 
## available in Isobe et al. 1990.
#############################################################################################################

# ##Examples of use of slopes.R 
# > source("/Users/yourdirectorypath/slopes.R")
# > slopes
#  ## let's use the slopes.R functions on a dataset with your X and Y variables of interest
# > btreg <- slopes(log10(dataset$X),log10(dataset$Y))
# > btreg
# > names(btreg)
# ## let's extract out the parameters for the intercept (A) and the slope (B) for each of the regression models
# > btreg$parameters
# ## the correlation 
# > btreg$correlation

##########################################

# slopes<-function(x,y) {
  #COMPUTE AVERAGES AND SUMS
  XAVG <- mean(x)
  YAVG <- mean(y)
  RN <- length(x)
  SXX <- sum((x-mean(x))^2)
  SYY <- sum((y-mean(y))^2)
  SXY <- sum((x-mean(x))*(y-mean(y)))
  
  #COMPUTE THE SLOPE COEFFICIENTS
  B <- 0
  B[1] <- SXY/SXX
  B[2] <- SYY/SXY
  B[3] <- (B[1]*B[2] - 1.0 + sqrt((1.0 + B[1]^2)*(1.0 + B[2]^2)))/(B[1] + B[2])
  B[4] <- 0.5*(B[2] - 1.0/B[1] + sign(SXY)*sqrt(4.0 + (B[2] - 1.0/B[1])^2))
  B[5] <- sign(SXY)*sqrt(B[1]*B[2])
  B[6] <- 0.5*(B[1] + B[2])
  
  #COMPUTE INTERCEPT COEFFICIENTS
  A <- 0
  for(i in 1:6) {
    A[i] <- YAVG - B[i]*XAVG
  }
  
  #PREPARE FOR COMPUTATION OF VARIANCES
  x <- x - XAVG
  y <- y - YAVG
  GAM1 <- B[3]/((B[1] + B[2])*sqrt((1.0 + B[1]^2)*(1.0 + B[2]^2)))
  GAM2 <- B[4]/(sqrt(4.0*B[1]^2 + (B[1]*B[2] - 1.0)^2))
  SUM1 <- sum((x*(y - B[1]*x))^2)
  SUM2 <- sum((y*(y - B[2]*x))^2)
  SUM3 <- sum(x*y*(y - B[1]*x)*(y - B[2]*x))
  COV <- SUM3/(B[1]*SXX^2)
  
  #COMPUTE VARIANCES OF THE SLOPE COEFFICIENTS
  SIGB <- 0
  SIGB[1] <- SUM1/(SXX^2)
  SIGB[2] <- SUM2/(SXY^2)
  SIGB[3] <- (GAM1^2)*(((1.0 + B[2]^2)^2)*SIGB[1] + 2.0*(1.0 + B[1]^2)*(1.0 + B[2]^2)*COV + ((1.0+B[1]^2)^2)*SIGB[2])
  SIGB[4] <- (GAM2^2)*(SIGB[1]/B[1]^2 + 2.0*COV + B[1]^2*SIGB[2])
  SIGB[5] <- 0.25*(B[2]*SIGB[1]/B[1] + 2.0*COV + B[1]*SIGB[2]/B[2])
  SIGB[6] <- 0.25*(SIGB[1] + 2.0*COV + SIGB[2])
  
  #COMPUTE STANDARD DEVIATIONS OF THE INTERCEPT COEFFICIENTS
  SIGA <- 0
  SIGA[1] <- sum(((y - B[1]*x)*(1.0 - RN*XAVG*x/SXX))^2)
  SIGA[2] <- sum(((y - B[2]*x)*(1.0 - RN*XAVG*y/SXY))^2)
  SIGA[3] <- sum(((x*(y - B[1]*x)*(1.0 + B[2]^2)/SXX + y*(y - B[2]*x)*(1.0 + B[1]^2)/SXY)*GAM1*XAVG*RN - y + B[3]*x)^2)
  SIGA[4] <- sum(((x*(y - B[1]*x)/SXX + y*(y - B[2]*x)*(B[1]^2)/SXY)*GAM2*XAVG*RN/sqrt(B[1]^2) - y + B[4]*x)^2)
  SIGA[5] <- sum(((x*(y - B[1]*x)*sqrt(B[2]/B[1])/SXX + y*(y - B[2]*x)*sqrt(B[1]/B[2])/SXY)*0.5*RN*XAVG - y + B[5]*x)^2)
  SIGA[6] <- sum(((x*(y - B[1]*x)/SXX + y*(y - B[2]*x)/SXY)*0.5*RN*XAVG - y + B[6]*x)^2)
  
  #generates R object
  correlation<-cor(x,y)
  parameters <- matrix(data=NA,nrow=6,ncol=4,dimnames=list(c('OLS(Y/X)','OLS(X/Y)','OLS BISECTOR','ORTHOGONAL','REDUCED MAJ AXIS','MEAN OLS'),
                                                           c('A','SD(A)','B','SD(B)')))
  parameters[,1] <- A
  parameters[,2] <- sqrt(SIGA)/RN
  parameters[,3] <- B
  parameters[,4] <- sqrt(SIGB)
  x <- x + XAVG
  y <- y + YAVG
  
  #return(x,y,correlation,parameters)
  #return(c(x,y,correlation,parameters))
  
  output<-list(x,y,correlation,parameters)
  
  names(output)<-c("x","y","correlation","parameters")
  
  return(output)
}