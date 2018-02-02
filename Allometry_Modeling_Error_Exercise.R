
########################################################################
##  Brian J. Enquist, Jan 2018
##  Assessing Regression models with differig leves of error
##  Building on exercise by Steven Holland http://strata.uga.edu/8370/lecturenotes/regression.html
##
########################################################################

# lets generate a distribution of x values with a given mean value but differing variance or standard deviation
# play with altering the variance in x with differing sd values

x <- rnorm(n=30, mean=15, sd=4)
hist(x, breaks=5)

# let's jazz up the plots and make plots comparable
hist(x, xlab="x variable", ylab="Frequency of x", 
     border="blue", col="green",
     xlim=c(1,30), ylim=c(0,20), # put some limits on the range of x and y variables
     breaks = 10   # number of breaks in the histogram
     )

#now let's create some error that we can model in our regression
errors <- rnorm(n=30, sd=1)
slope <- 2.15   # the true slope of the allometric relationship
intercept <- 4.64 # the true intercept of the allometric relationship

# let's now model our y variable of interest via the regression equation 

y <- slope*x + intercept + errors

plot(x, y, pch=16)

# let's fit the Model I Ordinary Least Squares Regression

myRegression <- lm(y ~ x)

# object contains the statistics for the y-intercept (called Intercept) and slope (called x here, but in general is called whatever the predictor variable is named).
myRegression
names(myRegression)
summary(myRegression)

# do the fitted intercept and slope differ from the 'true'

# Hypothesis tests - You can perform hypothesis tests on the regression with the summary() function.
summary(myRegression)

# Next shown are the residuals, and here you can evaluate the assumption of whether the residuals are normally distributed. 

confint(myRegression)

plot(myRegression)

#The first plot shows residuals (ε) against the fitted values (Ŷ). Two aspects are critical. First, the residuals should be centered around zero (shown by the gray dashed line) for all fitted values; the red line showing the trend of the residuals should be close to the gray dashed line. Second, the variance of the residuals should not change with fitted values: the cloud of points should not widen to the left or to the right.

#The second plot is a qqnorm() plot. The residuals are assumed to be normally distributed, so they should fall fairly close to the gray dashed line on this plot, a line that indicates perfect agreement between the theoretical and observed quantiles. This match should be especially good in the middle part of the range, as it is common for points near the ends of the range to lie off the dashed line. Strong systematic departures indicate non-normally distributed errors.

#The third plot focuses on trends in the sizes of the residuals. The square root of the absolute value of the residuals is shown, which allows one to see trends in the residuals with the fitted values. The residuals should show no trend (i.e., the red line should be relatively flat), and the residuals should not form a triangular region, with an increasing range in the size of residuals as the fitted values increase.

#The fourth plot is the standardized residuals versus leverage, giving insight as to which data points have the greatest influence on the estimates of slope and intercept. Note that this is a different plot than indicated in the Crawley text (p. 144), but the goal is the same. Standardized residuals should show no trend with leverage (i.e., the solid red line should lie close to the dotted gray line). Furthermore, no points should have large values of Cook’s distance, which would indicate that an outlier (a point with a large residual) also has strong leverage (ability to pull the regression line in some direction).

# Adding the regression to a plot
plot(x, y, pch=16)
abline(myRegression, col='red', lwd=2)

# Prediction using the regression - In some cases, you would like to use the regression for prediction, that is, you know the independent (predictor) variable and you would like to know the value of the dependent (result) variable. To do this, you need to supply the regression model and a list of values of the independent variable for which you would like predictions.
predict(myRegression, list(x=c(12.2, 13.7, 14.45)))


######   Model 2 regression  ######
# The SMA y-intercept is calculated as it is for the least-squares regression, that is, the line must pass through the centroid.
# Functions for the SMA slope and intercept are straightforward.
# here we can use our slopes.S routine or if you would like just calculate the model II regression by hand

smaSlope <- function(x,y) {
  b1 <- sd(y)/sd(x)
  if (cor(x,y)<0) b1 <- -b1 # change sign for negative correlation
  b1
}

smaIntercept <- function(x,y) {
  b1 <- sd(y)/sd(x)
  if (cor(x,y)<0) b1 <- -b1 # change sign for negative correlation
  b0 <- mean(y) - mean(x)*b1
  b0
}
smaSlope(x,y)
smaIntercept(x,y)



