#read data
data <- read.table("/Users/apple/Desktop/ID5059/Practicals/P01/auto-mpg.data", sep=" ", header=T) 

attach(data)

#remove rows with NAs
cleandata <- na.omit(data)

# convert to comma separated values for use in Excel
write.csv(cleandata,"auto-data.csv")

#####1 Linear Models
#The function linreg performs a linear regression
#Input Arguments:
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot

linreg <- function(x, y, output=1, ploto=1, opt=0) 
{
  #Scale data to [0,1] 
  x<-x-min(x)
  x<-x/max(x)
  #Create a Design Matrix DM 
  n <- length(x)
  q <- 2
  DM = matrix(1,n,q) 
  DM[,2] <- x
  
  #Perform regression
  # beta <- solve(t(DM)%*%DM)%*%t(DM)%*%y 
  reg <- lm(y~0+DM)
  
  #Calculate goodness of fit measures
  #Residual sum of squares
  rss <-  sum(sapply(residuals(reg), function(x) { x^2 }))
  #Coefficient of determination: R^2
  R2 <- 1 - (rss/ (t(y)%*%y-(mean(y)**2*n)))
  #Adjusted Coefficient of determination: R^2
  R2adj <- 1 - ( (n-1)/(n-q) ) * (1-R2)
  #AIC
  aic <- AIC(reg)
  
  if(output==1)
  {
    #Summary output  
    cat("RSS: ", rss, "\n")
    cat("TSS: ", t(y)%*%y-(mean(y)**2/n), "\n")
    cat("R-squared: ", R2, "\n")
    cat("Adjusted R-squared: ", R2adj, "\n")
    cat("AIC: ", aic, "\n")
    #cat("Coefficents: \n")
    #print(coef(reg))
    #print(summary(reg))
    #print(anova(reg))
    
    #Graphic
    xp <- 0:100/100
    n <- length(xp)
    DMp = matrix(1,n,q)
    DMp[,2] <- xp
    
    if(ploto==1) par(mfrow=c(1,2))
    if(ploto==1) matplot(xp, DMp, type="l", lwd=2, main="Individual functions")
    if(ploto==1) plot(x,y, main="Linear regression", pch=20, col="darkblue")
    lines(xp, DMp%*%coef(reg), col="orange", type="l", lwd=2)
  }
  
  if(opt==1) {return(c(R2adj,aic))}   
}

#show the graphs and RSS of linear regression and calculate the MSE
##1 displacement
linreg(cleandata$displacement, cleandata$mpg)
mse_displacement1 <- 8378.822/length(cleandata$mpg)

##2 horsepower
linreg(cleandata$horsepower, cleandata$mpg)
mse_horsepower1 <- 9385.916/length(cleandata$mpg)

##3 weight
linreg(cleandata$weight, cleandata$mpg) 
mse_weight1 <- 7321.234/length(cleandata$mpg) 

##4 acceleration
linreg(cleandata$acceleration, cleandata$mpg)
mse_acceleration1 <- 19550.46/length(cleandata$mpg)


#mse of all four attributes
mse_linear <- c(mse_displacement1, mse_horsepower1, mse_weight1, mse_acceleration1)
mse_linear



#####2 Bin Smooths
#The function binsmoothREG performs a binsmooth regression with a user defined binlength
#Input Arguments:
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       binlength - amount of x values per bin
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot



binsmoothREG <- function(x, y, binlength=20, knotsdef=NULL, output=1, ploto=1, opt=0)
{
  #Scale data to [0,1]
  x<-x-min(x)
  x<-x/max(x)
  #Sort x values in ascending order
  y <- y[order(x)]
  x <- sort(x)
  n <- length(x)
  #Devide data into bins
  if(is.vector(knotsdef)) bins = knotsdef
  else bins = ceiling(length(x) / binlength)
  #Create Design Matrix without intercept
  DM <- matrix(1, length(x), bins)
  #Set all elements not corresponding to region j equal 0
  for(i in 1:bins)
  {
    if(i==1) { xstart = 1 }
    if(i>1) { xstart = (i-1)*binlength+1 }
    xend = min(xstart + binlength-1, length(x))
    binelements <- xstart:xend
    elements <- 1:length(x)
    elements[binelements] <- 0
    DM[elements,i] <- 0
  }
  
  #Perform Linear Regreesion
  reg <- lm(y~0+DM)
  
  #Calculate goodness of fit measures
  q <- bins
  #Residual sum of squares
  rss <-  sum(sapply(residuals(reg), function(x) { x^2 }))
  #Coefficient of determination: R^2
  R2 <- 1 - (rss/ (t(y)%*%y-(mean(y)**2*n)))
  #Adjusted Coefficient of determination: R^2
  R2adj <- 1 - ( (n-1)/(n-q) ) * (1-R2)   
  #AIC
  aic <- AIC(reg)
  
  if(output==1)
  {
    #Summary output  
    cat("Elements per bin: ", binlength, "\n")
    cat("Number of bins: ", bins, "\n")
    cat("RSS: ", rss, "\n")
    cat("TSS: ", t(y)%*%y-(mean(y)**2/n), "\n")
    cat("R-squared: ", R2, "\n")
    cat("Adjusted R-squared: ", R2adj, "\n")
    cat("AIC: ", aic, "\n")
    #cat("Coefficents: \n")
    #print(coef(reg))
    #print(summary(reg))
    #print(anova(reg))
    
    #Graphic 
    if(ploto==1) plot(x,y, main="Binsmooth regression", pch=20, col="darkblue")
    j<-1
    for(i in 1:length(coef(reg)))
    {
      if(i>1) lines(c(x[xend],x[xend]), c(as.numeric(coef(reg)[i-1]), as.numeric(coef(reg)[i])), col="red", lwd=2)
      xstart = j
      if(i>1) lines(c(x[xend],x[xstart]), c(as.numeric(coef(reg)[i]), as.numeric(coef(reg)[i])), col="red", lwd=2)
      xend = min(j+binlength-1, length(x))
      lines(c(x[xstart],x[xend]), rep(as.numeric(coef(reg)[i]), 2), col="red", lwd=2)
      j<-j+binlength
      
    }
  }
  
  if(opt==1) return(c(R2adj,aic))    
}


#show the graphs and RSS of bin smooths model and calculate the MSE
##1 displacement
binsmoothREG(cleandata$displacement, cleandata$mpg)
mse_displacement2 <- 6547.693/length(cleandata$mpg)

##2 horsepower
binsmoothREG(cleandata$horsepower, cleandata$mpg)
mse_horsepower2 <- 7004.712/length(cleandata$mpg)

##3 weight
binsmoothREG(cleandata$weight, cleandata$mpg) 
mse_weight2 <- 6530.837/length(cleandata$mpg) 

##4 acceleration
binsmoothREG(cleandata$acceleration, cleandata$mpg)
mse_acceleration2 <- 17757.41/length(cleandata$mpg)


#mse of all four attributes
mse_binsmooths <- c(mse_displacement2, mse_horsepower2, mse_weight2, mse_acceleration2)
mse_binsmooths



#####3 b-spline bases
#The function bsplinereg performs a bspline regression with user defined knots
#Input Arguments:
#Input Arguments bsplinereg(...):
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       knots - number of knots in [0,1]
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot

#Calculate basis (rekursiv)
basis <- function(x, degree, i, knots) 
{ 
  
  if(degree == 0)
  { B <- ifelse((x >= knots[i]) & (x < knots[i+1]), 1, 0) 
  } else { 
    if((knots[degree+i] - knots[i]) == 0) 
    { alpha1 <- 0 
    } else { 
      alpha1 <- (x - knots[i])/(knots[degree+i] - knots[i]) } 
    
    if((knots[i+degree+1] - knots[i+1]) == 0) 
    { alpha2 <- 0 
    } else { alpha2 <- (knots[i+degree+1] - x)/(knots[i+degree+1] - knots[i+1]) } 
    B <- alpha1*basis(x, (degree-1), i, knots) + alpha2*basis(x, (degree-1), (i+1), knots) 
  } 
  
  return(B) 
}

#Create bspline Desin Matrix
bspline <- function(x, degree, knotpos) 
{ 
  #Number of basis
  K <- length(knotpos) + degree + 1 
  #Number of observations
  n <- length(x)
  #Set Boundary knots
  Boundary.knots = c(0,1)
  #create new vector with knot positons 
  knotpos <- c(rep(Boundary.knots[1], (degree+1)), knotpos, rep(Boundary.knots[2], (degree+1))) 
  
  
  #Create design matrix
  DM <- matrix(0,n,K) 
  for(j in 1:K) DM[,j] <- basis(x, degree, j, knotpos) 
  if(any(x == Boundary.knots[2])) DM[x == Boundary.knots[2], K] <- 1 
  #Return DM  
  return(DM) 
}


bsplinereg <- function(x, y, knots=0, knotsdef=NULL,  degree, output=1, ploto=1, opt=0)
{
  #Scale data to [0,1]
  x<-x-min(x)
  x<-x/max(x)
  #Sort x values in ascending order
  y <- y[order(x)]
  x <- sort(x)
  n <- length(x)
  
  #Calculate knot postions
  if(knots == 0) knotpos <- NULL
  if(knots != 0) knotpos <- 1:knots / (knots+1)
  if(length(knotsdef)>0) knotpos <- knotsdef
  #Create Design Matrix 
  DM <- bspline(x, degree, knotpos) 
  
  #Perform penalized regression
  reg <- lm(y ~ 0 + DM)
  print(summary(reg))
  
  #Calculate goodness of fit measures
  q <- length(knotpos) + degree + 1
  #Residual sum of squares
  rss <-  sum(sapply(residuals(reg), function(x) { x^2 }))
  #Coefficient of determination: R^2
  R2 <- 1 - (rss/ (t(y)%*%y-(mean(y)**2*n)))
  #Adjusted Coefficient of determination: R^2
  R2adj <- 1 - ( (n-1)/(n-q) ) * (1-R2)   
  #AIC
  aic <- AIC(reg)
  
  if(output==1)
  {
    #Summary output  
    cat("Number of knots = ", knots, "\n")
    cat("Knot positions = ", knotpos, "\n")
    cat("RSS: ", rss, "\n")
    cat("TSS: ", t(y)%*%y-(mean(y)**2/n), "\n")
    cat("R-squared: ", R2, "\n")
    cat("Adjusted R-squared: ", R2adj, "\n")
    cat("AIC: ", aic , "\n")
    #cat("Coefficents: \n")
    #print(coef(reg))
    #print(summary(reg))
    #print(anova(reg))
    
    #Graphics
    
    #Values for prediction
    xp <- 0:100/100
    DMp <- bspline(xp, degree, knotpos)
    
    if(ploto==1)par(mfrow=c(1,2))
    if(ploto==1) matplot(xp, (DMp), type="l", lwd=2, main="Individual spline functions")
    if(ploto==1) for(i in 1:length(knotpos)) abline(v=knotpos[i], col="red", lty=2) 
    if(ploto==1) plot(x,y, main="BSpline Regression", pch=20, col="darkblue")
    
    points(xp,DMp%*%coef(reg), type="l", lwd=2, col="brown")
    if(ploto==1) for(i in 1:length(knotpos)) abline(v=knotpos[i], col="red", lty=2) 
  }
  if(opt==1) return(c(R2adj, aic))       
}


#plot the relationship of mpg and x variables for suitable knots and degree
plot(cleandata$mpg~cleandata$displacement)
plot(cleandata$mpg~cleandata$horsepower)
plot(cleandata$mpg~cleandata$weight)
plot(cleandata$mpg~cleandata$acceleration)
#knots should be 2

#show the graphs and RSS of b-spline and calculate the MSE
##1 displacement
bsplinereg(cleandata$displacement, cleandata$mpg, degree = 8, knots = 2) 
mse_displacement3 <- 6608.903/length(cleandata$mpg)

##2 horsepower
bsplinereg(cleandata$horsepower, cleandata$mpg, degree = 8, knots = 2)
mse_horsepower3 <- 7062.692/length(cleandata$mpg)

##3 weight
bsplinereg(cleandata$weight, cleandata$mpg, degree = 8, knots = 2)
mse_weight3 <- 6669.5/length(cleandata$mpg)

##4 acceleration
bsplinereg(cleandata$acceleration, cleandata$mpg, degree = 5, knots = 2)
mse_acceleration3 <- 18626.84/length(cleandata$mpg)


#mse of all four attributes
mse_bspline <- c(mse_displacement3, mse_horsepower3, mse_weight3, mse_acceleration3)
mse_bspline
