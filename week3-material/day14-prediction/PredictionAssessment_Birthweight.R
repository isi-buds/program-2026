##
#####  Code for examples presented in ISI BUDS - Prediction Assessment
#####	 Author: D. Gillen
##
#####
#####	Source in needed libraries and functions
#####
##
library(MASS)
library(leaps)
##
#####
#####	Read in King county BW data
#####
##
path <- "/Users/dgillen/SynologyDrive/Teaching/ISI_BUDS/2022/Lecture4_Prediction2/"
weight <- read.table( paste(path,"KingCounty2001_data.txt",sep=""), header=TRUE )
names(weight)
## All records have 'plurality' = 1; so remove from dataset
weight <- weight[,-3]

##
#####
#####	Simple descriptives
#####
##
#####	Hist of birth weight
##
  hist(weight$bwt, xlab="", ylab="", main="", freq=FALSE, nclass=25, col="lightgrey", axes=FALSE)
  axis(1, at=seq(from=0, to=5000, by=1000))
  title(main="Birth weight, grams", cex.main=2, col.main="blue", font.main=1)

##
#####	Scatterplots
##
  par(mfrow=c(2,2), mar=c(3.1, 4.1, 3.1, 2.1))
  plot(weight$age, weight$bwt, xlab="", ylab="Birth Weight")
  title(main="Mother's age, years", cex.main=2, col.main="blue", font.main=1)
  plot(weight$smokeN, weight$bwt, xlab="", ylab="Birth Weight")
  title(main="Cigarettes smoked per day", cex.main=2, col.main="blue", font.main=1)
  plot(weight$drinkN, weight$bwt, xlab="", ylab="Birth Weight")
  title(main="Alchoholic drinks per week", cex.main=2, col.main="blue", font.main=1)
  plot(weight$educ, weight$bwt, xlab="", ylab="Birth Weight")
  title(main="Highest grade completed", cex.main=2, col.main="blue", font.main=1)



##
#####
#####	Best subsets regression (need to install 'leaps' library)
#####
##
## Model with only the intercept
##
fit0 <- lm(bwt ~ 1, data=weight)

## Perform best subsets analysis
##
maxModel <- as.formula(bwt ~ 	gender + age + race + parity + 
								married + smokeN + drinkN + 
								firstep + welfare + smoker + 
								drinker + wpre + educ)
fitF <- lm(maxModel, data=weight)
summary( fitF )

bestSubRSS  <- summary(regsubsets(maxModel, data=weight, nvmax=17, nbest=1))

## What is returned?
##
names(bestSubRSS)

## 'results' contains the subset size, k, and the residual sum of squares
##
results <- c(0, sum((weight$bwt- fitted(fit0))^2))
results <- rbind(results, cbind(apply(bestSubRSS$which, 1, sum)-1, bestSubRSS$rss))

##	##	Plot of RSS (naive) vs. subset size
  par(mfrow=c(1,1), mar=c(4.1, 4.1, 1, 1), cex=2)
  plot(results[,1], results[,2], type="n", axes=F, xlab="", ylab="", ylim=range(results[,2])*c(0.95, 1.05))
  points(results[,1], results[,2], pch=5)
  title(xlab="Subset size, k")
  mtext("Residual sum of squares", 2, cex=2, line=1)
  axis(1, at=c(0:17))
  lines(results[,1], results[,2], type="b", col="red", pch=20, lwd=3)


##
#####
#####	Now let's do best subsets with Cp as the criteria
#####
##
bestSubCp  <- leaps(x=model.matrix(fitF), y=weight$bwt, int=FALSE, nbest=1, method="Cp")
## 'results' contains the subset size, k, and the Cp value
##
results <- NULL
results <- rbind(results, cbind(apply(bestSubCp$which, 1, sum)-1, bestSubCp$Cp))

##	Plot of Mallow's Cp vs. subset size
  par(mar=c(4.1, 4.1, 1, 1), cex=2)
  plot(results[,1], results[,2], type="n", axes=F, xlab="", ylab="", ylim=range(results[,2])*c(0.95, 1.05))
  points(results[,1], results[,2], pch=5)
  title(xlab="Subset size, k")
  mtext("Mallow's Cp", 2, cex=2, line=1)
  axis(1, at=c(0:17))
  lines(results[,1], results[,2], type="b", col="red", pch=20, lwd=3)


##
#####
#####	 Which covariates were selected in the k=10 model...
#####
##
cbind( dimnames( model.matrix(fitF) )[[2]],bestSubCp$which[10,])

##
#####
#####	Stepwise selection using AIC
#####
##
fitStepAIC <- stepAIC( fit0, scope=maxModel, direction="forward" )
fitStepAIC
cbind( dimnames( model.matrix(fitF) )[[2]], bestSubCp$which[10,] )

##
#####
##### What is the cross-validated estimate of prediction error?
#####
##
##
#####		Computation of the CV or k-fold CV statistic for a lm fit (squared error loss)
##
cv.lm <- function( lmFit, data, K="n", GCV=FALSE ){
  y <- model.frame(lmFit)[,1]
  yhat <- lmFit$fitted
  n <- length(yhat)
  lmFormula <- formula( lmFit )
  
  Xmat <- model.matrix(lmFit)
  p <- dim( Xmat )[2]
  H <- Xmat %*% solve( t(Xmat)%*%Xmat ) %*% t(Xmat)	
  if( GCV==FALSE ) cv <- mean( ( (y-yhat) / (1-diag(H)) )^2 ) 
  else cv <- mean( ( (y-yhat) / (1-sum(diag(H))/n) )^2 ) 
  
  cv.k <- NULL
  if( K !="n" ) {
    ord <- sample(1:n, n)
    y <- y[ord]
    data <- data[ord,]
    rss.k <- rep(NA,n)
    for( i in 1:K ){
      keep <- 1:ceiling(n/K) + (i-1)*ceiling(n/K)
      if( max(keep) > n ) keep <- min(keep):n
      fit.k <- lm( lmFormula, data=data[!is.element( 1:n, keep ),] )
      yhat.k <- predict( fit.k, newdata=data[keep,] )
      rss.k[keep] <- (y[keep] - yhat.k)^2
    }
    cv.k <- mean( rss.k )
  }
  return( cbind(cv, cv.k) )
}

fit.10 <- lm( bwt ~ wpre + smoker + gender + married + race + parity + welfare + 
                    educ + smokeN, data=weight)
cv.lm(fit.10,weight,K=5)


