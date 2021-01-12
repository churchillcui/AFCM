## This script shows how to fit Additive Functional Cox Model (AFCM) (Cui et al. 2021) and
## Linear Functional Cox Model (LFCM) (Gellar et al. 2015) using mgcv package in R.

## Author: Erjia Cui
## Date: 01/12/2021

library(mgcv)
library(refund)

###########################################################################
# Simulate a dataset
###########################################################################
## The dataset "data_analysis" consists of:
##  event: a binary indicator of the event.
##  survtime: survival time (in year, for example)
##  X: functional predictor stored in a matrix, each row is the functional observation of one subject
##  Z: other predictor. we only include one other predictor in this example

N <- 2000 ## number of subjects
S <- 1000 ## number of functional observations per subject
event <- rbinom(N, 1, 0.3)
survtime <- runif(N, 0, 10)
X <- matrix(rnorm(N*S), nrow = N, ncol = S) ## in practice, transformations on X may be necessary for identifiability
X <- fpca.face(X)$Yhat ## smooth each curve using fast sandwich smoother (Xiao et al. 2013, 2016)
Z <- rnorm(N, 1, 1)

data_analysis <- data.frame(event, survtime, Z, X = I(X))
rm(event, survtime, X, Z)

###########################################################################
# Fit AFCM and LFCM using mgcv
###########################################################################
## Add variables related to numerical integration
### lmat: numeric integration of the functional term
data_analysis$lmat <- I(matrix(1/S, ncol=S, nrow=nrow(data_analysis)))
### tmat: time indices, we assume an equally-spaced grid of [0, 1] for example
data_analysis$tmat <- I(matrix(seq(0, 1, len=S), ncol=S, nrow=nrow(data_analysis), byrow=TRUE))

## Fit AFCM, using "ti" function to specify additional identifiability constraints in "mc" parameter
fit_afcm <- gam(survtime ~ Z + ti(tmat, X, by=lmat, bs=c("cr","cr"), k=c(10,10), mc=c(FALSE,TRUE)), 
                weights=event, data=data_analysis, family=cox.ph())

## Fit LFCM
fit_lfcm <- gam(survtime ~ Z + s(tmat, by=lmat*X, bs="cr", k=10), 
                weights=event, data=data_analysis, family=cox.ph())

###########################################################################
# Visualize the estimates
###########################################################################
par(mfrow = c(1,2))
vis.gam(fit_afcm, view = c("tmat", "X"), plot.type = "contour", color = "topo", main = "Estimates from AFCM")
vis.gam(fit_lfcm, view = c("tmat", "X"), plot.type = "contour", color = "topo", main = "Estimates from LFCM")

