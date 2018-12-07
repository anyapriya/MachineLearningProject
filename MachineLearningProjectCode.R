library(matrixcalc)
library(MASS)
library(ggplot2)





#covMatrix <- diag(corrMatrix) * corrMatrix * diag(corrMatrix)
#covMatrix
#cov2cor(covMatrix)



Simulation <- function(times, obs, test.valid.size, p, rho, errorSD, betaVals){
  Output <- NULL
  for(i in 1:times){
    Output <-rbind(Output, Testing(n = obs, sizes = test.valid.size, p = p, rho = rho, errorSD = errorSD, betaVals = betaVals))
  }
  return(Output)
}


#Called by Simulation, runs one round of the simulation, returns the number of non-zero coefficients and MSE from optimal lambda that minimizes MSE
Testing <- function(n, sizes, p, rho, errorSD, betaVals){
  SampleOut <- mySample(n, p, rho, errorSD, betaVals)
  x <- SampleOut[[1]]
  y <- SampleOut[[2]]
  Output <- Validation(x = x, y = y, sizes = sizes)
  NumbNonZeros.Lasso <- sum(Output[[1]] != 0)
  NumbNonZeros.ElNet <- sum(Output[[3]] != 0)
  MSE.Lasso <- Output[[2]]
  MSE.ElNet <- Output[[4]]
  return(c(NumbNonZeros.Lasso, MSE.Lasso, NumbNonZeros.ElNet, MSE.ElNet))
}

        #Called by Testing, randomly generates the data set
        mySample <- function(n, p, rho, errorSD, betaVals){
          SigmaMat <- covMatrix(rho, p)
          mu <- rep(0,p)
          x <- mvrnorm(n = n, mu = mu, Sigma = SigmaMat)
          error <- rnorm(n)
          y <- x %*% betaVals + errorSD*error
          return(list(x,y))
        }
        
                  #Called by mySample, creates variance-covariance matrix
                  covMatrix <- function(rho, p){
                    covMatrix <- matrix(rep(0,(p^2)), ncol = p)
                    for(i in 1:p){
                      for(j in 1:p){
                        covMatrix[i,j] <- rho^abs(i-j)
                      }
                    }
                    return(covMatrix)
                  }

#Called by Testing, finds the value of lambda to minimize the MSE and returns the BetaHats and MSE
Validation <- function(x,y,sizes){
  testSize <- sizes[1]
  
  MinMSE.Lasso <- optim(c(0.0001), MSEValid.Lasso,  x = x, y = y, sizes = sizes, method = "Brent", lower = 0, upper = 5)
  BetaHat.Lasso <- coordDescAlg.Lasso(MinMSE.Lasso$par, y[1:testSize], x[1:testSize,])
  
  MinMSE.ElasticNet <- optim(c(0.0001, 0.0001), MSEValid.ElasticNet,  x = x, y = y, sizes = sizes)
  BetaHat.ElasticNet <- coordDescAlg.ElNet(MinMSE.ElasticNet$par, y[1:testSize], x[1:testSize,])
  
  return(list(BetaHat.Lasso,MinMSE.Lasso$value, BetaHat.ElasticNet, MinMSE.ElasticNet$value))
}



#Called by Validation, calculates the value of MSE for a given lambda
MSEValid.Lasso <- function(lambda, x, y, sizes){
  testSize <- sizes[1]
  validSize <- sizes[2]
  ytrain = y[1:testSize]
  xtrain = x[1:testSize,]
  yvalid = y[(testSize+1):(testSize + validSize)]
  xvalid = x[(testSize+1):(testSize + validSize),]
  BetaHat <- coordDescAlg.Lasso(lambda, ytrain, xtrain)
  MSE <- mean((yvalid - xvalid %*% BetaHat)^2)
  return(MSE)
}


#Called by Validation, calculates the value of MSE for a given lambda
MSEValid.ElasticNet <- function(lambda, x, y, sizes){
  testSize <- sizes[1]
  validSize <- sizes[2]
  ytrain = y[1:testSize]
  xtrain = x[1:testSize,]
  yvalid = y[(testSize+1):(testSize + validSize)]
  xvalid = x[(testSize+1):(testSize + validSize),]
  BetaHat <- coordDescAlg.ElNet(lambda, ytrain, xtrain)
  MSE <- mean((yvalid - xvalid %*% BetaHat)^2)
  return(MSE)
}




#Called by MSEValid and Validation, uses the coordinate descent algorithm to calculate beta for the a given lambda
coordDescAlg.Lasso <- function(lambda, y, x){
  
  x <- as.matrix(x)
  converge <- FALSE
  n <- dim(x)[1]
  p <- dim(x)[2]
  BetaHat <- rep(0, p)
  BetaHat.star <- rep(0, p)
  newBetaHat <- rep(0,p)
  r <- matrix(rep(0,n*p), ncol = p)
  
  while(!converge){
    for(j in 1:p){
      r[,j] <- y - x[,-j]%*%BetaHat[-j]
      
      BetaHat.star[j] <- (t(x[,j])%*% r[,j]) / (t(x[,j])%*% x[,j])
      newBetaHat[j] <- sign(BetaHat.star[j])*pmax(BetaHat.star[j] - n*(lambda / (t(x[,j])%*% x[,j])), 0)
      
      if(all(round(BetaHat, 5) == round(newBetaHat,5))){
        converge <- TRUE
        break
      } else {
        BetaHat <- newBetaHat
      }
      
      
    }
    
    
  }
  
  return(BetaHat)
  
}























coordDescAlg.ElNet <- function(lambda, y, x){
  
  x <- as.matrix(x)
  converge <- FALSE
  n <- dim(x)[1]
  p <- dim(x)[2]
  BetaHat <- rep(0, p)
  BetaHat.star <- rep(0, p)
  newBetaHat <- rep(0,p)
  r <- matrix(rep(0,n*p), ncol = p)
  
  while(!converge){
    for(j in 1:p){
      r[,j] <- y - x[,-j]%*%BetaHat[-j]
      
      BetaHat.star[j] <- (t(x[,j])%*% r[,j]) / (t(x[,j])%*% x[,j])
      newBetaHat[j] <- (sign(BetaHat.star[j])*pmax(BetaHat.star[j] - (lambda[1] / (t(x[,j])%*% x[,j])), 0)) / ((1 + 2*lambda[2] ))
      
      
      if(all(round(BetaHat, 5) == round(newBetaHat,5))){
        converge <- TRUE
        break
      } else {
        BetaHat <- newBetaHat
      }
      
    }
  }
  
  return(BetaHat)
  
}









#Baseline
output <- Simulation(times = 50, obs = 240, test.valid.size = c(20,20), p = 8, rho = 0.5, errorSD = 3, betaVals = matrix(c(3,1.5,0,0,2,0,0,0), ncol = 1))


mysummary <- cbind(colMeans(output),apply(output, 2, var))
rownames(mysummary) <- c("NonZeroCoefs.Lasso", "MSE.Lasso", "NonZeroCoefs.ElNet", "MSE.ElNet")
colnames(mysummary) <- c("Mean", "Var")
mysummary




meanAndvar <- function(x){
  cbind(colMeans(x),apply(x, 2, var))
}

#Change number of simulations
NSims <- c(3:10)

output <- lapply(NSims, Simulation, obs = 240, test.valid.size = c(20,20), p = 8, rho = 0.5, errorSD = 3, betaVals = matrix(c(3,1.5,0,0,2,0,0,0), ncol = 1))
mysummary <- lapply(output, meanAndvar)

NonZeroCoefsSummary.Lasso <- data.frame("NSims" = NSims, Mod = "Lasso", "Mean" = sapply(mysummary, "[[", 1,1), "Var" = sapply(mysummary, "[[", 1,2))
NonZeroCoefsSummary.ElNet <- data.frame("NSims" = NSims, Mod = "ElNet", "Mean" = sapply(mysummary, "[[", 3,1), "Var" = sapply(mysummary, "[[", 3,2))
NonZeroCoefsSummary <- rbind(NonZeroCoefsSummary.Lasso, NonZeroCoefsSummary.ElNet)

MSEsummary.Lasso <- data.frame("NSims" = NSims, Mod = "Lasso", "Mean" = sapply(mysummary, "[[", 2,1), "Var" = sapply(mysummary, "[[", 2,2))
MSEsummary.ElNet <- data.frame("NSims" = NSims, Mod = "ElNet", "Mean" = sapply(mysummary, "[[", 4,1), "Var" = sapply(mysummary, "[[", 4,2))
MSEsummary <- rbind(MSEsummary.Lasso, MSEsummary.ElNet)

ggplot(NonZeroCoefsSummary) + ggtitle("Non-Zero Coefficients by Number of Simulations") + 
  geom_line(aes(x = NSims, y = Mean, color = Mod)) + 
  geom_ribbon(aes(x= NSims, ymin = Mean - sqrt(Var), ymax = Mean + sqrt(Var), fill = Mod), alpha = 0.3)

ggplot(MSEsummary) + ggtitle("MSE by Number of Simulations") + 
  geom_line(aes(x = NSims, y = Mean, color = Mod)) + 
  geom_ribbon(aes(x= NSims, ymin = Mean - sqrt(Var), ymax = Mean + sqrt(Var), fill = Mod), alpha = 0.3)



#Change observation number
Nobs <- c(50:2000)

output <- lapply(Nobs, Simulation, times = 50, test.valid.size = c(20,20), p = 8, rho = 0.5, errorSD = 3, betaVals = matrix(c(3,1.5,0,0,2,0,0,0), ncol = 1))
mysummary <- lapply(output, meanAndvar)

NonZeroCoefsSummary.Lasso <- data.frame("Nobs" = Nobs, Mod = "Lasso", "Mean" = sapply(mysummary, "[[", 1,1), "Var" = sapply(mysummary, "[[", 1,2))
NonZeroCoefsSummary.ElNet <- data.frame("Nobs" = Nobs, Mod = "ElNet", "Mean" = sapply(mysummary, "[[", 3,1), "Var" = sapply(mysummary, "[[", 3,2))
NonZeroCoefsSummary <- rbind(NonZeroCoefsSummary.Lasso, NonZeroCoefsSummary.ElNet)

MSEsummary.Lasso <- data.frame("Nobs" = Nobs, Mod = "Lasso", "Mean" = sapply(mysummary, "[[", 2,1), "Var" = sapply(mysummary, "[[", 2,2))
MSEsummary.ElNet <- data.frame("Nobs" = Nobs, Mod = "ElNet", "Mean" = sapply(mysummary, "[[", 4,1), "Var" = sapply(mysummary, "[[", 4,2))
MSEsummary <- rbind(MSEsummary.Lasso, MSEsummary.ElNet)

ggplot(NonZeroCoefsSummary) + ggtitle("Non-Zero Coefficients by Number of Observations") + 
  geom_line(aes(x = Nobs, y = Mean, color = Mod)) + 
  geom_ribbon(aes(x= Nobs, ymin = Mean - sqrt(Var), ymax = Mean + sqrt(Var), fill = Mod), alpha = 0.3)

ggplot(MSEsummary) + ggtitle("MSE by Number of Observations") + 
  geom_line(aes(x = Nobs, y = Mean, color = Mod)) + 
  geom_ribbon(aes(x= Nobs, ymin = Mean - sqrt(Var), ymax = Mean + sqrt(Var), fill = Mod), alpha = 0.3)


#Change train number


#Change validation number


#Change p


#Change sparsity levels of beta


#Change sigma

SigVals <- c(-2:2)

output <- lapply(SigVals, Simulation, times = 50, obs = 240, test.valid.size = c(20,20), p = 8, rho = 0.5, betaVals = matrix(c(3,1.5,0,0,2,0,0,0), ncol = 1))
mysummary <- lapply(output, meanAndvar)

NonZeroCoefsSummary.Lasso <- data.frame("SigVals" = SigVals, Mod = "Lasso", "Mean" = sapply(mysummary, "[[", 1,1), "Var" = sapply(mysummary, "[[", 1,2))
NonZeroCoefsSummary.ElNet <- data.frame("SigVals" = SigVals, Mod = "ElNet", "Mean" = sapply(mysummary, "[[", 3,1), "Var" = sapply(mysummary, "[[", 3,2))
NonZeroCoefsSummary <- rbind(NonZeroCoefsSummary.Lasso, NonZeroCoefsSummary.ElNet)

MSEsummary.Lasso <- data.frame("SigVals" = SigVals, Mod = "Lasso", "Mean" = sapply(mysummary, "[[", 2,1), "Var" = sapply(mysummary, "[[", 2,2))
MSEsummary.ElNet <- data.frame("SigVals" = SigVals, Mod = "ElNet", "Mean" = sapply(mysummary, "[[", 4,1), "Var" = sapply(mysummary, "[[", 4,2))
MSEsummary <- rbind(MSEsummary.Lasso, MSEsummary.ElNet)

ggplot(NonZeroCoefsSummary) + ggtitle("Non-Zero Coefficients by Standard Deviation of Error") + 
  geom_line(aes(x = SigVals, y = Mean, color = Mod)) + 
  geom_ribbon(aes(x= SigVals, ymin = Mean - sqrt(Var), ymax = Mean + sqrt(Var), fill = Mod), alpha = 0.3)

ggplot(MSEsummary) + ggtitle("MSE by Standard Deviation of Error") + 
  geom_line(aes(x = SigVals, y = Mean, color = Mod)) + 
  geom_ribbon(aes(x= SigVals, ymin = Mean - sqrt(Var), ymax = Mean + sqrt(Var), fill = Mod), alpha = 0.3)








n = 240; size = c(20,20); p = 8; rho = 0.5; errorSD = 3; betaVals = matrix(c(3,1.5,0,0,2,0,0,0), ncol = 1)

MSEValidResult <- sapply(seq(0, 1, by = 0.01), MSEValid.Lasso, x = x, y = y, sizes = sizes)

testerfunction.1 <- function(merp){
  myglmnet <- glmnet(x[1:20,], y[1:20], alpha = 1, lambda = merp)
  #intlasso <- lasso(x, y, merp)
  mylasso <- coordDescAlg.Lasso(merp, y, x)
  MSEglmnet <- mean((y[21:40] - x[21:40,] %*% myglmnet$beta)^2)
  #MSEIntLasso <- mean((y[21:40] - x[21:40,] %*% intlasso)^2)
  MSEmyLasso <- mean((y[21:40] - x[21:40,] %*% mylasso)^2)
  return(c(merp, MSEglmnet, MSEmyLasso))
}



testerfunction.2 <- function(merp){
  myglmnet <- glmnet(x[1:20,], y[1:20], alpha = .5, lambda = merp)
  myElNet <- coordDescAlg.ElNet(c(merp,0.5), y, x)
  MSEglmnet <- mean((y[21:40] - x[21:40,] %*% myglmnet$beta)^2)
  #MSEIntLasso <- mean((y[21:40] - x[21:40,] %*% intlasso)^2)
  MSEmyElNet <- mean((y[21:40] - x[21:40,] %*% myElNet)^2)
  return(c(merp, MSEglmnet, MSEmyElNet))
}


t(sapply(seq(0, 5, by = 0.05), testerfunction.2))
cbind(seq(0, 3, by = 0.01), MSEValidResult,  sapply(seq(0, 1, by = 0.01), testerfunction))



