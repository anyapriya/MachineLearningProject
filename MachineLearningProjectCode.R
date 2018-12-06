library(matrixcalc)
library(MASS)






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
  BetaHat.Lasso <- coordDescAlg.Lasso(MinMSE.Lasso$minimum, y[1:testSize], x[1:testSize,])
  
  MinMSE.ElasticNet <- optim(c(0.0001, 0.0001), MSEValid.ElasticNet,  x = x, y = y, sizes = sizes)
  BetaHat.ElasticNet <- coordDescAlg.ElasticNet(MinMSE.ElasticNet$par, y[1:testSize], x[1:testSize,])
  
  MSE <- MinMSE$objective
  return(list(BetaHat.Lasso,MinMSE.Lasso$objective, BetaHat.ElasticNet, MinMSE.ElasticNet$value))
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
  BetaHat <- coordDescAlg.ElasticNet(lambda, ytrain, xtrain)
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
    
      BetaHat.star[j] <- (1/n)*(t(x[,j])%*% r[,j]) / (t(x[,j])%*% x[,j])
      newBetaHat[j] <- sign(BetaHat.star[j])*pmax(BetaHat.star[j] - (lambda / (t(x[,j])%*% x[,j])), 0)
    
      if(all(round(BetaHat, 3) == round(newBetaHat,3))){
        converge <- TRUE
        break
      } else {
        BetaHat <- newBetaHat
      }

      
    }
    

  }
  
  return(BetaHat)
  
}



























coordDescAlg.ElasticNet <- function(lambda, y, x){
  
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
      
      BetaHat.star[j] <- (1/n)*(t(x[,j])%*% r[,j]) / (t(x[,j])%*% x[,j])
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



set.seed(3)

Simulation(times = 3, obs = 240, test.valid.size = c(20,20), p = 8, rho = 0.5, errorSD = 3, betaVals = matrix(c(3,1.5,0,0,2,0,0,0), ncol = 1))


