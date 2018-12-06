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
  NumbNonZeros <- sum(Output[[1]] != 0)
  MSE <- Output[[2]]
  return(c(NumbNonZeros, MSE))
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
  MinMSE <- optimise(MSEValid, c(0,100), x=x, y=y, sizes = sizes)
  BetaHat <- coordDescAlg.Lasso(MinMSE$minimum, y[1:testSize], x[1:testSize,])
  MSE <- MinMSE$objective
  return(list(BetaHat,MSE))
}



#Called by Validation, calculates the value of MSE for a given lambda
MSEValid <- function(lambda, x, y, sizes){
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



#Divide by variance of thing

#Called by MSEValid and Validation, uses the coordinate descent algorithm to calculate beta for the a given lambda
coordDescAlg.Lasso <- function(lambda, y, x){
  
  x <- as.matrix(scale(x))
  y <- scale(y)
  converge <- FALSE
  n <- dim(x)[1]
  p <- dim(x)[2]
  BetaHat <- rep(0, p)
  BetaHat.star <- rep(0, p)
  newBetaHat <- rep(0,p)
  r <- matrix(rep(0,n*p), ncol = p)
  
  while(!converge){
    # for(j in 1:p){
    #   r[,j] <- y - x[,-j]%*%BetaHat[-j]
    #   BetaHat.star[j] <- (1/n)*(t(x[,j])%*%r[,j]) / (t((x[,j])%*%x[,j]))
    #   newBetaHat[j] <- sign(BetaHat.star[j])*pmax((abs(BetaHat.star[j]) - lambda / (t((x[,j])%*%x[,j]))), 0)
    #   
    # }
    for(j in 1:p){
      r[,j] <- y - x[,-j]%*%BetaHat[-j]
    }
    BetaHat.star <- (1/n)*colSums(x*r)
    
    newBetaHat <- sign(BetaHat.star)*pmax((abs(BetaHat.star) - lambda), 0)
    

      if(all(round(BetaHat, 5) == round(newBetaHat,5))){
        converge <- TRUE
        break
      } else {
        BetaHat <- newBetaHat
      }
      
    
  }
  
  return(BetaHat)
  
}


set.seed(3)
placeholder <- 0

for(i in 1:100){
output <- mySample(240, 8, 0.5, 3, betaVals = matrix(c(3,1.5,0,0,2,0,0,0)))
coordDescAlg.Lasso(0.1, output[[2]][1:20], output[[1]][1:20,])
}
placeholder <- placeholder + 1



ErrorOutput <- output
x <- ErrorOutput[[1]][1:20,]
y <- ErrorOutput[[2]][1:20]
lambda <- 0






set.seed(3)

Simulation(times = 50, obs = 240, test.valid.size = c(20,20), p = 8, rho = 0.5, errorSD = 3, betaVals = matrix(c(3,1.5,0,0,2,0,0,0), ncol = 1))


#Don't need to standardize the data





coordDescAlg.ElasticNet <- function(lambda, y, x){
  
  x <- as.matrix(x)
  converge <- FALSE
  n <- dim(x)[1]
  p <- dim(x)[2]
  BetaHat <- rep(0, p)
  BetaHat.star <- rep(0, p)
  r <- matrix(rep(0,n*p), ncol = p)
  
  while(!converge){
    for(j in 1:p){
      r[,j] <- y - x[,-j]%*%BetaHat[-j]
    }
    BetaHat.star <- (1/n)*colSums(x*r)
    
    newBetaHat <- sign(BetaHat.star)*pmax((abs(BetaHat.star) - lambda), 0)
    
    if(all(round(BetaHat, 3) == round(newBetaHat,3))){
      converge <- TRUE
      break
    } else {
      BetaHat <- newBetaHat
    }
    
    
  }
  
  return(BetaHat)
  
}
