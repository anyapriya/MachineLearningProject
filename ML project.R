#simulation
beta <- matrix(c(3,1.5,0,0,2,0,0,0),ncol=1)
sigma <- 3
error <- rnorm(240)
sigmaMatrix <- matrix((rep(0,8*8)),ncol=8)
for(i in 1:8){
  for(j in 1:8){
    sigmaMatrix[i,j] <- 0.5^abs(i-j)
  }
}

library(MASS)

mysample <-mvrnorm(n = 240, mu = rep(0,8), Sigma = sigmaMatrix)
y <- mysample%*%beta + sigma*error

##


CA.Lasso <- function(y,x,lamda){
  
  converge <- FALSE
  n <- dim(x)[1]
  p <- dim(x)[2]
  beta.hat <- rep(0,p)
  r <- matrix(rep(NA,n*p),ncol=p)

  while(!converge){
   for (j in 1:p){
    r[,j] <- y - x[,-j]%*%beta.hat[-j]}
    
    beta.hat.star <- (1/n) * colSums(x*r)
    
    new.beta.hat <- sign(beta.hat.star)*pmax(abs(beta.hat.star)-lamda,0)
    
    if(all(round(beta.hat,6) == round(new.beta.hat,6))){
      converge <- TRUE
    }
    else{(beta.hat=new.beta.hat)
  }
    return(beta.hat)
}
}
CA.Lasso(y,mysample,1)

##########ELASTIC NET

CA.Elasticnet <- function(y,x,lamda1,lamda2){
  
  converge <- FALSE
  n <- dim(x)[1]
  p <- dim(x)[2]
  beta.hat <- rep(0,p)
  r <- matrix(rep(NA,n*p),ncol=p)
  
  while(!converge){
    for (j in 1:p){
      r[,j] <- y - x[,-j]%*%beta.hat[-j]}
    
    beta.hat.star <- (1/n) * colSums(x*r)
    
    new.beta.hat <- sign(beta.hat.star)*pmax(abs(beta.hat.star)-lamda1,0)*((beta.hat.star)**2 - lamda2)
    
    if(all(round(beta.hat,6) == round(new.beta.hat,6))){
      converge <- TRUE
    }
    else{(beta.hat=new.beta.hat)
    }
    return(beta.hat)
  }
}
CA.Elasticnet(y,mysample,0.5,0.5)









#lamda = seq(0,100)

#training = sample(lamda,size = 50,replace = F)
#testing = lamda[!training]

#first we fit the training data to our model to make predictions about testing error and then we use that to calculate the MSE 
#We want sth with the lowest MSE
#What model to fit it on? Linear??




  