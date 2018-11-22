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
      r[,j] <- y - rowSums(x[,-j]*beta.hat[-j])
    }
    beta.hat.star <- (1/n) * colSums(x*r)
    
    new.beta.hat <- sign(beta.hat.star)*pmax(abs(beta.hat.star)-lamda)
    
    if(all(beta.hat == new.beta.hat)){
      converge=TRUE
    }
    else{beta.hat=new.beta.hat}
  }
  return(beta.hat)
}
CA.Lasso(y,mysample,1)


CA.elasticnet <- function(x,y,lamda1,lamda2){
  
  converge <- FALSE
  n <- dim(x)[1]
  p <- dim(x)[2]
  beta.hat <- rep(0,p)
  r <- matrix(rep(NA,n*p),ncol=p)
  
  while(!converge){
    for (j in 1:p){
      r[,j] = y - rowSums(x[,-j]*beta.hat[-j])
      
      
      beta.hat.star <- (1/n) * colSums(x*r)
      
      new.beta.hat <- sign(beta.hat.star)*pmax(abs(beta.hat.star)-lamda1)*((beta.hat.star)**2 - lamda2)
      
      if(all(beta.hat == new.beta.hat)){
        converge=TRUE
      }
      else{beta.hat=new.beta.hat}
    }
    return(beta.hat)
  
  
}

  