library(matrixcalc)



set.seed(1234)


x <- data.frame("Int" = rep(1,100), "Var1" = rnorm(100,10,5), "Var2" = rnorm(100,40,1), "Var3" = rnorm(100,70,20))
beta <- matrix(c(4,5,-16,0.5), ncol = 1)

y <- as.matrix(x) %*% beta







coordDescAlg.Lasso <- function(y, x, lambda){
  
  converge <- FALSE
  n <- dim(x)[1]
  p <- dim(x)[2]
  BetaHat <- rep(0, p)
  BetaHat.star <- rep(0, p)
  r <- matrix(rep(NA,n*p), ncol = p)
  
  while(!converge){
    for(j in 1:p){
      r[,j] <- y - rowSums(x[,-j]*BetaHat[-j])
    }
    
    BetaHat.star <- (1/n)*colSums(hadamard.prod(x,r))
    
    newBetaHat <- sign(BetaHat.star)*pmax((abs(BetaHat.star) - lambda), 0)
    
    if(all(BetaHat == newBetaHat)){
      converge <- TRUE
    }
    else{
      BetaHat <- newBetaHat
    }
  }
  
  return(BetaHat)
  
}















coordDescAlg.ElasNet <- function(y, x, lambda){
  
  converge <- FALSE
  n <- dim(x)[1]
  p <- dim(x)[2]
  BetaHat <- rep(0, p)
  BetaHat.star <- rep(0, p)
  r <- matrix(rep(NA,n*p), ncol = p)
  
  while(!converge){
    for(j in 1:p){
      r[,j] <- y - rowSums(x[,-j]*BetaHat[-j])
    }
    
    BetaHat.star <- (1/n)*colSums(hadamard.prod(x,r))
    
    newBetaHat <- sign(BetaHat.star)*pmax((abs(BetaHat.star) - lambda), 0)
    
    if(all(BetaHat == newBetaHat)){
      converge <- TRUE
    }
    else{
      BetaHat <- newBetaHat
    }
  }
  
}