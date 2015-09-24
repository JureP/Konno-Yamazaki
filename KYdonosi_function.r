library(quantmod)
library(lpSolve)
library(lpSolveAPI)
library(tseries)

KYdonosi <- function(	donosi, ## matrika donosiv
						Edonosi, ## pricakovni donosi
						upperBound = NULL, ## zgornja omejitev deleza 
						lowerBound = NULL, ## spodnja omejitev deleza
						od=5*10^-4, ## od pricakovane donosnosti 
						do=0.002, ## do pricakovane donose
						po =1*10^-4 ## premik
						){
  ## KYdonosi izracuna KY portfelj pri danih donosih (za potrebe sensitivity analysis)
  CovMat <- cov(donosi)
  ## omejitve A
  ##1 
  omejitveMAD <- diag(2*nrow(donosi))
  omejitveMAD[row(omejitveMAD) - col(omejitveMAD) == -1] = -1
  omejitveMAD <- omejitveMAD[seq(1,ncol(omejitveMAD),by = 2),]
  ##2
  omejitveX <- -t(t(donosi) - Edonosi)
  omejitev1 <- cbind(omejitveMAD, omejitveX)
  ## omejitve donosov 
  vsota1 <- c(rep(0,ncol(omejitveMAD)), rep(1,ncol(donosi)))
  neenakost <- c(rep(0,ncol(omejitveMAD)), Edonosi)
  ##skupaj
  omejitve <- rbind(omejitev1, vsota1, neenakost)
  ## predznaki
  
  signs = c(rep('=', nrow(omejitve)-1),'>=')
  x <- NULL
  ret1 <- NULL
  sig1 <- NULL
  for(i in seq(od,do, by = po)){
    #for(i in seq(0.0005,0.002, by = 0.00005)){
    portfelj <- (2*nrow(donosi)+1):(2*nrow(donosi)+ncol(donosi))
    linProg <- make.lp(0, 2*nrow(donosi)+ncol(donosi))
    set.objfn(linProg, c(rep(1, 2*nrow(donosi)), rep(0,ncol(donosi))))
    set.bounds(linProg, upper = rep(upperBound, length(portfelj)), lower = rep(lowerBound, length(portfelj)), columns = portfelj)
    minR <- i
    #minR <- 0.0005
    b0 <- c(matrix(0, nrow(omejitveMAD),1), c(1,minR))
    
    for(j in 1:nrow(omejitve)){
      add.constraint(linProg, omejitve[j,], signs[j], b0[j])
    }
    
    solve(linProg)
    y <- get.variables(linProg)
    x <- cbind(x,y[portfelj])
    ret1 <- cbind(ret1,y[portfelj]%*% Edonosi)
    sig1 <- cbind(sig1,y[portfelj]%*%CovMat%*%y[portfelj])
  }
  #plot(sig1, ret1, type = 'o')
  rezultat <- list(delezi = x, donos = ret1,sigma =  sig1)
  return(rezultat)
}
