library(magic)
library(quantmod)
library(lpSolve)
library(lpSolveAPI)
library(tseries)

KYturnoverConstraint <- function(donosi,zacPort,tuCon = 10, upperBound = NULL, lowerBound = NULL,od=5*10^-4, do=0.002, po =1*10^-4 ){
  ## izracuna portfelj ce imamo turnover constraints tuCon in zacetni portfelj zacPort
  Edonosi <- apply(donosi,2, mean)
  CovMat <- cov(donosi)
  ## omejitve A
  dimObdobji <- nrow(donosi)
  dimPortfelj <- ncol(donosi)
  ## omejitve y,z,x
  omejitve <- matrix(0, dimObdobji+2+2*dimPortfelj+1+1,2*dimObdobji+ dimPortfelj+2*dimPortfelj)
  omejitveMAD <- diag(2*nrow(donosi))
  omejitveMAD[row(omejitveMAD) - col(omejitveMAD) == -1] = -1
  omejitveMAD <- omejitveMAD[seq(1,ncol(omejitveMAD),by = 2),]
  #B <- matrix(c(1,-1),1,2)
  #omejitve[1:dimObdobji,1:(2*dimObdobji)] <- do.call(adiag, replicate(dimObdobji, B, simplify = FALSE))
  omejitve[1:dimObdobji,1:(2*dimObdobji)] <- omejitveMAD
  omejitve[1:dimObdobji,(2*dimObdobji+1):(2*dimObdobji+ dimPortfelj)] <- -t(t(donosi) - Edonosi) 
  
  #omejitve Edx < R in xtx =1
  omejitve[(dimObdobji+1):(dimObdobji+2),(2*dimObdobji+1):(2*dimObdobji+ dimPortfelj)] <- rbind(rep(1,dimPortfelj), Edonosi)
  
  ## TURNOVER CONSTRAINTS
  
  omejitve[(dimObdobji+3):(dimObdobji+2+2*dimPortfelj),(2*dimObdobji+1):(2*dimObdobji+ dimPortfelj)] <- rbind(diag(dimPortfelj),-diag(dimPortfelj))
  omejitve[(dimObdobji+3):(dimObdobji+2+2*dimPortfelj),(2*dimObdobji+ dimPortfelj+1):(2*dimObdobji+ dimPortfelj+2*dimPortfelj)] <- -diag(2*dimPortfelj)
  omejitve[(dimObdobji+2+2*dimPortfelj+1),(2*dimObdobji+ dimPortfelj+1):(2*dimObdobji+ dimPortfelj+2*dimPortfelj)] <- 1
  omejitve[(dimObdobji+2+2*dimPortfelj+2),(2*dimObdobji+ dimPortfelj+1):(2*dimObdobji+ dimPortfelj+2*dimPortfelj)] <- c(rep(1, dimPortfelj), rep(-1, dimPortfelj))
  ## predznaki
  signs = c(rep('=', dimObdobji+1),'>=',rep('<=', dimPortfelj),rep('<=', dimPortfelj), '<=', '=')
  
  
  
  
  ## predznaki
  
  x <- NULL
  ret1 <- NULL
  sig1 <- NULL
    
  for(i in seq(od,do, by = po)){
    #i = 0.0005
    #for(i in seq(0.0005,0.002, by = 0.00005)){
    portfelj <- (2*nrow(donosi)+1):(2*nrow(donosi)+ncol(donosi))
    linProg <- make.lp(0, ncol(omejitve))
    set.objfn(linProg, c(rep(1, 2*nrow(donosi)), rep(0,ncol(donosi)+2*dimPortfelj)))
    set.bounds(linProg, upper = rep(upperBound, length(portfelj)), lower = rep(lowerBound, length(portfelj)), columns = portfelj)
    minR <- i
    #minR <- 0.0005
    b0 <- c(matrix(0, dimObdobji,1), c(1,minR), zacPort, -zacPort, tuCon,0)
    
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
