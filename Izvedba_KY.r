## okolje kjer so funkcije
## setwd()
source("KonnoYamazaki_function.r")
source("KYdonosi_function.r")       
source("KYtransactionCosts_function.r")
source("KYturnoverConstraint_function.r")

library(quantmod)
library(lpSolve)
library(lpSolveAPI)
library(tseries)

data <- new.env()
##IZBOR DELNIC######################################################################
####################################################################################
delnice <-'ORCL VRSN FOXA SNDK NKE MRK KSS CVX AMZN MSFT NOK TSM STJ HOG TEVA GE NFLX BA MON GS AIG ADBE CHL RDS-A CSCO AMX BUD WMT PTR GT MTU TTM CSIQ BX BRK-A AAPL PFE EBAY F FCAU'
stocks <- strsplit(delnice, ' ')[[1]]
getSymbols(stocks, src = 'yahoo', from = '2011-01-02', to = '2015-01-08', env = data, auto.assign = T)
priceAdj <- do.call(merge, (eapply(data, Ad)))
## donosi
donosi <- apply(priceAdj, 2, diff)/priceAdj[-nrow(priceAdj),]  
##PRICAKOVANI DONOSI, KOVARIANCNA MATRIKA###########################################
####################################################################################
## pricakovani donos 
Edonosi <- apply(donosi,2, mean)
## korelacijska matrika
CovMat <- cov(donosi)




##KONNO-YAMAZAKI (BREZ OMEJITEV KRATKE PRODAJE IN DELEZA POSAMEZNE DELNICE)#########
####################################################################################
BrezOmejitev <- KonnoYamazaki(donosi, upperBound = 10^4, lowerBound = -10^4, od=5*10^-4, do=0.0012, po =0.5*10^-4)


##KONNO-YAMAZAKI (Z OMEJITVIJO KRATKE PRODAJE IN MAKSIMALNEGA DELEZA ENE DELNICE)###
####################################################################################
omejen1 <- KonnoYamazaki(donosi, upperBound = 0.1, lowerBound = 0, od=5*10^-4, do=0.0012, po =0.5*10^-4 )
omejen2 <- KonnoYamazaki(donosi, upperBound = 10^4, lowerBound = 0, od=5*10^-4, do=0.0012, po =0.5*10^-4 )
omejen3 <- KonnoYamazaki(donosi, upperBound = 0.1, lowerBound = -10^4, od=5*10^-4, do=0.0012, po =0.5*10^-4 )

## primerjava L1 tvegana t volatilnostjo
plot(BrezOmejitev$sigma, BrezOmejitev$L1, type = 'o', xlab = 'volatilnost', ylab = 'L1-tveganje', col = 'blue')
plot(omejen1$sigma, omejen1$L1, type = 'o', xlab = 'volatilnost', ylab = 'L1-tveganje', col = 'blue')
plot(omejen2$sigma, omejen2$L1, type = 'o', xlab = 'volatilnost', ylab = 'L1-tveganje', col = 'blue')
plot(omejen3$sigma, omejen3$L1, type = 'o', xlab = 'volatilnost', ylab = 'L1-tveganje', col = 'blue')

## primerjava donosi(sigma) portfeljev
mejaGrafaU <- max(c(BrezOmejitev$sigma, omejen1$sigma))
mejaGrafaD <- min(c(BrezOmejitev$sigma, omejen1$sigma))
plot(BrezOmejitev$sigma, BrezOmejitev$donos, type = 'o', xlim = c(mejaGrafaD, mejaGrafaU), col = 'red',
		main = '', xlab = 'sigma', ylab = 'donos')
lines(omejen1$sigma[omejen1$sigma != 0], omejen1$donos[omejen1$sigma != 0], type = 'o', col = 'blue')
lines(omejen2$sigma[omejen2$sigma != 0], omejen2$donos[omejen2$sigma != 0], type = 'o', col = 'green')
lines(omejen3$sigma[omejen3$sigma != 0], omejen3$donos[omejen3$sigma != 0], type = 'o', col = 'black')
legend('bottomright', c('Brez omejitev', 'Z omejitvam 1', 'Z omejitvam 2', 'Z omejitvam 3'), lwd = c(2,2,2,2),
		col = c('red', 'blue', 'green', 'black'), bty="n")

		
		
		
##TESTIRANJE NA t+1 (dela dolgo cas)################################################
####################################################################################
# izbira t
n <- 300
t <- nrow(donosi)-n
## portfelj sestavljen s Konno-Yamazaki (brez omejitev) 
realizacija1 <- NULL
## portfelj sestavljen s Konno-Yamazaki (omejitev deleza na 0.1, brez short)
realizacija2 <- NULL
## portfelj sestavljen s portfioli.optim {iz paketa tseries} (brez omejitev)
realizacijaPO1 <- NULL
## portfelj sestavljen s portfioli.optim {iz paketa tseries}
realizacijaPO2 <- NULL
for(i in 1:(n-1)){
  t <- t+1
  if(i%%10 == 0 ){print(i)}
  kos <- donosi[1:t,]
  test <- donosi[t+1,]
  od=5*10^-4
  do=0.0012
  po =0.5*10^-4 
  zaporedjeN <- seq(od, do , po)
  potek1 <- KonnoYamazaki(kos, 10^4, -10^4, od, do, po)
  #potek2 <- KonnoYamazaki(kos, 0.3, -10^4, od, do, po)
  potek2 <- KonnoYamazaki(kos, 0.1, 0, od, do, po)
  PO1 <- NULL
  PO2 <- NULL
  for (j in zaporedjeN){
    p1 <- portfolio.optim(kos, pm = j, riskless = TRUE, skorts = TRUE, covmat = CovMat)
    PO1 <- cbind(PO1, p1$pw)
	#p2 <- portfolio.optim(kos, pm = j, reshigh = rep(0.3, ncol(kos)), shorts = FALSE)
	p2 <- portfolio.optim(kos, pm = j, reshigh = rep(0.3, ncol(kos)), shorts = FALSE, covmat = CovMat)
    PO2 <- cbind(PO2, p2$pw)
  }
  realizacija1 <- cbind(realizacija1,t(potek1$delezi)%*%t(test))
  realizacija2 <- cbind(realizacija2,t(potek2$delezi)%*%t(test))
  realizacijaPO1 <- cbind(realizacijaPO1, t(PO1)%*%t(test))
  realizacijaPO2 <- cbind(realizacijaPO2, t(PO2)%*%t(test))
}


donosi1 <- rowMeans(realizacija1)
donosi2 <- rowMeans(realizacija2)
donosiPO1 <- rowMeans(realizacijaPO1)
donosiPO2 <- rowMeans(realizacijaPO2)


## grafi donosov
slikaOd <- min(c(zaporedjeN, donosi1, donosi2, donosiPO1, donosiPO2))
slikaDo <- max(c(zaporedjeN, donosi1, donosi2, donosiPO1, donosiPO2))
plot(zaporedjeN, ylim =c(slikaOd-0.5*10^-4, slikaDo))
legend('bottomleft', c('predviden donos','ne omejen KY','omejen KY',
					'ne omejen p.opt' ,'omejen p.opt')
					,pch = c(1,17,17,19,19), col = c(1,2,3,2,3), bty = 'n')
# legend(x = 3.2, y = slikaOd+ slikaOd*0.15 , c('predviden donos'),
		# pch = c(1), col = c(1), bty="n")
# legend(x = 6.1, y = slikaOd+ slikaOd*0.3, c('ne omejen KY','omejen KY')
					# ,pch = c(17,17), col = c(2,3), bty="n")
# legend(x = 8.5, y = slikaOd+ slikaOd*0.3, c('ne omejen p.opt' ,'omejen p.opt')
					# ,pch = c(19,19), col = c(2,3), bty="n")
points(donosi1, col = 2, pch = 17)
points(donosi2, col = 3, pch = 17)
points(donosiPO1, col =2, pch = 19)
points(donosiPO2, col =3, pch = 19)


####### preostali grafi
var1 <- apply(realizacija1, 1, var)
var2 <- apply(realizacija2, 1, var)
varPO1 <- apply(realizacijaPO1, 1, var)
varPO2 <- apply(realizacijaPO2, 1, var)

## grafi donosnosti v odnvisnosti od  variance
plot(sqrt(var1), donosi1, type = 'o')
plot(sqrt(var2), donosi2, type = 'o')
plot(sqrt(varPO1), donosiPO1, type = 'o')
plot(sqrt(varPO2), donosiPO2, type = 'o')

## grafi  variance
yos <- c(min(sqrt(var1), sqrt(var2)), max(sqrt(var1),sqrt(var2)))
plot(sqrt(var1), type = 'o', xlab = 'zaporedni portfelj', ylab = 'volatilnost', col = 'blue',
	ylim = yos)
lines(sqrt(var2), type = 'o', col = 2)
legend('topleft', c('neomejen KY','omejen KY'),lwd = c(2,2), col = c('blue', 2), bty="n")


plot(sqrt(varPO1), type = 'o')
plot(sqrt(varPO2), type = 'o')





##OBCUTLJIVOST (dela dolgo casa)####################################################
####################################################################################
donosi <- apply(priceAdj, 2, diff)/priceAdj[-nrow(priceAdj),]  
Edonosi <- apply(donosi,2, mean)
## variance normalne porazdelitve
disturb <- c(2, 1, 0.7, 0.3, 0.15, 0.05, 0.0025, 0.001)
## ponovitev z vsako porazdelitvijo N(0, disturb)
rept = 1:100
## list spremembe portfelja za vsako zeljeno donosnost
od=5*10^-4
do=0.0012
po =0.5*10^-4
zDonosi <- seq(od, do , po)
sprememba <- list()
rez <- data.frame(matrix(NA, length(rept), length(disturb)))
for (k in zDonosi){
	sprememba[[paste0('zDonosi ',k)]] <- rez
}
## matrika v katero se shranjujejo vsota abs motenj
motnja <- data.frame(matrix(NA, length(rept), length(disturb)))
names(rez) <- disturb
## KY prtfelj s ocenjenimi parametri iz preteklih cen
original <- KYdonosi(donosi, Edonosi, upperBound = 10^4, lowerBound = -10^4, od=5*10^-4, do=0.0012, po =0.5*10^-4)
i <- 0
for (l in disturb){
	i <- i+1
	print(i)
	for (k in rept){
		print(k)
		Nak1 <- rnorm(length(Edonosi), mean = 0, sd = l*mean(abs(Edonosi)))
		EdonosiR1 <- Edonosi - Nak1
		motnja[k, i] <- mean(abs(Nak1))
		rand1 <- KYdonosi(donosi, EdonosiR1, upperBound = 10^4, lowerBound = -10^4, od=5*10^-4, do=0.0012, po =0.5*10^-4)
		for(j in 1:length(zDonosi)){
			sprememba[[j]][k, i] <- mean(abs(original$delezi-rand1$delezi)[ , j])
		}
	}
}


## graf spremembe portfelja v odvisnosti od motnje za razlicne pricakovane donose.
ymax <- max(colMeans(sprememba[[length(zDonosi)]]))
plot(colMeans(motnja), colMeans(sprememba[[1]]),
	type = 'o', col = 'blue', xlab = 'motnja', ylab = 'sprememba porfelja', ylim = c(0,ymax))
## barve na grafu
barve <- 1
## zeljeni donosi portfeljev
portfeljiI <- NULL
for (i in seq(1,length(zDonosi),3)){
	barve <- c(barve, tail(barve,1) + 1)
	lines(colMeans(motnja), colMeans(sprememba[[i]]),type = 'o', col = tail(barve,1))
	portfeljiI <- c(portfeljiI, paste0('donosi ', zDonosi[i]))
}
legend('topleft', portfeljiI,, lwd = rep(2,9), col = barve, bty="n")

	

##TURNOVER CONSTRAINTS##############################################################
####################################################################################

zacPort <- c(rep(0.2, 5), rep(0,ncol(donosi)-5))
ky <- KonnoYamazaki(donosi, upperBound = 10^4, lowerBound = -10^4, od=5*10^-4, do=0.0012, po =0.5*10^-4)
tc1 <- KYturnoverConstraint(donosi, zacPort, tuCon = 3, upperBound = 10^4, lowerBound = -10^4, od=5*10^-4, do=0.0012, po =0.5*10^-4)
tc2 <- KYturnoverConstraint(donosi, zacPort, tuCon = 2, upperBound = 10^4, lowerBound = -10^4, od=5*10^-4, do=0.0012, po =0.5*10^-4)
tc3 <- KYturnoverConstraint(donosi, zacPort, tuCon = 1.5, upperBound = 10^4, lowerBound = -10^4, od=5*10^-4, do=0.0012, po =0.5*10^-4)
tc4 <- KYturnoverConstraint(donosi, zacPort, tuCon = 1.3, upperBound = 10^4, lowerBound = -10^4,od=5*10^-4, do=0.0012, po =0.5*10^-4)
tc5 <- KYturnoverConstraint(donosi, zacPort, tuCon = 1, upperBound = 10^4, lowerBound = -10^4, od=5*10^-4, do=0.0012, po =0.5*10^-4)
tc6 <- KYturnoverConstraint(donosi, zacPort, tuCon = 0.7, upperBound = 10^4, lowerBound = -10^4, od=5*10^-4, do=0.0012, po =0.5*10^-4)


xos <- c(min(ky$sigma, tc6$sigma), max(ky$sigma, tc6$sigma))
yos <- c(min(ky$donos, tc6$donos), max(ky$donos, tc6$donos))
plot(ky$sigma, ky$donos, type = 'o', col = 'black', xlim = xos, ylim = yos, xlab = 'sigma', ylab = 'donos')
lines(tc1$sigma, tc1$donos, type = 'o', col = 'red')
lines(tc2$sigma, tc2$donos, type = 'o', col = 3)
lines(tc3$sigma[tc3$sigma != 0], tc3$donos[tc3$sigma != 0], type = 'o', col = 4)
lines(tc4$sigma[tc4$sigma != 0], tc4$donos[tc4$sigma != 0], type = 'o', col = 5)
lines(tc5$sigma[tc5$sigma != 0], tc5$donos[tc5$sigma != 0], type = 'o', col = 6)
lines(tc6$sigma[tc6$sigma != 0], tc6$donos[tc6$sigma != 0], type = 'o', col = 8)
legend('bottomright', c('h = 3', 'h = 2', 'h = 1,5', 'h = 1,3', 'h = 1', 'h = 0,7'),
		col = c("red","green3", "blue","cyan","magenta","gray" ), lwd = rep(2,6), bty="n")

##TRANSACTION COSTS#################################################################
####################################################################################


zacPort <- c(rep(0.2, 5), rep(0, ncol(donosi)-5))
costs1 <- KYtransactionCosts(donosi, zacPort, sellCosts = 0, buyCosts = 0, upperBound = 10^4, lowerBound = -10^4, od=5*10^-4, do=0.0012, po =0.5*10^-4)
ky <- KonnoYamazaki(donosi, upperBound = 10^4, lowerBound = -10^4, od=5*10^-4, do=0.0012, po =0.5*10^-4)

plot(costs1$sigma, costs1$donos, type = 'o', col = 'red')
legend('topleft', c('KYtranCost 0', 'KY'), lty = c(1,3), col = c('red', 'blue'))
lines(ky$sigma, ky$donos, lty = 3, col = 'blue', lwd = 2)

costs2 <- KYtransactionCosts(donosi, zacPort, sellCosts = 5*10^-5, buyCosts = 5*10^-5, upperBound = 10^4, lowerBound = -10^4, od=5*10^-4, do=0.0012, po =0.5*10^-4)
costs3 <- KYtransactionCosts(donosi, zacPort, sellCosts = 10^-4, buyCosts = 10^-4, upperBound = 10^4, lowerBound = -10^4, od=5*10^-4, do=0.0012, po =0.5*10^-4)
costs4 <- KYtransactionCosts(donosi, zacPort, sellCosts = 2*10^-4, buyCosts = 2*10^-4, upperBound = 10^4, lowerBound = -10^4, od=5*10^-4, do=0.0012, po =0.5*10^-4)
costs5 <- KYtransactionCosts(donosi, zacPort, sellCosts = 3*10^-4, buyCosts = 3*10^-4, upperBound = 10^4, lowerBound = -10^4, od=5*10^-4, do=0.0012, po =0.5*10^-4)

xos <- c(min(ky$sigma, costs5$sigma), max(ky$sigma, costs5$sigma))
yos <- c(min(ky$donos, costs5$donos), max(ky$donos, costs5$donos))

plot(costs1$sigma, costs1$donos, type = 'o', col = 'black', xlim = xos, ylim = yos)
lines(costs2$sigma, costs2$donos, type = 'o', col = 'red')
lines(costs3$sigma, costs3$donos, type = 'o', col = 3)
lines(costs4$sigma, costs4$donos, type = 'o', col = 4)
lines(costs5$sigma, costs5$donos, type = 'o', col = 5)


