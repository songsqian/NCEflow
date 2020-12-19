## Hierarchical Polr model for converting RichTol to Bioclass

## initializing

packages<-function(x, repos="http://cran.r-project.org", ...){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x, repos=repos,
                     dependencies=T, ...)
    require(x,character.only=TRUE)
  }
}

packages(arm)
packages(dplyr)
packages(lattice)
packages(rv)
packages(tikzDevice)
packages(rpart)
packages(rstan)
packages(MASS)

rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

nchains <-  min(c(parallel::detectCores(), 4))
niters <- 10000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

base <- getwd()
dataDIR <- paste(base, "Data", sep="/")
## Data
invCond <- read.csv(paste(dataDIR, "Inv_Condition.csv", sep="/"))
names(invCond)
plot(SEQ~RichTOL, data=invCond)

## polr model
POLR1 <- polr(ordered(SEQ) ~ RichTOL, data=invCond)
beta <- coef(POLR1)
kappa <- POLR1$zeta

c1.5 <- kappa[1]/beta
c2.5 <- kappa[2]/beta
c3.5 <- kappa[3]/beta
c4.5 <- kappa[4]/beta
sigma <- 1/beta

se.c <- summary(POLR1)$coef[2:5,2]/beta
IRichTOL <- seq(range(invCond$RichTOL)[1],range(invCond$RichTOL)[2],,50)
pE <- invlogit((c1.5 - IRichTOL)/sigma)
pG <- invlogit((c2.5 - IRichTOL)/sigma) - invlogit((c1.5 - IRichTOL)/sigma)
pGF <- invlogit((c3.5 - IRichTOL)/sigma) - invlogit((c2.5 - IRichTOL)/sigma)
pF <- invlogit((c4.5 - IRichTOL)/sigma) - invlogit((c3.5 - IRichTOL)/sigma)
pP <- 1.0 - invlogit((c4.5 - IRichTOL)/sigma)

## plot showing uncertainty in the cut-off points
pdf(file="polr0.pdf", height=3, width=3)
plot(range(invCond$RichTOL), c(0,1), xlab="RichTOL", ylab="Prob",
     type="n")
polygon(x=c(c1.5-se.c[1], c1.5+se.c[1], c1.5+se.c[1],c1.5-se.c[1]),
        y=c(0,0,1,1), col="gray", density=-1, border=NA)
polygon(x=c(c2.5-se.c[2], c2.5+se.c[2], c2.5+se.c[2],c2.5-se.c[2]),
        y=c(0,0,1,1), col="gray", density=-1, border=NA)
polygon(x=c(c3.5-se.c[3], c3.5+se.c[3], c3.5+se.c[3],c3.5-se.c[3]),
        y=c(0,0,1,1), col="gray", density=-1, border=NA)
polygon(x=c(c4.5-se.c[4], c4.5+se.c[4], c4.5+se.c[4],c4.5-se.c[4]),
        y=c(0,0,1,1), col="gray", density=-1, border=NA)
lines(IRichTOL, pE)
lines(IRichTOL, pG, lty=2)
lines(IRichTOL, pGF, lty=3)
lines(IRichTOL, pP, lty=4)
dev.off()

expected <- function(x, c1.5, c2.5, c3.5, c4.5, sigma=1/beta){
	p1.5 <- invlogit((x-c1.5)/sigma)
	p2.5 <- invlogit((x-c2.5)/sigma)
        p3.5 <- invlogit((x-c3.5)/sigma)
        p4.5 <- invlogit((x-c4.5)/sigma)
	return((1*(1-p1.5)+2*(p1.5-p2.5)+3*(p2.5-p3.5)+
                4*(p3.5-p4.5)+5*p4.5))
}

pdf(file="polr1.pdf", height=3, width=3)
par(mar=c(3,5,0.25,0.25), mgp=c(1.5,0.25,0), tck=-0.005)
plot(0, 0, xlim=range(invCond$RichTOL), ylim=c(1,5), xlab="RichTOL",
     ylab="", type="n", axes=F)
axis(1)
axis(2, at=1:5, labels=c("Excellent","Good","Good-Fair","Fair","Poor"), las=1)
lines(rep(c1.5, 2), c(1,2))
lines(rep(c2.5, 2), c(2,3))
lines(rep(c3.5, 2), c(3,4))
lines(rep(c4.5, 2), c(4,5))
curve(expected(x, c1.5, c2.5, c3.5, c4.5), add=TRUE)
points(invCond$RichTOL,
       jitter(invCond$SEQ, 0.25), pch=1,
            col="gray", cex=0.5)
dev.off()

## polr with interaction
invCond$Upland <- as.numeric(invCond$ERIII < 50)
POLR2 <- polr(ordered(SEQ) ~ RichTOL*Upland, data=invCond)
summary(POLR2)
beta2 <- coef(POLR2)
kappa2 <- POLR2$zeta

## lowland (`Upland=0`)
betaL <- beta2[1]
kappa2 <- POLR2$zeta

c1.5L <- kappa2[1]/betaL
c2.5L <- kappa2[2]/betaL
c3.5L <- kappa2[3]/betaL
c4.5L <- kappa2[4]/betaL
sigmaL <- 1/betaL

temp=invCond$Upland==0
pdf(file="polr1L.pdf", height=3, width=3)
par(mar=c(3,5,0.25,0.25), mgp=c(1.5,0.25,0), tck=-0.005)
plot(0, 0, xlim=range(invCond$RichTOL), ylim=c(1,5), xlab="RichTOL",
     ylab="", type="n", axes=F)
axis(1)
axis(2, at=1:5, labels=c("Excellent","Good","Good-Fair","Fair","Poor"), las=1)
lines(rep(c1.5L, 2), c(1,2))
lines(rep(c2.5L, 2), c(2,3))
lines(rep(c3.5L, 2), c(3,4))
lines(rep(c4.5L, 2), c(4,5))
curve(expected(x, c1.5L, c2.5L, c3.5L, c4.5L, sigma=sigmaL), add=TRUE)
points(invCond$RichTOL[temp],
       jitter(invCond$SEQ[temp], 0.25), pch=1,
            col="gray", cex=0.5)
dev.off()

## UpLand `Upland=1`
betaU <- beta2[1]+beta2[3]
alphaU <- beta2[2]

c1.5U <- (kappa2[1]-alphaU)/betaU
c2.5U <- (kappa2[2]-alphaU)/betaU
c3.5U <- (kappa2[3]-alphaU)/betaU
c4.5U <- (kappa2[4]-alphaU)/betaU
sigmaU <- 1/betaU

temp=invCond$Upland==1
pdf(file="polr1U.pdf", height=3, width=3)
par(mar=c(3,5,0.25,0.25), mgp=c(1.5,0.25,0), tck=-0.005)
plot(0, 0, xlim=range(invCond$RichTOL), ylim=c(1,5), xlab="RichTOL",
     ylab="", type="n", axes=F)
axis(1)
axis(2, at=1:5, labels=c("Excellent","Good","Good-Fair","Fair","Poor"), las=1)
lines(rep(c1.5U, 2), c(1,2))
lines(rep(c2.5U, 2), c(2,3))
lines(rep(c3.5U, 2), c(3,4))
lines(rep(c4.5U, 2), c(4,5))
curve(expected(x, c1.5U, c2.5U, c3.5U, c4.5U, sigma=sigmaU), add=TRUE)
points(invCond$RichTOL[temp],
       jitter(invCond$SEQ[temp], 0.25), pch=1,
            col="gray", cex=0.5)
dev.off()

POLR3 <- polr(ordered(SEQ) ~ RichTOL*factor(ERIII), data=invCond)
summary(POLR3)
beta3 <- coef(POLR3)
kappa3 <- POLR3$zeta
## ER 45
beta45 <- beta3[1]

c1.5ER45 <- (kappa3[1])/beta45
c2.5ER45 <- (kappa3[2])/beta45
c3.5ER45 <- (kappa3[3])/beta45
c4.5ER45 <- (kappa3[4])/beta45
sigma45 <- 1/beta45

temp=invCond$ERIII==45
pdf(file="polr1ER45.pdf", height=3, width=3)
par(mar=c(3,5,0.25,0.25), mgp=c(1.5,0.25,0), tck=-0.005)
plot(0, 0, xlim=range(invCond$RichTOL), ylim=c(1,5), xlab="RichTOL",
     ylab="", type="n", axes=F)
axis(1)
axis(2, at=1:5, labels=c("Excellent","Good","Good-Fair","Fair","Poor"), las=1)
lines(rep(c1.5ER45, 2), c(1,2))
lines(rep(c2.5ER45, 2), c(2,3))
lines(rep(c3.5ER45, 2), c(3,4))
lines(rep(c4.5ER45, 2), c(4,5))
curve(expected(x, c1.5ER45, c2.5ER45, c3.5ER45, c4.5ER45, sigma=sigma45),
      add=TRUE)
points(invCond$RichTOL[temp],
       jitter(invCond$SEQ[temp], 0.25), pch=1,
            col="gray", cex=0.5)
dev.off()

## ER 63
beta63 <- beta3[1]+beta3[4]
alpha63 <- beta3[2]

c1.5ER63 <- (kappa2[1]-alpha63)/beta63
c2.5ER63 <- (kappa2[2]-alpha63)/beta63
c3.5ER63 <- (kappa2[3]-alpha63)/beta63
c4.5ER63 <- (kappa2[4]-alpha63)/beta63
sigma63 <- 1/beta63


temp=invCond$ERIII==63
pdf(file="polr1ER63.pdf", height=3, width=3)
par(mar=c(3,5,0.25,0.25), mgp=c(1.5,0.25,0), tck=-0.005)
plot(0, 0, xlim=range(invCond$RichTOL), ylim=c(1,5), xlab="RichTOL",
     ylab="", type="n", axes=F)
axis(1)
axis(2, at=1:5, labels=c("Excellent","Good","Good-Fair","Fair","Poor"), las=1)
lines(rep(c1.5ER63, 2), c(1,2))
lines(rep(c2.5ER63, 2), c(2,3))
lines(rep(c3.5ER63, 2), c(3,4))
lines(rep(c4.5ER63, 2), c(4,5))
curve(expected(x, c1.5ER63, c2.5ER63, c3.5ER63, c4.5ER63, sigma=sigma63),
      add=TRUE)
points(invCond$RichTOL[temp],
       jitter(invCond$SEQ[temp], 0.25), pch=1,
            col="gray", cex=0.5)
dev.off()

## ER 65
beta65 <- beta3[1]+beta3[5]
alpha65 <- beta3[3]

c1.5ER65 <- (kappa2[1]-alpha65)/beta65
c2.5ER65 <- (kappa2[2]-alpha65)/beta65
c3.5ER65 <- (kappa2[3]-alpha65)/beta65
c4.5ER65 <- (kappa2[4]-alpha65)/beta65
sigma65 <- 1/beta65


temp=invCond$ERIII==65
pdf(file="polr1ER65.pdf", height=3, width=3)
par(mar=c(3,5,0.25,0.25), mgp=c(1.5,0.25,0), tck=-0.005)
plot(0, 0, xlim=range(invCond$RichTOL), ylim=c(1,5), xlab="RichTOL",
     ylab="", type="n", axes=F)
axis(1)
axis(2, at=1:5, labels=c("Excellent","Good","Good-Fair","Fair","Poor"), las=1)
lines(rep(c1.5ER65, 2), c(1,2))
lines(rep(c2.5ER65, 2), c(2,3))
lines(rep(c3.5ER65, 2), c(3,4))
lines(rep(c4.5ER65, 2), c(4,5))
curve(expected(x, c1.5ER65, c2.5ER65, c3.5ER65, c4.5ER65, sigma=sigma65),
      add=TRUE)
points(invCond$RichTOL[temp],
       jitter(invCond$SEQ[temp], 0.25), pch=1,
            col="gray", cex=0.5)
dev.off()
