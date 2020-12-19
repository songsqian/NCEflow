## Initializing
to3 <- function(x){
  ifelse (x < 10, paste("00", x, sep=""), ifelse(x < 100, paste("0", x, sep=""), as.character(x)))
}

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
packages(maptools)
packages(maps)
packages(mapproj)
packages(rpart)
packages(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

nchains <-  min(c(parallel::detectCores(), 4))
niters <- 10000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

packages(ggplot2)

base <- getwd()
##base <- "~/Dropbox/Flow-Ecology\ Modeling/R"
dataDIR <- paste(base, "Data", sep="/")##"~/Dropbox/Flow-Ecology\ Modeling/Data"
plotDIR <- paste(base, "Figs", sep="/")##"~/Dropbox/Flow-Ecology\ Modeling/Figs"

##source("https://github.com/songsqian/eesR/R/eesrfuns.r")
hockey <- function(x,alpha1,beta1,beta2,brk,eps=diff(range(x))/100,
                   delta=T) {
       ## alpha1 is the intercept of the left line segment
       ## beta1 is the slope of the left line segment
       ## beta2 is the slope of the right line segment
       ## brk is location of the break point
       ## 2*eps is the length of the connecting quadratic piece

       ## reference: Bacon & Watts "Estimating the Transition Between
       ## Two Intersecting Straight Lines", Biometrika, 1971
        x <- x-brk
        if (delta) beta2 <- beta1+beta2
        x1 <- -eps
        x2 <- +eps
        b <- (x2*beta1-x1*beta2)/(x2-x1)
        cc <- (beta2-b)/(2*x2)
        a <- alpha1+beta1*x1-b*x1-cc*x1^2
        alpha2 <- - beta2*x2 +(a + b*x2 + cc*x2^2)

        lebrk <- (x <= -eps)
        gebrk <- (x >= eps)
        eqbrk <- (x > -eps & x < eps)

        result <- rep(0,length(x))
        result[lebrk] <- alpha1 + beta1*x[lebrk]
        result[eqbrk] <- a + b*x[eqbrk] + cc*x[eqbrk]^2
        result[gebrk] <- alpha2 + beta2*x[gebrk]
        result
}

sim.nls <- function (object, n.sims=100){
# sim.nls:  get posterior simulations of sigma and beta from a nls object
#
# Arguments:
#
#     object:  the output of a call to "nls"
#              with n data points and k predictors
#     n.sims:  number of independent simulation draws to create
#
# Output is a list (sigma.sim, beta.sim):
#
#     sigma.sim:  vector of n.sims random draws of sigma
#       (for glm's, this just returns a vector of 1's or else of the
#       square root of the overdispersion parameter if that is in the model)
#     beta.sim:  matrix (dimensions n.sims x k) of n.sims random draws of beta
#

  object.class <- class(object)[[1]]
  if (object.class!="nls") stop("not a nls object")

    summ <- summary (object)
    coef <- summ$coef[,1:2,drop=FALSE]
    dimnames(coef)[[2]] <- c("coef.est","coef.sd")
    sigma.hat <- summ$sigma
    beta.hat <- coef[,1]
    V.beta <- summ$cov.unscaled
    n <- summ$df[1] + summ$df[2]
    k <- summ$df[1]
    sigma <- rep (NA, n.sims)
    beta <- array (NA, c(n.sims,k))
    dimnames(beta) <- list (NULL, names(beta.hat))
    for (s in 1:n.sims){
      sigma[s] <- sigma.hat*sqrt((n-k)/rchisq(1,n-k))
      beta[s,] <- mvrnorm (1, beta.hat, V.beta*sigma[s]^2)
    }
    return (list (beta=beta, sigma=sigma))
  }
#Preparation
##Reading data

## sites and samples
sites <- read.csv(paste(dataDIR, "Inv_Sites.csv", sep="/")) ## not needed, all infor are in the next file
sites_cf <- read.csv(paste(dataDIR, "CF_sites.csv", sep="/"))
sites_pd <- read.csv(paste(dataDIR, "PD_sites.csv", sep="/"))

sites <- rbind(sites_cf, sites_pd)
samples <- read.csv(paste(dataDIR, "Inv_Samples.csv", sep="/"))

## Inv data
cfRichTol <- read.csv(paste(dataDIR, "CF_NoAmbig_TOL_SWAT.csv", sep="/"))
ydRichTol <- read.csv(paste(dataDIR, "YAD_NCSC_NoAmbig_TOL_SWAT.csv", sep="/"))

cfRich <- read.csv(paste(dataDIR, "CF_NoAmbig_R_SWAT.csv", sep="/"))
ydRich <- read.csv(paste(dataDIR, "YAD_NCSC_NoAmbig_R_SWAT.csv", sep="/"))

cfRichtraits <- read.csv(paste(dataDIR, "CF_Rich_traits.csv", sep="/"))
ydRichtraits <- read.csv(paste(dataDIR, "YAD_Rich_traits.csv", sep="/"))
names(ydRichtraits)[22] <- "Sync_1"
names(ydRichtraits)[41] <- "Swim_1"

## LU data
cfLanduse <- read.csv(paste(dataDIR, "CF_Basin_LU_INV_samps.csv", sep="/"))
ydLanduse <- read.csv(paste(dataDIR, "YAD_Basin_LU_INV_Samp.csv", sep="/"))

cfLanduseP <- read.csv(paste(dataDIR, "CF_Basin_pLU_INV_Samp.csv", sep="/"))
ydLanduseP <- read.csv(paste(dataDIR, "YAD_Basin_pLU_INV_Samp.csv", sep="/"))

## Env data
envTemp <- read.csv(paste(dataDIR, "Inv_Air_Temp_Data_at_Site.csv", sep="/"))
envHumi <- read.csv(paste(dataDIR, "Inv_Humidity_Data_at_Site.csv", sep="/"))
envPrec <- read.csv(paste(dataDIR, "Inv_Precip_Data_at_Site.csv", sep="/"))

## Flow data
pdFlowyr <- read.csv(paste(dataDIR, "PeeDeeFlow", "PeeDee_EFlow_Stats_by_year.csv", sep="/"))
cfFlowyr <- read.csv(paste(dataDIR, "CapeFearFlow", "CapeFear_EFlow_Stats_by_year.csv", sep="/"))

pdFlowPOR <- read.csv(paste(dataDIR, "PeeDeeFlow", "PeeDee_Flow_Metrics_POR.csv", sep="/"))
cfFlowPOR <- read.csv(paste(dataDIR, "CapeFearFlow", "Cape_Fear_Flow_Metrics_POR.csv", sep="/"))
## Data processing
## merging by sample ID
### cf + yad
RichTol <- rbind(cfRichTol, ydRichTol[,-2])
Rich <- rbind(cfRich, ydRich[,-2])
Rtraits <- rbind(cfRichtraits, ydRichtraits[,-2])
landuse <- rbind(cfLanduse, ydLanduse)
landuseP <- rbind(cfLanduseP, ydLanduseP)

data <- merge(envTemp, envHumi[,c(5, 14:44)], by="SampleID")
data <- merge(data, envPrec[,c(5, 14:29)])
data <- merge(data, landuse, by="SampleID")
data <- merge(data, landuseP, by="SampleID")
data <- merge(data, Rich[,c(6, 15:45)])
data <- merge(data, RichTol[,c(6, 15:25)])
data$site_yr <- paste(data$EcositeID.x, data$YR.x)
cfFlowyr$site_yr <- paste(cfFlowyr$site_no,
                          format(as.Date(cfFlowyr$min_date, "%m/%d/%Y"), "%Y"))
pdFlowyr$site_yr <- paste(pdFlowyr$site_no,
                          format(as.Date(pdFlowyr$min_date, "%m/%d/%y"), "%Y"))

dataPOR <- merge(data, rbind(cfFlowPOR, pdFlowPOR), by.x="EcositeID.x",
                 by.y="site_no")
dataYR  <- merge(data, rbind(cfFlowyr, pdFlowyr), by="site_yr")
## something is wrong here, LU_xx disappeared
## resolved (EcositeID_.x_)
## not all years have matching flow data (297 of the 1371 obs have no flow)

temp <- data.frame(EcositeID=dataYR$EcositeID.x,
                   Lat=dataYR$Lat,
                   Long=dataYR$Long,
                   RICH=dataYR$RICH,
                   EPTR = dataYR$EPTR,
                   RICHtol = dataYR$RichTOL,
                   Dev=dataYR$Develop,
                   semiDev=dataYR$SemiDevel,
                   pDev=dataYR$pDevelop,
                   psemiDev=dataYR$pSemiDevel,
                   Ag=dataYR$LU_41+dataYR$LU_43+dataYR$LU_44+dataYR$LU_45,
                   pAg=dataYR$pLU_41+dataYR$pLU_43+dataYR$pLU_44+dataYR$pLU_45,
                   Natural=dataYR$LU_50+dataYR$LU_60,
                   pNatural=dataYR$pLU_50+dataYR$pLU_60,
                   spTemp=dataYR$SpMeanAvgTemp,
                   spSDTemp=dataYR$SpSDAvgTemp,
                   faTemp=dataYR$FallMeanAvgTemp,
                   faSDTemp=dataYR$FallSDAvgTemp,
                   suTemp=dataYR$SumMeanAvgTemp,
                   suSDTemp=dataYR$SumSDAvgTemp,
                   wnTemp=dataYR$WinMeanAvgTemp,
                   wnSDTemp=dataYR$WinSDAvgTemp,
                   spHumid=dataYR$SpMeanDailyAvgHumid,
                   spSDHumid=dataYR$SpSDDailyAvgHumid,
                   spPrecip=dataYR$SpMeanDailyPrecip,
                   spSDprecip = dataYR$SpSDDailyPrecip,
                   log_annual_mu= log(dataYR$ma1_mean_disc)-
                       log(1+(0.01*dataYR$ma3_mean_annual_var)^2),
                   log_annual_sigma=
                       sqrt(log(1+(0.01*dataYR$ma3_mean_annual_var)^2)))

temp<- merge(temp, sites[,c("EcositeID", "Elev", "DA", "Level.3.Ecoregion",
                            "Level.4.Ecoregion")], by="EcositeID")
tempFlow <- cbind(temp, dataYR[,210:361])

## Data Preparation 2
##Importing basin characteristics data files and meage with existing flow data
ncstreams <- read.csv(paste(dataDIR, "CFStreamStats.csv", sep="/"))
ncstreams$Name <- as.character(ncstreams$Name)
ncstreams$Name <- paste(substring(ncstreams$Name, 1, 2), to3(as.numeric(substring(ncstreams$Name, 3, 5))), sep="")

scstreams <- read.csv(paste(dataDIR, "SC_attributes.csv", sep="/"))
scstreams$Name <- as.character(scstreams$Name)
scstreams$Name <- paste(substring(scstreams$Name, 1, 2), to3(as.numeric(substring(scstreams$Name, 3, 5))), sep="")

BasinChara <- data.frame(Name=as.character(c(ncstreams$Name, scstreams$Name)),
                         DrainageA=c(ncstreams$DRNAREA, scstreams$ADJDA_SQMI),
                         elev_outlet=c(ncstreams$OUTLETELEV, scstreams$OUTLETELEV),
                         elev=c(ncstreams$ELEV, scstreams$ELEV),
                         elev_min=c(ncstreams$MINBELEV, scstreams$MINBELEV),
                         elev_max=c(ncstreams$ELEVMAX, scstreams$ELEVMAX),
                         slope_basin=c(ncstreams$BSLDEM30FT, scstreams$BSLDEM30FT),
                         slope_stream=c(ncstreams$CSL10_85FM, scstreams$CSL10_85FM),
                         perim_basin=c(ncstreams$BASINPERIM, scstreams$BASINPERIM),
                         precip=c(ncstreams$PRECIP, scstreams$PRECIP),
                         precip_24h50y=c(ncstreams$I24H50Y, scstreams$I24H50Y),
                         sgA = c(ncstreams$SSURGOA, scstreams$SSURGOA),
                         sgB = c(ncstreams$SSURGOB, scstreams$SSURGOB),
                         sgC = c(ncstreams$SSURGOC, scstreams$SSURGOC),
                         sgD = c(ncstreams$SSURGOD, scstreams$SSURGOD))

predFlow <- tempFlow[, 1:32]
predFlow$ml20 <- tempFlow$ml20
predFlow$dh20 <- tempFlow$dh20
predFlow$EcositeID <- as.character(predFlow$EcositeID)

dataAll <- merge(predFlow, BasinChara, by.x="EcositeID", by.y="Name")


