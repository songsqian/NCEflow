## NC Eflow Project
## Gibbs sampler version of the cBN model

## Front matters
##

source("EcFl1_data.R") ## R code for processing data
packages(reshape2)
packages(maps)
packages(maptools)
packages(mapproj)
packages(tikzDevice)
rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

## Data
wkdata_all <- data.frame(ye = dataAll$RICHtol,
                     y1 = dataAll$log_annual_sigma^2,
                     y2 = dataAll$log_annual_mu,
                     y3 = dataAll$ml20,
                     y4 = log(dataAll$dh20),
                     x1 = car::logit(dataAll$pNatural),
                     x2 = dataAll$sgC,
                     x3 = dataAll$precip,
                     x4 = log(dataAll$DA),
                     x5 = dataAll$suTemp,
                     x6 = dataAll$spTemp,
                     EcositeID = dataAll$EcositeID,
                     gr = dataAll$Level.3.Ecoregion,
                     Lat= dataAll$Lat,
                     Long=dataAll$Long)
## removel precip= 0
wkdata <- wkdata_all[wkdata_all$x3>0,]
wkdata$x3 <- log(wkdata$x3)
wkdata$gr2 <- as.numeric(wkdata$gr==63)
wkdata$gr3 <- as.numeric(wkdata$gr==63 | wkdata$gr==65 | wkdata$gr==75)

wkdata_melt <- melt(wkdata, id = "EcositeID",
                    measured=c("y1", "y2", "y3","y4",
                               "x1", "x2", "x3","x4", "x5","x6",
                               "gr","gr2","gr3"))
wkdata_site <- dcast(wkdata_melt, EcositeID~variable, mean)

wkdata_melt2 <- melt(wkdata_all, id = "EcositeID",
                    measured=c("y1", "y2", "y3","y4",
                               "x1", "x2", "x3","x4", "x5","x6",
                               "gr","gr2","gr3"))
wkdata_site2 <- dcast(wkdata_melt2, EcositeID~variable, mean)

## EcositeID mean

## Fitting individual models
## M1 <- glm(y1 ~ x1+x2+x3, data=wkdata, family=Gamma(link="log"))
M2 <- lm(y2 ~ x3+x4, data=wkdata)
M2_2 <- lm(y2 ~ x3*gr2+x4*gr2, data=wkdata)
M2_3 <- lm(y2 ~ x3*gr3+x4*gr3, data=wkdata)

M3 <- lm(y3 ~ x1+x2+x3, data=wkdata)
M3_2 <- lm(y3 ~ x1*gr2+x2*gr2+x3*gr2, data=wkdata)
M3_3 <- lm(y3 ~ x1*gr3+x2*gr3+x3*gr3, data=wkdata)

M4 <- lm(y4 ~ x1+x2, data=wkdata)
M4_2 <- lm(y4 ~ x1*gr2+x2*gr2, data=wkdata)
M4_3 <- lm(y4 ~ x1*gr3+x2*gr3, data=wkdata)

M5 <- lm(ye ~ x1+y2+y4+x5+y3, data=wkdata)
M5_2 <- lm(ye ~ x1*gr2+y2+y3*gr2+x5*gr2, data=wkdata)
M5_3 <- lm(ye ~ x1*gr3+y2+y3*gr3+x5*gr3, data=wkdata)

### Gibbs sampler version
nsims=10000
##simM1 <- sim(M1_2, n.sims=nsims)
simM2 <- sim(M2_2, n.sims=nsims)
simM3 <- sim(M3_2, n.sims=nsims)
simM4 <- sim(M4_2, n.sims=nsims)
coef5_sim <- matrix(0, nrow=nsims, ncol=length(coef(M5_2)))
sig5_sim <- numeric()
for (i in 1:nsims){
    mu2 <- model.matrix(~ x3 * gr2 + x4 * gr2, data = wkdata)%*%simM2@coef[i,]
    mu3 <- model.matrix(~ x1*gr2+x2*gr2+x3*gr2, data = wkdata)%*%simM3@coef[i,]
    lmData <- data.frame(ye=wkdata$ye, x1=wkdata$x1,
                         mu2=mu2, mu3=mu3,x3=wkdata$x3,
                         x5=wkdata$x5, gr2=wkdata$gr2)
    M <- lm(ye ~ x1*gr2+mu2+mu3*gr2+x5*gr2, data=lmData)
    coef5_sim[i,] <- coef(M)
    sig5_sim[i] <- summary(M)$sigma
}
coef5_rv <- rvsims(coef5_sim)
sig5_rv <- rvsims(sig5_sim)


## readin projected climate and developmental data
packages(readxl)
predDIR <- paste(base, "Future_clim", sep="/")
luDIR <- paste(base, "Future_LU", sep="/")

read_clim <- function(DIR=predDIR,
                      file="inmcm4_CF_pr.xlsx"){
    file <- paste(DIR, file, sep="/")
    shts <- excel_sheets(file)
    temp <- list()
    for (i in 1:length(shts))
        temp[[i]] <- read_excel(file, sheet=shts[i])
    return(temp)
}

read_LU <- function(DIR=luDIR,
                    file="CF/CF_sq_LU_time.xlsx"){
    file <- paste(DIR, file, sep="/")
    shts <- excel_sheets(file)
    temp <- list()
    for (i in 1:length(shts))
        temp[[i]] <- read_excel(file, sheet=shts[i])
    return(temp)
}

## data importing function climate_model plus LU status devST
climate_data <- function(climate_model="inmcm4"){
    CF_PR <- read_clim(file=paste(climate_model, "CF_pr.xlsx", sep="_"))
    CF_TM <- read_clim(file=paste(climate_model, "CF_temp.xlsx", sep="_"))
    PD_PR <- read_clim(file=paste(climate_model, "PD_pr.xlsx", sep="_"))
    PD_TM <- read_clim(file=paste(climate_model, "PD_temp.xlsx", sep="_"))
    x3_pred_cf <- CF_PR[[1]]
    x5_pred_cf <- CF_TM[[3]]
    x3_pred_pd <- PD_PR[[1]]
    x5_pred_pd <- PD_TM[[3]]

    ## adding leading 0(s) to Ecosit (cf files only)
    x3_pred_cf$EcositeID <- paste(substring(x3_pred_cf$Ecosite_Name, 1, 2),
                                  to3(as.numeric(substring(x3_pred_cf$Ecosite_Name,3, 5))), sep="")
    x3_pred_pd$EcositeID <- paste(substring(x3_pred_pd$Ecosite_Name, 1, 2),
                                  to3(as.numeric(substring(x3_pred_pd$Ecosite_Name,3, 5))), sep="")
    x5_pred_cf$EcositeID <- paste(substring(x5_pred_cf$Ecosite_Name, 1, 2),
                                  to3(as.numeric(substring(x5_pred_cf$Ecosite_Name,3, 5))),sep="")
    x5_pred_pd$EcositeID <- paste(substring(x5_pred_pd$Ecosite_Name, 1, 2),
                                  to3(as.numeric(substring(x5_pred_pd$Ecosite_Name,3, 5))),sep="")
    return(list(x3_pred=rbind(x3_pred_cf,x3_pred_pd),
                x5_pred=rbind(x5_pred_cf,x5_pred_pd)))
}

## Landuse data
LU_files <- function(data=CFsq_LU, EcoSite="Name"){
    forest <- data[[2]]
    forest$EcositeID <- paste(substring(forest[[EcoSite]], 1, 2),
                              to3(as.numeric(substring(forest[[EcoSite]],3, 5))),
                              sep="")
    oo <- order(forest$EcositeID)
    forest <- forest[oo,]
    forest1 <- forest[,names(forest)%in%as.character(seq(2015,2065,5))]

    shrub <- data[[3]]
    shrub$EcositeID <- paste(substring(shrub[[EcoSite]], 1, 2),
                             to3(as.numeric(substring(shrub[[EcoSite]],3, 5))),
                             sep="")
    oo <- order(shrub$EcositeID)
    shrub <- shrub[oo,]
    shrub1 <- shrub[,names(shrub)%in%as.character(seq(2015,2065,5))]

    wetland <- data[[5]]
    wetland$EcositeID <- paste(substring(wetland[[EcoSite]], 1, 2),
                        to3(as.numeric(substring(wetland[[EcoSite]],3, 5))),
                               sep="")
    oo <- order(wetland$EcositeID)
    wetland <- wetland[oo,]
    wetland1 <- wetland[,names(wetland)%in%as.character(seq(2015,2065,5))]

    pNature <- cbind(forest[["EcositeID"]], forest1 + shrub1 + wetland1)
    names(pNature)[1]  <- "EcositeID"
    return(na.omit(pNature))
}

landuse_data <- function(){
    CFsq_LU <- LU_files(data=read_LU(file="CF/CF_sq_LU_time.xlsx"),
                        EcoSite="Name")
    CFws_LU <- LU_files(data=read_LU(file="CF/CF_ws_LU_time.xlsx"),
                        EcoSite="Name")
    PDsq_LU <- LU_files(data=read_LU(file="PD/PD_sq_LU_time.xlsx"),
                        EcoSite="EcoSite")
    PDws_LU <- LU_files(data=read_LU(file="PD/PD_ws_LU_time.xlsx"),
                        EcoSite="EcoSite")
    return(list(stQ=rbind(CFsq_LU, PDsq_LU), wtS=rbind(CFws_LU, PDws_LU)))
}

LUdata_all  <- landuse_data()
pdf(file=paste(base, "Figs", "LandUseFuture.pdf", sep="/"),
    height=4, width=5)
par(mfrow=c(2,1), mar=c(2,2,1,2), oma=c(2,3,1,1),
    mgp=c(1.25,0.125,0), tck=0.01)
boxplot(LUdata_all$stQ[,-1])
mtext("sq",side=4, las=0)
boxplot(LUdata_all$wtS[,-1])
mtext("ws",side=4, las=0)
mtext("Year", side=1, outer=T)
mtext("% Natural", side=2, outer=T, las=0)
dev.off()

Pr_Temp1_all <- climate_data(climate_model="inmcm4")
Pr_Temp2_all <- climate_data(climate_model="MRI-CGCM3")
Pr_Temp3_all <- climate_data(climate_model="NorESM1-M")

##pdf(file=paste(plotDIR, "FuturePred.pdf", sep="/"),
##    height=7, width=6)
##png(file=paste(plotDIR, "FuturePred.png", sep="/"),
##    height=7*120, width=6*120)
tikz(file=paste(plotDIR, "FuturePred.tex", sep="/"),
     height=7, width=6, standAlone=F)
par(mfrow=c(3,2), mar=c(2.5,3.5,1,0.5), oma=c(3,3,1,1),
    mgp=c(1.25,0.125,0), tck=0.01)
ttt <- names(Pr_Temp1_all$x3_pred)%in%seq(2015,2065,5)
boxplot(Pr_Temp1_all$x3_pred[,ttt], ylim=c(880,2700))
mtext("Precipitation (mm)", side=3, outer=F)
mtext("inmcm4",side=2, las=0, line=1.75)
ttt <- names(Pr_Temp1_all$x5_pred)%in%seq(2015,2065,5)
boxplot(Pr_Temp1_all$x5_pred[,ttt], ylim=c(20,30))
mtext("Temperature (C)", side=3, outer=F)

ttt <- names(Pr_Temp2_all$x3_pred)%in%seq(2015,2065,5)
boxplot(Pr_Temp2_all$x3_pred[,ttt], ylim=c(880,2700))
mtext("MRI-CGCM3",side=2, las=0, line=1.75)
ttt <- names(Pr_Temp2_all$x5_pred)%in%seq(2015,2065,5)
boxplot(Pr_Temp2_all$x5_pred[,ttt], ylim=c(20,30))

ttt <- names(Pr_Temp3_all$x3_pred)%in%seq(2015,2065,5)
boxplot(Pr_Temp3_all$x3_pred[,ttt], ylim=c(880,2700))
mtext("NorESM1-M",side=2, las=0, line=1.75)
ttt <- names(Pr_Temp3_all$x5_pred)%in%seq(2015,2065,5)
boxplot(Pr_Temp3_all$x5_pred[,ttt], ylim=c(20,30))

mtext("Year", side=1, outer=T)
dev.off()

extract_data <- function(x1_pred=x1_all$stQ, x3_pred=x3x5_all$x3_pred,
                         x5_pred=x3x5_all$x5_pred,
                         fit_data=wkdata_site2, yr="2020"){
    if (yr=="Base") yr1 <- "2015"
    else yr1 <- yr  ## LU data do not have 2014
    x1_pred <- x1_pred[, c("EcositeID",yr1)]
    names(x1_pred) <- c("EcositeID","x1")
    x1_pred$x1 <- car::logit(x1_pred$x1)
    ##x1_pred <- x1_pred[!is.na(x1_pred$x1),]
    if (yr=="Base"){
        x3_pred <- cbind(x3_pred[,"EcositeID"],
                         apply(x3_pred[,as.character(2014:2019)], 1, mean))
        x5_pred <- cbind(x5_pred[,"EcositeID"],
                         apply(x5_pred[,as.character(2014:2019)], 1, mean))
    }else{
        x3_pred <- x3_pred[, c("EcositeID",yr)]
        x5_pred <- x5_pred[, c("EcositeID",yr)]
    }
    names(x3_pred) <- c("EcositeID", "x3")
    x3_pred$x3 <- log(x3_pred$x3/25.4)  ## log precipitation (in.)
    names(x5_pred) <- c("EcositeID", "x5")
    x3x5<- merge(x3_pred, x5_pred)
    x1x3x5 <- merge(x1_pred, x3x5)
    data <- merge(fit_data[, names(fit_data)!="x1" &names(fit_data)!="x3" &
                             names(fit_data)!="x5"],
                  x1x3x5, by="EcositeID")
    ## replacing precip & temp with projected values
    return(data)
}

## Function for making predictions
meansRV <- function(Fit_mu= simM2,
                    Fit_ml= simM3,
                    Fit_dh= simM4,
                    Fit_rt= coef5_sim, rt_sig=sig5_sim,
                    predDATA = pred_data) {
    coef2 <- rvsims(Fit_mu@coef)
    coef3 <- rvsims(Fit_ml@coef)
    coef4 <- rvsims(Fit_dh@coef)
    coef5 <- rvsims(Fit_rt)

    EcoR <- predDATA$gr2
    sigma2 <- rvsims(Fit_mu@sigma)
    sigma3 <- rvsims(Fit_ml@sigma)
    sigma4 <- rvsims(Fit_dh@sigma)
    sigma5 <- rvsims(rt_sig)

    mu2 <- coef2[1]+coef2[2]*predDATA$x3+
        coef2[3]*predDATA$gr2+coef2[4]*predDATA$x4+
        coef2[5]*predDATA$x3*predDATA$gr2+
        coef2[6]*predDATA$x4*predDATA$gr2
    mu3 <- coef3[1]+coef3[2]*predDATA$x1+coef3[3]*predDATA$gr2+
        coef3[4]*predDATA$x2 +coef3[5]*predDATA$x3+
        coef3[6]*predDATA$x1*predDATA$gr2+
        coef3[7]*predDATA$x2*predDATA$gr2+
        coef3[8]*predDATA$x3*predDATA$gr2
    mu4 <- coef4[1]+coef4[2]*predDATA$x1 +
        coef4[3]*predDATA$gr2 + coef4[4]*predDATA$x2 +
        coef4[5]*predDATA$gr2 * predDATA$x1+
        coef4[6]*predDATA$gr2 * predDATA$x2
    mu5 <- coef5[1]+coef5[2]*predDATA$x1 +
        coef5[3] * predDATA$gr2 + coef5[4]*mu2 + coef5[5]*mu3+
        coef5[6]*predDATA$x5 + coef5[7]*predDATA$x1*predDATA$gr2+
        coef5[8]*mu3*predDATA$gr2+coef5[9]*predDATA$x5*predDATA$gr2

    mus <- rvnorm(1, mu2, sigma2)
    y3 <- rvnorm(1, mu3, sigma3)
    y4 <- rvnorm(1, mu4, sigma4)
    ye <- rvnorm(1, mu5, sigma5)

    return(list(mu2=mu2, mu3=mu3, mu4=mu4, mu5=mu5, sigma2=sigma2,
                sigma3=sigma3, sigma4=sigma4, sigma5=sigma5,
                f_log_mus=mus, ml20=y3, log_dh20=y4, richtol=ye,
                EcositeID=predDATA$EcositeID, EcoR=predDATA$gr))
}

### six scenarios ###
### three climate models -- "inmcm4", "MRI-CGCM3", and "NorESM1-M"
### two development scenarios -- sq & ws

simulation_scen <- function(x1=x1_all, CM="inmcm4", devS="stQ"){
    return(list(
        base =meansRV(predDATA=extract_data(x1_pred=x1_all[[devS]],yr="Base")),
        f2020=meansRV(predDATA=extract_data(x1_pred=x1_all[[devS]],yr="2020")),
        f2025=meansRV(predDATA=extract_data(x1_pred=x1_all[[devS]],yr="2025")),
        f2030=meansRV(predDATA=extract_data(x1_pred=x1_all[[devS]],yr="2030")),
        f2035=meansRV(predDATA=extract_data(x1_pred=x1_all[[devS]],yr="2035")),
        f2040=meansRV(predDATA=extract_data(x1_pred=x1_all[[devS]],yr="2040")),
        f2045=meansRV(predDATA=extract_data(x1_pred=x1_all[[devS]],yr="2045")),
        f2050=meansRV(predDATA=extract_data(x1_pred=x1_all[[devS]],yr="2050")),
        f2055=meansRV(predDATA=extract_data(x1_pred=x1_all[[devS]],yr="2055")),
        f2060=meansRV(predDATA=extract_data(x1_pred=x1_all[[devS]],yr="2060")),
        f2065=meansRV(predDATA=extract_data(x1_pred=x1_all[[devS]],yr="2065")))
        )
}

x1_all  <- landuse_data()  ## includes both dev scenarios
### 1. "inmcm4" & sq: sc1
x3x5_all <- climate_data(climate_model="inmcm4")
scnr1 <- simulation_scen(CM="inmcm4", devS="stQ")
### 2. "inmcm4" & ws: sc2
scnr2 <- simulation_scen(CM="inmcm4", devS="wtS")
### 3. "MRI-CGCM3" & sq: sc3
x3x5_all <- climate_data(climate_model="MRI-CGCM3")
scnr3 <- simulation_scen(CM="MRI-CGCM3", devS="stQ")
### 4. "MRI-CGCM3" & ws: sc4
scnr4 <- simulation_scen(CM="MRI-CGCM3", devS="wtS")
### 5. "NorESM1-M" & sq: sc5
x3x5_all <- climate_data(climate_model="NorESM1-M")
scnr5 <- simulation_scen(CM="NorESM1-M", devS="stQ")
### 6. "NorESM1-M" & ws: sc6
scnr6 <- simulation_scen(CM="NorESM1-M", devS="wtS")

## Figures
##pdf(file=paste(plotDIR,"predictionMean.pdf",sep="/"),width=6,height=7)
png(file=paste(plotDIR,"predictionMean.png",sep="/"),width=6*120,height=7*120)
par(mfrow=c(3,2), oma=c(3,3,1,1), mar=c(2.5,3.5,1,0.5),
    mgp=c(1.25,0.125,0),las=1,tck=0.01)
boxplot(summary(scnr1$base$richtol)[["mean"]],
        summary(scnr1$f2020$richtol)[["mean"]],
        summary(scnr1$f2025$richtol)[["mean"]],
        summary(scnr1$f2030$richtol)[["mean"]],
        summary(scnr1$f2035$richtol)[["mean"]],
        summary(scnr1$f2040$richtol)[["mean"]],
        summary(scnr1$f2045$richtol)[["mean"]],
        summary(scnr1$f2050$richtol)[["mean"]],
        summary(scnr1$f2055$richtol)[["mean"]],
        summary(scnr1$f2060$richtol)[["mean"]],
        summary(scnr1$f2065$richtol)[["mean"]],
        names=c("baseline", as.character(seq(2020,2065,5))), ylab="inmcm4")
mtext(text="sq",side=3,outer=F)
abline(h=6)

boxplot(summary(scnr2$base$richtol)[["mean"]],
        summary(scnr2$f2020$richtol)[["mean"]],
        summary(scnr2$f2025$richtol)[["mean"]],
        summary(scnr2$f2030$richtol)[["mean"]],
        summary(scnr2$f2035$richtol)[["mean"]],
        summary(scnr2$f2040$richtol)[["mean"]],
        summary(scnr2$f2045$richtol)[["mean"]],
        summary(scnr2$f2050$richtol)[["mean"]],
        summary(scnr2$f2055$richtol)[["mean"]],
        summary(scnr2$f2060$richtol)[["mean"]],
        summary(scnr2$f2065$richtol)[["mean"]],
        names=c("baseline", as.character(seq(2020,2065,5))))
mtext(text="ws",side=3,outer=F)
abline(h=6)

boxplot(summary(scnr3$base$richtol)[["mean"]],
        summary(scnr3$f2020$richtol)[["mean"]],
        summary(scnr3$f2025$richtol)[["mean"]],
        summary(scnr3$f2030$richtol)[["mean"]],
        summary(scnr3$f2035$richtol)[["mean"]],
        summary(scnr3$f2040$richtol)[["mean"]],
        summary(scnr3$f2045$richtol)[["mean"]],
        summary(scnr3$f2050$richtol)[["mean"]],
        summary(scnr3$f2055$richtol)[["mean"]],
        summary(scnr3$f2060$richtol)[["mean"]],
        summary(scnr3$f2065$richtol)[["mean"]],
        names=c("baseline", as.character(seq(2020,2065,5))), ylab="MRI-CGCM3")
abline(h=6)

boxplot(summary(scnr4$base$richtol)[["mean"]],
        summary(scnr4$f2020$richtol)[["mean"]],
        summary(scnr4$f2025$richtol)[["mean"]],
        summary(scnr4$f2030$richtol)[["mean"]],
        summary(scnr4$f2035$richtol)[["mean"]],
        summary(scnr4$f2040$richtol)[["mean"]],
        summary(scnr4$f2045$richtol)[["mean"]],
        summary(scnr4$f2050$richtol)[["mean"]],
        summary(scnr4$f2055$richtol)[["mean"]],
        summary(scnr4$f2060$richtol)[["mean"]],
        summary(scnr4$f2065$richtol)[["mean"]],
        names=c("baseline", as.character(seq(2020,2065,5))))
abline(h=6)
boxplot(summary(scnr5$base$richtol)[["mean"]],
        summary(scnr5$f2020$richtol)[["mean"]],
        summary(scnr5$f2025$richtol)[["mean"]],
        summary(scnr5$f2030$richtol)[["mean"]],
        summary(scnr5$f2035$richtol)[["mean"]],
        summary(scnr5$f2040$richtol)[["mean"]],
        summary(scnr5$f2045$richtol)[["mean"]],
        summary(scnr5$f2050$richtol)[["mean"]],
        summary(scnr5$f2055$richtol)[["mean"]],
        summary(scnr5$f2060$richtol)[["mean"]],
        summary(scnr5$f2065$richtol)[["mean"]],
        names=c("baseline", as.character(seq(2020,2065,5))), ylab="NorESM1-M")
abline(h=6)
boxplot(summary(scnr6$base$richtol)[["mean"]],
        summary(scnr6$f2020$richtol)[["mean"]],
        summary(scnr6$f2025$richtol)[["mean"]],
        summary(scnr6$f2030$richtol)[["mean"]],
        summary(scnr6$f2035$richtol)[["mean"]],
        summary(scnr6$f2040$richtol)[["mean"]],
        summary(scnr6$f2045$richtol)[["mean"]],
        summary(scnr6$f2050$richtol)[["mean"]],
        summary(scnr6$f2055$richtol)[["mean"]],
        summary(scnr6$f2060$richtol)[["mean"]],
        summary(scnr6$f2065$richtol)[["mean"]],
        names=c("baseline", as.character(seq(2020,2065,5))))
abline(h=6)
mtext(text="Predicted RichTOL Mean", las=0,
      side=2, outer=T, line=0)
mtext(text="Year", las=0,
      side=1, outer=T, line=0)
dev.off()

tikz(file=paste(plotDIR,"pred_RichTOL.tex",sep="/"),width=6,height=2,
     standAlone=F)
par(mfrow=c(1,3), oma=c(1,1,1,1), mar=c(2,2.5,1,0.5),
    mgp=c(1.25,0.125,0),las=1,tck=0.01)
boxplot(summary(scnr1$base$richtol)[["mean"]],
        summary(scnr1$f2020$richtol)[["mean"]],
        summary(scnr1$f2025$richtol)[["mean"]],
        summary(scnr1$f2030$richtol)[["mean"]],
        summary(scnr1$f2035$richtol)[["mean"]],
        summary(scnr1$f2040$richtol)[["mean"]],
        summary(scnr1$f2045$richtol)[["mean"]],
        summary(scnr1$f2050$richtol)[["mean"]],
        summary(scnr1$f2055$richtol)[["mean"]],
        summary(scnr1$f2060$richtol)[["mean"]],
        summary(scnr1$f2065$richtol)[["mean"]],
        names=c("baseline", as.character(seq(2020,2065,5))), ylab="RichTOL")
mtext(text="inmcm4",side=3,outer=F)
abline(h=6)

boxplot(summary(scnr3$base$richtol)[["mean"]],
        summary(scnr3$f2020$richtol)[["mean"]],
        summary(scnr3$f2025$richtol)[["mean"]],
        summary(scnr3$f2030$richtol)[["mean"]],
        summary(scnr3$f2035$richtol)[["mean"]],
        summary(scnr3$f2040$richtol)[["mean"]],
        summary(scnr3$f2045$richtol)[["mean"]],
        summary(scnr3$f2050$richtol)[["mean"]],
        summary(scnr3$f2055$richtol)[["mean"]],
        summary(scnr3$f2060$richtol)[["mean"]],
        summary(scnr3$f2065$richtol)[["mean"]],
        names=c("baseline", as.character(seq(2020,2065,5))), ylab="")
abline(h=6)
mtext(text="MRI-CGCM3",side=3,outer=F)

boxplot(summary(scnr5$base$richtol)[["mean"]],
        summary(scnr5$f2020$richtol)[["mean"]],
        summary(scnr5$f2025$richtol)[["mean"]],
        summary(scnr5$f2030$richtol)[["mean"]],
        summary(scnr5$f2035$richtol)[["mean"]],
        summary(scnr5$f2040$richtol)[["mean"]],
        summary(scnr5$f2045$richtol)[["mean"]],
        summary(scnr5$f2050$richtol)[["mean"]],
        summary(scnr5$f2055$richtol)[["mean"]],
        summary(scnr5$f2060$richtol)[["mean"]],
        summary(scnr5$f2065$richtol)[["mean"]],
        names=c("baseline", as.character(seq(2020,2065,5))), ylab="")
abline(h=6)
mtext(text="NorESM1-M",side=3,outer=F)
dev.off()

##pdf(file=paste(plotDIR,"predictionProb.pdf",sep="/"),width=6,height=7)
png(file=paste(plotDIR,"predictionProb.png",sep="/"),width=6*120,height=7*120)
par(mfrow=c(3,2), oma=c(3,3,1,1), mar=c(2.5,3.5,1,0.5),
    mgp=c(1.25,0.125,0),las=1,tck=0.01)
boxplot(Pr(scnr1$f2020$richtol>scnr1$base$richtol),
        Pr(scnr1$f2025$richtol>scnr1$base$richtol),
        Pr(scnr1$f2030$richtol>scnr1$base$richtol),
        Pr(scnr1$f2035$richtol>scnr1$base$richtol),
        Pr(scnr1$f2040$richtol>scnr1$base$richtol),
        Pr(scnr1$f2045$richtol>scnr1$base$richtol),
        Pr(scnr1$f2050$richtol>scnr1$base$richtol),
        Pr(scnr1$f2055$richtol>scnr1$base$richtol),
        Pr(scnr1$f2060$richtol>scnr1$base$richtol),
        Pr(scnr1$f2065$richtol>scnr1$base$richtol),
        ylim=c(0,1), names=seq(2020,2065,5), ylab="inmcm4")
mtext(text="sq",side=3,outer=F)
abline(h=0.5)

boxplot(Pr(scnr2$f2020$richtol>scnr2$base$richtol),
        Pr(scnr2$f2025$richtol>scnr2$base$richtol),
        Pr(scnr2$f2030$richtol>scnr2$base$richtol),
        Pr(scnr2$f2035$richtol>scnr2$base$richtol),
        Pr(scnr2$f2040$richtol>scnr2$base$richtol),
        Pr(scnr2$f2045$richtol>scnr2$base$richtol),
        Pr(scnr2$f2050$richtol>scnr2$base$richtol),
        Pr(scnr2$f2055$richtol>scnr2$base$richtol),
        Pr(scnr2$f2060$richtol>scnr2$base$richtol),
        Pr(scnr2$f2065$richtol>scnr2$base$richtol),
        ylim=c(0,1), names=seq(2020,2065,5))
mtext(text="ws",side=3,outer=F)
abline(h=0.5)

boxplot(Pr(scnr3$f2020$richtol>scnr3$base$richtol),
        Pr(scnr3$f2025$richtol>scnr3$base$richtol),
        Pr(scnr3$f2030$richtol>scnr3$base$richtol),
        Pr(scnr3$f2035$richtol>scnr3$base$richtol),
        Pr(scnr3$f2040$richtol>scnr3$base$richtol),
        Pr(scnr3$f2045$richtol>scnr3$base$richtol),
        Pr(scnr3$f2050$richtol>scnr3$base$richtol),
        Pr(scnr3$f2055$richtol>scnr3$base$richtol),
        Pr(scnr3$f2060$richtol>scnr3$base$richtol),
        Pr(scnr3$f2065$richtol>scnr3$base$richtol),
        ylim=c(0,1), names=seq(2020,2065,5), ylab="MRI-CGCM3")
abline(h=0.5)

boxplot(Pr(scnr4$f2020$richtol>scnr4$base$richtol),
        Pr(scnr4$f2025$richtol>scnr4$base$richtol),
        Pr(scnr4$f2030$richtol>scnr4$base$richtol),
        Pr(scnr4$f2035$richtol>scnr4$base$richtol),
        Pr(scnr4$f2040$richtol>scnr4$base$richtol),
        Pr(scnr4$f2045$richtol>scnr4$base$richtol),
        Pr(scnr4$f2050$richtol>scnr4$base$richtol),
        Pr(scnr4$f2055$richtol>scnr4$base$richtol),
        Pr(scnr4$f2060$richtol>scnr4$base$richtol),
        Pr(scnr4$f2065$richtol>scnr4$base$richtol),
        ylim=c(0,1), names=seq(2020,2065,5))
abline(h=0.5)

boxplot(Pr(scnr5$f2020$richtol>scnr5$base$richtol),
        Pr(scnr5$f2025$richtol>scnr5$base$richtol),
        Pr(scnr5$f2030$richtol>scnr5$base$richtol),
        Pr(scnr5$f2035$richtol>scnr5$base$richtol),
        Pr(scnr5$f2040$richtol>scnr5$base$richtol),
        Pr(scnr5$f2045$richtol>scnr5$base$richtol),
        Pr(scnr5$f2050$richtol>scnr5$base$richtol),
        Pr(scnr5$f2055$richtol>scnr5$base$richtol),
        Pr(scnr5$f2060$richtol>scnr5$base$richtol),
        Pr(scnr5$f2065$richtol>scnr5$base$richtol),
        ylim=c(0,1), names=seq(2020,2065,5), ylab="NorESM1-M")
abline(h=0.5)

boxplot(Pr(scnr6$f2020$richtol>scnr6$base$richtol),
        Pr(scnr6$f2025$richtol>scnr6$base$richtol),
        Pr(scnr6$f2030$richtol>scnr6$base$richtol),
        Pr(scnr6$f2035$richtol>scnr6$base$richtol),
        Pr(scnr6$f2040$richtol>scnr6$base$richtol),
        Pr(scnr6$f2045$richtol>scnr6$base$richtol),
        Pr(scnr6$f2050$richtol>scnr6$base$richtol),
        Pr(scnr6$f2055$richtol>scnr6$base$richtol),
        Pr(scnr6$f2060$richtol>scnr6$base$richtol),
        Pr(scnr6$f2065$richtol>scnr6$base$richtol),
        ylim=c(0,1), names=seq(2020,2065,5))
abline(h=0.5)
mtext(text="Pr. RichTOL>baseline", las=0,
      side=2, outer=T, line=0)
mtext(text="Year", las=0,
      side=1, outer=T, line=0)
dev.off()

tikz(file=paste(plotDIR,"pred_ProbWorse.tex",sep="/"),width=6,height=2,
     standAlone=F)
par(mfrow=c(1,3), oma=c(1,1,1,1), mar=c(2,2.5,1,0.5),
    mgp=c(1.3,0.125,0),las=1, tck=0.01)
boxplot(Pr(scnr1$f2020$richtol>scnr1$base$richtol),
        Pr(scnr1$f2025$richtol>scnr1$base$richtol),
        Pr(scnr1$f2030$richtol>scnr1$base$richtol),
        Pr(scnr1$f2035$richtol>scnr1$base$richtol),
        Pr(scnr1$f2040$richtol>scnr1$base$richtol),
        Pr(scnr1$f2045$richtol>scnr1$base$richtol),
        Pr(scnr1$f2050$richtol>scnr1$base$richtol),
        Pr(scnr1$f2055$richtol>scnr1$base$richtol),
        Pr(scnr1$f2060$richtol>scnr1$base$richtol),
        Pr(scnr1$f2065$richtol>scnr1$base$richtol),
        ylim=c(0,1), names=seq(2020,2065,5), ylab="Pr. RichTOL$>$baseline")
mtext(text="inmcm4",side=3,outer=F)
abline(h=0.5)

boxplot(Pr(scnr3$f2020$richtol>scnr3$base$richtol),
        Pr(scnr3$f2025$richtol>scnr3$base$richtol),
        Pr(scnr3$f2030$richtol>scnr3$base$richtol),
        Pr(scnr3$f2035$richtol>scnr3$base$richtol),
        Pr(scnr3$f2040$richtol>scnr3$base$richtol),
        Pr(scnr3$f2045$richtol>scnr3$base$richtol),
        Pr(scnr3$f2050$richtol>scnr3$base$richtol),
        Pr(scnr3$f2055$richtol>scnr3$base$richtol),
        Pr(scnr3$f2060$richtol>scnr3$base$richtol),
        Pr(scnr3$f2065$richtol>scnr3$base$richtol),
        ylim=c(0,1), names=seq(2020,2065,5), ylab="")
abline(h=0.5)
mtext(text="MRI-CGCM3",side=3,outer=F)

boxplot(Pr(scnr5$f2020$richtol>scnr5$base$richtol),
        Pr(scnr5$f2025$richtol>scnr5$base$richtol),
        Pr(scnr5$f2030$richtol>scnr5$base$richtol),
        Pr(scnr5$f2035$richtol>scnr5$base$richtol),
        Pr(scnr5$f2040$richtol>scnr5$base$richtol),
        Pr(scnr5$f2045$richtol>scnr5$base$richtol),
        Pr(scnr5$f2050$richtol>scnr5$base$richtol),
        Pr(scnr5$f2055$richtol>scnr5$base$richtol),
        Pr(scnr5$f2060$richtol>scnr5$base$richtol),
        Pr(scnr5$f2065$richtol>scnr5$base$richtol),
        ylim=c(0,1), names=seq(2020,2065,5), ylab="")
abline(h=0.5)
mtext(text="NorESM1-M",side=3,outer=F)
mtext(text="Year", las=0, side=1, outer=T, line=0)
dev.off()

## More plots and saving data

##pdf(file=paste(plotDIR, "probGF6.pdf", sep="/"),
##    width=6,height=7)
png(file=paste(plotDIR, "probGF6.png", sep="/"),
    width=6*120,height=7*120)
par(mfrow=c(3,2), oma=c(3,3,1,1), mar=c(2.5,3.5,1,0.5),
    mgp=c(1.25,0.125,0),las=1,tck=0.01)
boxplot(Pr(scnr1$base$richtol>6),
        Pr(scnr1$f2020$richtol>6),
        Pr(scnr1$f2025$richtol>6),
        Pr(scnr1$f2030$richtol>6),
        Pr(scnr1$f2035$richtol>6),
        Pr(scnr1$f2040$richtol>6),
        Pr(scnr1$f2045$richtol>6),
        Pr(scnr1$f2050$richtol>6),
        Pr(scnr1$f2055$richtol>6),
        Pr(scnr1$f2060$richtol>6),
        Pr(scnr1$f2065$richtol>6),
        ylim=c(0,1), names=c("baseline", as.character(seq(2020,2065,5))),
        ylab="inmcm4")
mtext(text="sq",side=3,outer=F)
abline(h=0.5)

boxplot(Pr(scnr2$base$richtol>0),
        Pr(scnr2$f2020$richtol>6),
        Pr(scnr2$f2025$richtol>6),
        Pr(scnr2$f2030$richtol>6),
        Pr(scnr2$f2035$richtol>6),
        Pr(scnr2$f2040$richtol>6),
        Pr(scnr2$f2045$richtol>6),
        Pr(scnr2$f2050$richtol>6),
        Pr(scnr2$f2055$richtol>6),
        Pr(scnr2$f2060$richtol>6),
        Pr(scnr2$f2065$richtol>6),
        ylim=c(0,1), names=c("baseline", as.character(seq(2020,2065,5))))
mtext(text="ws",side=3,outer=F)
abline(h=0.5)

boxplot(Pr(scnr3$base$richtol>6),
        Pr(scnr3$f2020$richtol>6),
        Pr(scnr3$f2025$richtol>6),
        Pr(scnr3$f2030$richtol>6),
        Pr(scnr3$f2035$richtol>6),
        Pr(scnr3$f2040$richtol>6),
        Pr(scnr3$f2045$richtol>6),
        Pr(scnr3$f2050$richtol>6),
        Pr(scnr3$f2055$richtol>6),
        Pr(scnr3$f2060$richtol>6),
        Pr(scnr3$f2065$richtol>6),
        ylim=c(0,1), names=c("baseline", as.character(seq(2020,2065,5))),
        ylab="MRI-CGCM3")
abline(h=0.5)

boxplot(Pr(scnr4$base$richtol>6),
        Pr(scnr4$f2020$richtol>6),
        Pr(scnr4$f2025$richtol>6),
        Pr(scnr4$f2030$richtol>6),
        Pr(scnr4$f2035$richtol>6),
        Pr(scnr4$f2040$richtol>6),
        Pr(scnr4$f2045$richtol>6),
        Pr(scnr4$f2050$richtol>6),
        Pr(scnr4$f2055$richtol>6),
        Pr(scnr4$f2060$richtol>6),
        Pr(scnr4$f2065$richtol>6),
        ylim=c(0,1), names=c("baseline", as.character(seq(2020,2065,5))))
abline(h=0.5)

boxplot(Pr(scnr5$base$richtol>6),
        Pr(scnr5$f2020$richtol>6),
        Pr(scnr5$f2025$richtol>6),
        Pr(scnr5$f2030$richtol>6),
        Pr(scnr5$f2035$richtol>6),
        Pr(scnr5$f2040$richtol>6),
        Pr(scnr5$f2045$richtol>6),
        Pr(scnr5$f2050$richtol>6),
        Pr(scnr5$f2055$richtol>6),
        Pr(scnr5$f2060$richtol>6),
        Pr(scnr5$f2065$richtol>6),
        ylim=c(0,1), names=c("baseline", as.character(seq(2020,2065,5))),
        ylab="NorESM1-M")
abline(h=0.5)

boxplot(Pr(scnr6$base$richtol>6),
        Pr(scnr6$f2020$richtol>6),
        Pr(scnr6$f2025$richtol>6),
        Pr(scnr6$f2030$richtol>6),
        Pr(scnr6$f2035$richtol>6),
        Pr(scnr6$f2040$richtol>6),
        Pr(scnr6$f2045$richtol>6),
        Pr(scnr6$f2050$richtol>6),
        Pr(scnr6$f2055$richtol>6),
        Pr(scnr6$f2060$richtol>6),
        Pr(scnr6$f2065$richtol>6),
        ylim=c(0,1), names=c("baseline", as.character(seq(2020,2065,5))))
abline(h=0.5)
mtext(text="Pr. Bio. Cond. < Good-Fair", las=0,
      side=2, outer=T, line=0)

mtext(text="Year", las=0,
      side=1, outer=T, line=0)
dev.off()

tikz(file=paste(plotDIR,"pred_ProbGF6.tex",sep="/"),width=6,height=2,
     standAlone=F)
par(mfrow=c(1,3), oma=c(1,1,1,1), mar=c(2,2.5,1,0.5),
    mgp=c(1.3,0.125,0),las=1, tck=0.01)
boxplot(Pr(scnr1$base$richtol>6),
        Pr(scnr1$f2020$richtol>6),
        Pr(scnr1$f2025$richtol>6),
        Pr(scnr1$f2030$richtol>6),
        Pr(scnr1$f2035$richtol>6),
        Pr(scnr1$f2040$richtol>6),
        Pr(scnr1$f2045$richtol>6),
        Pr(scnr1$f2050$richtol>6),
        Pr(scnr1$f2055$richtol>6),
        Pr(scnr1$f2060$richtol>6),
        Pr(scnr1$f2065$richtol>6),
        ylim=c(0,1), names=c("baseline", as.character(seq(2020,2065,5))),
        ylab="Pr. BC $<$ Good-Fair")
mtext(text="inmcm4",side=3,outer=F)
abline(h=0.5)
boxplot(Pr(scnr3$base$richtol>6),
        Pr(scnr3$f2020$richtol>6),
        Pr(scnr3$f2025$richtol>6),
        Pr(scnr3$f2030$richtol>6),
        Pr(scnr3$f2035$richtol>6),
        Pr(scnr3$f2040$richtol>6),
        Pr(scnr3$f2045$richtol>6),
        Pr(scnr3$f2050$richtol>6),
        Pr(scnr3$f2055$richtol>6),
        Pr(scnr3$f2060$richtol>6),
        Pr(scnr3$f2065$richtol>6),
        ylim=c(0,1), names=c("baseline", as.character(seq(2020,2065,5))),
        ylab="")
abline(h=0.5)
mtext(text="MRI-CGCM3",side=3,outer=F)

boxplot(Pr(scnr5$base$richtol>6),
        Pr(scnr5$f2020$richtol>6),
        Pr(scnr5$f2025$richtol>6),
        Pr(scnr5$f2030$richtol>6),
        Pr(scnr5$f2035$richtol>6),
        Pr(scnr5$f2040$richtol>6),
        Pr(scnr5$f2045$richtol>6),
        Pr(scnr5$f2050$richtol>6),
        Pr(scnr5$f2055$richtol>6),
        Pr(scnr5$f2060$richtol>6),
        Pr(scnr5$f2065$richtol>6),
        ylim=c(0,1), names=c("baseline", as.character(seq(2020,2065,5))),
        ylab="")
abline(h=0.5)
mtext(text="NorESM1-M",side=3,outer=F)
mtext(text="Year", las=0,
      side=1, outer=T, line=0)
dev.off()




### Data by scenarios: predicted RichTOL, probability > baseline, prob > 6.5

outDIR <- paste(base, "output", sep="/")
### Plots: boxplots of predicted precip and temp and %Natural
###        boxplots of RichTOL, prob > 6.5

## Site Lat & Long
cf_sites <- read.csv(paste(dataDIR, "CF_sites.csv", sep="/"))
pd_sites <- read.csv(paste(dataDIR, "PD_sites.csv", sep="/"))
Sites_LatLong <- rbind(cf_sites[,c("EcositeID","HRU","Lat","Long",
                                   "Level.3.Ecoregion")],
                       pd_sites[,c("EcositeID","HRU","Lat","Long",
                                   "Level.3.Ecoregion")])

## Output data files
## 1. probability > 2020
output_prob <- function(output_file=scnr6,
                        LatLong=Sites_LatLong){
    prob2020 <- data.frame(
        sc_20=Pr(output_file$f2020$richtol>output_file$base$richtol),
        sc_25=Pr(output_file$f2025$richtol>output_file$base$richtol),
        sc_30=Pr(output_file$f2030$richtol>output_file$base$richtol),
        sc_35=Pr(output_file$f2035$richtol>output_file$base$richtol),
        sc_40=Pr(output_file$f2040$richtol>output_file$base$richtol),
        sc_45=Pr(output_file$f2045$richtol>output_file$base$richtol),
        sc_50=Pr(output_file$f2050$richtol>output_file$base$richtol),
        sc_55=Pr(output_file$f2055$richtol>output_file$base$richtol),
        sc_60=Pr(output_file$f2060$richtol>output_file$base$richtol),
        sc_65=Pr(output_file$f2065$richtol>output_file$base$richtol))
    prob2020$EcositeID <- output_file$base$EcositeID
    return(merge(prob2020, LatLong, by="EcositeID"))
}

## scenario 1
tmp <- output_prob(output_file=scnr1)
write.csv(tmp, file=paste(outDIR, "sc1_prob_base.csv", sep="/"))
## scenario 2
tmp <- output_prob(output_file=scnr2)
write.csv(tmp, file=paste(outDIR, "sc2_prob_base.csv", sep="/"))
## scenario 3
tmp <- output_prob(output_file=scnr3)
write.csv(tmp, file=paste(outDIR, "sc3_prob_base.csv", sep="/"))
## scenario 4
tmp <- output_prob(output_file=scnr4)
write.csv(tmp, file=paste(outDIR, "sc4_prob_base.csv", sep="/"))
## scenario 5
tmp <- output_prob(output_file=scnr5)
write.csv(tmp, file=paste(outDIR, "sc5_prob_base.csv", sep="/"))
## scenario 6
tmp <- output_prob(output_file=scnr6)
write.csv(tmp, file=paste(outDIR, "sc6_prob_base.csv", sep="/"))

## 2. Mean RichTOL

output_mean <- function(output_file=scnr6,
                        LatLong=Sites_LatLong){
    mean2020 <- data.frame(
        sc_base=summary(output_file$base$richtol)$mean,
        sc_20=summary(output_file$f2020$richtol)$mean,
        sc_25=summary(output_file$f2025$richtol)$mean,
        sc_30=summary(output_file$f2030$richtol)$mean,
        sc_35=summary(output_file$f2035$richtol)$mean,
        sc_40=summary(output_file$f2040$richtol)$mean,
        sc_45=summary(output_file$f2045$richtol)$mean,
        sc_50=summary(output_file$f2050$richtol)$mean,
        sc_55=summary(output_file$f2055$richtol)$mean,
        sc_60=summary(output_file$f2060$richtol)$mean,
        sc_65=summary(output_file$f2065$richtol)$mean)
    mean2020$EcositeID <- output_file$base$EcositeID
    return(merge(mean2020, LatLong, by="EcositeID"))
}

## scenario 1
tmp <- output_mean(output_file=scnr1)
write.csv(tmp, file=paste(outDIR, "sc1_mean.csv", sep="/"))
## scenario 2
tmp <- output_mean(output_file=scnr2)
write.csv(tmp, file=paste(outDIR, "sc2_mean.csv", sep="/"))
## scenario 3
tmp <- output_mean(output_file=scnr3)
write.csv(tmp, file=paste(outDIR, "sc3_mean.csv", sep="/"))
## scenario 4
tmp <- output_mean(output_file=scnr4)
write.csv(tmp, file=paste(outDIR, "sc4_mean.csv", sep="/"))
## scenario 5
tmp <- output_mean(output_file=scnr5)
write.csv(tmp, file=paste(outDIR, "sc5_mean.csv", sep="/"))
## scenario 6
tmp <- output_mean(output_file=scnr6)
write.csv(tmp, file=paste(outDIR, "sc6_mean.csv", sep="/"))

## 3. probability > 6 (worse than good-fair)
output_prob6 <- function(output_file=scnr6,
                        LatLong=Sites_LatLong){
    probGR6 <- data.frame(
        sc_base=Pr(output_file$base$richtol>6),
        sc_20=Pr(output_file$f2020$richtol>6),
        sc_25=Pr(output_file$f2025$richtol>6),
        sc_30=Pr(output_file$f2030$richtol>6),
        sc_35=Pr(output_file$f2035$richtol>6),
        sc_40=Pr(output_file$f2040$richtol>6),
        sc_45=Pr(output_file$f2045$richtol>6),
        sc_50=Pr(output_file$f2050$richtol>6),
        sc_55=Pr(output_file$f2055$richtol>6),
        sc_60=Pr(output_file$f2060$richtol>6),
        sc_65=Pr(output_file$f2065$richtol>6))
    probGR6$EcositeID <- output_file$base$EcositeID
    return(merge(probGR6, LatLong, by="EcositeID"))
}

## scenario 1
tmp <- output_prob6(output_file=scnr1)
write.csv(tmp, file=paste(outDIR, "sc1_probGF.csv", sep="/"))
probGF1 <- tmp
## scenario 2
tmp <- output_prob6(output_file=scnr2)
write.csv(tmp, file=paste(outDIR, "sc2_probGF.csv", sep="/"))
## scenario 3
tmp <- output_prob6(output_file=scnr3)
write.csv(tmp, file=paste(outDIR, "sc3_probGF.csv", sep="/"))
probGF3 <- tmp
## scenario 4
tmp <- output_prob6(output_file=scnr4)
write.csv(tmp, file=paste(outDIR, "sc4_probGF.csv", sep="/"))
## scenario 5
tmp <- output_prob6(output_file=scnr5)
write.csv(tmp, file=paste(outDIR, "sc5_probGF.csv", sep="/"))
probGF5 <- tmp
## scenario 6
tmp <- output_prob6(output_file=scnr6)
write.csv(tmp, file=paste(outDIR, "sc6_probGF.csv", sep="/"))

## mapping sc 1, 3, 5
##scenario 1
probGF1_2025 <- cut(probGF1$sc_25, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF1_2030 <- cut(probGF1$sc_30, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF1_2035 <- cut(probGF1$sc_35, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF1_2040 <- cut(probGF1$sc_40, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF1_2045 <- cut(probGF1$sc_45, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF1_2050 <- cut(probGF1$sc_50, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF1_2055 <- cut(probGF1$sc_55, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF1_2060 <- cut(probGF1$sc_60, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF1_2065 <- cut(probGF1$sc_65, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF1_base <- cut(probGF1$sc_base, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF1_sp <- data.frame(Lat=probGF1$Lat, Long=probGF1$Long, base=probGF1_base,
                         "2045"=probGF1_2030)

##scenario 3
probGF3_2025 <- cut(probGF3$sc_25, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF3_2030 <- cut(probGF3$sc_30, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF3_2035 <- cut(probGF3$sc_35, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF3_2040 <- cut(probGF3$sc_40, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF3_2045 <- cut(probGF3$sc_45, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF3_2050 <- cut(probGF3$sc_50, breaks=c(-0.00001,0.25,0.5,0.75,1),
                   labels=c("a","b","c","d"))
probGF3_2055 <- cut(probGF3$sc_55, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF3_2060 <- cut(probGF3$sc_60, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF3_2065 <- cut(probGF3$sc_65, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF3_base <- cut(probGF3$sc_base, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF3_sp <- data.frame(Lat=probGF3$Lat, Long=probGF3$Long, base=probGF3_base,
                         "2045"=probGF3_2030)
## scenario 5
probGF5_2025 <- cut(probGF5$sc_25, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF5_2030 <- cut(probGF5$sc_30, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF5_2035 <- cut(probGF5$sc_35, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF5_2040 <- cut(probGF5$sc_40, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF5_2045 <- cut(probGF5$sc_45, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF5_2050 <- cut(probGF5$sc_50, breaks=c(-0.00001,0.25,0.5,0.75,1),
                   labels=c("a","b","c","d"))
probGF5_2055 <- cut(probGF5$sc_55, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF5_2060 <- cut(probGF5$sc_60, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF5_2065 <- cut(probGF5$sc_65, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF5_base <- cut(probGF5$sc_base, breaks=c(-0.00001,0.25,0.5,0.75,1),
                    labels=c("a","b","c","d"))
probGF5_sp <- data.frame(Lat=probGF5$Lat, Long=probGF5$Long, base=probGF5_base,
                         "2045"=probGF5_2045)
## summary table


GF5_tbl <- as.data.frame(
    cbind(apply(table(probGF5_base, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF5_2025, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF5_2030, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF5_2035, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF5_2040, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF5_2045, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF5_2050, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF5_2055, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF5_2060, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF5_2065, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,]))
names(GF5_tbl) <- c("base", seq(2025,2065,5))

GF3_tbl <- as.data.frame(
    cbind(apply(table(probGF3_base, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF3_2025, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF3_2030, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF3_2035, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF3_2040, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF3_2045, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF3_2050, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF3_2055, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF3_2060, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF3_2065, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,]))
names(GF3_tbl) <- c("base", seq(2025,2065,5))

GF1_tbl <- as.data.frame(
    cbind(apply(table(probGF1_base, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF1_2025, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF1_2030, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF1_2035, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF1_2040, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF1_2045, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF1_2050, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF1_2055, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF1_2060, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,],
          apply(table(probGF1_2065, probGF5$Level.3.Ecoregion), 2,
                function(x)x/sum(x))[4,]))

names(GF1_tbl) <- c("base", seq(2025,2065,5))

write.table(round(rbind(GF1_tbl, GF3_tbl, GF5_tbl), 3),
            file="worstcategory.txt", sep="&")

cols <- c("blue", "cyan", "pink",  "red")
##pdf(file=paste(plotDIR, "pred_maps.pdf", sep="/"),
##    height=6, width=4)
##png(file=paste(plotDIR, "pred_maps.png", sep="/"),
##    height=6*120, width=4*120)
tikz(file=paste(plotDIR, "pred_maps.tex", sep="/"),
    height=6, width=4, standAlone=T)
par(mfrow=c(3,2), oma=c(2,2,2,0), mar=c(0.,0.,0.,0.))
map("state",region=c("north carolina","south carolina"))
points(probGF1_sp$Long, probGF1_sp$Lat, col=cols[as.numeric(probGF1_sp$base)],
       cex=0.35, pch=16)
mtext("Baseline", side=3, outer=F, line=1.5)
mtext("2014-2019", side=3, outer=F, line=0.)
mtext("inmcm4", side=2, outer=F, line=1.5, las=0)

map("state",region=c("north carolina","south carolina"))
points(probGF1_sp$Long, probGF1_sp$Lat, col=cols[as.numeric(probGF1_sp$X2045)],
       cex=0.35, pch=16)
mtext("Worst case", side=3, outer=F, line=1.5)
mtext("2035", side=3, outer=F, line=0)

map("state",region=c("north carolina","south carolina"))
points(probGF3_sp$Long, probGF3_sp$Lat, col=cols[as.numeric(probGF3_sp$base)],
       cex=0.35, pch=16)
mtext("MRI-CGCM3", side=2, outer=F, line=1.5, las=0)

map("state",region=c("north carolina","south carolina"))
points(probGF3_sp$Long, probGF3_sp$Lat, col=cols[as.numeric(probGF3_sp$X2045)],
       cex=0.35, pch=16)
mtext("2030", side=3, outer=F, line=0)

map("state",region=c("north carolina","south carolina"))
points(probGF5_sp$Long, probGF5_sp$Lat, col=cols[as.numeric(probGF5_sp$base)],
       cex=0.35, pch=16)
mtext("NorESM1-M", side=2, outer=F, line=1.5, las=0)

map("state",region=c("north carolina","south carolina"))
points(probGF5_sp$Long, probGF5_sp$Lat, col=cols[as.numeric(probGF5_sp$X2045)],
       cex=0.35, pch=16)
mtext("2045", side=3, outer=F, line=0)
dev.off()

## Graphic abstract
tikz(file=paste(plotDIR, "GrAb_pred_maps.tex", sep="/"),
    height=3, width=4, standAlone=T)
par(mfrow=c(1,2), oma=c(2,2,2,0), mar=c(0.,0.,0.,0.))
map("state",region=c("north carolina","south carolina"))
points(probGF1_sp$Long, probGF1_sp$Lat, col=cols[as.numeric(probGF1_sp$base)],
       cex=0.35, pch=16)
mtext("Baseline", side=3, outer=F, line=1.5)
mtext("2014-2019", side=3, outer=F, line=0.)
mtext("inmcm4", side=2, outer=F, line=1.5, las=0)

map("state",region=c("north carolina","south carolina"))
points(probGF1_sp$Long, probGF1_sp$Lat, col=cols[as.numeric(probGF1_sp$X2045)],
       cex=0.35, pch=16)
mtext("Worst case", side=3, outer=F, line=1.5)
mtext("2035", side=3, outer=F, line=0)
dev.off()
