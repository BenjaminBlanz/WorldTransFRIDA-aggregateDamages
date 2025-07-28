# read data
baselineFolder <- file.path('workOutput','determineAggDam_Baseline-S20000-policy_EMB-ClimateFeedback_On-ClimateSTAOverride_Off',
														'detectedParmSpace','PerVarFiles-RDS')
sta <- readRDS(file.path(baselineFolder,'energy_balance_model_surface_temperature_anomaly.RDS'))
sta <- sta[complete.cases(sta),]
gdp <- readRDS(file.path(baselineFolder,'gdp_real_gdp_in_2021c.RDS'))
gdp <- gdp[gdp$id %in% sta$id,]
gdppc <- readRDS(file.path(baselineFolder,'demographics_real_gdp_per_person.RDS'))
gdppc <- gdppc[gdppc$id %in% sta$id,]
nid <- nrow(sta)
nyear <- ncol(sta)-1
years <- as.numeric(colnames(sta)[-1])

# calc gro ###
gdpGro <- gdp
gdppcGro <- gdppc
gdpGro[,as.character(years[nyear])] <- NA
gdppcGro[,as.character(years[nyear])] <- NA
for(year in years[-nyear]){
	gdpGro[,as.character(year)] <- gdp[,as.character(year+1)] - gdp[,as.character(year)]
	gdppcGro[,as.character(year)] <- gdppc[,as.character(year+1)] - gdppc[,as.character(year)]
}

# restruction data for regresssion ####
regDFColnames <- c('id','year','sta','gdp','gdppc','gdpGro','gdpGroRt','gdppcGro','gdppcGroRt')
regDF <- matrix(NA,nrow=nid*nyear,ncol=length(regDFColnames))
colnames(regDF) <- regDFColnames
regDF[,'id'] <- rep(1:nid,nyear)
regDF[,'year'] <- rep(years,each=nid)
regDF[,'sta'] <- unname(unlist(sta[,-1]))
regDF[,'gdp'] <- unname(unlist(gdp[,-1]))
regDF[,'gdpGro'] <- unname(unlist(gdpGro[,-1]))
regDF[,'gdpGroRt'] <- regDF[,'gdpGro']/regDF[,'gdp']
regDF[,'gdppc'] <- unname(unlist(gdppc[,-1]))
regDF[,'gdppcGro'] <- unname(unlist(gdppcGro[,-1]))
regDF[,'gdppcGroRt'] <- regDF[,'gdppcGro']/regDF[,'gdppc']
regDF <- as.data.frame(regDF)
regDF <- regDF[complete.cases(regDF),]
nobs <- nrow(regDF)
# plot scatter ####
selection <- which(regDF$id==2) 

plot(regDF$year[selection],regDF$gdp[selection])
plot(regDF$year[selection],regDF$gdpGro[selection])
plot(regDF$year[selection],regDF$gdpGroRt[selection])

plot(regDF$gdp,regDF$gdpGro,pch='.')
plot(regDF$gdp,regDF$gdpGroRt,pch='.',ylim=c(-0.1,0.1), xlim=c(0,2e6))
plot(regDF$sta,regDF$gdpGro,pch='.',ylim=c(-2e4,4e4))
plot(regDF$sta,regDF$gdpGroRt,pch='.',ylim=c(-0.1,0.1))
plot(regDF$sta,regDF$gdp,pch='.',ylim=c(0,1e6))

subset <- which(regDF$gdp < 1e6 & abs(regDF$gdpGroRt)<0.1)
library(rgl)
plot3d(regDF$gdp[subset],regDF$sta[subset],regDF$gdpGroRt[subset],cex=0.5)

# model ####
lMod <- lm(gdpGro ~ -1 + gdp + sta + I(sta^2), data=regDF) 
summary(lMod)

lModGrRt <- lm(gdpGroRt ~ gdp + I(gdp^2) + I(gdp*sta) + sta + I(sta^2),data=regDF)
summary(lModGrRt)

library(plm)
feMod <- plm(gdpGro ~ gdp + I(gdp*sta) + sta + I(sta^2) + I(sta^3), index=('id'), model='within',data=regDF)
summary(feMod)
reMod <- plm(gdpGro ~ gdp + I(gdp*sta) + sta + I(sta^2) + I(sta^3), index=('id'), model='random',data=regDF)
summary(reMod)
phtest(feMod,reMod)

feModGrRt <- plm(gdpGroRt ~ gdp + I(gdp*sta) + sta + I(sta^2), index=('id'), model='within',data=regDF)
summary(feModGrRt)
reModGrRt <- plm(gdpGroRt ~ gdp + I(gdp*sta) + sta + I(sta^2), index=('id'), model='random',data=regDF)
summary(reModGrRt)


# pred ####
gdpSup <- seq(0,1e6,length.out=20)
staSup <- seq(0,7,0.5)
support <- data.frame(gdp=rep(gdpSup,length(staSup)),sta=rep(staSup,each=length(gdpSup)))
fePredGrRt <- predict(feModGrRt,newdata = support)
zvals <- matrix(fePredGrRt,nrow=length(gdpSup),ncol=length(staSup))

plot3d(regDF$gdp[subset],regDF$sta[subset],regDF$gdpGroRt[subset],cex=0.5)
surface3d(gdpSup,staSup,zvals,col='red')

# pred 0 deg gro ####
predDF <- data.frame(gdp=regDF$gdp, sta=rep(0,nobs), id=regDF$id)
regDF$gdpGroPred <- predict(lMod)
regDF$gdpGro0Deg <- predict(lMod,newdata=predDF)

predpDF <- pdata.frame(predDF, index='id')
regDF$gdpGroPredFE <- predict(feMod)
regDF$gdpGro0DegFE <- predict(feMod,newdata = predpDF)

staSup <- seq(0,7,0.1)

# plot GrRt ###
my1DPolyFun <- function(x,coefs){
	y = 0
	for(i in 1:length(coefs)){
		y = y + coefs[i] * x^(i-1)
	}
	return(y)
}
plot(regDF$sta,regDF$gdpGroRt,pch='.',ylim=c(-0.1,0.1))
points(staSup, my1DPolyFun(staSup, c( mean(fixef(feModGrRt)),feModGrRt$coefficients)),col='red',pch=20)

predGro <- my1DPolyFun(staSup, c( mean(fixef(feModGrRt)),feModGrRt$coefficients))
DamFac <- 1-predGro/predGro[1]

# plot ####
library(rgl)
randomSelectionForPlotting <- sample(plotConstrainedIdc,size = 10000)

open3d()
plotLim <- 4e4
plotConstrainedIdc <- which(abs(regDF$gdpGro)<=plotLim)
plot3d(regDF$gdp[plotConstrainedIdc],regDF$sta[plotConstrainedIdc],regDF$gdpGro[plotConstrainedIdc])
plotConstrainedIdc <- which(abs(regDF$gdpGroPred)<=plotLim)
points3d(regDF$gdp[plotConstrainedIdc],regDF$sta[plotConstrainedIdc],regDF$gdpGroPred[plotConstrainedIdc],col='red')
plotConstrainedIdc <- which(abs(regDF$gdpGroPredFE)<=plotLim & abs(regDF$gdp) <= 4e6)
points3d(regDF$gdp[plotConstrainedIdc],regDF$sta[plotConstrainedIdc],regDF$gdpGroPredFE[plotConstrainedIdc],col='green')

open3d()
plot3d(regDF$gdp,regDF$sta,regDF$gdpGroRt,zlim=c(-0.1,0.1),xlim=c(0,4e6))

plot(regDF$sta[randomSelectionForPlotting],regDF$gdpGro[randomSelectionForPlotting],pch='.',ylim=c(-2e4,4e4))
points(regDF$sta[randomSelectionForPlotting],regDF$gdpGroPred[randomSelectionForPlotting],pch='.',col='red')
points(regDF$sta,regDF$gdpGro0Deg,pch='.',col='red')

