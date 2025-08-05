my1DPolyFun <- function(x,coefs){
	y = 0
	for(i in 1:length(coefs)){
		y = y + coefs[i] * x^(i-1)
	}
	return(y)
}

# read data
baselineFolder <- file.path('workOutput','determineAggDam_Baseline-S20000-policy_EMB-ClimateFeedback_On-ClimateSTAOverride_Off',
														'detectedParmSpace','PerVarFiles-RDS')
sta <- readRDS(file.path(baselineFolder,'energy_balance_model_surface_temperature_anomaly.RDS'))
# sta <- sta[complete.cases(sta),] # drop incomplete runs
gdp <- readRDS(file.path(baselineFolder,'gdp_real_gdp_in_2021c.RDS'))
gdp <- gdp[gdp$id %in% sta$id,]
gdppc <- readRDS(file.path(baselineFolder,'demographics_real_gdp_per_person.RDS'))
gdppc <- gdppc[gdppc$id %in% sta$id,]
pop <- readRDS(file.path(baselineFolder,'demographics_population.RDS'))
pop <- pop[pop$id %in% sta$id,]
edem <- readRDS(file.path(baselineFolder,'energy_demand_demand_for_energy.RDS'))
edem <- edem[edem$id %in% sta$id,]
esup <- readRDS(file.path(baselineFolder,'energy_supply_total_energy_output.RDS'))
esup <- esup[esup$id %in% sta$id,]
esht <- edem - esup # energy shortage
inf <- readRDS(file.path(baselineFolder,'inflation_inflation_rate.RDS'))
inf <- inf[inf$id %in% sta$id,]
prinv <- readRDS(file.path(baselineFolder,'gdp_private_investment_in_in_2021c.RDS'))
prinv <- prinv[prinv$id %in% sta$id,]
puinv <- readRDS(file.path(baselineFolder,'gdp_public_investment_in_2021c.RDS'))
puinv <- puinv[puinv$id %in% sta$id,]

nid <- nrow(sta)
years <- as.numeric(colnames(sta)[-1])
nyear <- length(years)

# calc gro ####
gdpGro <- gdp
gdppcGro <- gdppc
gdpGro[,as.character(years[nyear])] <- NA
gdppcGro[,as.character(years[nyear])] <- NA
for(year in years[-nyear]){
	gdpGro[,as.character(year)] <- gdp[,as.character(year+1)] - gdp[,as.character(year)]
	gdppcGro[,as.character(year)] <- gdppc[,as.character(year+1)] - gdppc[,as.character(year)]
}

# calc relGDP ####
# gdp relative to gdp in the preceedig period ####
rgdp <- gdp
rgdp[,as.character(years[1])] <- NA
for(year in years[-nyear]){
	rgdp[,as.character(year+1)] <- gdp[,as.character(year+1)]/gdp[,as.character(year)]
}
# lagged GDP ####
lgdp <- gdp
lgdp[,as.character(years[1])] <- NA
for(year in years[-nyear]){
	lgdp[,as.character(year+1)] <- gdp[,as.character(year)]
}

# restructure data for regresssion ####
regDFColnames <- c('id','year','sta','gdp','lgdp','rgdp','gdppc','gdpGro','gdpGroRt','gdppcGro','gdppcGroRt',
									 'pop','edem','esup','esht','inf','prinv','puinv')
regDF <- matrix(NA,nrow=nid*nyear,ncol=length(regDFColnames))
colnames(regDF) <- regDFColnames
regDF[,'id'] <- rep(1:nid,nyear)
regDF[,'year'] <- rep(years,each=nid)
regDF[,'sta'] <- unname(unlist(sta[,-1]))
regDF[,'gdp'] <- unname(unlist(gdp[,-1]))
regDF[,'lgdp'] <- unname(unlist(lgdp[,-1]))
regDF[,'rgdp'] <- unname(unlist(rgdp[,-1]))
regDF[,'gdpGro'] <- unname(unlist(gdpGro[,-1]))
regDF[,'gdpGroRt'] <- regDF[,'gdpGro']/regDF[,'gdp']
regDF[,'gdppc'] <- unname(unlist(gdppc[,-1]))
regDF[,'gdppcGro'] <- unname(unlist(gdppcGro[,-1]))
regDF[,'gdppcGroRt'] <- regDF[,'gdppcGro']/regDF[,'gdppc']
regDF[,'pop'] <- unname(unlist(pop[,-1]))
regDF[,'edem'] <- unname(unlist(edem[,-1]))
regDF[,'esup'] <- unname(unlist(esup[,-1]))
regDF[,'esht'] <- unname(unlist(esht[,-1]))
regDF[,'inf'] <- unname(unlist(inf[,-1]))
regDF[,'prinv'] <- unname(unlist(prinv[,-1]))
regDF[,'puinv'] <- unname(unlist(puinv[,-1]))
regDF <- as.data.frame(regDF)
regDF <- regDF[complete.cases(regDF),]

nobs <- nrow(regDF)


library('datawizard')
regDF <- demean(regDF,c('gdp','lgdp','rgdp','gdpGro','gdpGroRt'),'id')

# model ####
# par(pch='.')
# plot(regDF$year[regDF$id==2],regDF$rgdp_within[regDF$id==2],type='l')
# plot(regDF$sta,regDF$rgdp_within,ylim=c(-0.1,0.1))
predDF <- regDF
predDF$sta <- 0
## OLS ####
# regular model
lMod <- lm(gdp ~ poly(lgdp,sta,degree=5), data=regDF)
summary(lMod)
par(pch='.')
plot(regDF$lgdp, resid(lMod),
		 xlim=c(0,5e6),ylim=c(-2e5,2e5),
		 main='prediction error')
abline(h=0,col='red')
plot(regDF$sta, resid(lMod),
		 ylim=c(-2e5,2e5),
		 main='prediction error')
abline(h=0,col='red')
lModPred <- predict(lMod)
lModPred0sta <- predict(lMod, newdata = predDF)
plot(regDF$lgdp, regDF$gdp, xlim=c(0,2e6), ylim=c(0,2e6))
points(regDF$lgdp, lModPred0sta,col='red')
lModLoss <- lModPred0sta - lModPred#regDF$gdp
lModLossRel <- lModLoss/lModPred#regDF$gdp
plot(regDF$sta,lModLossRel,ylim=c(-0.1,0.8),
		 main='year relative GDP loss')


# FE model
# discarded as errors are larger than in the corresponding standard model
# this is due to the fe only accounting for differences in the level between groups
# not for differences in trends between groups. In fact removing the average
# levels for each groups means that their origins no longer conincide for cases where
# data have strong trends, such as this.
# 
# lMod <- lm(gdp_within ~ poly(lgdp_within,sta,degree=5), data=regDF) 
# summary(lMod)
# par(pch='.')
# plot(regDF$sta,regDF$gdp_within,ylim=c(-1e6,2e6))
# points(regDF$sta,predict(lMod),col='red')
# plot(regDF$lgdp_within,regDF$gdp_within - predict(lMod),
# 		 xlim=c(-1e6,1e6),ylim=c(-1e5,1e5),
# 		 main='prediction error')
# 
# lModPred <- predict(lMod)
# lModPred0sta <- predict(lMod, newdata = predDF)
# lModLoss <- lModPred0sta - regDF$gdp_within
# lModLossRel <- lModLoss/regDF$gdp
# plot(regDF$sta,lModLoss,ylim=c(0,1e5),
# 		 main='GDP loss')
# plot(regDF$sta,lModLossRel,ylim=c(-0.1,0.8),
# 		 main='Relative GDP loss')

## OLS grRt ####
lModGrRt <- lm(gdpGroRt_within ~ gdp_within + sta + I(sta^2) + I(sta^3) + I(sta^4) + I(sta^5), data=regDF) 
summary(lModGrRt)
par(pch='.')
plot(regDF$sta,regDF$gdpGroRt_within,ylim=c(-0.1,0.1))
points(regDF$sta,predict(lModGrRt),col='red')


# 
# lModGrRt <- lm(gdpGroRt ~ gdp + I(gdp^2) + I(gdp*sta) + sta + I(sta^2),data=regDF)
# summary(lModGrRt)
# 
library(plm)
# feMod <- plm(gdpGro ~ gdp + I(gdp*sta) + sta + I(sta^2) + I(sta^3), index=('id'), model='within',data=regDF)
# summary(feMod)
# reMod <- plm(gdpGro ~ gdp + I(gdp*sta) + sta + I(sta^2) + I(sta^3), index=('id'), model='random',data=regDF)
# summary(reMod)
# phtest(feMod,reMod)
# 
# feModGrRt <- plm(gdpGroRt ~ gdp + I(gdp*sta) + sta + I(sta^2), index=('id'), model='within',data=regDF)
feModGrRt <- lm(gdpGroRt ~ gdp + I(gdp*sta) + sta + I(sta^2), data=regDF)
summary(feModGrRt)
feModGrRtPred <- predict(feModGrRt)
plot(regDF$gdp,regDF$gdpGroRt,pch='.', xlim=c(0,2e6), ylim=c(-0.1,0.1))
points(regDF$gdp,feModGrRtPred,pch='.',col='red')
plot(regDF$sta,regDF$gdpGroRt,pch='.',ylim=c(-0.1,0.1))
points(regDF$sta,feModGrRtPred,pch='.',col='red')

# reModGrRt <- plm(gdpGroRt ~ gdp + I(gdp*sta) + sta + I(sta^2), index=('id'), model='random',data=regDF)
# summary(reModGrRt)

feMod <- plm(gdp ~ lgdp + sta + I(sta^2), index='id', model='within',data=regDF)
# feMod <- plm(rgdp ~ sta + I(sta^2) + pop + esht, index='id', model='within',data=regDF)
summary(feMod)

# pred ####
predDF <- regDF
feModPred <- predict(feMod,newdata = predDF)
feModGrRtPred <- predict(feModGrRt,newdata = predDF)
predDF$sta <- 0
feModPred0sta <- predict(feMod,newdata = predDF)
feModGrRt0sta <- predict(feModGrRt,newdata = predDF)

# losses ####
gdpLoss <- feModPred0sta - feModPred
relGdpLoss <- gdpLoss/feModPred
plot(regDF$sta,relGdpLoss,pch='.',ylim=c(-0.1,1))
stas <- seq(0,7,1)
for(sta in stas){
	boxplot(relGdpLoss[regDF$sta < sta+0.5 & regDF$sta > sta-0.5],at=sta,add=T)
}

gdpGrRtLoss <- feModGrRt0sta - feModGrRtPred
relgdpGrRtLoss <- gdpGrRtLoss / (1+feModGrRt0sta)
plot(regDF$sta,regDF$gdpGroRt,pch='.',ylim=c(-0.1,0.1))
points(regDF$sta,feModGrRtPred,pch='.',col='blue')
points(regDF$sta,feModGrRt0sta,pch='.',col='red')

plot(regDF$sta,gdpGrRtLoss,pch='.',ylim=c(-0.1,0.1))
plot(regDF$sta,relgdpGrRtLoss,pch='.',ylim=c(-0.1,0.1))

# plot ####
if(F){
	par(pch='.')
	plot(regDF$sta,regDF$gdp,ylim=c(0,2e6))
	points(regDF$sta,predict(feMod),col='red')
	plot(regDF$lgdp,regDF$gdp - predict(feMod),xlim=c(0,1e6),ylim=c(-1e6,1e6))
	
	zvals <- matrix(fePredGrRt,nrow=length(gdpSup),ncol=length(staSup))
	plot3d(regDF$gdp[subset],regDF$sta[subset],regDF$gdpGroRt[subset],cex=0.5)
	surface3d(gdpSup,staSup,zvals,col='red')
	
	# plot GrRt 

	plot(regDF$sta,regDF$gdpGroRt,pch='.',ylim=c(-0.1,0.1))
	points(staSup, my1DPolyFun(staSup, c( mean(fixef(feModGrRt)),feModGrRt$coefficients)),col='red',pch=20)
	
	predGro <- my1DPolyFun(staSup, c( mean(fixef(feModGrRt)),feModGrRt$coefficients))
	DamFac <- 1-predGro/predGro[1]
	#plot scatter
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
}
