plw <- 600
plh <- 600

# read data
# on my machine
baselineFolder <- file.path('workOutput','determineAggDam_Baseline-S20000-policy_EMB-ClimateFeedback_On-ClimateSTAOverride_Off',
														'detectedParmSpace','PerVarFiles-RDS')
# on levante
# baselineFolder <- file.path('/work/mh0033/b383346/WorldTransFrida-Uncertainty/workOutput/UA_EMBv6Try2_nS100000')

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

rm(list=c('gdp','pdppc','pop','edem','esup','esht','inf','prinv','puinv','sta'))

# consider drop early years as burn in
# regDF <- regDF[regDF$year>=2050,]

nobs <- nrow(regDF)

# model ####
# par(pch='.')
# plot(regDF$year[regDF$id==2],regDF$rgdp_within[regDF$id==2],type='l')
# plot(regDF$sta,regDF$rgdp_within,ylim=c(-0.1,0.1))
## OLS ####
degree <- 5
# regular model
lMod <- lm(gdp ~ poly(lgdp,sta,degree=degree), data=regDF)
# lMod <- lm(gdp ~ poly(lgdp,sta,degree=degree) + poly(lgdp,pop,esht,inf,prinv,puinv,degree=1), data=regDF)

# 2stage model
# popMod <- lm(pop ~ poly(sta,degree=degree), data=regDF)
# regDF$popNoSta <- resid(popMod)
# eshtMod <- lm(esht ~ poly(sta,degree=degree), data=regDF)
# regDF$eshtNoSta <- resid(eshtMod)
# infMod <- lm(inf ~ poly(sta,degree=degree), data=regDF)
# regDF$infNoSta <- resid(eshtMod)
# prinvMod <- lm(prinv ~ poly(sta,degree=degree), data=regDF)
# regDF$prinvNoSta <- resid(prinvMod)
# puinvMod <- lm(puinv ~ poly(sta,degree=degree), data=regDF)
# regDF$puinvNoSta <- resid(puinvMod)
# predDF <- regDF
# predDF$sta <- 0 
# lMod <- lm(gdp ~ poly(lgdp,sta,degree=degree) + poly(lgdp,popNoSta,eshtNoSta,infNoSta,prinvNoSta,puinvNoSta,degree=1), data=regDF)

summary(lMod)
png(file.path('figures','residVsLgdp.png'),width = plw,height = plh)
par(pch='.')
plot(regDF$lgdp, resid(lMod),
		 xlim=c(0,5e6),ylim=c(-2e5,2e5),
		 main='prediction error')
abline(h=0,col='red')
dev.off()
png(file.path('figures','residVsSta.png'),width = plw,height = plh)
par(pch='.')
plot(regDF$sta, resid(lMod),
		 ylim=c(-2e5,2e5),
		 main='prediction error')
abline(h=0,col='red')
dev.off()
lModPred <- predict(lMod)
predDF <- regDF
predDF$sta <- 0
lModPred0sta <- predict(lMod, newdata = predDF)
png(file.path('figures','gdpVsLgdp.png'),width = plw,height = plh)
par(pch='.')
plot(regDF$lgdp, regDF$gdp, xlim=c(0,2e6), ylim=c(0,2e6),
		 main='gdp depending on gdp_t-1')
points(regDF$lgdp, lModPred0sta,col='red')
legend('topleft',
			 legend=c('FRIDA EMB','Pred. STA=0'),
			 col=c('black','red'),
			 pch=20)
dev.off()
lModLoss <- lModPred0sta - regDF$gdp
lModLossRel <- lModLoss/regDF$gdp
png(file.path('figures','relGDPloss.png'),width = plw,height = plh)
par(pch='.')
plot(regDF$sta,lModLossRel,
		 xlim=c(0,8),ylim=c(-0.1,0.8),
		 col=adjustcolor(1,alpha.f = 0.01),
		 # col=1,
		 main='year relative GDP loss')
dev.off()

# model2 ####
library(optimx)
predFitM <- function(par,sta){
	return(par[1] + par[2]*log(sta+par[3])+par[4]*sta)
}
fitRMSE <- function(par,sta){
	return(sqrt(mean((lModLossRel - predFitM(par,sta))^2)))
}
par0 <- c(0,0.04,0.5,0)
fittedPar <- optimx(par0,fitRMSE,method = 'BFGS',sta=regDF$sta)
par <- unlist(unname(fittedPar[which.min(fittedPar$value),1:length(par0)]))
staSup <- seq(-par[1],10,length.out=200)
png(file.path('figures','relGDPlossFunFit.png'),width = plw,height = plh)
par(pch='.')
plot(staSup,predFitM(par,staSup),col='red',type='l',ylim=c(-0.1,0.8),lwd=3)
dev.off()

# FE model
# discarded as errors are larger than in the corresponding standard model
# this is due to the fe only accounting for differences in the level between groups
# not for differences in trends between groups. In fact removing the average
# levels for each groups means that their origins no longer conincide for cases where
# data have strong trends, such as this.
# 
