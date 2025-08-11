# config

degree <- 3
covDegree <- 3
twoStage <- FALSE
plw <- 6+0.21
plh <- 6+0.51
plu <- 'in'
pld <- 2000

readData <- TRUE
makePredict <- TRUE
makeDFmod <- TRUE


# read data ####
if(readData){
	cat('reading data...')
	# on my machine
	# baselineFolder <- file.path('workOutput','determineAggDam_Baseline-S20000-policy_EMB-ClimateFeedback_On-ClimateSTAOverride_Off',
	# 														'detectedParmSpace','PerVarFiles-RDS')
	# on levante
	baselineFolder <- file.path('/work/mh0033/b383346/WorldTransFrida-Uncertainty/workOutput/UA_EMBv6Try2_nS100000/detectedParmSpace/PerVarFiles-RDS')
	
	gdp <- readRDS(file.path(baselineFolder,'gdp_real_gdp_in_2021c.RDS'))
	gdp <- gdp[complete.cases(gdp),] # drop incomplete runs
	gdp <- gdp[which(!is.infinite(gdp[,'2150'])),] # drop runs that end in inf
	sta <- readRDS(file.path(baselineFolder,'energy_balance_model_surface_temperature_anomaly.RDS'))
	sta <- sta[sta$id %in% gdp$id,]
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
	cat('done\n')
	
	cat('calculating relative and lagged vars.')
	# calc gro ####
	gdpGro <- gdp
	gdppcGro <- gdppc
	gdpGro[,as.character(years[nyear])] <- NA
	gdppcGro[,as.character(years[nyear])] <- NA
	for(year in years[-nyear]){
		gdpGro[,as.character(year)] <- gdp[,as.character(year+1)] - gdp[,as.character(year)]
		gdppcGro[,as.character(year)] <- gdppc[,as.character(year+1)] - gdppc[,as.character(year)]
	}
	# calc rel1gdp ####
	# gdp relative to gdp in the preceedig period ####
	cat('.')
	rgdp <- gdp
	rgdp[,as.character(years[1])] <- NA
	for(year in years[-nyear]){
		rgdp[,as.character(year+1)] <- gdp[,as.character(year+1)]/gdp[,as.character(year)]
	}
	# lagged GDP ####
	cat('.')
	l1gdp <- gdp
	l1gdp[,as.character(years[1])] <- NA
	for(year in years[-nyear]){
		l1gdp[,as.character(year+1)] <- gdp[,as.character(year)]
	}
	# twice lagged gdp ####
	cat('.')
	l2gdp <- gdp
	l2gdp[,as.character(years[1:2])] <- NA
	for(year in years[-((nyear-1):nyear)]){
		l2gdp[,as.character(year+2)] <- gdp[,as.character(year)]
	}
	# thrice lagged gdp ####
	cat('.')
	l3gdp <- gdp
	l3gdp[,as.character(years[1:3])] <- NA
	for(year in years[-((nyear-2):nyear)]){
		l3gdp[,as.character(year+3)] <- gdp[,as.character(year)]
	}
	# future GDP ####
	cat('.')
	fgdp <- gdp
	fgdp[,as.character(years[nyear])] <- NA
	for(year in years[-nyear]){
		fgdp[,as.character(year)] <- gdp[,as.character(year+1)]
	}
	# lagged STA ####
	cat('.')
	l1sta <- sta
	l1sta[,as.character(years[1])] <- NA
	for(year in years[-nyear]){
		l1sta[,as.character(year+1)] <- sta[,as.character(year)]
	}
	# twice lagged STA ####
	cat('.')
	l2sta <- sta
	l2sta[,as.character(years[1:2])] <- NA
	for(year in years[-((nyear-1):nyear)]){
		l2sta[,as.character(year+2)] <- sta[,as.character(year)]
	}
	# thrice lagged STA ####
	cat('.')
	l3sta <- sta
	l3sta[,as.character(years[1:3])] <- NA
	for(year in years[-((nyear-2):nyear)]){
		l3sta[,as.character(year+3)] <- sta[,as.character(year)]
	}
	# future STA ####
	cat('.')
	fsta <- sta
	fsta[,as.character(years[nyear])] <- NA
	for(year in years[-nyear]){
		fsta[,as.character(year)] <- sta[,as.character(year+1)]
	}
	cat('done\n')
	
	# restructure data for regresssion ####
	cat('restructuring data...')
	regDFColnames <- c('id','year','sta','fsta','l1sta','l2sta','l3sta',
										 'gdp','l1gdp','l2gdp','l3gdp','fgdp','rgdp','gdppc','gdpGro','gdpGroRt','gdppcGro','gdppcGroRt',
										 'pop','edem','esup','esht','inf','prinv','puinv')
	regDF <- matrix(NA,nrow=nid*nyear,ncol=length(regDFColnames))
	colnames(regDF) <- regDFColnames
	regDF[,'id'] <- rep(1:nid,nyear)
	regDF[,'year'] <- rep(years,each=nid)
	regDF[,'sta'] <- unname(unlist(sta[,-1]))
	regDF[,'fsta'] <- unname(unlist(fsta[,-1]))
	regDF[,'l1sta'] <- unname(unlist(l1sta[,-1]))
	regDF[,'l2sta'] <- unname(unlist(l2sta[,-1]))
	regDF[,'l3sta'] <- unname(unlist(l3sta[,-1]))
	regDF[,'gdp'] <- unname(unlist(gdp[,-1]))
	regDF[,'l1gdp'] <- unname(unlist(l1gdp[,-1]))
	regDF[,'l2gdp'] <- unname(unlist(l2gdp[,-1]))
	regDF[,'l3gdp'] <- unname(unlist(l3gdp[,-1]))
	regDF[,'fgdp'] <- unname(unlist(fgdp[,-1]))
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
	cat('done\n')
	cat('cleanup...')
	rm(list=c('gdp','fgdp','l1gdp','rgdp','gdpGro','gdppc','gdppcGro','pop','edem','esup','esht','inf','prinv','puinv',
						'sta','fsta','l1sta','l2sta','l3sta'))
	gc(verbose = F)
	# consider drop early years as burn in
	# regDF <- regDF[regDF$year>=2050,]
	nobs <- nrow(regDF)
	cat('done\n')
}

# data vis ####
cat('data vis...')
# png(file.path('figures','0-gdpVsl1gdp.png'),width = plw, height = plh, units = plu, res = pld)
# par(pch='.')
# plot(regDF$l1gdp,regDF$gdp)
# dev.off()

# png(file.path('figures','0-gdpVsSta.png'),width = plw, height = plh, units = plu, res = pld)
# par(pch='.')
# plot(regDF$sta,regDF$gdp)
# dev.off()

png(file.path('figures','0-gdpVsl1gdp-zoom.png'),width = plw, height = plh, units = plu, res = pld)
par(pch='.')
plot(regDF$l1gdp,regDF$gdp,xlim=c(0,10e6),ylim=c(0,10e6))
dev.off()

png(file.path('figures','0-gdpVsSta-zoom.png'),width = plw, height = plh, units = plu, res = pld)
par(pch='.')
plot(regDF$sta,regDF$gdp,xlim=c(0,10),ylim=c(0,5e6))
dev.off()

png(file.path('figures','0-gdpVsTime-zoom.png'),width = plw, height = plh, units = plu, res = pld)
par(pch='.')
plot(regDF$year,regDF$gdp,ylim=c(0,5e6))
dev.off()

png(file.path('figures','0-staVsTime-zoom.png'),width = plw, height = plh, units = plu, res = pld)
par(pch='.')
plot(regDF$year,regDF$sta,ylim=c(0,10))
dev.off()
cat('done\n')

# OLS ####
cat('empiric model...')
if(!twoStage){
	# regular model
	# lMod <- lm(gdp ~ poly(l1gdp,sta,degree=degree), data=regDF)
	
	# lMod <- lm(gdp ~ poly(l1gdp,l2gdp,l3gdp,sta,l1sta,l2sta,l3sta,degree=degree), data=regDF)
	
	# lMod <- lm(gdp ~ poly(l1gdp,sta,pop,esht,inf,prinv,puinv,degree=degree), data=regDF)
	 
	# lMod <- lm(fgdp ~ poly(gdp,l1gdp,fsta,sta,l1sta,pop,esht,inf,prinv,puinv,degree=degree), data=regDF)
	# lMod <- lm(fgdp ~ poly(gdp,l1gdp,l2gdp,fsta,sta,l1sta,l2sta,pop,esht,inf,prinv,puinv,degree=degree), data=regDF)
	#no interaction effects
	# lMod <- lm(fgdp ~ poly(gdp,degree=dgree)+
	# 					 	poly(l1gdp,degree=dgree)+
	# 					 	poly(l2gdp,degree=dgree)+
	# 					 	poly(fsta,degree=dgree)+
	# 					 	poly(sta,degree=dgree)+
	# 					 	poly(l1sta,degree=dgree)+
	# 					 	poly(l2sta,degree=dgree)+
	# 					 	poly(pop,degree=dgree)+
	# 					 	poly(esht,degree=dgree)+
	# 					 	poly(inf,degree=dgree)+
	# 					 	poly(prinv,degree=dgree)+
	# 					 	poly(puinv,degree=dgree),
	# 					 data=regDF)
	#no constant	
	lMod <- lm(fgdp ~ -1+poly(gdp,l1gdp,l2gdp,fsta,sta,l1sta,l2sta,pop,esht,inf,prinv,puinv,degree=degree), data=regDF)
	
	# lMod <- lm(gdp ~ poly(l1gdp,sta,degree=degree) + poly(l1gdp,pop,esht,inf,prinv,puinv,degree=covDegree), data=regDF)
} else{
	# 2 stage model
	# 1st stages remove the influence of sta from the other covariates
	l1gdpMod <- lm(l1gdp ~ poly(sta,degree=degree), data=regDF)
	regDF$l1gdpNoSta <- resid(l1gdpMod)
	rm(l1gdpMod)
	popMod <- lm(pop ~ poly(sta,degree=degree), data=regDF)
	regDF$popNoSta <- resid(popMod)
	rm(popMod)
	eshtMod <- lm(esht ~ poly(sta,degree=degree), data=regDF)
	regDF$eshtNoSta <- resid(eshtMod)
	rm(eshtMod)
	infMod <- lm(inf ~ poly(sta,degree=degree), data=regDF)
	regDF$infNoSta <- resid(infMod)
	rm(infMod)
	prinvMod <- lm(prinv ~ poly(sta,degree=degree), data=regDF)
	regDF$prinvNoSta <- resid(prinvMod)
	rm(prinvMod)
	puinvMod <- lm(puinv ~ poly(sta,degree=degree), data=regDF)
	regDF$puinvNoSta <- resid(puinvMod)
	rm(puinvMod)
	gc()
	# 2nd stage model with sta and covariates
	# lMod <- lm(gdp ~ poly(l1gdp,sta,degree=degree) + poly(l1gdp,popNoSta,eshtNoSta,infNoSta,prinvNoSta,puinvNoSta,degree=covDegree), data=regDF)
	lMod <- lm(gdp ~ poly(l1gdpNoSta,sta,popNoSta,eshtNoSta,infNoSta,prinvNoSta,puinvNoSta,degree=degree), data=regDF)
}
cat('done\n')
cat('diagnostic plots')
callString <- gsub(' ','',lMod$call)[2]
if(!twoStage){
	figFolder <- sprintf('OLS-N-%i-deg-%i-coDeg-%i-%s',nobs,degree,covDegree,callString)
} else{
	figFolder <- sprintf('2SLS-N-%i-deg-%i-coDeg-%i-%s',nobs,degree,covDegree,callString)
}
dir.create(file.path('figures',figFolder),recursive = T,showWarnings = F)
cat('.')
sink(file.path('figures',figFolder,'modelSummary.txt'))
summary(lMod)
sink()
cat('.')
png(file.path('figures',figFolder,'1-residVsl1gdp.png'),width = plw, height = plh, units = plu, res = pld)
par(pch='.')
plot(regDF$l1gdp, resid(lMod),
		 xlim=c(0,5e6),
		 # col=adjustcolor(1,alpha.f = 0.01),
		 main='prediction error')
abline(h=0,col='red')
mtext(figFolder,3,0.5,cex = 0.7)
dev.off()
cat('.')
png(file.path('figures',figFolder,'1-residVsSta.png'),width = plw, height = plh, units = plu, res = pld)
par(pch='.')
plot(regDF$sta, resid(lMod),
		 xlim=c(0,8),
		 # col=adjustcolor(1,alpha.f = 0.01),
		 main='prediction error')
abline(h=0,col='red')
mtext(figFolder,3,0.5,cex = 0.7)
dev.off()
cat('.')
residLims <- c(-2e6,2e6)
png(file.path('figures',figFolder,'1-residVsl1gdp-zoom.png'),width = plw, height = plh, units = plu, res = pld)
par(pch='.')
plot(regDF$l1gdp, resid(lMod),
		 xlim=c(0,5e6),ylim=residLims,
		 # col=adjustcolor(1,alpha.f = 0.01),
		 main='prediction error')
abline(h=0,col='red')
mtext(figFolder,3,0.5,cex = 0.7)
dev.off()
cat('.')
png(file.path('figures',figFolder,'1-residVsSta-zoom.png'),width = plw, height = plh, units = plu, res = pld)
par(pch='.')
plot(regDF$sta, resid(lMod),
		 xlim=c(0,8),
		 ylim=residLims,
		 # col=adjustcolor(1,alpha.f = 0.01),
		 main='prediction error')
abline(h=0,col='red')
mtext(figFolder,3,0.5,cex = 0.7)
dev.off()
cat('done\n')

cat('cleanup...')
gc(verbose = F)
cat('done\n')

# predict ####
if(makePredict){
	cat('Projections.')
	lModPred <- predict(lMod)
	cat('.')
	predDF <- regDF
	predDF$sta <- 0
	predDF$fsta <- 0
	predDF$l1sta <- 0
	predDF$l2sta <- 0
	predDF$l3sta <- 0
	lModPred0sta <- predict(lMod, newdata = predDF)
	cat('.')
	cat('done\n')
	cat('Damage Figures.')
	png(file.path('figures',figFolder,'2-gdpVsl1gdp-pred0.png'),width = plw, height = plh, units = plu, res = pld)
	par(pch='.')
	plot(regDF$l1gdp, regDF$gdp, xlim=c(0,2e6), ylim=c(-1e7,1e7),
			 main='gdp depending on gdp_t-1')
	points(regDF$l1gdp, lModPred,col='blue')
	points(regDF$l1gdp, lModPred0sta,col='red')
	legend('topleft',
				 legend=c('FRIDA EMB','Pred.','Pred. STA=0'),
				 col=c('black','blue','red'),
				 pch=20)
	mtext(figFolder,3,0.5,cex = 0.7)
	dev.off()
	cat('.')
	png(file.path('figures',figFolder,'2-gdpVsSta-pred0.png'),width = plw, height = plh, units = plu, res = pld)
	par(pch='.')
	plot(regDF$sta, regDF$gdp, xlim=c(0,8), ylim=c(-1e7,1e7),
			 main='gdp depending on gdp_t-1')
	points(regDF$sta, lModPred,col='blue')
	points(regDF$sta, lModPred0sta,col='red')
	legend('topleft',
				 legend=c('FRIDA EMB','Pred.','Pred. STA=0'),
				 col=c('black','blue','red'),
				 pch=20)
	mtext(figFolder,3,0.5,cex = 0.7)
	dev.off()
	cat('.')
		lModLoss <- lModPred0sta - regDF$gdp
	lModLossRel <- lModLoss/regDF$gdp
	cat('.')
	png(file.path('figures',figFolder,'2-relgdploss.png'),width = plw, height = plh, units = plu, res = pld)
	par(pch='.')
	plot(regDF$sta,lModLossRel,
			 xlim=c(0,8),
			 col=adjustcolor(1,alpha.f = 0.05),
			 # col=1,
			 main='year relative GDP loss')
	mtext(figFolder,3,0.5,cex = 0.7)
	dev.off()
	cat('.')
	png(file.path('figures',figFolder,'2-relgdploss-zoom.png'),width = plw, height = plh, units = plu, res = pld)
	par(pch='.')
	plot(0,0,
			 xlim=c(0,8),ylim=c(-20,80),
			 type='n',
			 # col=1,
			 xlab='STA',
			 ylab='relative yearly GDP loss in %',
			 main='year relative GDP loss')
	abline(h=seq(-20,80,10),lty=2,col='gray')
	abline(v=seq(0,8,1),lty=2,col='gray')
	points(regDF$sta,lModLossRel*100,
				 col=adjustcolor(1,alpha.f = 0.05))
	mtext(figFolder,3,0.5,cex = 0.7)
	dev.off()
	cat('done\n')
}

cat('cleanup...')
gc(verbose = F)
cat('done\n')

# model2 ####
if(makeDFmod){
	cat('Fitting damage function.')
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
	cat('done\n')
	cat('Damage function figure...')
	staSup <- seq(-par[1],10,length.out=200)
	png(file.path('figures',figFolder,'3-rel1gdplossFunFit.png'),width = plw, height = plh, units = plu, res = pld)
	par(pch='.')
	plot(staSup,predFitM(par,staSup),col='red',type='l',ylim=c(-0.1,0.8),lwd=3)
	dev.off()
	
	png(file.path('figures',figFolder,'2-relgdploss-zoom-funfit.png'),width = plw, height = plh, units = plu, res = pld)
	par(pch='.')
	plot(0,0,
			 xlim=c(0,8),ylim=c(-20,80),
			 type='n',
			 # col=1,
			 xlab='STA',
			 ylab='relative yearly GDP loss in %',
			 main='year relative GDP loss')
	abline(h=seq(-20,80,10),lty=2,col='gray')
	abline(v=seq(0,8,1),lty=2,col='gray')
	points(regDF$sta,lModLossRel*100,
				 col=adjustcolor(1,alpha.f = 0.05))
	lines(staSup,predFitM(par,staSup),col='red',type='l',lwd=3)
	mtext(figFolder,3,0.5,cex = 0.7)
	dev.off()
	
	cat('done\n')
}
# FE model
# discarded as errors are larger than in the corresponding standard model
# this is due to the fe only accounting for differences in the level between groups
# not for differences in trends between groups. In fact removing the average
# levels for each groups means that their origins no longer conincide for cases where
# data have strong trends, such as this.
# 
