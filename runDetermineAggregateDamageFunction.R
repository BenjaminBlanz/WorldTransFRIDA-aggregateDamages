aggDamWD <- getwd()
source('config-DetermineAggregateDamageFunction.R')

# function that writes the ClimateSTAOverride.csv file for use in 
# FRIDAforUncertaintyAnalysis
writeSTAForcing <- function(outputLocation,
														STAOverride,aggregateDamageFormationTime=1){
	filename <- paste0('ClimateSTAOverride_',gsub('\\.','_',STAOverride),'.csv')
	sink(file.path(outputLocation,filename),
			 append = F)
	cat('Energy Balance Model.SWITCH STA Override,1\n')
	cat(sprintf('Energy Balance Model.STA Override,%f\n',
							STAOverride))
	cat(sprintf('Energy Balance Model.Aggregate Damage Formation Time,%f\n',
							aggregateDamageFormationTime))
	sink()
	return(filename)
}

# enter the uncertainty analysis directory
setwd(location.fridaUncertaintyWD)

# baseline run ####
# send off the baseline run
expIDpreString <- 'determineAggDam'
baselineExpID <- paste0(expIDpreString,'_Baseline-S',numSample,'-policy_EMB-ClimateFeedback_On-ClimateSTAOverride_Off')
if(!file.exists(file.path(location.fridaUncertaintyWD,'workOutput',baselineExpID))){
	system(paste('./submit_UncertaintyAnalysisLevante.sh',
							 '-n',numSample,
							 '--pol','policy_EMB.csv',
							 '--cfb','ClimateFeedback_On.csv',
							 '--sta','ClimateSTAOverride_Off.csv',
							 '-s',paste0(expIDpreString,'_Baseline')))
	setwd(aggDamWD)
	stop('Baseline run has been submitted to SLURM, please restart this script once the baseline run has completed\n')
}
statusFile <- file.path(location.fridaUncertaintyWD,'workOutput',baselineExpID,'status')
if(!file.exists(statusFile)||
	 !readChar(statusFile,file.info(statusFile)$size-1)=='completed'){
	setwd(aggDamWD)
	stop('Baseline run has not completed yet. Run this script again once it has completed.\n')
}

# forced STA runs ####
# send off the forced STA runs
forcedRunsExpIDpreSring <- ''
forcedRuns <- data.frame(STA=STAs,expID=NA,staOverrideFileName=NA,status=NA)
for(STA.i in 1:length(STAs)){
	STA <- STAs[STA.i]
	staOverrideFileName <- writeSTAForcing(
		outputLocation=file.path(location.fridaUncertaintyWD,'FRIDA-configs'),
		STAOverride=STA)
	forcedRuns$staOverrideFileName[STA.i] <- staOverrideFileName
	expID <- paste0(expIDpreString,'_Baseline-S',numSample,'-policy_EMB-ClimateFeedback_On-',tools::file_path_sans_ext(staOverrideFileName))
	forcedRuns$expID[STA.i] <- expID
	statusFile <- file.path(location.fridaUncertaintyWD,'workOutput',expID)
	if(!file.exists(statusFile)){
		forcedRuns$status[STA.i] <- 'not present'
	} else if (readChar(statusFile,file.info(statusFile)$size-1)=='completed'){
		forcedRuns$status[STA.i] <- 'completed'
	} else {
		forcedRuns$status[STA.i] <- 'presumed running'
	}
}
if(sum(forcedRuns$status=='not present')>0){
	for(STA.i in which(forcedRuns$status=='not present')){
		staOverrideFileName <- forcedRuns$staOverrideFileName[STA.i]
		expID <- forcedRuns$expID[STA.i]
		expIDs$expID[STA.i] <- expID
		system(paste('./submit_UncertaintyAnalysisLevante.sh',
								 '-n',numSample,
								 '--pol','policy_EMB.csv',
								 '--cfb','ClimateFeedback_On.csv',
								 '--sta',staOverrideFileName,
								 '-s',expIDpreString,
								 '--cid',baselineExpID,
								 '--spps','true',
								 '--cpsp','true'))
	}
	setwd(aggDamWD)
	stop('Missing forced runs have been submitted to SLURM. Please restart this script once the forced runs have completed\n')
} else if(sum(forcedRuns$status=='completed')<nrow(forcedRuns)){
	setwd(aggDamWD)
	stop('Forced runs have not completed yet. Please restart this script when all of them have completed\n')
}

# return to the aggDam working directory
setwd(aggDamWD)

# read STA data ####
dataForDamFac <- data.frame(id=numeric(),y=numeric(),ylag=numeric(),yfut=numeric(),year=numeric(),STA=numeric())
dataForDamFacSTAs.lst <- list()
for(STAlevel.i in 1:length(STAs)){
	cat(sprintf('Processing STA level %f\n',STAs[STAlevel.i]))
	gdpForThisSTA <- readRDS(file.path(location.fridaUncertaintyWD,'workOutput',forcedRuns$expID[STAlevel.i],
																		 'detectedParmSpace','PerVarFiles-RDS','gdp_real_gdp_in_2021c.RDS'))
	# loglikeForThisSTA <- readRDS(file.path(location.fridaUncertaintyWD,'workOutput',forcedRuns$expID[STAlevel.i],
	# 																			 'detectedParmSpace','PerVarFiles-RDS','logLike.RDS'))
	
	# NA out any incomplete entry but do not remove the rows, so that the indexing is not broken
	# gdpForThisSTA[!complete.cases(gdpForThisSTA),-1] <- NA
	gdpForThisSTA.lag <- cbind(gdpForThisSTA[,1],matrix(NA,nrow=nrow(gdpForThisSTA),ncol=1),gdpForThisSTA[,c(2:(ncol(gdpForThisSTA)-1))])
	colnames(gdpForThisSTA.lag) <- colnames(gdpForThisSTA)
	gdpForThisSTA.fut <- cbind(gdpForThisSTA[,1],gdpForThisSTA[,c(3:(ncol(gdpForThisSTA)))],matrix(NA,nrow=nrow(gdpForThisSTA),ncol=1))
	colnames(gdpForThisSTA.fut) <- colnames(gdpForThisSTA)
	id <- rep(gdpForThisSTA$id,ncol(gdpForThisSTA)-1)
	year <- rep(as.numeric(colnames(gdpForThisSTA[,-1])),each=nrow(gdpForThisSTA))
	y <- unlist(gdpForThisSTA[,-1])
	ylag <- unlist(gdpForThisSTA.lag[,-1])
	yfut <- unlist(gdpForThisSTA.fut[,-1])
	ygro <- yfut - y
	dataForThisSTA <- data.frame(id=id,y=y,ylag=ylag,yfut=yfut,year=year,STA=rep(STAs[STAlevel.i],length(y)),ygro=ygro)
	dataForDamFacSTAs.lst[[STAlevel.i]] <- dataForThisSTA
	dataForDamFac <- rbind(dataForDamFac,dataForThisSTA)
	rm('y','ylag','yfut','ygro','dataForThisSTA','gdpForThisSTA','gdpForThisSTA.lag','gdpForThisSTA.fut')
}
numIDs <- length(unique(id))

# synth counterfact models ####
cat('calculating synthetic growth models...')
dataForSyntheticCounterfactual <- dataForDamFacSTAs.lst[[which(STAs==0)]]
# for testing
# summary(groMod0)
# plot(dataForSyntheticCounterfactualMemberI$y,dataForSyntheticCounterfactualMemberI$ygro,pch=20)
# points(dataForSyntheticCounterfactualMemberI$y[-length(dataForSyntheticCounterfactualMemberI$y)],predict(groMod0),col='red',pch=20)
parSynthPredGDP <- function(i){
	dataForSyntheticCounterfactualMemberI <- dataForSyntheticCounterfactual[dataForSyntheticCounterfactual$id==i,]
	if(is.na(dataForSyntheticCounterfactualMemberI$y[1])){
		return(rep(NA,2))
	} else {
		groMod0 <- lm(ygro~y,data=dataForSyntheticCounterfactualMemberI)
		return(groMod0$coefficients)
	}
}
library(parallel)
cl <- makePSOCKcluster(detectCores())
clusterExport(cl,'dataForSyntheticCounterfactual')
parRes <- parLapply(cl,unique(id),parSynthPredGDP)
stopCluster(cl)
cat('done\n')

# losses ####
cat('calculating climate related losses...')
groMod0Coefs <- matrix(unlist(parRes),nrow=numIDs,byrow=T)
colnames(groMod0Coefs) <- names(parRes[[1]])
# this relies on the dataFoorDamFac being ordered by sta, year, id. I.e.
# id year STA
# 1   1    1
# 2   1    1
# 3   1    1
# 1   2    1
# 2   2    1
numEntriesPerId <- length(unique(dataForDamFac$year))*length(unique(dataForDamFac$STA))
# predicted growth if STA were 0
dataForDamFac$predGro <- rep(groMod0Coefs[,'(Intercept)'],numEntriesPerId) + dataForDamFac$y * rep(groMod0Coefs[,'y'],numEntriesPerId)
# predicted gdp if STA were 0
dataForDamFac$predy <- dataForDamFac$ylag+
	rep(groMod0Coefs[,'(Intercept)'],numEntriesPerId) + dataForDamFac$ylag * rep(groMod0Coefs[,'y'],numEntriesPerId)
# loss from climate is difference between predicted (if STA were 0) - actual gdp
dataForDamFac$yloss <- dataForDamFac$predy - dataForDamFac$y
# as relative
dataForDamFac$yRelLoss <- dataForDamFac$yloss/dataForDamFac$predy
# mean losses per ensemble member
dataForDamFacAgg <- matrix(NA,nrow=numIDs,ncol=1+length(STAs))
colnames(dataForDamFacAgg) <- c('id',paste('STA',STAs))
dataForDamFacAgg <- as.data.frame(dataForDamFacAgg)
dataForDamFacAgg$id <- unique(id)
for(STAlevel.i in 1:length(STAs)){
	dataForThisSTA <- dataForDamFac[dataForDamFac$STA==STAs[STAlevel.i]&!is.na(dataForDamFac$yRelLoss),c('id','yRelLoss')]
	dataForDamFacAgg[,STAlevel.i+1] <- tapply(dataForThisSTA$yRelLoss,dataForThisSTA$id,mean,na.rm=T)
}
cat('done\n')

# plots ####
cat('plotting...')
# 
# plot(0,0,type='n',
# 		 ylim=c(-4e4,1.5e5),xlim=c(0,4e6),
# 		 xlab='real GDP',ylab='real GDP growth')
# grid()
# for(STAlevel.i in 1:length(STAs)){
# 	# sampleIdc <- sample(1:nrow(dataForDamFacSTAs.lst[[STAlevel.i]]),
# 	# 													 min(nrow(dataForDamFacSTAs.lst[[STAlevel.i]]), 4000))
# 	sampleIdc <- seq(1,nrow(dataForDamFacSTAs.lst[[STAlevel.i]]),1e3)
# 	points(dataForDamFacSTAs.lst[[STAlevel.i]]$y[sampleIdc],
# 				 dataForDamFacSTAs.lst[[STAlevel.i]]$ygro[sampleIdc],col=adjustcolor(STAs[STAlevel.i]+1,alpha.f=0.2),
# 				 pch=20)
# }
# points(dataForSyntheticCounterfactual$y,dataForSyntheticCounterfactual$ygro,pch=20,col='red')
# 
# legend('topleft',legend=paste0('STA ',STAs,'°'),pch=20,col=STAs+1)

source('myRudiVioPlot.R')
plot(0,0,
		 xlab='STA degC',ylab='Annual percentage loss of GDP',
		 xlim=c(0, 7),
		 ylim=c(-10,80),
		 type='n',yaxs='i')
grid()	
# points(dat$STA,dat$yRelLoss*100,col=1,pch=20)
for(STA in STAs){
	myRudiViolinPlot(dataForDamFac$yRelLoss[dataForDamFac$STA==STA]*100,
					at=STA,col='white',add=T,width=0.3,border='black',lwd=3,
					equiprobspacing = T,n=100)
	boxplot(dataForDamFac$yRelLoss[dataForDamFac$STA==STA]*100,at=STA,col='white',add=T,boxwex=0.4,axes=F,range = 0,border='black',lwd=2)
	myRudiViolinPlot(as.vector(dataForDamFacAgg[[paste('STA',STA)]])*100,
									 at=STA,col=NA,border.col = 'gray',add=T,width=0.3,lwd=2,
									 equiprobspacing = T,n=100)
	boxplot(dataForDamFacAgg[[paste('STA',STA)]]*100,at=STA,col='white',border='gray',add=T,boxwex=0.15,axes=F,range = 0,lwd=2)
}
figDir <- file.path('figures','ipccDamages')
dir.create(figDir,F,T)
dev.copy(png,file.path(figDir,'ipccDamgeFunction'),width='960',height='1380')
dev.off()

cat('done\n')
