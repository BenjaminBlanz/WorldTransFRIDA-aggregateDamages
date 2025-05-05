aggDamWD <- getwd()
source('config-DetermineAggregateDamageFunction.R')
configStr <- paste0('S-',numSample,'-nSTAs-',length(STAs),'-nTW-',length(timeWindows),'-',paste0(format(Sys.time(), "%Y%m%d-%H%M%S")))

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
} else {
	cat('Baseline completed continuing\n')
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
	expID <- paste0(expIDpreString,'-S',numSample,'-policy_EMB-ClimateFeedback_On-',tools::file_path_sans_ext(staOverrideFileName))
	forcedRuns$expID[STA.i] <- expID
	statusFile <- file.path(location.fridaUncertaintyWD,'workOutput',expID,'status')
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
		system(paste('./submit_UncertaintyAnalysisLevante.sh',
								 '-n',numSample,
								 '--pol','policy_EMB.csv',
								 '--cfb','ClimateFeedback_On.csv',
								 '--sta',staOverrideFileName,
								 '-s',expIDpreString,
								 '--cid',baselineExpID,
								 '--cpps','true',
								 '--cpsp','true'))
	}
	setwd(aggDamWD)
	stop('Missing forced runs have been submitted to SLURM. Please restart this script once the forced runs have completed\n')
} else if(sum(forcedRuns$status=='completed')<nrow(forcedRuns)){
	setwd(aggDamWD)
	stop('Forced runs have not completed yet. Please restart this script when all of them have completed\n')
} else {
	cat('Forced runs completed continuing\n')
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

# save ####
cat('writing dataForDamFac...')
dir.create('outputData',F,T)
saveRDS(dataForDamFac,file = file.path('outputData',paste0('dataForDamFac',configStr,'.RDS')))
cat('done\n')

# agg losses ####
cat('cacluating aggregate damFac data...\n')
dataForDamFacTemplate <- matrix(NA,nrow=numIDs,ncol=1+length(STAs))
colnames(dataForDamFacTemplate) <- c('id',paste('STA',STAs))
dataForDamFacTemplate <- as.data.frame(dataForDamFacTemplate)
dataForDamFacTemplate$id <- unique(id)
dataForDamFacAgg <- list()
cat('  starting cluster...')
getFreeMemoryKB <- function() {
	osName <- Sys.info()[["sysname"]]
	if (osName == "Windows") {
		x <- system2("wmic", args =  "OS get FreePhysicalMemory /Value", stdout = TRUE)
		x <- x[grepl("FreePhysicalMemory", x)]
		x <- gsub("FreePhysicalMemory=", "", x, fixed = TRUE)
		x <- gsub("\r", "", x, fixed = TRUE)
		return(as.integer(x))
	} else if (osName == 'Linux') {
		x <- system2('free', args='-k', stdout=TRUE)
		x <- strsplit(x[2], " +")[[1]][4]
		return(as.integer(x))
	} else {
		stop("Only supported on Windows and Linux")
	}
}
cl <- makePSOCKcluster(min(length(STAs),
													 detectCores(),
													 max(1,
													 		unname(floor(getFreeMemoryKB()/(ps::ps_memory_info()[1]/1024)))
													 		)
													 )
											 )
cat('done\n  workers reading dataForDamFac...')
gobble <- clusterCall(cl,function(){
	dataForDamFac <<- readRDS(file.path('outputData',paste0('dataForDamFac',configStr,'.RDS')))
	return()})
gobble <- clusterEvalQ(cl,source('config-DetermineAggregateDamageFunction.R'))
cat('done\n  calculating per STA...')
parDamAggPerSTA <- function(STAlevel.i){
	dataForThisSTA <- dataForDamFac[dataForDamFac$STA==STAs[STAlevel.i],c('id','yRelLoss')]
	retlist <- list()
	retlist$Mean <- tapply(dataForThisSTA$yRelLoss,dataForThisSTA$id,mean,na.rm=T)
	retlist$Median <- tapply(dataForThisSTA$yRelLoss,dataForThisSTA$id,median,na.rm=T)
	return(retlist)
}
parRes <- parLapply(cl,1:length(STAs),parDamAggPerSTA)
dataForDamFacAgg$Mean <- dataForDamFacTemplate
dataForDamFacAgg$Median <- dataForDamFacTemplate
for(STAlevel.i in 1:length(STAs)){
	dataForDamFacAgg$Mean[,STAlevel.i+1] <- parRes[[STAlevel.i]]$Mean
	dataForDamFacAgg$Median[,STAlevel.i+1] <- parRes[[STAlevel.i]]$Median
}
cat('done\n  calculating per TW and STA')
dataForDamFacAgg$timewindow <- list()
parDamAggPerSTAandTW <- function(STAlevel.i,tw.i){
	dataForThisSTAandTW <- dataForDamFac[dataForDamFac$year>=timeWindows[[tw.i]][1] &
																			 	dataForDamFac$year<=timeWindows[[tw.i]][2] &
																			 	dataForDamFac$STA==STAs[STAlevel.i],
																			 c('id','yRelLoss')]
	retlist <- list()
	retlist$Mean <- tapply(dataForThisSTAandTW$yRelLoss,dataForThisSTAandTW$id,mean,na.rm=T)
	retlist$Median <- tapply(dataForThisSTAandTW$yRelLoss,dataForThisSTAandTW$id,median,na.rm=T)
	return(retlist)
}
for(tw.i in 1:length(timeWindows)){
	cat('.')
	dataForDamFacAgg$timewindow[[tw.i]] <- list()
	dataForDamFacAgg$timewindow[[tw.i]][['startYear']] <- timeWindows[[tw.i]][1]
	dataForDamFacAgg$timewindow[[tw.i]][['endYear']] <- timeWindows[[tw.i]][2]
	dataForDamFacAgg$timewindow[[tw.i]][['Mean']] <- dataForDamFacTemplate
	dataForDamFacAgg$timewindow[[tw.i]][['Median']] <- dataForDamFacTemplate
	parRes <- parLapply(cl,1:length(STAs),parDamAggPerSTAandTW,tw.i=tw.i)
	for(STAlevel.i in 1:length(STAs)){
		dataForDamFacAgg$timewindow[[tw.i]][['Mean']][,STAlevel.i+1] <- parRes[[STAlevel.i]]$Mean
		dataForDamFacAgg$timewindow[[tw.i]][['Median']][,STAlevel.i+1] <- parRes[[STAlevel.i]]$Median
	}
}
cat('done\n  stopping cluster...')
stopCluster(cl)
cat('done\n')
cat('saving dataForDamFacAgg...')
saveRDS(dataForDamFacAgg,file = file.path('outputData',paste0('dataForDamFacAgg',configStr,'.RDS')))
cat('done\n')

