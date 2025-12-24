aggDamWD <- getwd()
source('config-DetermineAggregateDamageFunction.R')
configStr <- paste0('S-',numSample,'-nSTAts-',numSTAts,'-',identifier)

# function that writes the ClimateSTAOverride.csv file for use in 
# FRIDAforUncertaintyAnalysis
writeSTAForcingTS <- function(outputLocation,staTimeseriesID){
	filename <- paste0('ClimateSTAOverrideTS_',gsub('\\.','_',staTimeseriesID),'.csv')
	sink(file.path(outputLocation,filename),
			 append = F)
	cat('year,Energy Balance Model.STA Override TS\n')
	for(c.i in 2:ncol(embSta)){
		year <- colnames(embSta)[c.i]
		cat(paste(year,embSta[staTimeseriesID,c.i],sep=','))
		cat('\n')
	}
	sink()
	return(filename)
}

dataForDamFacFile <- file.path('outputData',paste('dataForDamFac',configStr,'.RDS',sep='-'))

if(exists(dataForDamFacFile)){
	cat(paste0('dataForDamFac exists: ',dataForDamFacFile,'\nreading...'))
	dataForDamFac <- readRDS(dataForDamFacFile)
	cat('done\n')
} else {
	# enter the uncertainty analysis directory
	setwd(location.fridaUncertaintyWD)
	
	# baseline run ####
	# send off the baseline run
	expIDpreString <- 'determineAggDam'
	baselineExpID <- paste0(expIDpreString,'_Baseline-S',numSample,'-policy_EMB-ClimateFeedback_On-ClimateSTAOverride_Off')
	if(!file.exists(file.path(location.fridaUncertaintyWD,'workOutput',baselineExpID))){
		system(paste('./submit_UncertaintyAnalysisLevante.sh',
								 '-n',numSample,
								 '-h',7,
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
	
	# read STA from baseline ####
	embSta <- readRDS(file.path(location.fridaUncertaintyWD,'workOutput',baselineExpID,
													 'detectedParmSpace','PerVarFiles-RDS',
													 'energy_balance_model_surface_temperature_anomaly.RDS'))
	years <- as.numeric(colnames(embSta[-1]))
	numYears <- length(years)
	# embSta <- readRDS('workOutput/determineAggDam_Baseline-S20000-policy_EMB-ClimateFeedback_On-ClimateSTAOverride_Off/detectedParmSpace/PerVarFiles-RDS/energy_balance_model_surface_temperature_anomaly.RDS')
	embStaTempOrd <- embSta[order(embSta[['2150']]),]
	idsToSample <- embStaTempOrd$id[round(1+(nrow(embStaTempOrd)-1)*(0:(numSTAts-1))/(numSTAts-1))]
	rm(embStaTempOrd)
	
	# forced STA runs ####
	# send off the forced STA runs
	forcedRunsExpIDpreSring <- 'forcedSTAts'
	forcedRuns <- data.frame(STAid=idsToSample,expID=NA,staOverrideFileName=NA,status=NA)
	for(STA.i in 1:numSTAts){
		staOverrideFileName <- writeSTAForcingTS(
			outputLocation=file.path(location.fridaUncertaintyWD,'FRIDA-configs'),
			staTimeseriesID=idsToSample[STA.i])
		forcedRuns$staOverrideFileName[STA.i] <- staOverrideFileName
		expID <- paste0(expIDpreString,'-',tools::file_path_sans_ext(staOverrideFileName),
										'-S',numSample,'-policy_EMB-ClimateFeedback_On-ClimateSTAOverride_TS')
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
			system(paste('./submit_UncertaintyAnalysisLevante.sh',
									 '-n',numSample,
									 '--pol','policy_EMB.csv',
									 '--cfb','ClimateFeedback_On.csv',
									 '--sta','ClimateSTAOverride_TS.csv',
									 '--stats',staOverrideFileName,
									 '-s',paste0(expIDpreString,'-',tools::file_path_sans_ext(staOverrideFileName)),
									 '--cid',baselineExpID,
									 '--cpps','true',
									 '--cpsp','true',
									 '--plt','false'))
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
	
	# timeshift function ####
	timeshiftData <- function(data, amount){
		if(amount > 0){
			resData <- cbind(data[,1],data[,(amount+2):ncol(data)],matrix(NA,nrow=nrow(data),ncol=amount))
		} else if (amount < 0){
			resData <- cbind(data[,1],matrix(NA,nrow=nrow(data),ncol=-amount),data[,c(2:(ncol(data)+amount))])
		} else {
			resData <- data
		}
		colnames(resData) <- colnames(data)
		return(resData)
	}
	
	# read STA data ####
	timeshifts <- c(1,-1:-5,-10,-20)
	nrowOutput <- numSample*numYears*numSTAts
	dataForDamFac <- data.frame(id=numeric(nrowOutput),staID=numeric(nrowOutput),
															year=numeric(nrowOutput),STA=numeric(nrowOutput),
															GDP=numeric(nrowOutput))
	timeshiftColNamesSTA <- timeshiftColNamesGDP <- c()
	for(ts.i in 1:length(timeshifts)){
		if(timeshifts[ts.i]<0){
			eval(parse(text=paste0('newcolgdp <- data.frame(gdplag',-timeshifts[ts.i],'=numeric(nrowOutput))')))
			eval(parse(text=paste0('newcolsta <- data.frame(stalag',-timeshifts[ts.i],'=numeric(nrowOutput))')))
		} else {
			eval(parse(text=paste0('newcolgdp <- data.frame(gdpfut',timeshifts[ts.i],'=numeric(nrowOutput))')))
			eval(parse(text=paste0('newcolsta <- data.frame(stafut',timeshifts[ts.i],'=numeric(nrowOutput))')))
		}
		timeshiftColNamesGDP[ts.i] <- colnames(newcolgdp)
		timeshiftColNamesSTA[ts.i] <- colnames(newcolsta)
		dataForDamFac<- cbind(dataForDamFac, newcolgdp, newcolsta)
	}
	cat('Processing STAts\n')
	for(STA.i in 1:numSTAts){
		cat(sprintf('\rProcessing STA.i %i of %i',STA.i,numSTAts))
		gdpForThisSTAts <- readRDS(file.path(location.fridaUncertaintyWD,'workOutput',forcedRuns$expID[STA.i],
																			 'detectedParmSpace','PerVarFiles-RDS','gdp_real_gdp_in_2021c.RDS'))
		staForThisSTAts <- readRDS(file.path(location.fridaUncertaintyWD,'workOutput',forcedRuns$expID[STA.i],
																			 'detectedParmSpace','PerVarFiles-RDS','energy_balance_model_surface_temperature_anomaly.RDS'))
		loglikeForThisSTAts <- readRDS(file.path(location.fridaUncertaintyWD,'workOutput',forcedRuns$expID[STA.i],
																					 'detectedParmSpace','PerVarFiles-RDS','logLike.RDS'))
		# NA out any incomplete entry but do not remove the rows, so that the indexing is not broken
		gdpForThisSTAts[loglikeForThisSTAts$logLike < -1e200, -1] <- NA
		
		gdpForThisSTAts.timeshifts <- list()
		staForThisSTAts.timeshifts <- list()
		for(ts.i in 1:length(timeshifts)){
			gdpForThisSTAts.timeshifts[[ts.i]] <- timeshiftData(gdpForThisSTAts,timeshifts[ts.i])
			staForThisSTAts.timeshifts[[ts.i]] <- timeshiftData(staForThisSTAts,timeshifts[ts.i])
		}
		
		idxForThisSTAts <- 1:(numSample*numYears)+(numSample*numYears)*(STA.i-1)
		dataForDamFac$id[idxForThisSTAts]     <- rep(gdpForThisSTAts$id,ncol(gdpForThisSTAts)-1)
		dataForDamFac$staID[idxForThisSTAts]  <- STA.i
		dataForDamFac$year[idxForThisSTAts]   <- rep(as.numeric(colnames(gdpForThisSTAts[,-1])),each=nrow(gdpForThisSTAts))
		dataForDamFac$STA[idxForThisSTAts]    <- unlist(staForThisSTAts[,-1])
		dataForDamFac$GDP[idxForThisSTAts]    <- unlist(gdpForThisSTAts[,-1])
		for(ts.i in 1:length(timeshifts)){
			dataForDamFac[[timeshiftColNamesGDP[ts.i]]]   <- unlist(gdpForThisSTAts.timeshifts[[ts.i]][-1])
			dataForDamFac[[timeshiftColNamesSTA[ts.i]]] <- unlist(staForThisSTAts.timeshifts[[ts.i]][-1])
		}
	}
	numIDs <- length(unique(dataForDamFac$id))
	cat('\n')
	
	# save data frame ####
	cat('Saving data frame...RDS...')
	saveRDS(dataForDamFac,dataForDamFacFile)
	cat('csv...')
	write.table(dataForDamFac,file.path('outputData',paste('dataForDamFac',configStr,'.csv',sep='-')),
							sep = ',',row.names = F)
	cat('done\n')
}

# Data plot ####
fig.dir <- 'figures/forcedSTAts'
fig.w <- 15
fig.h <- 15
fig.res <- 300
fig.u <- 'cm'

hist2d <- function(x,y=NULL,xbreaks=NULL,ybreaks=NULL,
									 col='black',
									 countAlpha = seq(1,0,length.out=10),
									 countbreaks=NULL,
									 z=NULL,zbreaks=NULL,zcols=NULL,
									 plot = T,
									 returnCounts = F,
									 dens=F,verbose=T,){
	if(is.null(xbreaks)){
		xbreaks <- hist(x,plot=F)$breaks
	}
	if(is.null(ybreaks)){
		ybreaks <- hist(y,plot=F)$breaks
	}
	if(is.data.frame(x)){
		if(ncol(x)>=2){
			xy <- x[,1:2]
			colnames(xy) <- c('x','y')
		} 
		if(ncol(x)==3){
			z <- x[,3]
		}
		if(ncol(x)>3){
			if(verbose){
				cat('Columns >3 ignored.\n')
			}
		}
	} else {
		xy <- data.frame(x=x, y=y)
	}
	if(!is.null(z)){
		if(is.null(zbreaks)|is.null(zcols)){
			stop('If z is provided both zbreaks and zcols have to be specified. zcols should have one fewer entry than zbreaks.\n')
		}
		zbreaks[length(zbreaks)] <- Inf
		counts <- array(NA, dim=c(length(xbreaks)-1,length(ybreaks)-1,length(zbreaks)-1))
		for(zb.i in 1:(length(zbreaks)-1)){
			xy.sel <- xy[which(z>=zbreaks[zb.i]&z<zbreaks[zb.i+1]),]
			counts[,,zb.i] <- hist2d(xy.sel,xbreaks=xbreaks,ybreaks=ybreaks,plot=F,
															 returnCounts=T,dens=dens)
		}
	} else {
		counts <- array(NA,dim=c(length(xbreaks)-1,length(ybreaks)-1))
		for(x.i in 1:(length(xbreaks)-1)){
			histDat <- hist(xy$y[xy$x>=xbreaks[x.i] & xy$x<xbreaks[x.i+1]],
											breaks=c(-Inf,ybreaks,Inf),plot=F)
			counts[x.i,] <- histDat$counts[c(-1,-length(ybreaks+2))]
		}
		if(dens){
			counts <- counts/nrow(xy)
		}
	}
	if(plot){
		if(is.null(countbreaks)){
			if(length(dim(counts)==3)){
				maxCounts <- max(colSums(counts,na.rm=T,dims=2))
				zcols <- col2rgb(zcols)
			} else {
				maxCounts <- max(counts)
			}
			countbreaks <- c(0,sqrt(seq(1,maxCounts^2,length.out=length(countAlpha)+1)))
		} else if(length(countbreaks)!=(length(countcols)+1)){
			stop('Length of countcols has to be length of countbreaks minus one.\n')
		}
		for(x.i in 1:(length(xbreaks)-1)){
			for(y.i in 1:(length(ybreaks)-1)){
				if(!is.null(z)){
					mixedCol <- c()
					for(component in 1:3){
						mixedCol[component] <- weighted.mean(x = zcols[component,], w = counts[x.i,y.i,])
					}
					mixedCol[is.nan(mixedCol)] <- 0
					mixedCol <- rgb(mixedCol[1],mixedCol[2],mixedCol[3],maxColorValue=255)
					count <- sum(counts[x.i,y.i,])
				} else {
					mixedCol <- col
					count <- counts[x.i,y.i]
				}
				if(count==0){
					alpha <- 0	
				} else {
					alpha <- countAlpha[which(diff(count<=countbreaks)==1)]
				}
				rect(xbreaks[x.i],ybreaks[y.i],xbreaks[x.i+1],ybreaks[y.i+1],
						 density=-1,
						 col = adjustcolor(mixedCol,alpha),
						 border = NA)
			}
		}
	}
	if(returnCounts){
		return(counts)
	}
}

dir.create(fig.dir,F,T)
staCols <- rainbow(numSTAts*2+1)[1:numSTAts]
densAlpha <- (seq(0,1,length.out=20))^(1/2)
yearBreaks <- seq(1979.5,2150.5,1)
gdpBreaks <- seq(0,2e6,length.out=1000)

cat('plotting GDP...')
png(file.path(fig.dir,'1-GDP.png'),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
plot(0,0,type='n',
		 xlim=range(dataForDamFac$year),ylim=c(0,2e6),
		 xlab='Year',ylab='billion constant 2021 $',
		 main='Real GDP')
hist2d(x = dataForDamFac$year, y = dataForDamFac$GDP,
			 xbreaks = yearBreaks, ybreaks = gdpBreaks)
dev.off()
cat('done\n')

staCounts <- array(NA,dim=c(length(yearBreaks)-1,length(gdpBreaks)-1,numSTAts))
for(sta.i in 1:numSTAts){
	staCounts[,,sta.i] <- hist2d(x = dataForDamFac$year[dataForDamFac$staid==sta.i],
															 y = dataForDamFac$GDP[dataForDamFac$staid==sta.i],
															 xbreaks = seq(1979.5,2150.5,1), ybreaks = seq(0,2e6,length.out=1000),
															 returnCounts= T,plot = F)
}

cat('plotting STA...')
png(file.path(fig.dir,'2-STA.png'),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
plot(dataForDamFac$year[which(dataForDamFac$id==1)],
		 dataForDamFac$STA[which(dataForDamFac$id==1)],
		 col=staCols[dataForDamFac$staID[which(dataForDamFac$id==1)]],
		 pch=20,
		 xlab='Year',ylab='K',
		 main='Surface Temperature Anomaly',
		 pch=20)
dev.off()
cat('done\n')

cat('plotting STA vs GDP...')
dir.create(fig.dir,F,T)
png(file.path(fig.dir,'3-GDPvsSTA.png'),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
plot(dataForDamFac$year,dataForDamFac$GDP,
		 col=staCols[dataForDamFac$staID],
		 xlab='Year',ylab='billion 2021 $',
		 pch=20)
dev.off()
cat('done\n')
