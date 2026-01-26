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

dataFile <- file.path('outputData',paste0('data-',configStr,'.RDS'))

if(file.exists(dataFile)){
	cat(paste0('data exists: ',dataFile,'\nreading...'))
	data <- readRDS(dataFile)
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
	data <- data.frame(id=numeric(nrowOutput),staID=numeric(nrowOutput),
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
		data<- cbind(data, newcolgdp, newcolsta)
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
		data$id[idxForThisSTAts]     <- rep(gdpForThisSTAts$id,ncol(gdpForThisSTAts)-1)
		data$staID[idxForThisSTAts]  <- STA.i
		data$year[idxForThisSTAts]   <- rep(as.numeric(colnames(gdpForThisSTAts[,-1])),each=nrow(gdpForThisSTAts))
		data$STA[idxForThisSTAts]    <- unlist(staForThisSTAts[,-1])
		data$GDP[idxForThisSTAts]    <- unlist(gdpForThisSTAts[,-1])
		for(ts.i in 1:length(timeshifts)){
			data[[timeshiftColNamesGDP[ts.i]]]   <- unlist(gdpForThisSTAts.timeshifts[[ts.i]][-1])
			data[[timeshiftColNamesSTA[ts.i]]] <- unlist(staForThisSTAts.timeshifts[[ts.i]][-1])
		}
	}
	numIDs <- length(unique(data$id))
	cat('\n')
	
	# save data frame ####
	cat('Saving data frame...RDS...')
	saveRDS(data,dataFile)
	cat('csv...')
	write.table(data,file.path('outputData',paste('data-',configStr,'.csv')),
							sep = ',',row.names = F)
	cat('done\n')
}

# Data plot ####
fig.dir <- 'figures/forcedSTAts'
fig.w <- 15
fig.h <- 15
fig.res <- 300
fig.u <- 'cm'

staCols <- rainbow(numSTAts*3+1)[(numSTAts*2+1):(numSTAts*3)]
countAlpha <- seq(0.3,1,length.out=20)
countBreaks <- c(0,1e-10,exp(1:(length(countAlpha)))[-1]/exp(length(countAlpha)))
yearBreaks <- seq(1979.5,2150.5,1)
gdpBreaks <- seq(0,2e6,length.out=500)
staIDbreaks <- 0.5:(numSTAts+0.5)
staBreaks <- seq(0,max(data$STA,na.rm=T),length.out=200)
staBreaksCol <- seq(0,max(data$STA,na.rm=T),length.out=length(staCols)+1)


hist2d <- function(x,y=NULL,xbreaks=NULL,ybreaks=NULL,
									 col='black',
									 countAlpha = seq(0,1,length.out=10),
									 countBreaks=NULL,
									 z=NULL,zbreaks=NULL,zcols=NULL,
									 plot = T,
									 returnCounts = F,
									 dens=F,verbose=T,counts=NULL){
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
	redoCounts <- F
	if(is.null(counts)){
		redoCounts <- T
	}	else if(!is.null(z) && sum(dim(counts)==c(length(xbreaks)-1,length(ybreaks)-1,length(zbreaks)-1))<3){
		redoCounts <- T
	} else if(sum(dim(counts)==c(length(xbreaks)-1,length(ybreaks)-1))<2){
		redoCounts <- T
	}
	if(redoCounts){
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
				counts[x.i,] <- histDat$counts[c(-1,-(length(ybreaks)+1))]
			}
			if(dens){
				counts <- counts/nrow(xy)
			}
		}
	}
	if(plot){
		if(length(dim(counts))==3){
			maxCounts <- max(colSums(counts,na.rm=T,dims=2))
			zcols <- col2rgb(zcols)
		} else {
			maxCounts <- max(counts)
		}
		if(is.null(countBreaks)){
			countBreaks <- c(0,1,exp(1:(length(countAlpha-1)))[-1]/exp(length(countAlpha-1))*maxCounts)
		} else if(length(countBreaks)!=(length(countAlpha)+1)){
			stop('Length of countAlpha has to be length of countbreaks minus one.\n')
		} else if(max(countBreaks)<=maxCounts){
			countBreaks <- countBreaks*maxCounts
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
					alpha <- countAlpha[which(diff(count<=countBreaks)==1)]
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

addCountLegend <- function(xleft,ybottom,xright,ytop,countBreaks,countAlpa,
													 maxCount,
													 col='black',xpd=T){
	countBreaks <- maxCount*countBreaks/max(countBreaks)
	oldXpd <- par('xpd')
	ystep <- (ytop-ybottom)/(1+length(countAlpha))
	vertBuff <- ystep/2
	xwidth <- xright-xleft
	par(xpd=xpd)
	rect(xleft,ybottom,xright,ytop,xpd=T,col='white')
	for(i in 1:length(countAlpa)){
		rect(xleft,
				 ybottom+vertBuff+ystep*(i-1),
				 xleft+xwidth*0.1,
				 ybottom+vertBuff+ystep*i,
				 col=adjustcolor(col,countAlpa[i]),
				 border = 'white')
	}
	for(i in 1:length(countBreaks)){
		text(xleft+xwidth*0.1,ybottom+vertBuff+ystep*(i-1),
				 sprintf(' %1.2e',countBreaks[i]),
				 adj=0)
	}
	rect(xleft,ybottom,xright,ytop,xpd=T)
	par(xpd=oldXpd)
}

dir.create(fig.dir,F,T)

## count alpha ####
cat('plotting countAlpha...')
png(file.path(fig.dir,'0-countVsAlpha.png'),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
par(mar=c(4.1,4.1,4.1,6))
plot(countBreaks[-length(countBreaks)],countAlpha,pch=20,
		 xlim=c(0,1),ylim=c(0,1))
grid()
addCountLegend(xleft = 1.05,ybottom = 0,xright = 1.35,ytop = 1,
							 countBreaks = countBreaks,countAlpa = countAlpha,
							 maxCount = 1)
dev.off()
cat('done\n')

## GDP ####
if(!file.exists(file.path(fig.dir,'1-GDP.png'))){
	cat('plotting GDP...')
	png(file.path(fig.dir,'1-GDP.png'),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
	plot(0,0,type='n',
			 xlim=range(data$year),ylim=c(0,2e6),
			 xlab='Year',ylab='billion constant 2021 $',
			 main='Real GDP')
	counts <- hist2d(x = data$year, y = data$GDP,
				 xbreaks = yearBreaks, ybreaks = gdpBreaks,
				 countAlpha = countAlpha, countBreaks = countBreaks,
				 returnCounts = T)
	addCountLegend(xleft = 2140,ybottom = 0,xright = 2200,ytop = 2e6,
								 countBreaks = countBreaks,countAlpa = countAlpha,
								 maxCount = max(counts))
	dev.off()
	cat('done\n')
}

## GDP per STAts ####
cat('plotting GDP per STAts')
for(sta.i in 1:numSTAts){
	if(!file.exists(file.path(fig.dir,paste0('1-GDP-STAtsID-',sta.i,'.png')))){
		cat('.')
		png(file.path(fig.dir,paste0('1-GDP-STAtsID-',sta.i,'.png')),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
		plot(0,0,type='n',
				 xlim=range(data$year),ylim=c(0,2e6),
				 xlab='Year',ylab='billion constant 2021 $',
				 main=paste0('Real GDP forced with STAts ',sta.i))
		hist2d(x = data$year[data$staID==sta.i],
					 y = data$GDP[data$staID==sta.i],
					 xbreaks = yearBreaks, ybreaks = gdpBreaks)
		dev.off()
	}
}
cat('done\n')

## GDP colored STAts ####
if(!file.exists(file.path(fig.dir,'1-GDP-STAtsID-color.png'))){
	cat('plotting GDP colored STAts...')
	png(file.path(fig.dir,'1-GDP-STAtsID-color.png'),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
	plot(0,0,type='n',
			 xlim=range(data$year),ylim=c(0,2e6),
			 xlab='Year',ylab='billion constant 2021 $',
			 main='Real GDP')
	hist2d(x = data$year, y = data$GDP,
				 xbreaks = yearBreaks, ybreaks = gdpBreaks,
				 z = data$staID,zbreaks = staIDbreaks,zcols = staCols)
	dev.off()
	cat('done\n')
}

## GDP colored STA ####
if(!file.exists(file.path(fig.dir,'1-GDP-STA-color.png'))){
	cat('plotting GDP colored STA...')
	png(file.path(fig.dir,'1-GDP-STA-color.png'),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
	plot(0,0,type='n',
			 xlim=range(data$year),ylim=c(0,2e6),
			 xlab='Year',ylab='billion constant 2021 $',
			 main='Real GDP')
	hist2d(x = data$year, y = data$GDP,
				 xbreaks = yearBreaks, ybreaks = gdpBreaks,
				 z = data$STA,zbreaks = staBreaksCol,zcols = staCols)
	dev.off()
	cat('done\n')
}

## STA ####
if(!file.exists(file.path(fig.dir,'2-STA.png'))){
	cat('plotting STA...')
	png(file.path(fig.dir,'2-STA.png'),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
	plot(data$year[which(data$id==1)],
			 data$STA[which(data$id==1)],
			 col=staCols[data$staID[which(data$id==1)]],
			 xlab='Year',ylab='K',
			 main='Surface Temperature Anomaly',
			 pch=20)
	dev.off()
	cat('done\n')
}

## STA vs GDP ####
if(!file.exists(file.path(fig.dir,'3-GDPvsSTA.png'))){
	cat('plotting STA vs GDP...')
	dir.create(fig.dir,F,T)
	png(file.path(fig.dir,'3-GDPvsSTA.png'),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
	plot(0,0,type='n',
			 xlim=c(0,max(data$STA,na.rm=T)),ylim=c(0,2e6),
			 xlab='Surface Temperature Anomaly K',ylab='GDP billion constant 2021 $',
			 main='GDP vs STA')
	hist2d(x = data$STA, y = data$GDP,
				 xbreaks = staBreaks, ybreaks = gdpBreaks)
	dev.off()
	cat('done\n')
}

## STA vs GDP per STAts ####
cat('plotting STA vs GDP per STAts')
for(sta.i in 1:numSTAts){
	if(!file.exists(file.path(fig.dir,paste0('1-GDPvsSTA-STAtsID-',sta.i,'.png')))){
		cat('.')
		png(file.path(fig.dir,paste0('1-GDPvsSTA-STAtsID-',sta.i,'.png')),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
		plot(0,0,type='n',
				 xlim=c(0,max(data$STA,na.rm=T)),ylim=c(0,2e6),
				 xlab='Surface Temperature Anomaly K',ylab='GDP billion constant 2021 $',
				 main=paste0('GDP vs STA for staID ',sta.i))
		hist2d(x = data$STA[data$staID==sta.i],
					 y = data$GDP[data$staID==sta.i],
					 xbreaks = staBreaks, ybreaks = gdpBreaks)
		dev.off()
	}
	cat('done\n')
}

## STA vs GDP ####
if(!file.exists(file.path(fig.dir,'3-GDPvsSTAcolByStaID.png'))){
	cat('plotting STA vs GDP...')
	dir.create(fig.dir,F,T)
	png(file.path(fig.dir,'3-GDPvsSTAcolByStaID.png'),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
	plot(0,0,type='n',
			 xlim=c(0,max(data$STA,na.rm=T)),ylim=c(0,2e6),
			 xlab='Surface Temperature Anomaly K',ylab='GDP billion constant 2021 $',
			 main='GDP vs STA')
	hist2d(x = data$STA, y = data$GDP,
				 xbreaks = staBreaks, ybreaks = gdpBreaks,
				 z = data$staID,zbreaks = staIDbreaks,zcols = staCols)
	dev.off()
	cat('done\n')
}

# fit model ####
dataComplete <- data[complete.cases(data),]
# counter factual data with zero sta
dataCF <- dataComplete
dataCF[,which(grepl('sta',colnames(dataCF),ignore.case = T))] <- 0

## prediction ####
predictDamFac <- function(){
	retList <- list()
	cat(sprintf('Predicting for %s\n', modelName))
	retList$pred <- predict(fitModel)
	cat(sprintf('Resid for %s\n', modelName))
	retList$resid <- resid(fitModel)
	cat(sprintf('Predicting CF for %s\n', modelName))
	retList$predCF <- predict(fitModel,dataCF)
	cat(sprintf('Calculting DF for %s\n', modelName))
	retList$DF <- 1-((retList$predCF - retList$pred)/retList$predCF)
	return(retList)
}


## plot ####
histResX <- 500
histResY <- 500
numColors <- 100
countAlpha <- seq(0.1,1,length.out=numColors)
countBreaks <- c(0,pnorm(countAlpha,mean=1.2,sd=0.2,log.p=F))
countBreaks[length(countBreaks)] <- 1
### resid ####
residCounts <- list() # cache for replotting
plotResid <- function(){
	cat(sprintf('Plotting obs vs resid for %s\n',modelName))
	png(file.path(fig.dir,paste0('4-',i,'-',modelName,'-GDPVsResid.png')),
			width = fig.w,height = fig.h,units = fig.u,res = fig.res)
	xrange <- range(dataComplete$GDP,na.rm=T)
	yrange <- range(resid,na.rm=T)
	plot(0,0,type='n',
			 xlim=xrange,ylim=yrange,
			 xlab='GDP billion constant 2021 $',ylab='resid',
			 main='GDP vs Resid')
	mtext(modelName,3,0.25)
	residCounts[[modelName]] <- hist2d(x=dataComplete$GDP,y=resid,
																		 xbreaks = seq(xrange[1],xrange[2],length.out=histResX),
																		 ybreaks = seq(yrange[1],yrange[2],length.out=histResY),
																		 countAlpha = countAlpha, countBreaks = countBreaks,
																		 counts = residCounts[[modelName]])
	dev.off()
}
### DF ####
DFCounts <- list() # cache for replotting
plotDF <- function(){
	cat(sprintf('Plotting DF for %s\n',modelName))
	png(file.path(fig.dir,paste0('5-',i,'-',modelName,'-STAVsDF.png')),
			width = fig.w,height = fig.h,units = fig.u,res = fig.res)
	xrange <- c(0,7)
	yrange <- c(0,0.8)
	plot(0,0,type='n',
			 xlim=xrange,ylim=yrange,
			 xlab='GDP billion constant 2021 $',ylab='resid',
			 main='GDP vs Resid')
	mtext(modelName,3,0.25)
	DFCounts[[modelName]] <- hist2d(x=dataComplete$STA,y=damFacRes$DF,
																	xbreaks = seq(xrange[1],xrange[2],length.out=histResX),
																	ybreaks = seq(yrange[1],yrange[2],length.out=histResY),
																	countAlpha = countAlpha, countBreaks = countBreaks,
																	counts = DFCounts[[modelName]])
	dev.off()
}



# OLS ####
cat('fitting ols model...')
fitModel <- lm(GDP~poly(stalag1,2)+poly(gdplag1,2)+I(gdplag1*stalag1),data=dataComplete)
modelName <- 'fitOLS'
eval(parse(text=fitStr))
sink(file.path(fig.dir,paste0(modelName,'-summary.txt')),append = F)
summary(fitModel)
sink()
cat('plotting Resid...')
plotResid()
cat('predict dam fac...')
damFacRes <- predictDamFac()
cat('plot dam fac...')
plotDF()
cat('done\n')

# FeOLS ####
degrees <- c(1,2,3,4)
lagss <- list(c(1),
							c(1,2),
							c(1,2,3),
							c(1,2,3,4,5,10,20))
## polynomial with interactions, id fixed effects
for(degree in degrees){
	for(lags.i in 1:length(lagss)){
		lags <- lagss[[lags.i]]
		modelName <- paste0('fit.feols.d',degree,'.l',paste(lags,collapse=','))
		cat(sprintf('fitting %s...',modelName))
		fitStr <- paste0('fitModel <- feols(GDP~')
		for(gdp.l in lags){
			for(gdp.d in 1:degree){
				fitStr <- paste0(fitStr,'+gdplag',gdp.l,'^',gdp.d)
			}
		}
		for(sta.l in lags){
			for(sta.d in 1:degree){
				fitStr <- paste0(fitStr,'+stalag',sta.l,'^',sta.d)
			}
		}
		for(gdp.l in lags){
			for(gdp.d in 1:degree){
				for(sta.l in lags){
					for(sta.d in 1:degree){
						if(gdp.d+sta.d<=degree){
							fitStr <- paste0(fitStr,'+stalag',sta.l,'^',sta.d,'*gdplag',gdp.l,'^',gdp.d)
						}
					}
				}
			}
		}
		fitStr <- paste0(fitStr,' | id, data=dataComplete)')
		library(fixest)
		eval(parse(text=fitStr))
		sink(file.path(fig.dir,paste0(modelName,'-summary.txt')),append = F)
		summary(fitModel)
		sink()
		cat('plotting Resid...')
		plotResid()
		cat('predict dam fac...')
		damFacRes <- predictDamFac()
		cat('plot dam fac...')
		plotDF()
		cat('done\n')
	}
}

## polynomial with interactions, id fixed effects interacting with gdp
for(degree in degrees){
	for(lags.i in 1:length(lagss)){
		modelName <- paste0('fit.feols.intFE.d',degree,'.l',paste(lags,collapse=','))
		cat(sprintf('fitting %s...',modelName))
		fitStr <- paste0('fitModel <- feols(GDP~')
		for(gdp.l in lags){
			for(gdp.d in 1:degree){
				fitStr <- paste0(fitStr,'+gdplag',gdp.l,'^',gdp.d)
			}
		}
		for(sta.l in lags){
			for(sta.d in 1:degree){
				fitStr <- paste0(fitStr,'+stalag',sta.l,'^',sta.d)
			}
		}
		for(gdp.l in lags){
			for(gdp.d in 1:degree){
				for(sta.l in lags){
					for(sta.d in 1:degree){
						if(gdp.d+sta.d<=degree){
							fitStr <- paste0(fitStr,'+stalag',sta.l,'^',sta.d,'*gdplag',gdp.l,'^',gdp.d)
						}
					}
				}
			}
		}
		fitStr <- paste0(fitStr,' | id[gdplag1], data=dataComplete)')
		library(fixest)
		eval(parse(text=fitStr))
		sink(file.path(fig.dir,paste0(modelName,'-summary.txt')),append = F)
		summary(fitModel)
		sink()
		cat('plotting Resid...')
		plotResid()
		cat('predict dam fac...')
		damFacRes <- predictDamFac()
		cat('plot dam fac...')
		plotDF()
		cat('done\n')
	}
}
