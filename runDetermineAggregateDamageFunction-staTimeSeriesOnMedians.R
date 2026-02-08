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
medianDataFile <- file.path('outputData',paste0('medianData-',configStr,'.RDS'))
if(file.exists(medianDataFile)){
	mData <- readRDS(medianDataFile)
} else {
	if(file.exists(dataFile)){
		cat(paste0('data exists: ',dataFile,'\nreading...'))
		data <- readRDS(dataFile)
		cat('done\n')
	} else {
		stop('run config-DetermineAggregateDamageFunction.R to generate the data')
	}
	# check for completeness /  drop IDs that have an NA at any time or STAtsID
	GDPagg <- aggregate(data,by=list(data$id),mean)
	idsToRemove <- GDPagg$id[which(is.na(GDPagg$GDP))]
	data <- data[!data$id%in%idsToRemove,]
	# determine median data
	# - 1,2,3 removes staID, year, and id columns from data to avoid duplicate columns
	mData <- aggregate(data[,-c(1,2,3)], by=list(staID=data$staID,year=data$year), median)
	saveRDS(mData,medianDataFile)
}

# manually drop staID 10 because it looks broken
mData <- mData[mData$staID!=10,]

# calculate GDP growth
if(is.null(mData$GDPd1)){
	mData$GDPd1 <- mData$GDP - mData$gdplag1
} 
if(is.null(mData$gdpGrRt)){
	mData$gdpGrRt <- mData$GDPd1/mData$gdplag1
}

# data subsets ####
mDataFull <- mData
dataSubsets <- list()
dataSubsets$allYears <- mDataFull
dataSubsets$fiftyOnwards <- mDataFull[mDataFull$year>2050,]
dataSubsets$seventyOnwards <- mDataFull[mDataFull$year>2070,]
dataSubsets$twentyOnwards <- mDataFull[mDataFull$year>2020,]

stats.ids <- unique(mData$staID)
numSTAts <- length(stats.ids)

# plot data ####
fig.dir <- fig.dir.orig <- 'figures/forcedSTAtsMedians'
dir.create(fig.dir,F,T)
fig.w <- 15
fig.h <- 15
fig.res <- 600
fig.u <- 'cm'

staCols <- rainbow(numSTAts*3+1)[(numSTAts*2+1):(numSTAts*3)]

addSTAtsColLegend <- function(xleft,ybottom,xright,ytop){
	oldXpd <- par('xpd')
	ystep <- (ytop-ybottom)/(1+numSTAts)
	vertBuff <- ystep/2
	xwidth <- xright-xleft
	par(xpd=T)
	rect(xleft,ybottom,xright,ytop,xpd=T,col='white')
	for(i in 1:numSTAts){
		points(xleft+xwidth*0.2,
					 ybottom+vertBuff+ystep*(i-1),
					 pch=20,col=staCols[i])
	}
	for(i in 1:numSTAts){
		text(xleft+xwidth*0.3,ybottom+vertBuff+ystep*(i-1),
				 sprintf(' %i',stats.ids[i]),
				 adj=0)
	}
	text(xleft+xwidth*0.1,ytop-vertBuff,'STAts',adj=0)
	rect(xleft,ybottom,xright,ytop,xpd=T)
	par(xpd=oldXpd)
}
for(dataSubsetName in names(dataSubsets)){
	cat(sprintf('\nPlotting dataset %s:\n',dataSubsetName))
	fig.dir <- file.path(fig.dir.orig,dataSubsetName)
	dir.create(fig.dir,F,T)
	mData <- dataSubsets[[dataSubsetName]]
	## GDP ####
	cat('plotting GDP...')
	png(file.path(fig.dir,'1-GDP.png'),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
	plot(0,0,type='n',
			 xlim=range(mData$year),ylim=c(0,2e6),
			 xlab='Year',ylab='billion constant 2021 $',
			 main='Real GDP')
	# points(medianData$year,medianData$GDP, pch=20,
	# 			 col=staCols[medianData$staID])
	for(stats.id in stats.ids){
		lines(mData$year[mData$staID==stats.id],
					mData$GDP[mData$staID==stats.id],
					col=staCols[stats.id],
					lwd=4)
	}
	legend('topleft',legend = paste('STAts',stats.ids),
				 pch=20,col=staCols[stats.ids])
	dev.off()
	cat('done\n')
	
	## GDP Growth ####
	cat('plotting GDP Growth...')
	png(file.path(fig.dir,'1-GDPGrowth.png'),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
	plot(0,0,type='n',
			 xlim=range(mData$year),ylim=c(-1e3,1.2e4),
			 xlab='Year',ylab='billion constant 2021 $',
			 main='Real GDP Growth')
	# points(medianData$year,medianData$GDP, pch=20,
	# 			 col=staCols[medianData$staID])
	for(stats.id in stats.ids){
		lines(mData$year[mData$staID==stats.id],
					mData$GDPd1[mData$staID==stats.id],
					col=staCols[stats.id],
					lwd=4)
	}
	legend('topleft',legend = paste('STAts',stats.ids),
				 pch=20,col=staCols[stats.ids])
	dev.off()
	cat('done\n')
	
	## GDP Growth Rate ####
	cat('plotting GDP Growth Rate...')
	png(file.path(fig.dir,'1-GDPGrowthRate.png'),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
	plot(0,0,type='n',
			 xlim=range(mData$year),ylim=c(-0.01,0.06),
			 xlab='Year',ylab='Rate',
			 main='Real GDP Growth Rate')
	# points(medianData$year,medianData$GDP, pch=20,
	# 			 col=staCols[medianData$staID])
	for(stats.id in stats.ids){
		lines(mData$year[mData$staID==stats.id],
					mData$gdpGrRt[mData$staID==stats.id],
					col=staCols[stats.id],
					lwd=4)
	}
	legend('topright',legend = paste('STAts',stats.ids),
				 pch=20,col=staCols[stats.ids])
	dev.off()
	cat('done\n')
	
	## STA ####
	cat('plotting STA...')
	png(file.path(fig.dir,'2-STA.png'),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
	plot(0,0,type='n',
			 xlim=range(mData$year),ylim=c(0,7),
			 xlab='Year',ylab='K',
			 main='Surface Temperature Anomaly')
	# points(medianData$year,medianData$STA, pch=20,
	# 			 col=staCols[medianData$staID])
	for(stats.id in stats.ids){
		lines(mData$year[mData$staID==stats.id],
					mData$STA[mData$staID==stats.id],
					col=staCols[stats.id],
					lwd=4)
	}
	legend('topleft',legend = paste('STAts',stats.ids),
				 pch=20,col=staCols[stats.ids])
	dev.off()
	cat('done\n')
	
	## STA vs GDP ####
	cat('plotting STA vs GDP...')
	png(file.path(fig.dir,'3-GDPvsSTA.png'),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
	plot(mData$STA,mData$GDP,type='p',pch=20,
			 col=staCols[mData$staID],
			 xlim=c(0,max(mData$STA,na.rm=T)),ylim=c(0,2e6),
			 xlab='Surface Temperature Anomaly K',ylab='GDP billion constant 2021 $',
			 main='GDP vs STA')
	legend('topleft',legend = paste('STAts',stats.ids),
				 pch=20,col=staCols[stats.ids])
	dev.off()
	cat('done\n')
	
	## STA vs GDP Growth ####
	cat('plotting STA vs GDP Growth...')
	png(file.path(fig.dir,'3-GDPd1vsSTA.png'),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
	plot(mData$STA,mData$GDPd1,type='p',pch=20,
			 col=staCols[mData$staID],
			 xlim=c(0,max(mData$STA,na.rm=T)),ylim=c(-1e3,1.2e4),
			 xlab='Surface Temperature Anomaly K',ylab='GDP billion constant 2021 $',
			 main='GDPd1 vs STA')
	legend('topleft',legend = paste('STAts',stats.ids),
				 pch=20,col=staCols[stats.ids])
	dev.off()
	cat('done\n')
	
	## STA vs GDP Growth Rate ####
	cat('plotting STA vs GDP Growth...')
	png(file.path(fig.dir,'3-GDPGrRtvsSTA.png'),width = fig.w,height = fig.h,units = fig.u,res = fig.res)
	plot(mData$STA,mData$gdpGrRt,type='p',pch=20,
			 col=staCols[mData$staID],
			 xlim=c(0,max(mData$STA,na.rm=T)),ylim=c(-0.01,0.06),
			 xlab='Surface Temperature Anomaly K',ylab='Rate',
			 main='GDP Growth Rate vs STA')
	legend('topright',legend = paste('STAts',stats.ids),
				 pch=20,col=staCols[stats.ids])
	dev.off()
	cat('done\n')
}
fig.dir <- fig.dir.orig

## STA vs GDP vs Time ####
library(rgl)
plot3d(mData$year,mData$STA,mData$GDP,
			 col=staCols[mData$staID],
			 xlab='Year',ylab='STA (K)',
			 zlab='GDP billion constant 2021 $')

# model fitting ####
mDataCF <- mDataFull
mDataCF[,which(grepl('sta',colnames(mDataCF),ignore.case = T))] <- 0*mDataCF[,which(grepl('sta',colnames(mDataCF),ignore.case = T))] 
fits <- list()


# fit functions ####
makePredAndDF <- function(fit,data=NULL){
	if(is.null(data)){
		data <- fit$model
	}
	dataCF <- data
	dataCF[,which(grepl('sta',colnames(dataCF),ignore.case = T))] <- 0*dataCF[,which(grepl('sta',colnames(dataCF),ignore.case = T))] 
	retList <- list()
	retList$fit <- fit
	retList$pred <- predict(retList$fit,newdata=data)
	retList$predCF <- predict(retList$fit,newdata = dataCF)
	if(colnames(fit$model)[1]=='GDPd1') {
		retList$df <- (retList$predCF - retList$pred)/(retList$predCF+data$gdplag1)
		retList$df2 <- (retList$predCF - retList$pred)/retList$predCF
		retList$df3 <- (retList$predCF - retList$pred)
	} else if (colnames(fit$model)[1]=='gdpGrRt') {
		retList$df <- (fit$model$gdplag1*(retList$predCF - retList$pred))/(fit$model$gdplag1*(1+retList$predCF))
		retList$df2 <- (retList$predCF - retList$pred)/retList$predCF
		retList$df3 <- (retList$predCF - retList$pred)
	} else {
		retList$df <-  (retList$predCF - retList$pred)/retList$predCF
		retList$df2 <-  (retList$predCF - retList$pred)/retList$predCF
		retList$df3 <-  (retList$predCF - retList$pred)
	}
	retList$df1 <- retList$df 
	retList$data <- data
	retList$dataCF <- dataCF
	retList$resid <- retList$pred-data[[colnames(fit$model)[1]]]
	return(retList)
}

surfVals <- list()
gdpRange <- c(0,1.5e6)
gdpSurfCount <- 200
surfVals$gdplag1 <- surfVals$GDP <- seq(gdpRange[1],gdpRange[2],length.out=gdpSurfCount)
gdpGrRtRange <- c(-0.01,0.06)
surfVals$gdpGrRt <- seq(gdpGrRtRange[1],gdpGrRtRange[2],length.out=gdpSurfCount)
staRange <- c(0,6)
staSurfCount <- 200
surfVals$STA <- seq(staRange[1],staRange[2],length.out=staSurfCount)
dataPredSurf <- data.frame(STA=rep(surfVals$STA,each=gdpSurfCount),
													 gdplag1=rep(surfVals$GDP,staSurfCount))
dataPredSurfCF <- data.frame(STA=0,
													 gdplag1=rep(surfVals$GDP,staSurfCount))


predSurface <- function(fit,newdata){
	surfVals <- predict(fit,newdata=newdata)
	surfValsMat <- matrix(surfVals,nrow=gdpSurfCount,byrow = T)
	# surfValsMat.testSTA <- matrix(dataPredSurf$STA,nrow=gdpSurfCount,byrow = T)
	# surfValsMat.testGDP <- matrix(dataPredSurf$gdplag1,nrow=gdpSurfCount,byrow = T)
	return(surfValsMat)
}

fitPlot3d <- function(fit,yvar,zvar,zvarSurfVals){
	plot3d(mData$STA,mData[[yvar]],mData[[zvar]],
				 col=staCols[mData$staID],
				 xlab='STA',ylab=yvar,
				 zlab=zvar)
	# points3d(medianData$year,medianData$STA,fits$ols2$pred)
	surface3d(x=surfVals$STA,y=surfVals[[yvar]],z=predSurface(fit,newdata=dataPredSurf),
						front='lines',back='lines')
	surface3d(x=surfVals$STA,y=surfVals[[yvar]],z=predSurface(fit,newdata=dataPredSurfCF),
						col='red',
						front='lines',back='lines')
}

fitPlot <- function(fitAndDataList,pch=20,dfToUse=1){
	# par(mfrow=c(3,3),pch=pch,mar=c(4,4,1,1),mgp=c(2.5,1,0))
	layout(mat = matrix(c(1,1,1,1,
												2,3,4,5,
												6,7,8,9,
												10,11,12,13,
												14,15,16,17),nrow=5,byrow=T),
				 widths = c(0.1,1,1,1),heights = c(0.1,1,1,1,1))
	par(mar=c(0,0,0,0))
	plot(0,0,type='n',axes=F,xlab='',ylab='')
	text(0,0,paste('Empirical model: ',
								 gsub('~',' explained by ',as.character(fitAndDataList$fit$call)[2]),
								 'using dataset',mdataSubsetName),
			 adj=c(0.5,0.5),
			 cex=1.8)
	
	plot(0,0,type='n',axes=F,xlab='',ylab='')
	text(0,0,'FRIDA',srt=90,adj=c(0.5,0.5))
	par(pch=pch,mar=c(4,4,1,1),mgp=c(2.5,1,0))
	
	if(colnames(fitAndDataList$fit$model)[1]=='GDP'){
		plot(fitAndDataList$data$year,fitAndDataList$data$GDP,
				 xlab='Year',ylab='GDP')
	} else if(colnames(fitAndDataList$fit$model)[1]=='gdpGrRt') {
		plot(fitAndDataList$data$year,fitAndDataList$data$gdpGrRt,
		 xlab='Year',ylab='GDP growth rate')
	} else if(colnames(fitAndDataList$fit$model)[1]=='GDPd1') {
		plot(fitAndDataList$data$year,fitAndDataList$data$GDPd1,
				 xlab='Year',ylab='GDP growth')
	}
	plot(fitAndDataList$data$year,fitAndDataList$data$gdplag1,
			 xlab='Year',ylab='GDP lag 1')
	plot(fitAndDataList$data$year,fitAndDataList$data$STA,
			 xlab='Year',ylab='STA')
	
	par(mar=c(0,0,0,0))
	plot(0,0,type='n',axes=F,xlab='',ylab='')
	text(0,0,'RESIDUALS',srt=90,adj=c(0.5,0.5))
	par(pch=pch,mar=c(4,4,1,1),mgp=c(2.5,1,0))
	plot(fitAndDataList$data$year,fitAndDataList$resid,
			 xlab='Year',ylab='resid',col='red')
	abline(h=0, col='gray')
	plot(fitAndDataList$data$gdplag1,fitAndDataList$resid,
			 xlab='gdp lag 1',ylab='resid',col='red')
	abline(h=0, col='gray')
	plot(fitAndDataList$data$STA,fitAndDataList$resid,
			 xlab='STA',ylab='resid',col='red')
	abline(h=0, col='gray')
	
	par(mar=c(0,0,0,0))
	plot(0,0,type='n',axes=F,xlab='',ylab='')
	text(0,0,'Comparing FRIDA to Empirical Mod.',srt=90,adj=c(0.5,0.5))
	par(pch=pch,mar=c(4,4,1,1),mgp=c(2.5,1,0))
	if(colnames(fitAndDataList$fit$model)[1]=='GDP'){
		plot(fitAndDataList$data$gdplag1,fitAndDataList$data$GDP,
				 xlab='GDP lag 1',ylab='GDP')
		points(fitAndDataList$data$gdplag1,fitAndDataList$pred,col='red')
		points(fitAndDataList$data$gdplag1,fitAndDataList$predCF,col='blue')
		plot(fitAndDataList$data$STA,fitAndDataList$data$GDP,
				 xlab='STA',ylab='GDP')
		points(fitAndDataList$data$STA,fitAndDataList$pred,col='red')
		points(fitAndDataList$data$STA,fitAndDataList$predCF,col='blue')
		
	} else if(colnames(fitAndDataList$fit$model)[1]=='gdpGrRt') {
		plot(fitAndDataList$data$gdplag1,fitAndDataList$data$gdpGrRt,
				 xlab='GDP lag 1',ylab='GDP Growth Rate')
		points(fitAndDataList$data$gdplag1,fitAndDataList$pred,col='red')
		points(fitAndDataList$data$gdplag1,fitAndDataList$predCF,col='blue')
		plot(fitAndDataList$data$STA,fitAndDataList$data$gdpGrRt,
				 xlab='STA',ylab='GDP Growth Rate')
		points(fitAndDataList$data$STA,fitAndDataList$pred,col='red')
		points(fitAndDataList$data$STA,fitAndDataList$predCF,col='blue')
	} else if(colnames(fitAndDataList$fit$model)[1]=='GDPd1') {
		plot(fitAndDataList$data$gdplag1,fitAndDataList$data$GDPd1,
				 xlab='GDP lag 1',ylab='GDP Growth')
		points(fitAndDataList$data$gdplag1,fitAndDataList$pred,col='red')
		points(fitAndDataList$data$gdplag1,fitAndDataList$predCF,col='blue')
		plot(fitAndDataList$data$STA,fitAndDataList$data$GDPd1,
				 xlab='STA',ylab='GDP Growth')
		points(fitAndDataList$data$STA,fitAndDataList$pred,col='red')
		points(fitAndDataList$data$STA,fitAndDataList$predCF,col='blue')
	}
	if(colnames(fitAndDataList$fit$model)[1]=='GDP'){
		plot(fitAndDataList$data$year,fitAndDataList$data$GDP,
				 xlab='Year',ylab='GDP')
	} else if(colnames(fitAndDataList$fit$model)[1]=='gdpGrRt') {
		plot(fitAndDataList$data$year,fitAndDataList$data$gdpGrRt,
				 xlab='Year',ylab='GDP growth rate')
	} else if(colnames(fitAndDataList$fit$model)[1]=='GDPd1') {
		plot(fitAndDataList$data$year,fitAndDataList$data$GDPd1,
				 xlab='Year',ylab='GDP growth')
	}
	points(fitAndDataList$data$year,fitAndDataList$pred,
				 col='red')
	points(fitAndDataList$data$year,fitAndDataList$predCF,
				 col='blue')
	
	
	par(mar=c(0,0,0,0))
	plot(0,0,type='n',axes=F,xlab='',ylab='')
	text(0,0,'DAMAGE FACTORS',srt=90,adj=c(0.5,0.5))
	par(pch=pch,mar=c(4,4,1,1),mgp=c(2.5,1,0))
	if(!is.null(fitAndDataList[[paste0('df',dfToUse)]])){
	plot(fitAndDataList$data$STA,fitAndDataList[[paste0('df',dfToUse)]],ylim=c(-0.05,1),
			 xlab='STA',ylab='damage factor',col='purple',xlim=c(0,7))	
	plot(fitAndDataList$data$STA,fitAndDataList[[paste0('df',dfToUse)]],
			 xlab='STA',ylab='damage factor',col='purple')	
	} else {
		plot(0,0)
		plot(0,0)
	}
	plot(0,0,type='n',xlab='',ylab='',axes=F)
	legend('center',legend=c('frida','pred','pred CF'),col=c('black','red','blue'),
				 pch=20,cex=2)
}

## OLS ####
mdataSubsetName <- 'fiftyOnwards'
mdataSubset <- dataSubsets[[mdataSubsetName]]
dfsToUse <- 1:3
### GDP 1 ####
fits$ols1 <- makePredAndDF(lm(GDP~STA+I(STA^2)+gdplag1,data=mdataSubset),mdataSubset)
sink(file.path(fig.dir,'OLS-GDP1.txt'))
summary(fits$ols1$fit)
sink()
png(file.path(fig.dir,'OLS-GDP1.png'),width = fig.w*3*0.8,height = fig.h*4*0.8,units = fig.u,res = fig.res*0.25)
fitPlot(fits$ols1)
dev.off()
# fitPlot3d(fits$ols1$fit,'gdplag1','GDP',gdpSurfVals)

### GDP 4 ####
fits$ols4 <- makePredAndDF(lm(GDP~STA+I(STA^2)+gdplag1+gdplag2,data=mdataSubset),mdataSubset)
sink(file.path(fig.dir,'OLS-GDP4.txt'))
summary(fits$ols4$fit)
sink()
png(file.path(fig.dir,'OLS-GDP4.png'),width = fig.w*3*0.8,height = fig.h*4*0.8,units = fig.u,res = fig.res*0.25)
fitPlot(fits$ols4)
dev.off()

### GDP 2 ####
fits$ols2 <- makePredAndDF(lm(GDP~STA+I(STA^2)+I(STA*gdplag1)+gdplag1+I(gdplag1^2),data=mdataSubset),mdataSubset)
sink(file.path(fig.dir,'OLS-GDP2.txt'))
summary(fits$ols2$fit)
sink()
png(file.path(fig.dir,'OLS-GDP2.png'),width = fig.w*3*0.8,height = fig.h*4*0.8,units = fig.u,res = fig.res*0.25)
fitPlot(fits$ols2)
dev.off()
# fitPlot3d(fits$ols2$fit,'gdplag1','GDP',gdpSurfVals)

### GDP 3 time ####
fits$ols3 <- makePredAndDF(lm(GDP~STA+I(STA^2)+year+I(year^2),data=mdataSubset),mdataSubset)
sink(file.path(fig.dir,'OLS-GDP3.txt'))
summary(fits$ols3$fit)
sink()
png(file.path(fig.dir,'OLS-GDP3.png'),width = fig.w*3*0.8,height = fig.h*4*0.8,units = fig.u,res = fig.res*0.25)
fitPlot(fits$ols3)
dev.off()
# fitPlot3d(fits$ols2$fit,'gdplag1','GDP',gdpSurfVals)

### GDP d1 1 ####
fits$olsGDPd11 <- makePredAndDF(lm(GDPd1~STA+I(STA^2)+gdplag1,data=mdataSubset),mdataSubset)
sink(file.path(fig.dir,'OLS-GDPd11.txt'))
summary(fits$olsGDPd11$fit)
sink()
for(dfToUse in dfsToUse){
	png(file.path(fig.dir,paste0('OLS-GDPd11-df',dfToUse,'.png')),width = fig.w*3*0.8,height = fig.h*4*0.8,units = fig.u,res = fig.res*0.25)
	fitPlot(fits$olsGDPd11,dfToUse=dfToUse)
	dev.off()
}

# fitPlot3d(fits$olsGDPd11$fit,'gdplag1','GDP',gdpSurfVals)
### GDP d1 2 ####
fits$olsGDPd12 <- makePredAndDF(lm(GDPd1~STA+I(STA^2)+year+I(year^2),
																	 # +I(sin(year*2*pi/13))+I(cos(year*2*pi/13)),
																	 data=mdataSubset),mdataSubset)
sink(file.path(fig.dir,'OLS-GDPd12.txt'))
summary(fits$olsGDPd12$fit)
sink()
for(dfToUse in dfsToUse){
	png(file.path(fig.dir,paste0('OLS-GDPd12-df',dfToUse,'.png')),width = fig.w*3*0.8,height = fig.h*4*0.8,units = fig.u,res = fig.res*0.25)
	fitPlot(fits$olsGDPd12,dfToUse = dfToUse)
	dev.off()
}
# fitPlot3d(fits$olsGDPd11$fit,'gdplag1','GDP',gdpSurfVals)

### GDP Gr Rt 1 ####
fits$olsGrRt1 <- makePredAndDF(lm(gdpGrRt~
																 	STA+I(STA^2)+gdplag1,data=mdataSubset),mdataSubset)
sink(file.path(fig.dir,'OLS-GDPGrRt1.txt'))
summary(fits$olsGrRt1$fit)
sink()
for(dfToUse in dfsToUse){
	png(file.path(fig.dir,paste0('OLS-GDPGrRt1-df',dfToUse,'.png')),width = fig.w*3*0.8,height = fig.h*4*0.8,units = fig.u,res = fig.res*0.25)
	fitPlot(fits$olsGrRt1,dfToUse = dfToUse)
	dev.off()
}
# fitPlot3d(fits$olsGrRt1$fit,'gdplag1','gdpGrRt',gdpGrRtSurfVals)

### GDP Gr Rt 4 ####
fits$olsGrRt4 <- makePredAndDF(lm(gdpGrRt~
																 	STA+I(STA^2)+gdplag1+year+I(year^2),data=mdataSubset),mdataSubset)
sink(file.path(fig.dir,'OLS-GDPGrRt4.txt'))
summary(fits$olsGrRt4$fit)
sink()
for(dfToUse in dfsToUse){
	png(file.path(fig.dir,paste0('OLS-GDPGrRt4-df',dfToUse,'.png')),width = fig.w*3*0.8,height = fig.h*4*0.8,units = fig.u,res = fig.res*0.25)
	fitPlot(fits$olsGrRt4,dfToUse = dfToUse)
	dev.off()
}
# fitPlot3d(fits$olsGrRt4$fit,'gdplag1','gdpGrRt',gdpGrRtSurfVals)

### GDP Gr Rt 2 ####
fits$olsGrRt <- makePredAndDF(lm(gdpGrRt~
																 	STA+I(STA^2)+gdplag1+I(gdplag1^2)+
																 	stalag10+I(stalag10^2)+
																 	stalag20+I(stalag20^2),
																 data=mdataSubset),mdataSubset)
summary(fits$olsGrRt$fit)
fitPlot(fits$olsGrRt)
fitPlot3d(fits$olsGrRt$fit,'gdplag1','gdpGrRt',gdpGrRtSurfVals)
### GDP Gr Rt 3 time ####
fits$olsGrRt3 <- makePredAndDF(lm(gdpGrRt~
																 	STA+I(STA^2)+year+I(year^2),
																 data=mdataSubset),mdataSubset)
sink(file.path(fig.dir,'OLS-GDPGrRt3.txt'))
summary(fits$olsGrRt3$fit)
sink()
png(file.path(fig.dir,'OLS-GDPGrRt3-df1.png'),width = fig.w*3*0.8,height = fig.h*4*0.8,units = fig.u,res = fig.res*0.25)
fitPlot(fits$olsGrRt3)
dev.off()
png(file.path(fig.dir,'OLS-GDPGrRt3-df2.png'),width = fig.w*3*0.8,height = fig.h*4*0.8,units = fig.u,res = fig.res*0.25)
fitPlot(fits$olsGrRt3,dfToUse = 2)
dev.off()
png(file.path(fig.dir,'OLS-GDPGrRt3-df3.png'),width = fig.w*3*0.8,height = fig.h*4*0.8,units = fig.u,res = fig.res*0.25)
fitPlot(fits$olsGrRt3,dfToUse = 3)
dev.off()
# fitPlot3d(fits$olsGrRt3$fit,'gdplag1','gdpGrRt',gdpGrRtSurfVals)



## polynomial with interactions, id ####
degree <-6
lags <- c(1,2,3,4,5)#,2,3)#,4,5,10,20) # max c(1,2,3,4,5,10,20)
modelName <- paste0('fit.ols.d',degree,'.l',paste(lags,collapse=','))
fitStr <- paste0('fitModel <- lm(GDP~')
for(year.d in 1:degree){
	fitStr <- paste0(fitStr,'+I(year')
	if(year.d>1){
		fitStr <- paste0(fitStr,'^',year.d)
	}
	fitStr <- paste0(fitStr,')')
}
# for(gdp.l in lags){
# 	for(gdp.d in 1:degree){
# 		fitStr <- paste0(fitStr,'+I(gdplag',gdp.l)
# 		if(gdp.d>1){
# 			fitStr <- paste0(fitStr,'^',gdp.d)
# 		}
# 		fitStr <- paste0(fitStr,')')
# 	}
# }
for(sta.l in lags){
	for(sta.d in 1:degree){
		fitStr <- paste0(fitStr,'+I(stalag',sta.l)
		if(sta.d>1){
			fitStr <- paste0(fitStr,'^',sta.d)
		}
		fitStr <- paste0(fitStr,')')
	}
}
# for(gdp.l in lags){
# 	for(gdp.d in 1:degree){
# 		for(sta.l in lags){
# 			for(sta.d in 1:degree){
# 				if(gdp.d+sta.d<=degree){
# 					fitStr <- paste0(fitStr,'+I(stalag',sta.l)
# 					if(sta.d>1){
# 						fitStr <- paste0(fitStr,'^',sta.d)
# 					}
# 					fitStr <- paste0(fitStr,'*gdplag',gdp.l)
# 					if(gdp.d>1){
# 						fitStr <- paste0(fitStr,'^',gdp.d)
# 					}
# 					fitStr <- paste0(fitStr,')')
# 				}
# 			}
# 		}
# 	}
# }
fitStr <- paste0(fitStr,', data=medianData)')
eval(parse(text=fitStr))
# sink(file.path(fig.dir,paste0(modelName,'-summary.txt')),append = F)
summary(fitModel)
# sink()
fits[[modelName]] <- makePredAndDF(fitModel)
plot3d(mData$year,mData$STA,fits[[modelName]]$fit$residuals)
par(mfrow=c(2,2))
par(mar=c(3,3,1,1))
par(pch='.')
plot(mData$STA,fits[[modelName]]$df,ylim=c(0,1))
plot(mData$STA,fits[[modelName]]$df)
plot(mData$year,fits[[modelName]]$pred,ylim=c(0,2e6))
points(mData$year,fits[[modelName]]$predCF,col='red')
legend('topleft',c('pred','predCF'),col=c('black','red'),pch=par()$pch)
plot(mData$year,fits[[modelName]]$predCF,col='red')
points(mData$year,fits[[modelName]]$pred)

