source('config-DetermineAggregateDamageFunction.R')
configStr <- paste0('S-',numSample,'-nSTAs-',length(STAs),'-nTW-',length(timeWindows),'-',identifier)
if(!exists('dataForDamFac')){
	cat('Reading dataForDamFac...')
	dataForDamFac <- readRDS(file.path('outputData',paste0('dataForDamFac-',configStr,'.RDS')))
	cat('done\n')
}
if(!exists('dataForDamFacAgg')){
	cat('Reading dataForDamFacAgg...')
	dataForDamFacAgg <- readRDS(file.path('outputData',paste0('dataForDamFacAgg-',configStr,'.RDS')))
	cat('done\n')
}

# plots ####
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


figDir <- file.path('figures','ipccDamages',configStr)
dir.create(figDir,F,T)
source('myRudiVioPlot.R')
plotTypes <- list()
plotTypes[['boxplotVioplot']] <- c('boxplot','vioplot')
plotTypes[['boxplot']] <- c('boxplot')
plotTypes[['vioplot']] <- c('vioplot')
timeWindows[[length(timeWindows)+1]] <- 'all'

# per year per member ####
for(tw.i in 1:length(timeWindows)){
	for(pt.i  in 1:length(plotTypes)){
		cat(sprintf('Plotting per year per ensemble member for time window %s %s',
								paste(timeWindows[[tw.i]],collapse = '-'),
								paste(plotTypes[[pt.i]],collapse=' and ')))
		png(file.path(figDir,sprintf('ipccDamgeFunction-PerYearPerEnsembleMember-%s-%s.png',
																 paste(plotTypes[[pt.i]],collapse='-'),
																 paste(timeWindows[[tw.i]],collapse = '-'))),width='960',height='1380')
		plot(0,0,
				 xlab='STA degC',ylab='Annual percentage loss of GDP',
				 xlim=c(0, 7),
				 ylim=c(-10,80),
				 type='n',yaxs='i',
				 main=sprintf('Per year damages per ensemble member years %s',paste(timeWindows[[tw.i]],collapse = '-')))
		cat('.')
		grid()	
		cat('.')
		# points(dat$STA,dat$yRelLoss*100,col=1,pch=20)
		for(STA in STAs){
			if(length(timeWindows[[tw.i]])==1 && timeWindows[[tw.i]]=='all'){
				if('vioplot' %in% plotTypes[[pt.i]]){
					myRudiViolinPlot(dataForDamFac$yRelLoss[dataForDamFac$STA==STA]*100,
													 at=STA,col='white',add=T,area=vioplot.area,border='black',lwd=3,
													 equiprobspacing = T,n=100)
					cat('.')
				}
				if('boxplot' %in% plotTypes[[pt.i]]){
					boxplot(dataForDamFac$yRelLoss[dataForDamFac$STA==STA]*100,
									at=STA,col='white',add=T,boxwex=boxplot.boxwex,axes=F,range = 0,border='black',lwd=2)
					cat('.')
				}
			} else {
				if('vioplot' %in% plotTypes[[pt.i]]){
					myRudiViolinPlot(dataForDamFac$yRelLoss[dataForDamFac$STA==STA & 
																										dataForDamFac$year>=timeWindows[[tw.i]][1] &
																										dataForDamFac$year<=timeWindows[[tw.i]][2]]*100,
													 at=STA,col='white',add=T,area=vioplot.area,border='black',lwd=3,
													 equiprobspacing = T,n=100)
					cat('.')
				}
				if('boxplot' %in% plotTypes[[pt.i]]){
					boxplot(dataForDamFac$yRelLoss[dataForDamFac$STA==STA & 
																				 	dataForDamFac$year>=timeWindows[[tw.i]][1] &
																				 	dataForDamFac$year<=timeWindows[[tw.i]][2]]*100,
									at=STA,col='white',add=T,boxwex=boxplot.boxwex,axes=F,range = 0,border='black',lwd=2)
					cat('.')
				}
			}
		}
		dev.off()
		cat('done\n')
	}
}

# per ensemble member measures ####
measures <- c('Mean','Median')
for(tw.i in 1:length(timeWindows)){
	for(measure in measures){
		for(pt.i  in 1:length(plotTypes)){
			if(measure %in% names(dataForDamFacAgg)){
				cat(sprintf('Plotting %s yearly per ensemble member for time window %s %s',
										measure,
										paste(timeWindows[[tw.i]],collapse = '-'),
										paste(plotTypes[[pt.i]],collapse=' and ')))
				png(file.path(figDir,sprintf('ipccDamgeFunction-%sPerEnsembleMember-%s-%s.png',
																		 measure,
																		 paste(plotTypes[[pt.i]],collapse='-'),
																		 paste(timeWindows[[tw.i]],collapse = '-'))),width='960',height='1380')
				plot(0,0,
						 xlab='STA degC',ylab='Annual percentage loss of GDP',
						 xlim=c(0, 7),
						 ylim=c(-10,80),
						 type='n',yaxs='i',
						 main=sprintf('%s yearly damages per ensemble member years %s',measure,paste(timeWindows[[tw.i]],collapse = '-')))
				cat('.')
				grid()	
				cat('.')
				# points(dat$STA,dat$yRelLoss*100,col=1,pch=20)
				for(STA in STAs){
					if(length(timeWindows[[tw.i]])==1 && timeWindows[[tw.i]]=='all'){
						if('vioplot' %in% plotTypes[[pt.i]]){
							myRudiViolinPlot(as.vector(dataForDamFacAgg[[measure]][[paste('STA',STA)]])*100,
															 at=STA,col=NA,border.col = 'black',add=T,area=vioplot.area,lwd=2,
															 equiprobspacing = T,n=100)
							cat('.')
						}
						if('boxplot' %in% plotTypes[[pt.i]]){
							boxplot(dataForDamFacAgg[[measure]][[paste('STA',STA)]]*100,
											at=STA,col='white',border='black',add=T,boxwex=boxplot.boxwex,axes=F,range = 0,lwd=2)
							cat('.')
						}
					} else {
						if('vioplot' %in% plotTypes[[pt.i]]){
							myRudiViolinPlot(as.vector(dataForDamFacAgg$timewindow[[tw.i]][[measure]][[paste('STA',STA)]])*100,
															 at=STA,col=NA,border.col = 'black',add=T,area=vioplot.area,lwd=2,
															 equiprobspacing = T,n=100)
							cat('.')
						}
						if('boxplot' %in% plotTypes[[pt.i]]){
							boxplot(dataForDamFacAgg$timewindow[[tw.i]][[measure]][[paste('STA',STA)]]*100,
											at=STA,col='white',border='black',add=T,boxwex=boxplot.boxwex,axes=F,range = 0,lwd=2)
							cat('.')
						}
					}
				}
				dev.off()
				cat('done\n')
			}
		}
	}
}
