source('config-DetermineAggregateDamageFunction.R')
dataForDamFac <- readRDS(file.path('outputData',paste0('dataForDamFac',configStr,'.RDS')))
dataForDamFacAgg <- readRDS(file.path('outputData',paste0('dataForDamFacAgg',configStr,'.RDS')))

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

vioplot.area <- 0.5
boxplot.boxwex <- 0.1

figDir <- file.path('figures','ipccDamages')
dir.create(figDir,F,T)
source('myRudiVioPlot.R')
png(file.path(figDir,paste0(format(Sys.time(), "%Y%m%d-%H%M%S"),'-ipccDamgeFunction.png')),width='960',height='1380')
plot(0,0,
		 xlab='STA degC',ylab='Annual percentage loss of GDP',
		 xlim=c(0, 7),
		 ylim=c(-10,80),
		 type='n',yaxs='i',
		 main='Damages per year per ensemble member')
grid()
# points(dat$STA,dat$yRelLoss*100,col=1,pch=20)
for(STA in STAs){
	myRudiViolinPlot(dataForDamFac$yRelLoss[dataForDamFac$STA==STA]*100,
									 at=STA,col='white',add=T,area=vioplot.area,border='black',lwd=3,
									 equiprobspacing = T,n=100)
	boxplot(dataForDamFac$yRelLoss[dataForDamFac$STA==STA]*100,at=STA,col='white',add=T,boxwex=boxplot.boxwex,axes=F,range = 0,border='black',lwd=2)
}
dev.off()
png(file.path(figDir,paste0(format(Sys.time(), "%Y%m%d-%H%M%S"),'-ipccDamgeFunction-boxplot.png')),width='960',height='1380')
plot(0,0,
		 xlab='STA degC',ylab='Annual percentage loss of GDP',
		 xlim=c(0, 7),
		 ylim=c(-10,80),
		 type='n',yaxs='i',
		 main='Damages per year per ensemble member')
grid()
# points(dat$STA,dat$yRelLoss*100,col=1,pch=20)
for(STA in STAs){
	boxplot(dataForDamFac$yRelLoss[dataForDamFac$STA==STA]*100,at=STA,col='white',add=T,boxwex=boxplot.boxwex,axes=F,range = 0,border='black',lwd=2)
}
dev.off()
png(file.path(figDir,paste0(format(Sys.time(), "%Y%m%d-%H%M%S"),'-ipccDamgeFunction-vioplot.png')),width='960',height='1380')
plot(0,0,
		 xlab='STA degC',ylab='Annual percentage loss of GDP',
		 xlim=c(0, 7),
		 ylim=c(-10,80),
		 type='n',yaxs='i',
		 main='Damages per year per ensemble member')
grid()
# points(dat$STA,dat$yRelLoss*100,col=1,pch=20)
for(STA in STAs){
	myRudiViolinPlot(dataForDamFac$yRelLoss[dataForDamFac$STA==STA]*100,
									 at=STA,col='white',add=T,area=vioplot.area,border='black',lwd=3,
									 equiprobspacing = T,n=100)
}
dev.off()
png(file.path(figDir,paste0(format(Sys.time(), "%Y%m%d-%H%M%S"),'-ipccDamgeFunction-meanPerEnsembleMember.png')),width='960',height='1380')
plot(0,0,
		 xlab='STA degC',ylab='Annual percentage loss of GDP',
		 xlim=c(0, 7),
		 ylim=c(-10,80),
		 type='n',yaxs='i',
		 main='Mean yearly damages per ensemble member')
grid()	
# points(dat$STA,dat$yRelLoss*100,col=1,pch=20)
for(STA in STAs){
	myRudiViolinPlot(as.vector(dataForDamFacAgg[[paste('STA',STA)]])*100,
									 at=STA,col=NA,border.col = 'black',add=T,area=vioplot.area,lwd=2,
									 equiprobspacing = T,n=100)
	boxplot(dataForDamFacAgg[[paste('STA',STA)]]*100,at=STA,col='white',border='black',add=T,boxwex=boxplot.boxwex,axes=F,range = 0,lwd=2)
}
dev.off()
png(file.path(figDir,paste0(format(Sys.time(), "%Y%m%d-%H%M%S"),'-ipccDamgeFunction-meanPerEnsembleMember-boxplot.png')),width='960',height='1380')
plot(0,0,
		 xlab='STA degC',ylab='Annual percentage loss of GDP',
		 xlim=c(0, 7),
		 ylim=c(-10,80),
		 type='n',yaxs='i',
		 main='Mean yearly damages per ensemble member')
grid()	
# points(dat$STA,dat$yRelLoss*100,col=1,pch=20)
for(STA in STAs){
	boxplot(dataForDamFacAgg[[paste('STA',STA)]]*100,at=STA,col='white',border='black',add=T,boxwex=boxplot.boxwex,axes=F,range = 0,lwd=2)
}
dev.off()

png(file.path(figDir,paste0(format(Sys.time(), "%Y%m%d-%H%M%S"),'-ipccDamgeFunction-meanPerEnsembleMember-vioplot.png')),width='960',height='1380')
plot(0,0,
		 xlab='STA degC',ylab='Annual percentage loss of GDP',
		 xlim=c(0, 7),
		 ylim=c(-10,80),
		 type='n',yaxs='i',
		 main='Mean yearly damages per ensemble member')
grid()	
# points(dat$STA,dat$yRelLoss*100,col=1,pch=20)
for(STA in STAs){
	myRudiViolinPlot(as.vector(dataForDamFacAgg[[paste('STA',STA)]])*100,
									 at=STA,col=NA,border.col = 'black',add=T,area=vioplot.area,lwd=2,
									 equiprobspacing = T,n=100)
}
dev.off()

cat('done\n')
