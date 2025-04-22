source('funRunFRIDA.R')
location.frida <- 'FRIDAforUncertaintyAnalysis'

writeSTAForcing <- function(location.frida,
														STAOverride,aggregateDamageFormationTime=1){
	sink(file.path(location.frida,'Data','ClimateSTAOverride.csv'),
			 append = F)
	cat('Energy Balance Model.SWITCH STA Override,1\n')
	cat(sprintf('Energy Balance Model.STA Override,%f\n',
							 STAOverride))
	cat(sprintf('Energy Balance Model.Aggregate Damage Formation Time,%f\n',
							 aggregateDamageFormationTime))
	sink()
}

# The first STA has to be 0 for the model code below to work
STAs <- c(0,1,2,3)

# run forced ensemble runs for all the STAs and then proceed
# TODO: automate this step


# read STA data ####
dataForDamFac <- data.frame(id=numeric(),y=numeric(),ylag=numeric(),yfut=numeric(),year=numeric(),STA=numeric())
dataForDamFacSTAs.lst <- list()
for(STAlevel.i in 1:length(STAs)){
	cat(sprintf('Processing STA level %f\n',STAs[STAlevel.i]))
	gdpForThisSTA <- readRDS(file.path('workOutput',
																		 paste0('N-20000-ChS-100-LCR-1000-IgB-FALSE-FrB-FALSE-KcE-FALSE-Sym-Min-AAZ-FALSE-CFB-TRUE-Pol-policy_EMB-CTO-',STAs[STAlevel.i]),
																		 'detectedParmSpace','PerVarFiles-RDS','gdp_real_gdp_in_2021c.RDS'))
	# loglikeForThisSTA <- readRDS(file.path('workOutput',
	# 																			 paste0('N-20000-ChS-100-LCR-1000-IgB-FALSE-FrB-FALSE-KcE-FALSE-Sym-Min-AAZ-FALSE-CFB-TRUE-Pol-policy_EMB-CTO-',STAs[STAlevel.i]),
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
