source('funRunFRIDA.R')

writeSTAForcing <- function(STAOverride,stopTime,aggregateDamageFormationTime=1){
	sink('FRIDAforAggregateDamageFunction/Data/ClimateSTAOverride.csv',
			 append = F)
	cat('Energy Balance Model.SWITCH STA Override,1\n')
	cat(sprintf('Energy Balance Model.STA Override,%f\n',
							 STAOverride))
	cat(sprintf('Energy Balance Model.Aggregate Damage Formation Time,%f\n',
							 aggregateDamageFormationTime))
	cat(sprintf('stop time,%f',stopTime))
	sink()
}

justRunFRIDA <- function(location.frida='FRIDAforAggregateDamageFunction',
												 name.fridaOutputFile='uncertainty_analysis_exported_variables.csv'){
	system(paste(file.path('Stella_Simulator_Linux','stella_simulator'),'-i','-x','-q',#'-s', #to output isdb
							 file.path(location.frida,'FRIDA.stmx')),
				 ignore.stdout = T,ignore.stderr = T,wait = T)
	runDat <- read.csv(file.path(location.frida,'Data',name.fridaOutputFile))
	colnames(runDat) <- cleanNames(colnames(runDat))
	rownames(runDat) <- runDat$year
	runDat <- runDat[,-1]
	return(runDat)
}
times <- 1982:2130
STAs <- 0:7
dmFData <- array(NA,dim=c(length(times),length(STAs)))
for(time.i in 1:length(times)){
	for(STAlevel.i in 1:length(STAs)){
		cat(sprintf('\rrunning time %i STA %.2f',
				times[time.i],STAs[STAlevel.i]))
		writeSTAForcing(STAs[STAlevel.i],times[time.i],aggregateDamageFormationTime = 2)
		runDat <- justRunFRIDA()
		dmFData[time.i,STAlevel.i] <- runDat[as.character(times[time.i]),'gdp_real_gdp_in_2021c']
	}
}

dmgFac <- array(NA,dim=c(length(times),length(STAs)-1))
for(STA.i in 1:(length(STAs)-1)){
	dmgFac[,STA.i] <- (dmFData[,1]-dmFData[,STA.i+1])/dmFData[,1]
}

par(pch=20)
plot(0,0,
		 xlab='STA degC',ylab='Annual percentage loss of GDP',ylim=c(-10,30),
		 type='n',xlim=range(STAs))
grid()
for(STA.i in 1:(length(STAs)-1)){
	points(rep(STAs[STA.i],length(times)),dmgFac[,STA.i]*100,col=1)
}
boxplot(dmgFac*100,at=0:7,add=T,boxwex=0.4,range=0,axes=F)

