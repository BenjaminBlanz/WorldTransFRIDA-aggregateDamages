# this file is used to configure the data preparation for 
# runDetermineAggregateDamageFunction which will then determine the aggregate damage 
# function.
location.fridaUncertaintyWD <- '../WorldTransFrida-Uncertainty'

# the number of samples to use for the damage function
numSample <- 2e4

# The number of sta timeseries to pull from emb for running damages
numSTAts <- 10

# The first STA has to be 0 
STAs <- c(0,1,2,3)

# time windows for calculating mean/median damages
timeWindows <- list() 
timeWindows[[1]] <- c(1981,1990)
timeWindows[[2]] <- c(1991,2000)
timeWindows[[3]] <- c(2001,2010)
timeWindows[[4]] <- c(2011,2020)
timeWindows[[5]] <- c(2011,2030)
timeWindows[[6]] <- c(2011,2040)
timeWindows[[7]] <- c(2011,2050)
timeWindows[[8]] <- c(2011,2060)
timeWindows[[9]] <- c(2011,2070)
timeWindows[[10]] <- c(2011,2080)
timeWindows[[11]] <- c(2011,2090)
timeWindows[[12]] <- c(2011,2100)

# identifier for this config
identifier <- 'default'

# plot properties
vioplot.area <- 0.5
boxplot.boxwex <- 0.1

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
staBreaks <- seq(0,max(dataForDamFac$STA,na.rm=T),length.out=200)
staBreaksCol <- seq(0,max(dataForDamFac$STA,na.rm=T),length.out=length(staCols)+1)
