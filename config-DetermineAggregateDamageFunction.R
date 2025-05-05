# this file is used to configure the data preparation for 
# runDetermineAggregateDamageFunction which will then determine the aggregate damage 
# function.
location.fridaUncertaintyWD <- '../WorldTransFrida-Uncertainty'

# the number of samples to use for the damage function
numSample <- 2e4

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
