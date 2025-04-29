# this file is used to configure the data preparation for 
# runDetermineAggregateDamageFunction which will then determine the aggregate damage 
# function.
location.fridaUncertaintyWD <- '../WorldTransFrida-Uncertainty'

# the number of samples to use for the damage function
numSample <- 2e4

# The first STA has to be 0 
STAs <- c(0,1,2,3)
