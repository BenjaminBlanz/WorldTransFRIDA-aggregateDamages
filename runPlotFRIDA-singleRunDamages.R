# a file that has time in the rows and fixed temperatures in columns
y <- read.csv('realGDP.csv',row.names = 1)
# the lagged value of gdp
ylag <- y
ylag[2:nrow(ylag),]<- y[1:(nrow(y)-1),]
ylag[1,] <- NA
# the value of gdp one step ahead
yfut <- y
yfut[1:(nrow(y)-1),] <- y[2:nrow(y),]
yfut[nrow(y),] <- NA
# data frame for the models below
dat <- data.frame(y=unname(unlist(y)),ylag=unname(unlist(ylag)),yfut=unname(unlist(yfut)),
									year=rep(as.numeric(rownames(y)),ncol(y)),STA=rep(0:(ncol(y)-1),each=nrow(y)))
# gdp growth is future value - current
dat$ygro <- dat$yfut - dat$y
plot(dat$y,dat$ygro,col=dat$STA+1,pch=20)

# linear model explaining GDP growth by current gdp, but using only data from STA==0
groMod0 <- lm(ygro~y,data=dat[dat$STA==0,])
summary(groMod0)
# predicted growth if STA were 0
dat$predGro <- predict(groMod0,newdata = data.frame(y=dat$y))
# predicted gdp if STA were 0
dat$predy <- dat$ylag+predict(groMod0,newdata = data.frame(y=dat$ylag))

plot(dat$y[dat$STA==0],dat$ygro[dat$STA==0],pch=20,
		 xlab='GDP',ylab='GDP Growth')
points(dat$y,dat$predGro,pch=20,col='red')
legend('topleft',legend=c('FRIDA 0 deg','Fit'),col=c('black','red'),pch=20)


# loss from climate is difference between predicted (if STA were 0) - actual gdp
dat$yloss <- dat$predy -dat$y
# as relative
dat$yRelLoss <- dat$yloss/dat$predy
yRelLoss <- matrix(dat$yRelLoss,ncol=ncol(y))

plot(dat$STA,dat$yRelLoss*100,col=dat$STA+1,pch=20,
		 xlab='STA degC',ylab='Annual percentage loss of GDP',ylim=c(-10,30),
		 type='n')
grid()	
# points(dat$STA,dat$yRelLoss*100,col=1,pch=20)
boxplot(yRelLoss*100,at=0:(ncol(y)-1),add=T,boxwex=0.15,range=0,axes=F)
