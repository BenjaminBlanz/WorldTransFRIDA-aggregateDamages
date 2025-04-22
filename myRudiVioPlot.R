library(spatstat)
# rudimentary violinplot
# x the data to be plotted
# w the weights of the data
# n number of segments
# area area filled by the polygon
# at the center line of the polygon
myRudiViolinPlot <- function(x,w=rep(1/length(x),length(x)),n=11,area=0.1,width=NULL,at=0,horz=F,
														 equiprobspacing=F,add=T,lim=NULL,col=gray(0.1),border.col=NA,
														 whiskerExt = 1e-4,...){
	x <- x[!is.na(x)]
	w <- w[!is.na(x)]
	if(!add){
		plot(rep(0,length(x)),x,type='n',xlim=c(at-1,at+1))
	}
	if(is.null(lim)){
		lim <- c(par('usr')[3]+0.01,par('usr')[4]-0.01)
	}
	w <- w[x>lim[1]&x<lim[2]]
	x <- x[x>lim[1]&x<lim[2]]
	x.ewcdf <- ewcdf(x,w)
	if(equiprobspacing){
		qs <- unique(quantile(x.ewcdf,seq(whiskerExt,1-whiskerExt,length.out=n)))
		n <- length(qs)
		qwidths <- area/n/diff(qs)
	} else { #equidistant spacing
		minval <- floor(quantile(x.ewcdf,whiskerExt))
		maxval <- ceiling(quantile(x.ewcdf,1-whiskerExt))
		qs <- seq(min(x),max(x),l=n)
		names(qs) <- seq(0,1,l=n)
		percentiles <- quantile(x.ewcdf,seq(0,1,length.out=10001))
		qwidths <- rep(NA,n-1)
		for(i in 1:(length(qs)-1)){
			qwidths[i] <- n * area *
				((max(which(percentiles<=qs[i+1]))-1) - (max(0,which(percentiles<qs[i]))-1))/10000
		}
	}
	if(!is.null(width)){
		qwidths <- width*qwidths/max(qwidths)
	} else {
		width <- area*n
	}
	qs <- pmax(qs,lim[1])
	qs <- pmin(qs,lim[2])
	qwidths <- pmax(qwidths,0)
	qwidths <- pmin(qwidths,width)
	qRange <- weighted.quantile(x,w,probs = c(whiskerExt,1-whiskerExt))
	if(horz){
		qwidths <- qs
		qs <-  area/n/diff(qwidths)
		lines(qRange,rep(at,2),col=col,lwd=0.5)
	} else {
		lines(rep(at,2),qRange,col=col,lwd=0.5)
	}
	polygon(x=c(rep(at+0.5*qwidths,each=2),rev(rep(at-0.5*qwidths,each=2))),
					y=c(qs[1],rep(qs[2:n],each=2),rev(rep(qs[c(-1,-n)],each=2)),qs[1]),
					border=border.col,col=col,
					...)
	if(horz){
		qs <- qwidths
	}
	return(qs)
}
