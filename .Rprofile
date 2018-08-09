library(colorout)
setOutputColors256(202, 214, 209, 184, 172, 179, verbose=FALSE)
options(width=160)
library(txtplot)
library(NostalgiR)

mode=function(x){
	dst=density(x)
	return(dst$x[which(dst$y==max(dst$y))])
}

HPD=function(x,P){
	#x=vecteur
	#P=90 si on veut 90% de la distribution
	return(x[which(x>quantile(x,(100-P)/200) & x<quantile(x, 1-(100-P)/200))])
}	

source("~/Programmes/R/twobar.R")

#fonction pour mieux visualiser une relation avec un très grand nombre de points
plotbin=function(a, b, quantiles=5, nameA=NULL, nameB=NULL, main="", span=0.5, cex=1, cex.lab=1, ylim=c(0,1), xlim=c(min(a), max(a)), lwd=1, cex.axis=1){
	#a=expression ; b=Fop
	tmp.a=na.omit(a)
	tmp.a=na.action(tmp.a)
	
	tmp.b=na.omit(b)
	tmp.b=na.action(tmp.b)

	tmp=union(as.vector(tmp.a), as.vector(tmp.b))
	if(length(tmp)>0){
		a=a[-tmp]
		b=b[-tmp]
	}
	
	L=length(a)
	binsize=round(L*quantiles/100)
	val=L
	while(val%%binsize!=0){
		val=val-1
	}
	bins=seq(from=L-val, to=L, by=binsize)

	cov.ordered=order(a)

	res.fop=res.cov=sd.fop=NULL
	for(i in 1:(length(bins)-1)){
			sub.tmp=cbind(a,b)[cov.ordered,][(bins[i]:(bins[i+1]-1)),]
			res.cov=c(res.cov, mean(as.numeric(sub.tmp[,1]), na.rm=T))
			res.fop=c(res.fop, mean(as.numeric(sub.tmp[,2]), na.rm=T))
			
			sd.fop=c(sd.fop, sd(as.numeric(sub.tmp[,2]), na.rm=T))
			}
#	plot(res.cov, res.fop,xlab="expression", ylab=expression(F["op"]), pch=16, col="black")
#	plot(res.cov, res.fop,xlab=nameA, ylab=nameB, main=main, pch=16, col="black", ylim=c(c(min(res.fop, sd.fop)), max(c(res.fop, sd.fop))))
#	plot(res.cov, res.fop,xlab=nameA, ylab=nameB, cex.lab=cex.lab, main=main, pch=16, col="black", ylim=c(min(res.fop)-max(sd.fop),max(res.fop)+max(sd.fop)), cex=cex)	#THE GOOD
	plot(res.cov, res.fop,xlab=nameA, ylab=nameB, main=main, pch=16, col="black", ylim=ylim, xlim=xlim, cex=cex, cex.lab=cex.lab, cex.axis=cex.axis)

	for(i in 1:(length(bins)-1)){
			segments(res.cov[i], res.fop[i]+sd.fop[i], res.cov[i], res.fop[i]-sd.fop[i])
	}

	mod=loess(res.fop~res.cov, span=span)
	lines(res.cov, fitted(mod), col="red", lwd=lwd)
}


#fonction image avec couleur spéciale pour les cases "nan"
image.nan=function(z,  zlim, col, na.color='gray', outside.below.color='black', outside.above.color='white',...)
{
      zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
      newz.below.outside <- zlim[1] - 2 * zstep # new z for values below zlim
      newz.above.outside <- zlim[2] + zstep # new z for values above zlim
      newz.na <- zlim[2] + 2 * zstep # new z for NA

      z[which(z<zlim[1])] <- newz.below.outside # we affect newz.below.outside
      z[which(z>zlim[2])] <- newz.above.outside # we affect newz.above.outside
      z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na

      zlim[1] <- zlim[1] - 2 * zstep # extend lower limit to include below value
      zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na

      col <- c(outside.below.color, col[1], col, outside.above.color, na.color) #correct by including col[1] at bottom of range

      image(z=z,  zlim=zlim, col=col, ...) # we finally call image(...)
}

pvalue=function(obs=x, ntrl=a){
	#calcule la pvalue de la valeur observée "x" par rapport à la distribution "a"
	if(obs<median(ntrl)){
		return(length(which(ntrl<obs))/length(ntrl))
		}else{
		return(length(which(ntrl>obs))/length(ntrl))
		}
}

newSample=function(data, size, weights){
	#échantillonne dans un vecteur "data" un nombre n=size éléments (avec remise), en tenant compte des poids contenus dans weights
	#exemple: newSample(data=letters[1:10], size=5, weights=1:10)
	tirage=rmultinom(1, size, weights)
	return(rep(data, tirage))
}

fst=function(nDemes, nIndPerDem, extinction, nRecolonizer, migration){
	numQ = 1/(2*nIndPerDem) + extinction/(2 * nRecolonizer) - extinction/(2 * nRecolonizer * 2 * nIndPerDem)
	denomQ = 1- (1-1/(2*nIndPerDem)) * ((1-migration)^2 * (1-extinction) + extinction * (1-1/(2 * nRecolonizer)) * 1/(2*nRecolonizer - 1))
	# phi =  1/(2*nRecolonizer - 1)

	#numQ = (1 - extinction) * (1/(2*nIndPerDem)) + extinction / (2 * nRecolonizer) 
	#denomQ = 1 - (1 - extinction) * ( 1 - migration)^2 * (1 - 1/(2 * nIndPerDem)) - extinction * 1/(2*nRecolonizer - 1) * (1-1/(2 * nRecolonizer)) 
	QR = numQ/denomQ

	return((QR-1/(2*nIndPerDem))*(2*nIndPerDem)/(2*nIndPerDem-1))
}


testFst=function(x, e){
	y = x[which(x$extRate == e), ]
	maxMig = max(y$migRate)
	plot(y$migRate, y$Fwc_ST, main=paste("ext.rate = ", e, sep=""), xlab = "N.m", ylab = expression(F["ST"]), pch=16, cex=1.25, ylim = c(0, 1))
	lines((0:(10*maxMig))/10, fst(y$nDemes[1], y$nInd[1], e, y$recolonization[1], (0:(10*maxMig))/10/maxMig), col="red")
}


twoden = function(vec1,vec2,xbks="auto",ybks="auto",space=1,L=0,xl="",yl="",mn="",lowx="",hix="",lowy="",hiy="",limy=0,hpoint=F,hcol="white",leg=T,return.data=F){

	total=length(vec1);

	if(xbks=="auto"){
		xs=hist(vec1,plot=F)$breaks
		bklongx=space*length(xs)
		xs=hist(vec1,plot=F,breaks=bklongx)$breaks
	}
	else {
		xs=xbks;
	}

	if(ybks=="auto"){
		ys=hist(vec2,plot=F)$breaks
		bklongy=space*length(ys)
		ys=hist(vec2,plot=F,breaks=bklongy)$breaks
	}
	else {
		ys=ybks;
	}

	comb=data.frame(vec1,vec2);

	library(grDevices);

	c=matrix(0,(length(xs)-1),(length(ys)-1))

	for( i in 1:(length(xs)-1)){
		for( j in 1:(length(ys)-1)){

	ay=subset(comb,comb$vec1>=xs[i]);
	bee=subset(ay,ay$vec1<xs[i+1]);
	cee=subset(bee,bee$vec2>=ys[j]);
	d=subset(cee,cee$vec2<ys[j+1])
	c[i,j]=length(d$vec2)/total;

		}
	}

	if(leg==T){
		layout(matrix(c(1,2),2,2,byrow=TRUE), c(3.75,.5), TRUE); 
		par(mar=c(6,5,1.5,1.5)); 
	}

	if(lowy==""){ lowy=ys[1] }
	if(lowx==""){ lowx=xs[1] }
	if(hix==""){ hix=xs[length(xs)] }
	if(hiy==""){ hiy=ys[length(ys)] }

	plot(vec1[1]~vec2[1],col="white",xlim=c(lowx,hix),ylim=c(lowy,hiy),xlab=xl,ylab=yl,main=mn)

	for( i in 1:(length(xs)-1)){
		for( j in 1:(length(ys)-1)){


		den=c[i,j]/max(c);
		rect(xs[i],ys[j],xs[i+1],ys[j+1],border=rgb(red=0,blue=0,green=0,alpha=L),col=(rgb(red=0,blue=0,green=0,alpha=den)))
		if(den==1 && hpoint==T){ points((xs[i+1]+xs[i])/2,(ys[j+1]+ys[j])/2,pch=19,col=hcol); }
		}
	}

	if(leg==T){
		#empty plot
		par(mar=c(8,0,6,4));
		plot(0:10/10,0:10/10,ylim=c(0,1),xlab="",xaxt="n",cex.main=0.5,yaxt="n",ylab="", cex=0,cex.lab=1);  

		#draw gradient
		for(i in 1:99){ rect(0,(i-1)/100,1,(i+2)/100,lwd=0,lty=0, col=rgb(red=0,blue=0,green=0,alpha=i/100)) };
		axis(side=4,at=0:10/10,cex.axis=0.8);  
		text(4,0.5, srt = 90, labels = "relative density", xpd = TRUE) ;
		rect(0,0.99,1,1.05,col=rgb(red=0,blue=0,green=0,alpha=1),lty=0, lwd=0)
	}
	if(return.data==T){ return(xs,ys,c); }
}

plot_quantiSex=function(x, y){
	# x = path to the quantiSex output file
	# y = number of generations to display
	x=read.table(x,h=T)
	layout(matrix(1:2, ncol=2), width=c(2/3, 1/3))
	par(mar=c(5,4,4,0), las=1)	
	if(x$sexSystem[1] == 0){sex="\nherma only\n"}
	if(x$sexSystem[1] == 1){sex="\nandrodioecy\n"}
	if(x$sexSystem[1] == 2){sex="\ngynodioecy\n"}
	
	plot(x$atGeneration, x$meanFemAllocCosexual, xlim=c(0, y), type="l", ylim=c(0,1), lwd=2, xlab="Generation", ylab="frequencies", cex.lab=1.1, main=paste("M=", x$migRate[1], "\tE=", x$extRate[1], "\tk=", x$recolonization[1], sex, "unisex advantage=", x$sexAvantage[1], sep=""))
	lines(x$atGeneration, x$cosexualProportion, col="red", lwd=1.5)
	abline(h=0.5, lty=2)
	plot.new()
	par(mar=c(5,0,4,2))
	legend("left", col=c("black", "red"), c("% Fem. alloc.\nin cosexuals" ,"% of cosexuals"), lty=1, bty="n", lwd=c(2,2), inset=-0.5)
}


image.scale <- function(z, zlim, col = heat.colors(12), breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
	if(!missing(breaks)){
		if(length(breaks) != (length(col)+1)){
			stop("must have one more break than colour")
		}
	}
	if(missing(breaks) & !missing(zlim)){
		breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
	}
	if(missing(breaks) & missing(zlim)){
		zlim <- range(z, na.rm=TRUE)
		zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
		zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
		breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
	}
	poly <- vector(mode="list", length(col))
	for(i in seq(poly)){
		poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
	}
	xaxt <- ifelse(horiz, "s", "n")
	yaxt <- ifelse(horiz, "n", "s")
	if(horiz){
		YLIM<-c(0,1)
		XLIM<-range(breaks)
	}
	if(!horiz){
		YLIM<-range(breaks)
		XLIM<-c(0,1)
	}
	if(missing(xlim)) xlim=XLIM
	if(missing(ylim)) ylim=YLIM
	plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
	for(i in seq(poly)){
		if(horiz){
			polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
		}
		if(!horiz){
			polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
		}
	}
}

plot_time=function(x){
	n = 12 # nDemes to plot
	nGen = 50
	x=read.table(x, h=F)
	toPlot = seq(nrow(x)-nGen, nrow(x), 1)
	par(mfrow=c(3,4), mar=c(4, 3.5, 1, 0.5))
	for(i in 1:n){
#		plot(1:nrow(x), x[,i], type="l", xlab="Time (generation)", ylab="nInd", ylim=c(0, max(x)), lwd=5)
#		polygon(c(0,1:nrow(x),nrow(x)), c(0,x[,i],0), col=rainbow(12)[i], border=NA)
		plot(toPlot, x[toPlot, i], type="l", xlab="", ylab="", ylim=c(0, max(x)), lwd=5, axes=F)
		polygon(c(toPlot[1], toPlot, nrow(x), toPlot[1]), c(x[toPlot[1], i], x[toPlot, i], par("usr")[3], par("usr")[3]), col=rainbow(n)[i], border=NA)
		axis(side=1, at = seq(min(toPlot)-1, max(toPlot)-1, 10), labels = c(seq(min(toPlot)-nrow(x)+nGen, max(toPlot)-nrow(x)+nGen, 10)), cex.axis=1.2)
		mtext("Generations", side = 1, line = 2.2)
		axis(side=2, at = round(seq(0, max(x), length.out=2), 0), cex.axis=1.2)
		mtext("#Individus", side = 2, line = 2.2)
	}
}


compare_matrix=function(x, y, xlab="", ylab="", zlab="", zlim = c(min(c(x,y)), max(c(x,y))), zlimResiduals=c(min(x-y), max(x-y)), cex.lab=1, couleurs=c("green", "white", "red"), watermark=F){
	# plot 4 graphes to show matrices x, y and x-y, as well as the distribution of x-y values
	# x = matrix of expected values
	# y = matrix of observed values
	gradient = colorRampPalette(couleurs)
	
	dev.new(width=8, height=7)
	layout(matrix(c(1,2,3,4,5,6,7,7), byrow=T, nrow=2), width=c(4/5, 1/5, 4/5, 1/5, 4/5, 1/5, 1/2,1/2))

	par(mar=c(4.5, 4, 4, 1), las=1)

	# matrice x
	if(is.null(colnames(x))){
		plot_axes = T
	}else{
		plot_axes = F
	}

	image(x, xlab="", ylab="", col=gradient(100), cex.axis=cex.lab, axes=plot_axes, zlim=zlim)
	mtext(side=3, text="expected", line=0.75, cex=cex.lab)

	if(watermark){watermark()}

	if(is.null(colnames(x))){
		mtext(side=1, text=xlab, line=2.5, cex=cex.lab)
		par(las=3)
		mtext(side=2, text=ylab, line=2.75, cex=cex.lab)
	}else{
		# axe des x
		migRates = rownames(x)
		posX = c((seq(1, length(migRates), 2)), length(migRates))
		axis(1, at=(posX-1)/(length(migRates)-1), labels = migRates[posX])
		mtext(xlab, 1, line=2.5, cex=cex.lab)

		# axe des y
		extRates = colnames(x)
		posY = c((seq(1, length(extRates), 2)), length(extRates))
		axis(2, at=(posY-1)/(length(extRates)-1), labels = extRates[posY])
		par(las=0)
		mtext(ylab, 2, line=2.75, cex=cex.lab)
	}

	par(las=1)

	image.scale(x, horiz=F, col=gradient(100), xlab="", ylab="", cex.lab=cex.lab, cex.axis=cex.lab, zlim=zlim)
	par(las=3)
	mtext(side=2, text=zlab, line=2.5, cex=cex.lab)


	# matrice y
	par(las=1)
	if(is.null(colnames(x))){
		plot_axes = T
	}else{
		plot_axes = F
	}

	image(y, xlab="", ylab="", col=gradient(100), cex.axis=cex.lab, axes=plot_axes, zlim=zlim)
	mtext(side=3, text="simulated", line=0.75, cex=cex.lab)

	if(watermark){watermark()}
	
	if(is.null(colnames(y))){
		mtext(side=1, text=xlab, line=2.5, cex=cex.lab)
		par(las=3)
		mtext(side=2, text=ylab, line=2.75, cex=cex.lab)
	}else{
		# axe des x
		migRates = rownames(y)
		posX = c((seq(1, length(migRates), 2)), length(migRates))
		axis(1, at=(posX-1)/(length(migRates)-1), labels = migRates[posX])
		mtext(xlab, 1, line=2.5, cex=cex.lab)

		# axe des y
		extRates = colnames(y)
		posY = c((seq(1, length(extRates), 2)), length(extRates))
		axis(2, at=(posY-1)/(length(extRates)-1), labels = extRates[posY])
		par(las=0)
		mtext(ylab, 2, line=2.75, cex=cex.lab)
	}

	par(las=1)

	image.scale(y, horiz=F, col=gradient(100), xlab="", ylab="", cex.lab=cex.lab, cex.axis=cex.lab, zlim=zlim)
	par(las=3)
	mtext(side=2, text=zlab, line=2.5, cex=cex.lab)


	# residuals = x - y

	par(las=1)
	if(is.null(colnames(x))){
		plot_axes = T
	}else{
		plot_axes = F
	}

	image(x-y, xlab="", ylab="", col=gradient(100), cex.axis=cex.lab, axes=plot_axes, zlim=zlimResiduals)
	mtext(side=3, text="expected - simulated", line=0.75, cex=cex.lab)

	if(watermark){watermark()}
	
	if(is.null(colnames(y))){
		mtext(side=1, text=xlab, line=2.5, cex=cex.lab)
		par(las=3)
		mtext(side=2, text=ylab, line=2.75, cex=cex.lab)
	}else{
		# axe des x
		migRates = rownames(y)
		posX = c((seq(1, length(migRates), 2)), length(migRates))
		axis(1, at=(posX-1)/(length(migRates)-1), labels = migRates[posX])
		mtext(xlab, 1, line=2.5, cex=cex.lab)

		# axe des y
		extRates = colnames(y)
		posY = c((seq(1, length(extRates), 2)), length(extRates))
		axis(2, at=(posY-1)/(length(extRates)-1), labels = extRates[posY])
		par(las=0)
		mtext(ylab, 2, line=2.75, cex=cex.lab)
	}

	par(las=1)

	image.scale(x-y, horiz=F, col=gradient(100), xlab="", ylab="", cex.lab=cex.lab, cex.axis=cex.lab, zlim=zlimResiduals)
	par(las=3)
	mtext(side=2, text="residuals", line=2.75, cex=cex.lab)

	z=c(x,y)
	par(mar=c(4.5, 4, 4, 3), las=1)
	hist(x-y, xlab="", ylab="", main="", cex.lab=cex.lab, cex.axis=cex.lab, xlim=c(-max(z), max(z)), n=20)
	mtext(side=1, text="residuals", line=2.5, cex=cex.lab)

	if(watermark){watermark()}
}

compare_matrix_fst=function(x, y, xlab="", ylab="", zlab="", cex.lab=1, couleurs=c("green", "white", "red"), watermark=F){
	# plot 4 graphes to show matrices x, y and x-y, as well as the distribution of x-y values
	# x = matrix of expected values
	# y = matrix of observed values
	gradient = colorRampPalette(couleurs)
	
	dev.new(width=8, height=7)
	layout(matrix(c(1,2,3,4,5,6,7,7), byrow=T, nrow=2), width=c(4/5, 1/5, 4/5, 1/5, 4/5, 1/5, 1/2,1/2))

	par(mar=c(4.5, 4, 4, 1), las=1)

	# matrice x
	if(is.null(colnames(x))){
		plot_axes = T
	}else{
		plot_axes = F
	}

	image(x, xlab="", ylab="", col=gradient(100), cex.axis=cex.lab, axes=plot_axes)
	if(watermark==T){
		watermark()
	}
	mtext(side=3, text="expected", line=0.75, cex=cex.lab)

	if(is.null(colnames(x))){
		mtext(side=1, text=xlab, line=2.5, cex=cex.lab)
		par(las=3)
		mtext(side=2, text=ylab, line=2.75, cex=cex.lab)
	}else{
		# axe des x
		migRates = rownames(x)
		posX = c((seq(1, length(migRates), 2)), length(migRates))
		axis(1, at=(posX-1)/(length(migRates)-1), labels = migRates[posX])
		mtext(xlab, 1, line=2.5, cex=cex.lab)

		# axe des y
		extRates = colnames(x)
		posY = c((seq(1, length(extRates), 2)), length(extRates))
		axis(2, at=(posY-1)/(length(extRates)-1), labels = extRates[posY])
		par(las=0)
		mtext(ylab, 2, line=2.75, cex=cex.lab)
	}

	par(las=1)

	image.scale(x, horiz=F, col=gradient(100), xlab="", ylab="", cex.lab=cex.lab, cex.axis=cex.lab)
	par(las=3)
	mtext(side=2, text=zlab, line=2.5, cex=cex.lab)


	# matrice y
	par(las=1)
	if(is.null(colnames(x))){
		plot_axes = T
	}else{
		plot_axes = F
	}

	image(y, xlab="", ylab="", col=gradient(100), cex.axis=cex.lab, axes=plot_axes)
	mtext(side=3, text="simulated", line=0.75, cex=cex.lab)
	if(watermark==T){
		watermark()
	}

	if(is.null(colnames(y))){
		mtext(side=1, text=xlab, line=2.5, cex=cex.lab)
		par(las=3)
		mtext(side=2, text=ylab, line=2.75, cex=cex.lab)
	}else{
		# axe des x
		migRates = rownames(y)
		posX = c((seq(1, length(migRates), 2)), length(migRates))
		axis(1, at=(posX-1)/(length(migRates)-1), labels = migRates[posX])
		mtext(xlab, 1, line=2.5, cex=cex.lab)

		# axe des y
		extRates = colnames(y)
		posY = c((seq(1, length(extRates), 2)), length(extRates))
		axis(2, at=(posY-1)/(length(extRates)-1), labels = extRates[posY])
		par(las=0)
		mtext(ylab, 2, line=2.75, cex=cex.lab)
	}

	par(las=1)

	image.scale(y, horiz=F, col=gradient(100), xlab="", ylab="", cex.lab=cex.lab, cex.axis=cex.lab)
	par(las=3)
	mtext(side=2, text=zlab, line=2.5, cex=cex.lab)


	# residuals = x - y

	par(las=1)
	if(is.null(colnames(x))){
		plot_axes = T
	}else{
		plot_axes = F
	}

	image(x-y, xlab="", ylab="", col=gradient(100), cex.axis=cex.lab, axes=plot_axes, zlim = c(-1, 1))
	mtext(side=3, text="expected - simulated", line=0.75, cex=cex.lab)
	if(watermark==T){
		watermark()
	}

	if(is.null(colnames(y))){
		mtext(side=1, text=xlab, line=2.5, cex=cex.lab)
		par(las=3)
		mtext(side=2, text=ylab, line=2.75, cex=cex.lab)
	}else{
		# axe des x
		migRates = rownames(y)
		posX = c((seq(1, length(migRates), 2)), length(migRates))
		axis(1, at=(posX-1)/(length(migRates)-1), labels = migRates[posX])
		mtext(xlab, 1, line=2.5, cex=cex.lab)

		# axe des y
		extRates = colnames(y)
		posY = c((seq(1, length(extRates), 2)), length(extRates))
		axis(2, at=(posY-1)/(length(extRates)-1), labels = extRates[posY])
		par(las=0)
		mtext(ylab, 2, line=2.75, cex=cex.lab)
	}

	par(las=1)

	image.scale(x-y, horiz=F, col=gradient(100), xlab="", ylab="", cex.lab=cex.lab, cex.axis=cex.lab, zlim = c(-1, 1))
	par(las=3)
	mtext(side=2, text="residuals", line=2.75, cex=cex.lab)

	z=c(x,y)
	par(mar=c(4.5, 4, 4, 3), las=1)
	hist(x-y, xlab="", ylab="", main="", cex.lab=cex.lab, cex.axis=cex.lab, xlim=c(-max(z), max(z)), n=20)
	mtext(side=1, text="residuals", line=2.5, cex=cex.lab)
	if(watermark==T){
		watermark()
	}

}


watermark = function(){
	tag1 = "DRAFT FIGURE"
	tag2 = ""
	#tag2 = "camille.roux.1@unil.ch"
	run.date <- format(Sys.Date(), "%m-%d-%Y")
	text(x = grconvertX(0.5, from = "npc"),  # aligner au centre des X 
	y = grconvertY(0.5, from = "npc"), # aligner au centre des Y 
	labels = tag1, # filigrane central
	cex = 5, font = 2, # en gros et gras
	col = rgb(1, 0, 0, .15), # transparent 
	srt = 45) # angle du texte = 45° 
	texte = paste(tag2, run.date)
	mtext(texte, side = 1, line = -1, adj = 1, col = rgb(1, 0, 0, .15), cex = 1.5)
}

cleanQuantiSexTable = function(x){
	collage = function(a){
		return(paste(a, collapse="_"))
	}
	y = cbind(x$nDemes, x$nIndMaxPerDeme, x$nQuantiLoci, x$selfingRate, x$fecundity, x$migRate, x$extRate, x$recolonization, x$sexSystem, x$sexAvantage)
	z = apply(y, MARGIN=1, FUN="collage")

	res = NULL

	for(i in names(table(z))){
		tmp = which(z==i)
		tmp = x[tmp,]
		res = rbind(res, apply(tmp, MARGIN=2, FUN="mean"))
	}

	return(res)
}

triVariableTable = function(x, y, z){
	# produces a table 'res' with x as rows, y as columns: res[x,y] = z
	# x = migration rates
	# y = extinction rates
	# z = Fst
	xValues = as.numeric(names(table(x)))
	yValues = as.numeric(names(table(y)))
	res = matrix(NA, ncol = length(xValues), nrow = length(yValues))
	colonne = 0
	for(i in xValues){
		ligne = 0
		colonne = colonne + 1
		for(j in yValues){
			ligne = ligne + 1
			tmp = mean(z[which(x == i & y == j)], na.rm = T)
			res[ligne, colonne] = tmp
		}
	}
	colnames(res) = xValues
	rownames(res) = yValues
	return(res)
}


plot3var = function(x, xlab="", ylab="", zlab="", main="", cex.lab=1, couleurs=c("green", "white", "red"), zlim = c(min(x), max(x)), watermark=F, nlevels=10){
	# x = table with values z as a function of 'rows' and 'columns'
	# plot z as a function of 2 variables
	# best continuous gradient: c("#ffffd9", "#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58")
	gradient = colorRampPalette(couleurs)
	
	dev.new(width=8, height=7)
	layout(matrix(c(1,2), byrow=T, ncol=2), width=c(4/5, 1/5))

	par(mar=c(4.5, 4, 4, 1), las=1)

	# matrice x
	if(is.null(colnames(x))){
		plot_axes = T
	}else{
		plot_axes = F
	}

	image(x, xlab="", ylab="", col=gradient(nlevels), cex.axis=cex.lab, axes=plot_axes, zlim=zlim)
	mtext(side=3, text=main, line=0.75, cex=cex.lab)

	if(is.null(colnames(x))){
		mtext(side=1, text=xlab, line=2.5, cex=cex.lab)
		par(las=3)
		mtext(side=2, text=ylab, line=2.75, cex=cex.lab)
	}else{
		# axe des x
		migRates = rownames(x)
		posX = c((seq(1, length(migRates), 2)), length(migRates))
		axis(1, at=(posX-1)/(length(migRates)-1), labels = migRates[posX])
		mtext(xlab, 1, line=2.5, cex=cex.lab)

		# axe des y
		extRates = colnames(x)
		posY = c((seq(1, length(extRates), 2)), length(extRates))
		axis(2, at=(posY-1)/(length(extRates)-1), labels = extRates[posY])
		par(las=0)
		mtext(ylab, 2, line=2.75, cex=cex.lab)
	}
	if(watermark){watermark()}

	par(las=1)

	image.scale(x, horiz=F, col=gradient(nlevels), xlab="", ylab="", cex.lab=cex.lab, cex.axis=cex.lab, zlim=zlim)
	par(las=3)
	mtext(side=2, text=zlab, line=2.5, cex=cex.lab)
}


correlation = function(x, y){
	# function that performs cor.test with the 3 methods
	# and returns the $estimate and the $p.value
	res = NULL
	for(i in c("pearson", "spearman", "kendall")){
		tmp = cor.test(x, y, method=i)
		res = rbind(res, c(tmp$estimate, tmp$p.value))
	}
	rownames(res) = c("pearson", "spearman", "kendall")
	colnames(res) = c("coefficient", "pvalue")
	return(round(res, 5))
}

plotpi = function (pi, pos, par, i, plot_std = T, xlab="", ylab=""){
	maxPi = max(as.numeric(unlist(pi[, grep("avg", colnames(pi))])))
	avg = as.numeric(pi[i, grep("avg", colnames(pi))])
	std = as.numeric(pi[i, grep("std", colnames(pi))])
	pos = as.numeric(pos[i, ])
	if(plot_std==T){
		ylim = c(0, max(avg + std) * 1.05)
	}else{
		ylim = c(0, maxPi*1.01)
	}
	titre1 = paste(colnames(par), collapse = "\t")
	titre2 = paste(par[i, ], collapse = "\t")
	titre = paste(titre1, titre2, sep = "\n")
	plot(pos, avg, xlim=c(0,1), ylim = ylim, col = "white", main = titre, xlab=xlab, ylab=ylab)
	if(plot_std==T){
		polygon(c(pos, rev(pos)), c(avg + std, rev(avg - std)), col = rgb(1, 0, 0, 0.25))
	}
	points(pos, avg, pch = 16, cex = 1.1)
}

plotABCsweep = function(i){
	# i = Id of replicate
	par(mfrow=c(3,2), mar=c(4.5, 3.95, 3, 1.95))
	plot(as.numeric(pos[i,]), as.numeric(tajD[i,]), xlim=c(0,1), xlab="position", ylab="tajD", type="l", cex.lab=1.2, lwd=2); abline(v=par$Sp[i], col="red")
	plot(as.numeric(pos[i,]), as.numeric(achazY[i,]), xlim=c(0,1), xlab="position", ylab="achazY", type="l", cex.lab=1.2, lwd=2); abline(v=par$Sp[i], col="red")
	plot(as.numeric(pos[i,]), as.numeric(pi_avg[i,]), xlim=c(0,1), xlab="position", ylab="pi avg", type="l", cex.lab=1.2, lwd=2); abline(v=par$Sp[i], col="red")
	plot(as.numeric(pos[i,]), as.numeric(pi_std[i,]), xlim=c(0,1), xlab="position", ylab="pi std", type="l", cex.lab=1.2, lwd=2); abline(v=par$Sp[i], col="red")
	plot(as.numeric(pos[i,]), as.numeric(pearsonR[i,]), xlim=c(0,1), xlab="position", ylab="pearsonR", type="l", cex.lab=1.2, lwd=2); abline(v=par$Sp[i], col="red")
	plot(as.numeric(pos[i,]), as.numeric(pearsonP[i,]), xlim=c(0,1), xlab="position", ylab="pearsonP", type="l", cex.lab=1.2, lwd=2); abline(v=par$Sp[i], col="red")
}

plot3var_v2 = function (x, y, z, xlab = "", ylab = "", zlab = "", main = "", cex.lab = 1, couleurs = c("green", "white", "red"), zlim = c(min(z), max(z)), watermark = F, nlevels = 10){
    mat = matrix(NA, length(table(y)), length(table(x)))
    colnames(mat) = names(table(x)) 
    rownames(mat) = names(table(y)) 
    # convert the 3 vectors x, y and z in a matrix mat[y, x] = z 
    ligne = 0
    colonne = 0 
    for( x_i in as.numeric(names(table(x))) ){ # variable 'x' in column
        colonne = colonne + 1
        ligne = 0
        for( y_i in as.numeric(names(table(y))) ){ # variable 'y' in raw
            ligne = ligne + 1
            mat[ligne, colonne] = z[which(x==x_i & y==y_i)]
        }
    }
        
    # create a gradient of colors and a graphic window
    gradient = colorRampPalette(couleurs)
    dev.new(width = 8, height = 7)
    layout(matrix(c(1, 2), byrow = T, ncol = 2), width = c(4/5, 1/5))
    par(mar = c(4.5, 4, 4, 1), las = 1)
    
    # plot
    image(t(mat), xlab = "", ylab = "", col = gradient(nlevels), cex.axis = cex.lab, axes = F, zlim = zlim)
    mtext(side = 3, text = main, line = 0.75, cex = cex.lab)
    
    if (is.null(colnames(mat))) {
        mtext(side = 1, text = xlab, line = 2.5, cex = cex.lab)
        par(las = 3)
        mtext(side = 2, text = ylab, line = 2.75, cex = cex.lab)
    }
    else {
        migRates = rownames(mat)
        posX = c((seq(1, length(migRates), 2)), length(migRates))
        axis(1, at=0:(length(table(x))-1)/(length(table(x))-1), labels=names(table(x)))

        mtext(xlab, 1, line = 2.5, cex = cex.lab)
        extRates = colnames(mat)
        posY = c((seq(1, length(extRates), 2)), length(extRates))
        axis(2, at=0:(length(table(y))-1)/(length(table(y))-1), labels=names(table(y)))
        par(las = 0)
        mtext(ylab, 2, line = 2.75, cex = cex.lab)
    }
    if (watermark) {
        watermark()
    }
    par(las = 1)
    image.scale(mat, horiz = F, col = gradient(nlevels), xlab = "",
        ylab = "", cex.lab = cex.lab, cex.axis = cex.lab, zlim = zlim)
    par(las = 3)
    mtext(side = 2, text = zlab, line = 2.5, cex = cex.lab)
}

