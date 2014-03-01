TNplot <- function(x, disease, consensus=NULL, view="Global", palette="jet")
{
	require("matlab")
	def.par <- par(no.readonly=TRUE)
	
	if(!is(x, "SimpleList")){
		stop("Data should be a SimpleList object \n")
	}
	if(missing(disease)){
		stop("Please provide the name of the Tumor type being visualized \n")
	}

	datT <- x$CLUSTER$Tumor.Clustered
	datN <- x$CLUSTER$Normal.Clustered
	fit.tumor.sample <- x$FIT$Fit.Tumor.Sample
	probeOrder <- x$FIT$Fit.Tumor.Probe$order
	sampleOrder <- fit.tumor.sample$order

	if(match.arg(palette) == "jet"){
		colors <- jet.colors(100)
	} else {
		colors <- palette
	}

	# Generate CpG Island annotation
	data(probes.cgi)
	cgi <- rep(0L, nrow(datT))
	cgi[which(rownames(datT) %in% probes.cgi)] <- 1L
	cgi <- cgi[probeOrder]

	#par(mar = c(5, 2, 0, 3))
	# Set the plot layout
	par(oma = c(0, 0, 4, 0))
	par(mar = c(5, 0.5, 0.25, 0))
	#layout(matrix(c(5,0,4,0,3,2,1,6), nrow=2, ncol=4, byrow=T), widths=c(1.25,0.25,3.5,1.0), heights=c(1.5,4))
	lmat <- rbind(c(5, NA, 4, NA), c(3, 2, 1, 6))
	lhei <- c(1, 4)
	lwid <- c(1.25, 0.25, 3.5, 1)
	if(!is.null(consensus)){
		consensus <- consensus[, sampleOrder]
		clustColors <- apply(consensus, 1, function(x){levels(as.factor(x))})
		for(i in 1:nrow(consensus)){
			lmat <- rbind(lmat[1, ], c(NA, NA, max(lmat, na.rm=T)+1, NA), lmat[-1,])
			lhei <- c(lhei[1], 0.1, lhei[-1])
		}
	}
	#layout(matrix(c(5,0,4,0,0,0,7,0,3,2,1,6), nrow=3, ncol=4, byrow=T),
	#       widths=c(1.25,0.25,3.5,1.0), heights=c(1.4,0.1,4))
	lmat[is.na(lmat)] <- 0
	on.exit(par(def.par))
	layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
	image(1:ncol(datT), 1:nrow(datT), t(datT), axes=F, col=colors, xlab="", ylab="")
	#box()
	#axis(side=1, at=1:ncol(datT), labels=substr(colnames(datT), 6, 12), tick=FALSE, las=3, cex.axis=0.75)
	mtext(paste(ncol(datT), "Tumors", sep=" "), side=1, line=1, cex=0.75)
	par(mar = c(5, 1.5, 0.25, 0))
	image(t(cgi), col=c("gray", "black"), axes=F)
	mtext("CpG Island", side=1, line=0.5, cex=0.5, las=2)
	mtext(paste(nrow(datT), "Probes", sep=" "), side=2, line=1, cex=0.75)
	par(mar = c(5, 2, 0.25, 5))
	image(1:ncol(datN), 1:nrow(datN), t(datN), axes=F, col=colors, xlab="", ylab="")
	#box()
	mtext(paste(ncol(datN), "Normals", sep=" "), side=1, line=1, cex=0.75)
	par(mar = c(0, 0.5, 0, 0))
	plot(as.dendrogram(fit.tumor.sample), axes=F, xaxs="i", leaflab="none")
	par(mar=c(5, 2, 0, 0))
	colorstrip(colors, description=expression(paste(beta, " value 0 -> 1", sep="")), cex=0.5)
	par(mar=c(5, 0.5, 0.25, 0))
	plot.new()	
	legend("topleft", legend=c("Island", "Non-Island"), fill=c("black", "gray"), title="UCSC CpG Island")
	if(!is.null(consensus)){
		for(i in nrow(consensus):1){
			par(mar = c(0.1, 0.5, 0, 0))
			image(as.matrix(as.numeric(as.factor(consensus[i,]))), col=clustColors[[i]], axes=F)
			mtext(paste("k", i+1, sep="="), side=2, line=0.5, cex=0.5, las=2)
		}
	}
	title(paste(view, "overview of", disease, "Methylation data :", ncol(datT), "Tumors", ncol(datN), "Normals", sep=" "), outer=TRUE)
	#par(def.par)
}


colorstrip <- function(colors, description, ...)
{
	count <- length(colors)
	m <- matrix(1:count, count, 1)
	image(m, col=colors, ylab="", axes=FALSE)
	#box()
	mtext(description, 1, adj=0.5, line=0.5, ...)
}
