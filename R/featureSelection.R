getGlobalOverview <- function(x, index.T, index.N)
{
	if(missing(x)) stop("No data provided \n")

	if(!is(x, "matrix"))
		stop("Expect a valid matrix of beta values as input \n")

	if(missing(index.T) || missing(index.N))
	{
		warning("Index of tumor or normal samples not provided. Inferring from sample names \n")
		designation <- rep("TUMOR", dim(x)[[2]])
		designation[which(substr(dimnames(x)[[2]], 14, 15) == "11")] <- "NORMAL"
		designation[which(substr(dimnames(x)[[2]], 14, 15) == "20")] <- "CONTROL"
		cat("\n", "Table of number of Tumor, Normal and Control Samples found", "\n\n")
		print(table(designation))
		index.T <- which(designation == "TUMOR")
		index.N <- which(designation == "NORMAL")
	}

	tumorSD <- apply(x[, index.T], 1, sd, na.rm=TRUE)
	clusterProbes <- names(tumorSD)[tumorSD > quantile(tumorSD, .95, na.rm=TRUE)]
	clusterProbes <- sample(clusterProbes, size=as.integer(length(clusterProbes)/5))

	#Clustering
	d.tumor.sample <- dist(t(x[clusterProbes, index.T]), method="euclidean")
	cluster.tumor.sample <- hclust(d.tumor.sample, method="ward")
	d.tumor.probe <- dist(x[clusterProbes, index.T], method="euclidean")
	cluster.tumor.probe <- hclust(d.tumor.probe, method="ward")

	d.normal.sample <- dist(t(x[clusterProbes, index.N]), method="euclidean")
	cluster.normal.sample <- hclust(d.normal.sample, method="ward")
	d.normal.probe <- dist(x[clusterProbes, index.N], method="euclidean")
	cluster.normal.probe <- hclust(d.normal.probe, method="ward")

	datT.clustered <- x[clusterProbes[cluster.tumor.probe$order], index.T[cluster.tumor.sample$order], drop=FALSE]
	datN.clustered <- x[clusterProbes[cluster.normal.probe$order], index.N[cluster.normal.sample$order], drop=FALSE]

	
}

TNplot <- function(x, disease, palette="jet")
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
	fit.sample <- x$FIT$Fit.Sample

	# Restrict total no. of data points to 1 million for plotting purposes
	#n <- ceiling(1000000 / ncol(datT))
	#if(nrow(datT) > n){
#		datT <- datT[sample(rownames(datT), n), , drop=FALSE]
#		datN <- datN[rownames(datT), , drop=FALSE]
#	}

	if(match.arg(palette) == "jet"){
		colors <- jet.colors(100)
	} else {
		colors <- palette
	}

	#par(mar = c(5, 2, 0, 3))
	#layout(matrix(c(4,3,2,1), nrow=2, ncol=2, byrow=T), widths=c(1.5,4), heights=c(1.5,4))
	par(mar = c(5, 0.5, 2, 3))
	layout(matrix(c(5,0,4,3,2,1), nrow=2, ncol=3, byrow=T), widths=c(1.5,0.25,3.75), heights=c(1.5,4))
	image(1:ncol(datT), 1:nrow(datT), t(datT), axes=F, col=colors, xlab="", ylab="")
	box()
	axis(side=1, at=1:ncol(datT), labels=substr(colnames(datT), 6, 12), tick=FALSE, las=3, cex.axis=0.75)
	#title(paste("Global overview of Methylation of", disease, "data", sep=" "))
	mtext(paste(ncol(datT), disease, "Tumor Samples", sep=" "), side=3, line=1, cex=1.0)
	par(mar = c(5,1.5,2,0))
	image(t(cgi2), col=c("white", "black"), axes=F)
	#mtext(paste(ncol(datT), disease, "Tumor Samples", sep=" "), side=1, line=1)
	mtext(paste(nrow(datT), "Probes", sep=" "), side=2, line=1)
	par(mar = c(5,2,2,5))
	image(1:ncol(datN), 1:nrow(datN), t(datN), axes=F, col=colors, xlab="", ylab="")
	box()
	mtext(paste(ncol(datN), "Normals", sep=" "), side=3, line=1)
	par(mar = c(2, 0.5, 0, 3))
	plot(as.dendrogram(fit.sample), axes=F, xaxs="i", leaflab="none")
	par(mar=c(7,2,5,5), cex=0.75)
	colorstrip(colors, description=expression(paste(beta, " value 0 -> 1", sep="")))
	par(def.par)
}


colorstrip <- function(colors, description)
{
	count <- length(colors)
	m <- matrix(1:count, count, 1)
	image(m, col=colors, ylab="", axes=FALSE)
	box()
	mtext(description, 1, adj=0.5, line=0.5)
}
