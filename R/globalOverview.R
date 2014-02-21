getGlobalOverview <- function(x, disease, index.T, index.N, outdir=".", type="pdf")
{
	if(missing(x)) stop("No data provided \n")

	if(!is(x, "matrix"))
		stop("Expect a valid matrix of beta values as input \n")

	if(missing(index.T) || missing(index.N))
	{
		message("Index of tumor or normal samples not provided. Inferring from sample names \n")
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
	set.seed(12345)
	clusterProbes <- sample(clusterProbes, size=5000)

	#Clustering
	cluster.tumor <- clusterData(x[clusterProbes, index.T, drop=FALSE], dist.method="euclidean")
	cluster.normal <- clusterData(x[clusterProbes, index.N, drop=FALSE], dist.method="euclidean")
	
	betaT.clustered <- x[clusterProbes[cluster.tumor$fit.probe$order], index.T[cluster.tumor$fit.sample$order], drop=FALSE]
	betaN.clustered <- x[clusterProbes[cluster.tumor$fit.probe$order], index.N[cluster.normal$fit.sample$order], drop=FALSE]

	retval <- SimpleList()
	retval$CLUSTER <- SimpleList("Tumor.Clustered" = betaT.clustered, "Normal.Clustered" = betaN.clustered)
	retval$FIT <- SimpleList("Fit.Tumor.Sample" = cluster.tumor$fit.sample, "Fit.Tumor.Probe" = cluster.tumor$fit.probe,
				 "Fit.Normal.Sample" = cluster.normal$fit.sample, "Fit.Normal.Probe" = cluster.normal$fit.probe)
	retval$DIST <- SimpleList("Dist.Tumor.Sample" = cluster.tumor$d.sample, "Dist.Tumor.Probe" = cluster.tumor$d.probe,
				 "Dist.Normal.Sample" = cluster.normal$d.sample, "Dist.Normal.Probe" = cluster.normal$d.probe)

	save(retval, file=file.path(outdir, paste(gsub("-","", Sys.Date()), paste(disease, "clustering", "output", "rda", sep="."), sep="_")))

	# Open a device if type is specified
	if(type=="pdf")
	{
		pdf(file=file.path(outdir, paste(disease, "global", "overview", "pdf", sep=".")))
		TNplot(retval, type="Global", disease=disease)
		dev.off()
	} else
	{
		TNplot(retval, disease=disease)
	}
}

