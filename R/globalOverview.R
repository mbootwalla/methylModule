getGlobalOverview <- function(x, disease, index.T, index.N, outdir=".")
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

	betaT <- x[clusterProbes, index.T, drop=FALSE]
	betaN <- x[clusterProbes, index.N, drop=FALSE]
	
	#Clustering
	cluster.tumor <- clusterData(betaT, dist.method="euclidean")
	cluster.normal <- clusterData(betaN, dist.method="euclidean")
	
	betaT.clustered <- betaT[cluster.tumor$fit.probe$order, cluster.tumor$fit.sample$order, drop=FALSE]
	betaN.clustered <- betaN[cluster.tumor$fit.probe$order, cluster.normal$fit.sample$order, drop=FALSE]

	retval <- SimpleList()
	retval$CLUSTER <- SimpleList("Tumor.Clustered" = betaT.clustered, "Normal.Clustered" = betaN.clustered)
	retval$FIT <- SimpleList("Fit.Tumor.Sample" = cluster.tumor$fit.sample, "Fit.Tumor.Probe" = cluster.tumor$fit.probe,
				 "Fit.Normal.Sample" = cluster.normal$fit.sample, "Fit.Normal.Probe" = cluster.normal$fit.probe)
	retval$DIST <- SimpleList("Dist.Tumor.Sample" = cluster.tumor$d.sample, "Dist.Tumor.Probe" = cluster.tumor$d.probe,
				 "Dist.Normal.Sample" = cluster.normal$d.sample, "Dist.Normal.Probe" = cluster.normal$d.probe)

	save(retval, file=file.path(outdir, paste(gsub("-","", Sys.Date()), paste(disease, "clustering", "output", "rda", sep="."), sep="_")))

	#cluster assessment using ConsensusClusterPlus
	message("\nPerforming Consensus Clustering\n\n")
	#consensusT <- clusterConsensus(betaT, view="Global")
	title <- file.path(".", paste(gsub("-","", Sys.Date()), "Global", "ConsensusClusterResults", sep="_"))
	consensusT <- ConsensusClusterPlus(betaT, maxK=7, clusterAlg="pam", distance="euclidean",
					   plot="pdf", seed=1234, reps=1000, writeTable=TRUE,
					   title=title)
	consensus <- consensusT[[1]]
	colnames(consensus) <- colnames(betaT)

	# Write out a table containing the sample cluster number for each k
	consensus.table <- sampleConsensusTable(consensusT)
	write.table(consensus.table,
		    file=file.path(outdir, paste(gsub("-","", Sys.Date()),
		    paste(disease, "global", "sample", "consensus", "table", "txt", sep="."), sep="_")),
		    sep="\t", quote=FALSE)
	
	# generate low res png file
	message("\nGenerating low resolution image of Global Overview\n\n")
	png(height=640, width=480, file=file.path(outdir, paste(gsub("-","", Sys.Date()),
	    paste(disease, "global", "overview", "lowres", "png", sep="."), sep="_")), res=72)
	TNplot(retval, consensus=consensus, view="Global", disease=disease)
	dev.off()
	
	# generate high res pdf file
	message("\nGenerating high resolution image of Global Overview\n\n")
	pdf(file=file.path(outdir, paste(gsub("-","", Sys.Date()),
	    paste(disease, "global", "overview", "highres", "pdf", sep="."), sep="_")))
	TNplot(retval, consensus=consensus, view="Global", disease=disease)
	dev.off()
}

