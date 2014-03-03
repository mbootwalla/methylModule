getHypermethylatedOverview <- function(x, disease, index.T, index.N, outdir=".")
{
	if(missing(x)) stop("No data provided \n")

	if(!is(x, "matrix"))
		stop("Expect a valid matrix of beta values as input \n")

	if(missing(index.T) || missing(index.N))
	{
		index <- inferSampleType(colnames(x))
		index.T <- index$Index.Tumor
		index.N <- index$Index.Normal
		index.C <- index$Index.Control
	}

	if(length(index.N) == 0) stop("No normal samples were found, please provide a file containing location of Normal samples\n\n")

	betaT.raw <- x[, index.T, drop=FALSE]
	betaN.raw <- x[, index.N, drop=FALSE]

	beta <- selectFeatures(T=betaT.raw, N=betaN.raw, disease=disease, type="hyper")

	betaT <- beta$T.hyper
	betaN <- beta$N.hyper

	betaT.dichot <- betaT > 0.3
	storage.mode(betaT.dichot) <- "numeric"
	betaN.dichot <- betaN > 0.3
	storage.mode(betaN.dichot) <- "numeric"

	# clustering using ward's method on the jaccard distance
	cluster.tumor <- clusterData(betaT.dichot)
	cluster.normal <- clusterData(betaN.dichot)

	betaT.clustered <- betaT[cluster.tumor$fit.probe$order, cluster.tumor$fit.sample$order, drop=FALSE]
	betaN.clustered <- betaN[cluster.tumor$fit.probe$order, cluster.normal$fit.sample$order, drop=FALSE]

	retval <- SimpleList()
	retval$CLUSTER <- SimpleList("Tumor.Clustered" = betaT.clustered, "Normal.Clustered" = betaN.clustered)
	retval$FIT <- SimpleList("Fit.Tumor.Sample" = cluster.tumor$fit.sample, "Fit.Tumor.Probe" = cluster.tumor$fit.probe,
				       "Fit.Normal.Sample" = cluster.normal$fit.sample, "Fit.Normal.Probe" = cluster.normal$fit.probe)				 
	retval$DIST <- SimpleList("Dist.Tumor.Sample" = cluster.tumor$d.sample, "Dist.Tumor.Probe" = cluster.tumor$d.probe,
					"Dist.Normal.Sample" = cluster.normal$d.sample, "Dist.Normal.Probe" = cluster.normal$d.probe)

	save(retval, file=file.path(outdir, paste(gsub("-","", Sys.Date()),
			   paste(disease, "dichotomized", "hypermethylated", "clustering", "output", "rda", sep="."), sep="_")))

	# cluster assessment using ConsensusClusterPlus
	message("\nPerforming Consensus Clustering\n\n")
	#consensusT <- clusterConsensus(betaT.dichot, view="Hypermethylated")
	title <- file.path(".", paste(gsub("-","", Sys.Date()), "Hypermethylated", "ConsensusClusterResults", sep="_"))
	consensusT <- ConsensusClusterPlus(betaT.dichot, maxK=7, clusterAlg="pam", distance="euclidean",
					   plot="pdf", seed=1234, reps=1000, writeTable=TRUE,
					   title=title)
	consensus <- consensusT[[1]]
	colnames(consensus) <- colnames(betaT.dichot)

	# Write out a table containing the sample cluster number for each k
	consensus.table <- sampleConsensusTable(consensusT)
	write.table(consensus.table,
		    file=file.path(outdir, paste(gsub("-","", Sys.Date()),
		    paste(disease, "dichotomized", "hypermethylated", "sample", "consensus", "table", "txt", sep="."), sep="_")),
		    sep="\t", quote=FALSE, row.names=FALSE)
	
	# generate low res png file
	message("\nGenerating low resolution image of Hypermethylated Overview\n\n")
	png(height=640, width=480, file=file.path(outdir, paste(gsub("-","", Sys.Date()),
	    paste(disease, "dichotomized", "hypermethylated", "overview", "lowres", "png", sep="."), sep="_")), res=72)
	TNplot(retval, consensus=consensus, view="Hypermethylated", disease=disease)
	dev.off()
	
	# generate high res pdf file
	message("\nGenerating high resolution image of Hypermethylated Overview\n\n")
	pdf(file=file.path(outdir, paste(gsub("-","", Sys.Date()),
	    paste(disease, "dichotomized", "hypermethylated", "overview", "highres", "pdf", sep="."), sep="_")))
	TNplot(retval, consensus=consensus, view="Hypermethylated", disease=disease)
	dev.off()
}
