getHypomethylatedOverview <- function(x, disease, index.T, index.N, outdir=".")
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
	
	beta.hypo <- selectFeatures(T=betaT.raw, N=betaN.raw, disease=disease, type="hypo")

	betaT.hypo <- beta.hypo$T.hypo
	betaN.hypo <- beta.hypo$N.hypo

	betaT.dichot.hypo <- betaT.hypo > 0.3
	storage.mode(betaT.dichot.hypo) <- "numeric"
	betaN.dichot.hypo <- betaN.hypo > 0.3
	storage.mode(betaN.dichot.hypo) <- "numeric"

	# clustering using ward's method on the jaccard distance
	cluster.tumor.hypo <- clusterData(betaT.dichot.hypo)
	cluster.normal.hypo <- clusterData(betaN.dichot.hypo)

	betaT.clustered.hypo <- betaT.hypo[cluster.tumor.hypo$fit.probe$order, cluster.tumor.hypo$fit.sample$order, drop=FALSE]
	betaN.clustered.hypo <- betaN.hypo[cluster.tumor.hypo$fit.probe$order, cluster.normal.hypo$fit.sample$order, drop=FALSE]

	retval.hypo <- SimpleList()
	retval.hypo$CLUSTER <- SimpleList("Tumor.Clustered" = betaT.clustered.hypo, "Normal.Clustered" = betaN.clustered.hypo)
	retval.hypo$FIT <- SimpleList("Fit.Tumor.Sample" = cluster.tumor.hypo$fit.sample, "Fit.Tumor.Probe" = cluster.tumor.hypo$fit.probe,
				       "Fit.Normal.Sample" = cluster.normal.hypo$fit.sample, "Fit.Normal.Probe" = cluster.normal.hypo$fit.probe)				 
	retval.hypo$DIST <- SimpleList("Dist.Tumor.Sample" = cluster.tumor.hypo$d.sample, "Dist.Tumor.Probe" = cluster.tumor.hypo$d.probe,
					"Dist.Normal.Sample" = cluster.normal.hypo$d.sample, "Dist.Normal.Probe" = cluster.normal.hypo$d.probe)

	save(retval.hypo, file=file.path(outdir, paste(gsub("-","", Sys.Date()),
			   paste(disease, "dichotomized", "hypomethylated", "clustering", "output", "rda", sep="."), sep="_")))

	# generate low res png file
	png(file=file.path(outdir, paste(disease, "dichotomized", "hypomethylated", "overview", "lowres", "png", sep=".")), res=72)
	TNplot(retval.hypo, type="Hypomethylated", disease=disease)
	dev.off()
	# generate high res pdf file
	pdf(file=file.path(outdir, paste(disease, "dichotomized", "hypomethylated", "overview", "highres", "pdf", sep=".")))
	TNplot(retval.hyper, type="Hypomethylated", disease=disease)
	dev.off()
}
