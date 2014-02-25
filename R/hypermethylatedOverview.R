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

	beta.hyper <- selectFeatures(T=betaT.raw, N=betaN.raw, disease=disease, type="hyper")

	betaT.hyper <- beta.hyper$T.hyper
	betaN.hyper <- beta.hyper$N.hyper

	betaT.dichot.hyper <- betaT.hyper > 0.3
	storage.mode(betaT.dichot.hyper) <- "numeric"
	betaN.dichot.hyper <- betaN.hyper > 0.3
	storage.mode(betaN.dichot.hyper) <- "numeric"

	# clustering using ward's method on the jaccard distance
	cluster.tumor.hyper <- clusterData(betaT.dichot.hyper)
	cluster.normal.hyper <- clusterData(betaN.dichot.hyper)

	betaT.clustered.hyper <- betaT.hyper[cluster.tumor.hyper$fit.probe$order, cluster.tumor.hyper$fit.sample$order, drop=FALSE]
	betaN.clustered.hyper <- betaN.hyper[cluster.tumor.hyper$fit.probe$order, cluster.normal.hyper$fit.sample$order, drop=FALSE]

	retval.hyper <- SimpleList()
	retval.hyper$CLUSTER <- SimpleList("Tumor.Clustered" = betaT.clustered.hyper, "Normal.Clustered" = betaN.clustered.hyper)
	retval.hyper$FIT <- SimpleList("Fit.Tumor.Sample" = cluster.tumor.hyper$fit.sample, "Fit.Tumor.Probe" = cluster.tumor.hyper$fit.probe,
				       "Fit.Normal.Sample" = cluster.normal.hyper$fit.sample, "Fit.Normal.Probe" = cluster.normal.hyper$fit.probe)				 
	retval.hyper$DIST <- SimpleList("Dist.Tumor.Sample" = cluster.tumor.hyper$d.sample, "Dist.Tumor.Probe" = cluster.tumor.hyper$d.probe,
					"Dist.Normal.Sample" = cluster.normal.hyper$d.sample, "Dist.Normal.Probe" = cluster.normal.hyper$d.probe)

	save(retval.hyper, file=file.path(outdir, paste(gsub("-","", Sys.Date()),
			   paste(disease, "dichotomized", "hypermethylated", "clustering", "output", "rda", sep="."), sep="_")))

	# generate low res png file
	png(file=file.path(outdir, paste(disease, "dichotomized", "hypermethylated", "overview", "lowres", "png", sep=".")), res=72)
	TNplot(retval.hypo, type="Hypermethylated", disease=disease)
	dev.off()
	# generate high res pdf file
	pdf(file=file.path(outdir, paste(disease, "dichotomized", "hypermethylated", "overview", "highres", "pdf", sep=".")))
	TNplot(retval.hyper, type="Hypermethylated", disease=disease)
	dev.off()
}
