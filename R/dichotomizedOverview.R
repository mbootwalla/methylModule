getDichotomizedOverview <- function(x, disease, index.T, index.N, outdir=".", type="pdf")
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
	beta.hypo <- selectFeatures(T=betaT.raw, N=betaN.raw, disease=disease, type="hypo")

	betaT.hyper <- beta.hyper$T.hyper
	betaN.hyper <- beta.hyper$N.hyper
	betaT.hypo <- beta.hypo$T.hypo
	betaN.hypo <- beta.hypo$N.hypo
	
	betaT.dichot.hyper <- betaT.hyper > 0.3
	storage.mode(betaT.dichot.hyper) <- "numeric"
	betaN.dichot.hyper <- betaN.hyper > 0.3
	storage.mode(betaN.dichot.hyper) <- "numeric"

	betaT.dichot.hypo <- betaT.hypo > 0.3
	storage.mode(betaT.dichot.hypo) <- "numeric"
	betaN.dichot.hypo <- betaN.hypo > 0.3
	storage.mode(betaN.dichot.hypo) <- "numeric"

	# clustering using ward's method on the jaccard distance
	cluster.tumor.hyper <- clusterData(betaT.dichot.hyper)
	cluster.normal.hyper <- clusterData(betaN.dichot.hyper)
	cluster.tumor.hypo <- clusterData(betaT.dichot.hypo)
	cluster.normal.hypo <- clusterData(betaN.dichot.hypo)

	betaT.clustered.hyper <- betaT.hyper[cluster.tumor.hyper$fit.probe$order, cluster.tumor.hyper$fit.sample$order, drop=FALSE]
	betaN.clustered.hyper <- betaN.hyper[cluster.tumor.hyper$fit.probe$order, cluster.normal.hyper$fit.sample$order, drop=FALSE]
	betaT.clustered.hypo <- betaT.hypo[cluster.tumor.hypo$fit.probe$order, cluster.tumor.hypo$fit.sample$order, drop=FALSE]
	betaN.clustered.hypo <- betaN.hypo[cluster.tumor.hypo$fit.probe$order, cluster.normal.hypo$fit.sample$order, drop=FALSE]

	retval.hyper <- SimpleList()
	retval.hyper$CLUSTER <- SimpleList("Tumor.Clustered" = betaT.clustered.hyper, "Normal.Clustered" = betaN.clustered.hyper)
	retval.hyper$FIT <- SimpleList("Fit.Tumor.Sample" = cluster.tumor.hyper$fit.sample, "Fit.Tumor.Probe" = cluster.tumor.hyper$fit.probe,
				       "Fit.Normal.Sample" = cluster.normal.hyper$fit.sample, "Fit.Normal.Probe" = cluster.normal.hyper$fit.probe)				 
	retval.hyper$DIST <- SimpleList("Dist.Tumor.Sample" = cluster.tumor.hyper$d.sample, "Dist.Tumor.Probe" = cluster.tumor.hyper$d.probe,
					"Dist.Normal.Sample" = cluster.normal.hyper$d.sample, "Dist.Normal.Probe" = cluster.normal.hyper$d.probe)


	retval.hypo <- SimpleList()
	retval.hypo$CLUSTER <- SimpleList("Tumor.Clustered" = betaT.clustered.hypo, "Normal.Clustered" = betaN.clustered.hypo)
	retval.hypo$FIT <- SimpleList("Fit.Tumor.Sample" = cluster.tumor.hypo$fit.sample, "Fit.Tumor.Probe" = cluster.tumor.hypo$fit.probe,
				       "Fit.Normal.Sample" = cluster.normal.hypo$fit.sample, "Fit.Normal.Probe" = cluster.normal.hypo$fit.probe)				 
	retval.hypo$DIST <- SimpleList("Dist.Tumor.Sample" = cluster.tumor.hypo$d.sample, "Dist.Tumor.Probe" = cluster.tumor.hypo$d.probe,
					"Dist.Normal.Sample" = cluster.normal.hypo$d.sample, "Dist.Normal.Probe" = cluster.normal.hypo$d.probe)

	save(retval.hyper, file=file.path(outdir, paste(gsub("-","", Sys.Date()),
			   paste(disease, "dichotomized", "hypermethylated", "clustering", "output", "rda", sep="."), sep="_")))
	save(retval.hypo, file=file.path(outdir, paste(gsub("-","", Sys.Date()),
			   paste(disease, "dichotomized", "hypomethylated", "clustering", "output", "rda", sep="."), sep="_")))

	if(type=="pdf")
	{
		pdf(file=file.path(outdir, paste(disease, "dichotomized", "hypermethylated", "overview", "pdf", sep=".")))
		TNplot(retval.hyper, type="Hypermethylated", disease=disease)
		dev.off()

		pdf(file=file.path(outdir, paste(disease, "dichotomized", "hypomethylated", "overview", "pdf", sep=".")))
		TNplot(retval.hypo, type="Hypermethylated", disease=disease)
		dev.off()
	} else
	{
		TNplot(retval.hyper, type="Hypomethylated", disease=disease)
		dev.new()
		TNplot(retval.hypo, type="Hypomethylated", disease=disease)
	}
}
