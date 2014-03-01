clusterData <- function(x, dist.method="binary", hclust.method="ward")
{
	#dist.method <- match.arg(dist.method)
	#hclust.method <- match.arg(hclust.method)
	d.sample <- dist(t(x), method=dist.method)
	d.probe <- dist(x, method=dist.method)
	fit.sample <- hclust(d.sample, method=hclust.method)
	fit.probe <- hclust(d.probe, method=hclust.method)
	retval <- list("d.sample"=d.sample, "d.probe"=d.probe,
		       "fit.sample"=fit.sample, "fit.probe"=fit.probe)
	return(retval)
}

clusterConsensus <- function(x, alg="pam", dist.method="euclidean", type="pdf", view="Global", outdir="ConsensusClusterResults")
{
	suppressPackageStartupMessages(require('ConsensusClusterPlus'))
	suppressPackageStartupMessages(require('cluster'))
	#alg <- match.arg(alg)
	#dist.method <- match.arg(dist.method)
	title <- paste(gsub("-","", Sys.Date()), view, outdir, sep="_")
	result <- ConsensusClusterPlus(x, maxK=7, clusterAlg=alg, distance=dist.method,
				    plot=type, seed=1234, reps=1000, writeTable=TRUE,
				    title=title)
	return(result)
}
