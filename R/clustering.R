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

clusterConsensus <- function(x, alg="pam", dist.method="euclidean", type="png", outdir="ConsensusClusterResults")
{
	require('ConsensusClusterPlus')
	require('cluster')
	#alg <- match.arg(alg)
	#dist.method <- match.arg(dist.method)
	title <- file.path(getwd(), outdir)
	return(ConsensusClusterPlus(x, maxk=7, clusterAlg=alg, distance=dist.method,
				    plot=type, seed=1234, reps=1000, writeTable=TRUE,
				    title=title))
}
