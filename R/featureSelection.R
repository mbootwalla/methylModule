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
