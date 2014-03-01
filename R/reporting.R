sampleConsensusTable <- function(consensus)
{
	if(!is(consensus, "list"))
		stop("Input should be a list that is output by the ConsensusClusterPlus function\n\n")

	retval <- do.call(rbind, lapply(consensus[-1],
						function(x){
							with(x, consensusClass)}))
	rownames(retval) <- paste("k", 1:nrow(retval), sep="=")
	return(retval)
}

generateReport <- function(disease, outdir=".")
{
	require("Nozzle.R1")

	report <- newReport("Clustering of Methylation Data")
	report <- setCollectionDate(report, format(Sys.Date(), "%b %d %Y"))
	report <- addToIntroduction(report,
				    newParagraph("This pipeline performs clustering of Level 3 methylation data obtained from the TCGA Project. The following steps make up the pipeline:"),
				    newList(newParagraph("Feature Selection"),
					    newParagraph("Data Dichotomization"),
					    newParagraph("Cluster Assessment using Consensus Clustering"),
					    newParagraph("Clustering the data using an agglomerative hierarchical clustering algorithm and visualizing it as a Heatmap")),
				    newParagraph("The analysis pipeline provides three overviews of the methylation data:"),
				    newList(isNumbered=TRUE,
					    newParagraph("Global Overview: This overview provides a top-level view of the differences in the methylation profiles between the Normal and Tumor Samples"),
					    newParagraph("Hypermethylated Overview: This overview zooms in on cancer-specific increase in the methylation levels of the Tumor samples when compared to the Normal Samples"),
					    newParagraph("Hypomethylated Overview: This overview zooms in on cancer-specific decrease in the methylation levels of the Tumor samples when compared to the Normal Samples")))

	report <- addToSummary(report,
			       newParagraph("We picked 5,000 features using the feature selection criteria outlined in the ", asStrong("Methods"), " section below. Both the samples and features were clustered using an unsupervised clustering algorithm. The clustered samples was visualized using an annotated heatmap. Cluster assessment was done using consensus clustering using the partition around medoids (PAM) algorithm for increasing number of clusters from ", asEmph("k=2"), " to ", asEmph("k=7"), ". Various metrics associated with consensus clustering are provided to assist the user in making a decision regarding the number of stable cluster subtypes associated with ", disease))

	writeReport(report, filename=file.path(outdir, "nozzle"))
}
					    
