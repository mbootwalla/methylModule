sampleConsensusTable <- function(consensus)
{
	if(!is(consensus, "list"))
		stop("Input should be a list that is output by the ConsensusClusterPlus function\n\n")

	df <- do.call(rbind, lapply(consensus[-1],
					function(x){
						with(x, consensusClass)}))
	rownames(df) <- paste("k", 1:nrow(df), sep="=")
	retval <- data.frame("SampleName"=colnames(df), stringsAsFactors=FALSE)
	retval <- cbind(retval, t(df))
	rownames(retval) <- NULL
	return(retval)
}

generateReport <- function(disease, outdir=".")
{
	require("Nozzle.R1")

	oldwd <- getwd()
	setwd(outdir)
	figures.lowres <- list.files(pattern="png")
	figures.highres <- list.files(pattern="pdf")
	tables <- list.files(pattern="table.txt")

	# Setup the Report
	report <- newReport("Clustering of Methylation Data")
	report <- setCollectionDate(report, format(Sys.Date(), "%d %b %Y"))

	# Add the Introduction
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

	# Add the Summary Section
	report <- addToSummary(report,
			       newParagraph("We picked 5,000 features using the feature selection criteria outlined in the ", asStrong("Methods"), " section below. Both the samples and features wre clustered using an unsupervised clustering algorithm. The clustered samples were visualized using an annotated heatmap. Cluster assessment was done using consensus clustering using the partition around medoids (PAM) algorithm for increasing number of clusters from ", asEmph("k=2"), " to ", asEmph("k=7"), ". Various metrics associated with consensus clustering are provided to assist the user in making a decision regarding the number of stable cluster subtypes associated with ", disease))

	# Setup and add the Methods
	method1 <- addTo(newSubSection("Feature Selection"),
			 newParagraph("The following features are removed from the data for all three overviews:"),
			 newList(newParagraph("Features that interrogate non-CpG sites"),
				 newParagraph("Features that interrogate sites on the sex chromosomes"),
				 newParagraph("Any features that are masked (represented as NA) in the Level 3 data")),
			 newParagraph("Additional feature selection is done for the hypo and hypermethylated overviews. In order to zoom in on cancer specific increase/decrease in methylation in the Tumor samples we select the corresponding features that are unmethylated (beta value < 0.2) / highly methylated (beta value > 0.7) in Normal samples."),
			 newParagraph("Next, we select features that exhibit the most variation across the Tumor samples. This step is performed for all three overviews. For the Global overview we select the top 5% most variable features whereas for the other two overviews we select the top 20% most variable features. The different proportions of features selected reflect the different number of features remaining after the additional selection step performed above for the hypo/hypermethylated overviews."),
			 newParagraph("Finally we randomly sample 5,000 features from the above cohorts and use these for the next steps in the pipeline. We use a ", asEmph("seed"), " of ", asEmph("1234"), " for the random sampling procedure in order to make it reproducible."),
			 newParagraph("Note: Features in the above context refer to the probes on the Illumina Infinium HumanMethylation450k array."))

	method2 <- addTo(newSubSection("Data Dichotomization"),
			 newParagraph("Data dichotomization greatly ameliorates the effect of tumor sample purity on the clustering and also removes a great portion of residual batch/platofrm effects that are mostly reflected in small variations near the two ends. Dichotomization results in the continuous data being converted to a binary form with ", asEmph("1 = methylated"), " and ", asEmph("0 = unmethylated"), "."),
			 newParagraph("Data dichotomization is done at a beta value of 0.3. Features having a beta value > than 0.3 are designated as methylated and those with a beta value <= 0.3 are designated as unmethylated. Data dichotomization is performed only for the hypo and hypermethylated overviews."))

	method3 <- addTo(newSubSection("Cluster Assessment: Consensus Clustering"),
			 newParagraph("Consensus clustering is a method that provides quantitative evidence for determining the number and membership of possible clusters within a dataset. It involves subsampling from a set of items and determines clusterings of specified counts (k). Then, pairwise consensus values, the proportion that two items occupied the same cluster out of the number of times they ouccrred in the same subsample, are calculated and stored in a consensus matrix. The consensus matrix is summarized in several graphical displays that enable a user to decide upon a reasonable cluster number and membership. This approach can help the user in identifying novel disease subtypes that are represented as distinct clusters within the data."),
			 newParagraph("Consensus clustering is done using the ", asEmph("R/Bioconductor"), " package ", asEmph("ConsensusClusterPlus"), ". Consensus clustering is performed using the Partitioning Around Medoids (PAM) algorithm to cluster the samples for increasing values of ", asEmph("k"), " from ", asEmph("k=2"), " to ", asEmph("k=7"), ". Consensus clustering is performed on the dichotomized data for the hypo and hypermethylated overviews while it is done on the continuous data for the global overview. We use the ", asEmph("ConsensusClusterPlus"), " function with 1000 repetitions and default settings for all other parameters. A ", asEmph("seed"), " of ", asEmph("1234"), " is used for the purpose of reproducibility. All the output from consensus cluster is made available below."))

	method4 <- addTo(newSubSection("Clustering and Data Visualization"),
			 newParagraph("Samples are clustered using an unsupervised agglomerative hierarchical clustering algorithm such as Ward's Minimum Variance method. Ward's method aims at finding compact, spherical clusters. Clustering using Ward's method is performed on the Jaccard distance, a distance measure that best suits binary data i.e. the hypo and hypermethylated overviews. For the Global overview clustering is done using Ward's method on the Euclidean distance. Hierarchical clustering is performed using the R function ", asEmph("hclust"), ". For the purpose of visualization the features are clustered as well using the same parameters as the samples."),
			 newParagraph("The continuous data is visualized as a heatmap. The beta values (0 to 1) are represented on a color gradient from blue to red with blue indicating low to no methylation and red indicating high methylation. The color scale is displayed on the figure. Both the Normal samples and the Tumor samples are visualized. The same set of features selected for the Tumor samples are used to visualize the Normal samples. The features on both heatmaps i.e. Normal and Tumor have the same ordering which is determined based on the clustering of the features in the Tumor samples. We also include a dendrogram representing the cluster hierarchy of the Tumor samples. There is also a panel idicating the relationship of a feature with known UCSC CpG Islands in the human genome version 19. We also include a panel which tracks a sample across different clusters for diferent values of ", asEmph("k"), " from ", asEmph("k=2"), " to ", asEmph("k=7"), ". This plot is similar to the Sample Tracking Plot generated by ConsensusClusterPlus. The cluster membership is inferred from consensus clustering."))

	report <- addToMethods(report, method1, method2, method3, method4)
	
	# Setup the Figures and Tables for the Results Section
	figure1 <- newFigure(grep("global", figures.lowres, value=TRUE),
			     fileHighRes = grep("global", figures.highres, value=TRUE),
			     exportId = "FIGURE_1",
			     "Heatmap visualizing both Normal and Tumor samples. For a detailed explanation of the figure refer to the Methods and Data section")

	figure2 <- newFigure(grep("hypo", figures.lowres, value=TRUE),
			     fileHighRes = grep("hypo", figures.highres, value=TRUE),
			     exportId = "FIGURE_2",
			     "Heatmap visualizing both Normal and Tumor samples. For a detailed explanation of the figure refer to the Methods and Data section")

	figure3 <- newFigure(grep("hyper", figures.lowres, value=TRUE),
			     fileHighRes = grep("hyper", figures.highres, value=TRUE),
			     exportId = "FIGURE_3",
			     "Heatmap visualizing both Normal and Tumor samples. For a detailed explanation of the figure refer to the Methods and Data section")

	tableData1 <- read.table(grep("global", tables, value=T), stringsAsFactors=F, header=T)[1:5, ]
	tableData2 <- read.table(grep("hypo", tables, value=T), stringsAsFactors=F, header=T)[1:5, ]
	tableData3 <- read.table(grep("hyper", tables, value=T), stringsAsFactors=F, header=T)[1:5, ]

	table1 <- newTable(tableData1, file=grep("global", tables, value=T), exportId="TABLE_1",
			   "Table containing Sample cluster membership for different k from k=2 to k=7 from consensus clustering")
	table2 <- newTable(tableData2, file=grep("hypo", tables, value=T), exportId="TABLE_2",
			   "Table containing Sample cluster membership for different k from k=2 to k=7 from consensus clustering")
	table3 <- newTable(tableData3, file=grep("hyper", tables, value=T), exportId="TABLE_3",
			   "Table containing Sample cluster membership for different k from k=2 to k=7 from consensus clustering")

	report <- addToResults(report,
			       addTo(newSubSection(paste("Global Overview of", disease, "Methylation Data", sep=" ")), figure1),
			       addTo(newSubSection(paste("Hypomethylated Overview of", disease, "Methylation Data", sep=" ")), figure2),
			       addTo(newSubSection(paste("Hypermethylated Overview of", disease, "Methylation Data", sep=" ")), figure3),
			       addTo(newSubSection("Sample Tracking Information for Global Overview"), table1),
			       addTo(newSubSection("Sample Tracking Information for Hypomethylated Overview"), table2),
			       addTo(newSubSection("Sample Tracking Information for Hypermethylated Overview"), table3))

	writeReport(report, filename=file.path(outdir, "nozzle"))
}
					    
