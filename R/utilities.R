inferSampleType <- function(subjects)
{
	message("Index of tumor or normal samples not provided. Inferring from sample names \n")
	designation <- rep("TUMOR", length(subjects))
	designation[which(substr(subjects, 14, 15) == "11")] <- "NORMAL"
	designation[which(substr(subjects, 14, 15) == "20")] <- "CONTROL"
	cat("\n", "Table of number of Tumor, Normal and Control Samples found", "\n\n")
	print(table(designation))
	index.T <- which(designation == "TUMOR")
	index.N <- which(designation == "NORMAL")
	index.C <- which(designation == "CONTROL")
	retval <- list("Index.Tumor"=index.T, "Index.Normal"=index.N, "Index.Control"=index.C)
	return(retval)
}
