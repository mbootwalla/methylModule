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

