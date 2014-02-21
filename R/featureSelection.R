selectFeatures <- function(T, N=NULL, disease, type="hyper")
{
	if(missing(disease)) stop("Please provide a valid TCGA disease type\n")
	if(missing(T)) stop("No data for tumor samples provided\n")
	if(is.null(N)) {
		data(TTmappings)
		mappings <- mappings[which(map$diseaseabr == disease), ]
		# TODO: Add in logic for checking tissue specific normals
	}

	if(!(type %in% c("hyper", "hypo"))) stop("Invalid feature type\n\n")

	if(type == "hyper") {
		T.hyper <- T[pickFeatures(N, type="hyper"), , drop=FALSE]
		sd.hyper <- apply(T.hyper, 1, sd, na.rm=TRUE)
		probes.hyper <- names(sd.hyper)[sd.hyper > quantile(sd.hyper, .80, na.rm=TRUE)]
		set.seed(1234)
		probes.hyper <- sample(probes.hyper, size=5000)
		T.hyper <- T.hyper[probes.hyper, , drop=FALSE]
		N.hyper <- N[probes.hyper, , drop=FALSE]
		retval <- list("T.hyper"=T.hyper, "N.hyper"=N.hyper)
		return(retval)
	}

	if(type == "hypo") {
		T.hypo <- T[pickFeatures(N, type="hypo"), , drop=FALSE]
		sd.hypo <- apply(T.hypo, 1, sd, na.rm=TRUE)
		probes.hypo <- names(sd.hypo)[sd.hypo > quantile(sd.hypo, .80, na.rm=TRUE)]
		set.seed(1234)
		probes.hypo <- sample(probes.hypo, size=5000)
		T.hypo <- T.hypo[probes.hypo, , drop=FALSE]
		N.hypo <- N[probes.hypo, , drop=FALSE]
		retval <- list("T.hypo"=T.hypo, "N.hypo"=N.hypo)
		return(retval)
	}
}

pickFeatures <- function(x, type="hyper")
{
	medianN <- apply(x, 1, median, na.rm = T)
	if(type=="hyper")
	{
		probes.hyper <- names(subset(medianN, medianN < 0.2))
		return(probes.hyper)
	} else
	{
		probes.hypo <- names(subset(medianN, medianN > 0.7))
		return(probes.hypo)
	}
}

