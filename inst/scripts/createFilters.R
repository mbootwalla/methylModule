require(FDb.InfiniumMethylation.hg19)
require(GenomicRanges)
# Retrieve mapping information for the Illumina 450k chip based on hg19
hm450k <- get450k()

# Create the cg filter which retains the names of only CpG probes
filter.noncg <- names(hm450k[-which(values(hm450k)$probeType == "cg")])
filter.noncg <- filter.noncg[order(filter.noncg)]

# Create the sex filter which retains all the probes not on sex chromosomes
filter.sex <- names(hm450k[which(seqnames(hm450k) %in% c("chrX", "chrY"))])
filter.sex <- filter.sex[order(filter.sex)]

# Combine all filters to create the all filter
filter.all <- unique(c(filter.noncg, filter.sex))
filter.all <- filter.all[order(filter.all)]

# Store all the filters into a named list
filters <- list("filter.noncg" = filter.noncg, "filter.sex" = filter.sex, "filter.all" = filter.all)
