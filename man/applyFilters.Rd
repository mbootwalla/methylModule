\name{applyFilters}
\alias{applyFilters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Applies the specified set of filters to the data
}
\description{
The function filters out the data for the probes that are specified by the user
}
\usage{
applyFilters(x, filter = "all", custom.filter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
            A matrix
}
  \item{filter}{
            A character vector specifying the types of filters to apply. Can specify multiple filters. See details
}
  \item{custom.filter}{
            A character vector containing the names of the probes to be filtered
}
}
\details{
The function currently supports filtering of probes based only on the Illumina Methylation 450k platform. The following filters have been implemented:
\describe{
  \item{none}{No filtering is done}
  \item{noncg}{All non-CpG probes are filtered out. Non-CpG probes are probes that do not begin with \sQuote{cg}}
  \item{sex}{All probes that are associated with the sex chromosomes are filtered out}
  \item{all}{All the filters are applied to the data i.e. noncg, sex, snp and rpt}
  \item{custom}{This option is provided so that the user can apply their own customized filter instead of using the filters provided with the package. A character vector containing the names of the probes to be filtered out is expected}
}
The script that was used to create the standard filters that have been provided with this package can be found in \file{inst/scripts/createFilters.R}
}
\value{
A matrix
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Moiz Bootwalla
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
