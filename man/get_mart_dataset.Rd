\name{get_mart_dataset}
\alias{get_mart_dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns a BioMart connection adapted to the expression data
}
\description{
This function returns a connection to the BioMart dataset of the species matching
the expression dataset. The first gene identifier in the dataset is used to guess the 
corresponding species from the prefix of the identifier.
The automatic guessing can be overriden by the user providing a non-empty BioMart dataset name.
}
\usage{
get_mart_dataset(expr_data, biomart_dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr_data}{
Gene expression matrix. Row names are ensembl gene identifiers. Column names are sample identifiers.
  }
  \item{biomart_dataset}{
A valid BioMart ensembl dataset name (e.g. "btaurus_gene_ensembl")
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A BioMart connection as returned by the getBM(...) method of package biomaRt.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Kevin Rue-Albrecht
}
\note{
%%  ~~further notes~~
}
\seealso{
Method \code{\link[biomaRt:getBM]{biomaRt::getBM}}.
}
\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
}
\keyword{ GO }
\keyword{ gene }
\keyword{ ontology }
\keyword{expression}


