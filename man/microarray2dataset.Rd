\name{microarray2dataset}
\alias{microarray2dataset}
\docType{data}
\title{
Table mapping probeset identifier prefixes to BioMart datasets.
}
\description{
The species corresponding to an probeset identifier can typically be identified
from the prefix of the identifier (e.g. Bt.457.1.S1_at corresponds to Bos taurus).
This table maps each known prefix to the corresponding species.
}
\usage{data(microarray2dataset)}
\format{
  A data frame with 105 observations on the following 5 variables.
  \describe{
    \item{\code{dataset}}{BioMart dataset.}
    \item{\code{microarray}}{BioMart microarray name for each species.}
    \item{\code{sample}}{An example probeset of the microarray.}
    \item{\code{prefix}}{The prefix pattern.}
    \item{\code{unique}}{Whether the prefix pattern is unique to the microarray.}
  }
}
\details{
All Agilent microarray share the same prefix pattern, making it very difficult to differentiate.
Many affymetrix microarrays also share the same prefix pattern.
}
\source{
A separate script was used to query the BioMart database and build this table.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(microarray2dataset)
## maybe str(microarray2dataset) ; plot(microarray2dataset) ...
}
\keyword{datasets}
