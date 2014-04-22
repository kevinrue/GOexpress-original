\name{prefix2dataset}
\alias{prefix2dataset}
\docType{data}
\title{
Table mapping ensembl gene identifier prefixes to BioMart datasets.
}
\description{
The species corresponding to an ensembl gene identifier can typically be identified
from the prefix of the identifier (e.g. ENSBTAG corresponds to Bos taurus).
This table maps each known prefix to the corresponding species.
}
\usage{data(prefix2dataset)}
\format{
  A data frame with 66 rows.
  \describe{
    \item{\code{prefix}}{Unique prefix.}
    \item{\code{dataset}}{BioMart dataset.}
  }
}
\details{
C. elegans, D. melanogaster, and S. cerevisiae have atypical identifier pattern and prefixes.
Prefixes for those are included in the table, but checked independently of the other species.
}
\source{
A separate script was used to query the BioMart database and build this table.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(prefix2dataset)
str(prefix2dataset)
}
\keyword{datasets}
