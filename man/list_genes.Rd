\name{list_genes}
\alias{list_genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns the genes associated with a Gene Ontology
}
\description{
Given a Gene Ontology (GO) identifier represented in the dataset,
returns a character vector listing the gene identifiers annotated to it.
}
\usage{
list_genes(result, go_id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
The output of \code{GO_anova()} or a subset of it obtained from
\code{subset_scores()}.
}
  \item{go_id}{
A Gene Ontology (GO) identifier represented in the dataset.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A character vector listing the gene identifiers of the genes associated
with the GO term.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--  or do  help(data=index)  for the standard data sets.
}

\keyword{ GO }
\keyword{ gene }
\keyword{ ontology }
\keyword{expression}
