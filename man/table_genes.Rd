\name{table_genes}
\alias{table_genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns a table listing the genes associated with a given Gene Ontology
}
\description{
Given a Gene Ontology (GO) identifier represented in the dataset and the output
variable of a \code{GO_analyse()} function, \code{table_genes()} returns a table
listing the genes associated with that go_id, their associated gene name,
and description.
}
\usage{
table_genes(result, go_id)
}
\arguments{
  \item{result}{
The output of the \code{GO_analyse()} function.
}
  \item{go_id}{
A Gene Ontology (GO) identifier.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame listing the statistics and BioMart annotations for the genes 
present in the expression dataset and associated with the GO term.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  Method \code{\link[GOexpress:GO_analyse]{GOexpress::GO_analyse}}.
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