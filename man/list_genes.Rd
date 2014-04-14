\name{list_genes}
\alias{list_genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns the genes associated with a given Gene Ontology
}
\description{
Given a Gene Ontology (GO) identifier (go_id) represented in the dataset,
returns a character vector listing the ensembl_gene_ids associated with it.
}
\usage{
list_genes(result, go_id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
The output of \code{GO_anova()} or a subset of it obtained from \code{subset_scores()}.
}
  \item{go_id}{
A Gene Ontology (GO) identifier represented in the dataset.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A character vector listing the ensembl identifiers of the genes associated
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GO }
\keyword{ gene }
\keyword{ ontology }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
