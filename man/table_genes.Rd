\name{table_genes}
\alias{table_genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns a table listing the genes associated with a given Gene Ontology
}
\description{
Given a Gene Ontology (GO) identifier \code{go_id} and a result from an anovaGO analysis,
return an annotated table listing the genes associated with that go_id, their 
associated gene name and description.
}
\usage{
table_genes(result, go_id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
The output of the \code{GO_anova()} function.
}
  \item{go_id}{
A Gene Ontology (GO) identifier.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame with three columns for ensembl_gene_id, associated_gene_name,
and description.
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
\code{\link[anovaGO:GO_anova]{anovaGO::GO_anova}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GO }
\keyword{ gene }
\keyword{ ontology }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line