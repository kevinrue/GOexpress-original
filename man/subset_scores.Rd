\name{subset_scores}
\alias{subset_scores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns a filtered GO_anova result object.
}
\description{
Given an output variable from a GO_anova analysis and a set of valid filters and thresholds,
returns an identically formatted list keeping only the rows of the score table passing the filters.
}
\usage{
subset_scores(result, gene_count=5, sig_count=1, ...)
}

\arguments{
  \item{result}{
The output of the \code{GO_anova()} function.
}
  \item{gene_count}{
Filter on the number of genes annotated to the GO term.
}
  \item{sig_count}{
Filter on the number of significant genes (FDR < 0.05) annotated to the GO term.
}
  \item{\dots}{
Pairs of filter and threshold values in the format "filter=threshold". 
Filters must be valid names from colnames(result$scores).
}
}
\details{
Currently, threshold are all considered as the minimal value allowed in the filtered table
for the corresponding filter. No maximal threshold has been implemented yet.
}
\value{
A list of two elements containing the results of the analysis.
  \item{scores }{A filtered table ranked in the same order as the original table.}
  \item{mapping }{The table mapping genes present in the dataset to GO terms.}
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
