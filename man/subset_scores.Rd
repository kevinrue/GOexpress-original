\name{subset_scores}
\alias{subset_scores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns a filtered list from GO_analyse results.
}
\description{
Given an output variable from a \code{GO_analyse} analysis and a set of valid filters
and thresholds, returns an identically formatted list keeping only the rows
of the score table passing all the filters.
}
\usage{
subset_scores(result, ...)

# Suggested use:
# subset_scores(result, total_count=5, namespace_1003="BP")
}
\arguments{
  \item{result}{
The output of the \code{GO_analyse()} function.
}
  \item{\dots}{
Additional pairs of filter and threshold values in the format "filter=threshold". 
Filters must be valid names from colnames(result$scores).
}
}
\details{
It is critical to filter out GO terms with very few genes (e.g. less than 5 genes),
as the scoring function is biased for those GO terms (see Wiki and FAQ). 

In addition, it is useful to retain only GO terms of one type (i.e. namespace) among 
the three possible: "biological process", "molecular function", and "cellular
component".

Suggested filters:
 \tabular{ll}{
 total_count: \tab Filter on the number of genes annotated to the GO term.\cr
 namespace_1003: \tab Filter keeping only the GO terms of a given type. Valid values are
 "biological_process", "molecular_function", and "cellular_component".
 Abbreviations "BP", "MF", and "CC" are also accepted.\cr
 }
 
Other filters:
 \tabular{ll}{
 data_count: \tab Filter on the number of genes in the dataset annotated to the GO term.\cr
 sig_count: \tab Filter on the number of genes in the dataset, annotated to the GO terms,
 and having a significant ANOVA result on the factor tested.\cr
 ave.F.score: \tab Average of the F.scores of all genes annotated to the GO term. Genes
 annotated to the GO term but absent from the dataset are assigned a F.value of 0, similarly
 to the genes with a non-significant ANOVA result on the factor of interest.\cr
 }

}
\value{
A list of four elements containing the results of the analysis.
  \item{ scores }{A table ranking the remaining GO terms after filtering.}
  \item{ mapping }{A table mapping the remaining GO terms to gene identifiers
  in the dataset.}
  \item{ anova }{The table ranking the remaining genes present in the dataset and
  mapped to the remaining GO terms.}
  \item{ factor }{The predefined grouping factor.}
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
