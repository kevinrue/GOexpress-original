\name{subset_scores}
\alias{subset_scores}
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
# subset_scores(result, total_count=5, namespace_1003="biological_process")
# synonym to 
# subset_scores(result, total=5, namespace="BP")
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
as the scoring function is biased for those GO terms (see Wiki, link below). 

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
\author{
Kevin Rue-Albrecht
}
\seealso{
  \href{https://github.com/kevinrue/GOexpress/wiki/Typical-analysis-pipeline}{Tutorial} on GitHub..
}
\examples{
# load the sample output data
data(raw_result)

# have an overview of the result variable
str(raw_result)

# filter for Biological Processes associated with 10+ genes
filtered_result = GOexpress::subset_scores(result=raw_result, total_count=5, namespace="BP")

# have an overview of the filtered result variable
str(filtered_result)
}

\keyword{ GOexpress }
